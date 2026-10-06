# =============================================================================
#  hh_stage2_gradient.jl — Stage 2: gradient-based multi-start closure fits
#  and a misspecified-class test.
#
#  Hidden current in the data-generating model (set HIDDEN below):
#    :M    I = gM·w·(V−EK),            ẇ = (σ((V−Vh)/k) − w)/τw      (in-class for :struct)
#    :KCa  I = gKCa·c/(c+Kd)·(V−EK),   ċ = −c/τc + β·σ((V+20)/5)     (out-of-class: saturating
#                                                                       conductance, spike-driven latent)
#  Closures (all share HH(gNa,gK) and a scalar latent z):
#    :struct  g·z·(V−E), g = exp(θ) ≥ 0, z∞ = σ((V−Vh)/k), τ const          (7 params)
#    :poly    a0 + a1 V + a2 z + a3 Vz + a4 z², same latent kinetics, no sign constraint
#    :nn      softplus(net_g(z))·(V−E) with the SAME bounded sigmoid kinetics as :struct —
#             kinetics structured, conductance nonlinearity free (6 + 3NH+1 params).
#             (A fully free-kinetics variant died in a dead-latent minimum: z∞ → 0.)
#    :sat     g0·z/(z+Kz)·(V−E), same kinetics — the TRUE class for :KCa up to latent rescaling;
#             a reference that separates model-class failure from optimizer failure.
#  Extrapolation test: burst / hyperpolarise / release protocol
#    200 ms at DC = 30 (latent pushed past the training range), 300 ms at DC = −15
#    (quiet AHP phase: closure evaluated at decaying z and hyperpolarised V, where a
#    leak-like polynomial term is wrong), 200 ms at DC = 8 (recovery firing).
#  Optimizer: ForwardDiff gradients + L-BFGS, anneal K = 3 → 0. Multi-start over a
#  structured grid of kinetics inits (sharp spike-counter prior: Vh ∈ {−30, −10},
#  τ ∈ {50, 150}, k = 3) plus the generic warm start; best training loss wins.
#  Gaussian jitter around the plain-HH warm start mostly lands in the dead-latent
#  basin (train loss = plain-HH loss), which is the dominant failure mode.
#  Training: four 400-ms sweeps at DC = 4, 8, 12, 16; test: unseen sweep DC = 10.
#
#  Deps: Optim, ForwardDiff, CairoMakie, Random, Statistics, Printf
# =============================================================================
using Optim, ForwardDiff, CairoMakie, Random, Statistics, Printf, LinearAlgebra

const HIDDEN   = :KCa          # :M or :KCa
const NRESTART = 4
const NH       = 4             # hidden units per sub-network in :nn

const C = 1.0; const ENa, EK, EL = 50.0, -77.0, -54.4; const gL = 0.3
const TRUE = [120.0, 36.0]
const gM_T, Vh_T, k_T, τ_T = 1.5, -35.0, 10.0, 100.0
const gKCa_T, Kd_T, τc_T, β_T = 3.0, 0.5, 150.0, 0.05

αm(V) = abs(V+40) < 1e-6 ? one(V) : 0.1*(V+40)/(1-exp(-(V+40)/10))
βm(V) = 4*exp(-(V+65)/18)
αh(V) = 0.07*exp(-(V+65)/20)
βh(V) = 1/(1+exp(-(V+35)/10))
αn(V) = abs(V+55) < 1e-6 ? 0.1*one(V) : 0.01*(V+55)/(1-exp(-(V+55)/10))
βn(V) = 0.125*exp(-(V+65)/80)
σ(x) = 1/(1+exp(-x)); softplus(x) = x > 30 ? x : log1p(exp(x))
Vh_of(u) = -80 + 80*σ(u);  u_of_Vh(Vh) = log((Vh+80)/(-Vh))
k_of(u)  =   2 + 28*σ(u);  u_of_k(k)   = log((k-2)/(30-k))

# tiny 1-hidden-layer MLP  R -> R with NH units; params laid out [W1(NH) b1(NH) W2(NH) b2]
nparams_mlp() = 3NH + 1
@inline function mlp(q, off, x)           # q::AbstractVector, off = index of first param
    s = zero(x) * zero(eltype(q))
    @inbounds for i in 1:NH
        s += q[off+2NH+i-1] * tanh(q[off+i-1]*x + q[off+NH+i-1])
    end
    return s + q[off+3NH]
end
# :nn parameter layout: q = [uVh, uk, log τ, E, net_g(3NH+1)]
const OFF_G = 5
nparams(kind) = kind === :struct ? 5 : kind === :poly ? 8 : kind === :nn ? 4 + nparams_mlp() : kind === :sat ? 6 : 0
# :sat layout: q = [uVh, uk, log τ, E, log g0, log Kz]

@inline function Ic(kind, q, V, z)
    kind === :none   && return zero(V)
    kind === :M      && return gM_T*z*(V-EK)
    kind === :KCa    && return gKCa_T*z/(z+Kd_T)*(V-EK)
    kind === :struct && return exp(q[1])*z*(V-q[5])
    kind === :poly   && return q[4] + q[5]*V + q[6]*z + q[7]*V*z + q[8]*z^2
    kind === :sat    && return exp(q[5])*z/(z+exp(q[6]))*(V-q[4])
    return softplus(mlp(q, OFF_G, z))*(V-q[4])                       # :nn
end
@inline function zkin(kind, q, V, z)
    kind === :none && return zero(V)
    kind === :M    && return (σ((V-Vh_T)/k_T) - z)/τ_T
    kind === :KCa  && return -z/τc_T + β_T*σ((V+20)/5)
    uVh, uk, τ = kind === :struct ? (q[2], q[3], exp(q[4])) : (q[1], q[2], exp(q[3]))   # :poly, :nn, :sat share layout here
    return (σ((V-Vh_of(uVh))/k_of(uk)) - z)/τ
end
@inline function rhs(x, I, p, kind, q, K, Vm)
    V, m, h, n, z = x
    Iion = p[1]*m^3*h*(V-ENa) + p[2]*n^4*(V-EK) + gL*(V-EL) + Ic(kind, q, V, z)
    return ( (I-Iion)/C + K*(Vm-V),
             αm(V)*(1-m)-βm(V)*m, αh(V)*(1-h)-βh(V)*h, αn(V)*(1-n)-βn(V)*n,
             zkin(kind, q, V, z) )
end
# generic in the parameter element type so ForwardDiff duals propagate
function rk4g(x0, I, dt, p, kind=:none, q=Float64[]; K=0.0, Vmeas=nothing)
    # element type from an actual RHS evaluation, so Duals in p or q propagate
    Tp = typeof(sum(rhs(x0, I[1], p, kind, q, K, 0.0)) + 0.0)
    N = length(I); X = Matrix{Tp}(undef, N, 5); x = ntuple(j -> Tp(x0[j]), 5); X[1,:] .= x
    @inbounds for i in 1:N-1
        Vm = Vmeas === nothing ? 0.0 : Vmeas[i]
        k1 = rhs(x, I[i], p, kind, q, K, Vm)
        k2 = rhs(x .+ 0.5dt .* k1, I[i], p, kind, q, K, Vm)
        k3 = rhs(x .+ 0.5dt .* k2, I[i], p, kind, q, K, Vm)
        k4 = rhs(x .+ dt .* k3,    I[i], p, kind, q, K, Vm)
        x = x .+ (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4); X[i+1,:] .= x
    end
    return X
end
function stim(N, dt, seed; dc=10.0, amp=6.0, τ=3.0)
    r = MersenneTwister(seed); I = zeros(N)
    for i in 1:N-1; I[i+1] = I[i] + dt*(-I[i]/τ) + sqrt(2dt/τ)*randn(r); end
    return dc .+ amp .* I
end
spikes(V) = [i for i in 1:length(V)-1 if V[i+1] > 0 && V[i] <= 0]

# ------------------------------ twin data -------------------------------------
rng = MersenneTwister(7); dt = 0.025
x0 = (-65.0, 0.05, 0.6, 0.32, 0.0)
DCs = [4.0, 8.0, 12.0, 16.0]; T_tr = 400.0; N_tr = Int(T_tr/dt)
sweeps = map(enumerate(DCs)) do (j, dc)
    I = stim(N_tr, dt, 20+j; dc=dc); X = rk4g(x0, I, dt, TRUE, HIDDEN)
    (I=I, X=X, V=X[:,1] .+ randn(rng, N_tr))
end
T, N = 300.0, Int(300.0/dt); t = (0:N-1) .* dt
I_te = stim(N, dt, 12; dc=10.0); Xte = rk4g(x0, I_te, dt, TRUE, HIDDEN); Vte = Xte[:,1] .+ randn(rng, N)
I_ex = vcat(stim(Int(200/dt), dt, 31; dc=30.0), stim(Int(300/dt), dt, 32; dc=-15.0, amp=2.0), stim(Int(200/dt), dt, 33; dc=8.0))
N_ex = length(I_ex); t_ex = (0:N_ex-1) .* dt
Xex = rk4g(x0, I_ex, dt, TRUE, HIDDEN); Vex = Xex[:,1] .+ randn(rng, N_ex)
PH2 = Int(200/dt)+1 : Int(500/dt); PH3 = Int(500/dt)+1 : N_ex
println("hidden current: ", HIDDEN)
for (j, sw) in enumerate(sweeps)
    @printf("train sweep DC = %4.1f: %2d spikes, latent range %.3f–%.3f\n", DCs[j], length(spikes(sw.X[:,1])), extrema(sw.X[:,5])...)
end
@printf("test  sweep DC = 10.0: %2d spikes\n", length(spikes(Xte[:,1])))
@printf("protocol sweep: %2d burst spikes, %2d recovery spikes, latent max %.3f (training max %.3f), quiet-phase AHP min %.1f mV\n",
        length(spikes(Xex[1:PH2[1],1])), length(spikes(Xex[PH3,1])), maximum(Xex[:,5]), maximum(maximum(sw.X[:,5]) for sw in sweeps), minimum(Xex[PH2,1]))

# ------------------------------ fitting ---------------------------------------
splitθ(kind, θ) = kind === :none ? (θ, eltype(θ)[]) : (view(θ, 1:2), view(θ, 3:length(θ)))
train_loss(kind) = (θ, K) -> begin p, q = splitθ(kind, θ)
    mean(mean((rk4g(x0, sw.I, dt, p, kind, q; K=K, Vmeas=sw.V)[:,1] .- sw.V).^2) for sw in sweeps) end
heldout(kind, θ) = begin p, q = splitθ(kind, θ)
    mean((rk4g(x0, I_te, dt, p, kind, q)[:,1] .- Vte).^2) end
function extrap(kind, θ)          # protocol sweep: quiet-phase MSE (subthreshold → meaningful), recovery spikes, 1-ms match, sign check
    p, q = splitθ(kind, θ); X = rk4g(x0, I_ex, dt, p, kind, q)
    Ivals = [Ic(kind, q, X[i,1], X[i,5]) for i in 1:N_ex]
    neg = kind === :none ? 0.0 : mean((Ivals .* (X[:,1] .- EK)) .< -1e-9)
    s1 = spikes(Xex[:,1]); s2 = spikes(X[:,1])
    match = isempty(s2) ? 0.0 : mean(any(abs.(s2 .- i) .* dt .< 1.0) for i in s1)
    (quiet = mean((X[PH2,1] .- Vex[PH2]).^2), ahp = minimum(X[PH2,1]), nrec = length(spikes(X[PH3,1])),
     match = match, zmax = maximum(X[:,5]), neg = neg)
end
function report_extrap(kind, θ)
    e = extrap(kind, θ)
    @printf("             protocol: quiet-phase MSE = %7.2f  AHP min = %.1f (true %.1f)  recovery spikes = %2d (true %2d)  1-ms match = %.2f  z_max = %.3f  anti-dissip = %.3f\n",
            e.quiet, e.ahp, minimum(Xex[PH2,1]), e.nrec, length(spikes(Xex[PH3,1])), e.match, e.zmax, e.neg)
end
function spike_match(kind, θ; tol=1.0)
    p, q = splitθ(kind, θ); s1 = spikes(Xte[:,1]); s2 = spikes(rk4g(x0, I_te, dt, p, kind, q)[:,1])
    isempty(s2) ? 0.0 : mean(any(abs.(s2 .- i) .* dt .< tol) for i in s1)
end
function lbfgs(f, θ0; maxit=200)
    g!(G, θ) = (G .= ForwardDiff.gradient(f, θ))
    optimize(f, g!, θ0, LBFGS(), Optim.Options(iterations=maxit, g_tol=1e-6, show_trace=false))
end
function anneal_grad(lossK, θ0; maxit=200)
    r = lbfgs(θ -> lossK(θ, 3.0), θ0; maxit=maxit)
    lbfgs(θ -> lossK(θ, 0.0), Optim.minimizer(r); maxit=maxit)
end
# kinetics parameters (uVh, uk, log τ) live at indices kin_idx(kind) within θ
kin_idx(kind) = kind === :struct ? (4, 5, 6) : (3, 4, 5)
const KIN_GRID = [(Vh, 3.0, τ) for Vh in (-30.0, -10.0) for τ in (50.0, 150.0)]
function multistart(kind, θ0, scale; n=NRESTART)
    best = nothing; r0 = MersenneTwister(99); iv, ik, iτ = kin_idx(kind)
    starts = Any[("warm", θ0)]
    for (Vh, k, τ) in KIN_GRID
        θs = copy(θ0); θs[iv] = u_of_Vh(Vh); θs[ik] = u_of_k(k); θs[iτ] = log(τ)
        push!(starts, (@sprintf("Vh=%.0f τ=%.0f", Vh, τ), θs))
    end
    for j in 1:max(0, n - length(KIN_GRID))
        push!(starts, ("jitter $j", θ0 .+ scale .* randn(r0, length(θ0))))
    end
    for (label, θs) in starts
        r = try anneal_grad(train_loss(kind), θs) catch e; @warn "start $label failed: $e"; nothing end
        r === nothing && continue
        dead = abs(Optim.minimum(r) - PLAIN_LOSS) < 1e-2 ? "  [dead latent]" : ""
        @printf("    %-6s start %-14s: train = %9.3f%s\n", kind, label, Optim.minimum(r), dead)
        best = (best === nothing || Optim.minimum(r) < Optim.minimum(best)) ? r : best
    end
    return best
end

# ------------------- E2: plain HH (misspecified) warm start -------------------
rp = anneal_grad(train_loss(:none), 0.7 .* TRUE); θp = Optim.minimizer(rp); const PLAIN_LOSS = Optim.minimum(rp)
@printf("\nE2 plain HH: gNa = %.1f gK = %.1f  train = %.2f  held-out = %.2f  spike match = %.2f\n",
        θp..., Optim.minimum(rp), heldout(:none, θp), spike_match(:none, θp)); report_extrap(:none, θp)

# ------------------- E4': closures, gradient multi-start ----------------------
println("\nE4' closure fits (noise floor 1.0)")
results = Dict{Symbol,Any}()
# :struct
θ0 = vcat(θp, log(0.5), u_of_Vh(-40.0), u_of_k(8.0), log(50.0), -70.0)
rs = multistart(:struct, θ0, vcat(5.0, 5.0, 0.5, 0.5, 0.5, 0.5, 5.0)); θs = Optim.minimizer(rs); results[:struct] = θs
@printf("  struct : g = %.2f Vh = %.1f k = %.1f τ = %.0f E = %.1f | train = %.3f held-out = %.2f spike match = %.2f\n",
        exp(θs[3]), Vh_of(θs[4]), k_of(θs[5]), exp(θs[6]), θs[7], Optim.minimum(rs), heldout(:struct, θs), spike_match(:struct, θs)); report_extrap(:struct, θs)
# :poly
θ0 = vcat(θp, u_of_Vh(-40.0), u_of_k(8.0), log(50.0), zeros(5))
rq = multistart(:poly, θ0, vcat(5.0, 5.0, 0.5, 0.5, 0.5, 0.1 .* ones(5))); θq = Optim.minimizer(rq); results[:poly] = θq
@printf("  poly   : a = %s Vh = %.1f k = %.1f τ = %.0f | train = %.3f held-out = %.2f spike match = %.2f\n",
        string(round.(θq[6:end], digits=3)), Vh_of(θq[3]), k_of(θq[4]), exp(θq[5]), Optim.minimum(rq), heldout(:poly, θq), spike_match(:poly, θq)); report_extrap(:poly, θq)
# :sat  — true-class reference (saturating conductance), same kinetics
θ0 = vcat(θp, u_of_Vh(-40.0), u_of_k(8.0), log(50.0), -70.0, log(1.0), log(0.1))
rt = multistart(:sat, θ0, vcat(5.0, 5.0, 0.5, 0.5, 0.5, 5.0, 0.5, 0.5)); θt = Optim.minimizer(rt); results[:sat] = θt
@printf("  sat    : Vh = %.1f k = %.1f τ = %.0f E = %.1f g0 = %.2f Kz = %.4f | train = %.3f held-out = %.2f spike match = %.2f\n",
        Vh_of(θt[3]), k_of(θt[4]), exp(θt[5]), θt[6], exp(θt[7]), exp(θt[8]), Optim.minimum(rt), heldout(:sat, θt), spike_match(:sat, θt)); report_extrap(:sat, θt)
# :nn   — bounded sigmoid kinetics (as :struct) + MLP conductance g(z) = softplus(net(z)),
#          net initialised to a positive ramp (W1 = 10, W2 = 1) so z carries gradient from the start
netg0 = vcat(10.0 .* ones(NH) .+ randn(rng, NH), randn(rng, NH), ones(NH) .+ 0.3 .* randn(rng, NH), -2.0)
θ0 = vcat(θp, u_of_Vh(-40.0), u_of_k(8.0), log(50.0), -70.0, netg0)
rn = multistart(:nn, θ0, vcat(5.0, 5.0, 0.5, 0.5, 0.5, 5.0, 0.5 .* ones(nparams_mlp()))); θn = Optim.minimizer(rn); results[:nn] = θn
@printf("  nn     : Vh = %.1f k = %.1f τ = %.0f E = %.1f | train = %.3f held-out = %.2f spike match = %.2f\n",
        Vh_of(θn[3]), k_of(θn[4]), exp(θn[5]), θn[6], Optim.minimum(rn), heldout(:nn, θn), spike_match(:nn, θn)); report_extrap(:nn, θn)
println("  learned conductance g(z) (nn) vs truth:")
for z in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)
    @printf("    z = %.2f  g_nn = %.3f", z, softplus(mlp(θn[3:end], OFF_G, z)))
    HIDDEN === :KCa && @printf("   true gKCa·c/(c+Kd) = %.3f", gKCa_T*z/(z+Kd_T))
    HIDDEN === :M   && @printf("   true gM·z = %.3f", gM_T*z); println()
end
@printf("  truth  : held-out = %.2f\n", heldout(HIDDEN, TRUE)); report_extrap(HIDDEN, TRUE)

# ------------------------------- figure ---------------------------------------
okabe = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000"]
set_theme!(fonts=(; regular="CMU Serif", bold="CMU Serif Bold"), fontsize=14)
outdir = joinpath(@__DIR__, "RUN_stage2_$(HIDDEN)"); mkpath(outdir)
f = Figure(size=(720,720))
for (row, (Iw, Vw, tw, ttl)) in enumerate(((I_te, Vte, t, "held-out, DC = 10"), (I_ex, Vex, t_ex, "burst / hyperpolarise / release protocol")))
    ax = Axis(f[row,1], xlabel="t [ms]", ylabel="V [mV]", title="$(ttl), hidden = $(HIDDEN)")
    lines!(ax, tw, Vw, color=(:gray,0.5), label="measured")
    lines!(ax, tw, rk4g(x0, Iw, dt, θs[1:2], :struct, θs[3:end])[:,1], color=okabe[5], label="struct")
    lines!(ax, tw, rk4g(x0, Iw, dt, θq[1:2], :poly,   θq[3:end])[:,1], color=okabe[6], label="poly")
    lines!(ax, tw, rk4g(x0, Iw, dt, θt[1:2], :sat,    θt[3:end])[:,1], color=okabe[1], label="sat (true class)")
    lines!(ax, tw, rk4g(x0, Iw, dt, θn[1:2], :nn,     θn[3:end])[:,1], color=okabe[3], label="nn (g(z) free)")
    row == 1 && axislegend(ax, position=:rt)
end
save(joinpath(outdir, "fig_heldout.pdf"), f)
println("\nfigure written to ", outdir)