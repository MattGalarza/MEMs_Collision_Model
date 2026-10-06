# =============================================================================
#  hh_stage3_allen.jl — Stage 3: real data. Voltage-only closure discovery on
#  Allen Cell Types whole-cell recordings (output of allen_pull.py).
#
#  Base model: Pospischil et al. 2008 cortical HH (I_Na, I_Kd with threshold
#  shift V_T, leak) WITHOUT the M-current. Fitted per cell by the nudging anneal.
#  Diagnostic: spike-triggered innovation (STA) of the base fit at K = 3.
#  Closures (shared scalar latent z, bounded sigmoid kinetics):
#    :struct  g·z·(V−E)                       linear conductance
#    :nn      softplus(net(z))·(V−E)           free sign-constrained conductance
#  If the cell is a regular-spiking pyramidal cell, the expectation is an
#  M-like outward current: E ≈ E_K, τ ~ 100–1000 ms, V_h ≈ −35 mV.  Allen's own
#  biophysical fits include Im, so this is checkable.
#
#  Units: per-area (mV, ms, μA/cm², mS/cm², μF/cm²). Injected pA → μA/cm² via
#  the area in the HDF5 (C_m / 1 μF/cm²); C is a free parameter to absorb error.
#
#  Deps: HDF5, Optim, ForwardDiff, CairoMakie, Statistics, Printf, Random
# =============================================================================
using HDF5, Optim, ForwardDiff, CairoMakie, Statistics, Printf, Random, LinearAlgebra

const FILE      = length(ARGS) > 0 ? ARGS[1] : "allen_485909730.h5"
const DT_TARGET = 0.05          # ms; sweeps are decimated to the nearest native multiple
const N_TRAIN   = 3             # long-square sweeps used for fitting
const K_FINAL   = 0.3           # last anneal stage; pure K = 0 trajectory matching is hopeless on real spikes
const NH        = 4

const ENa, EK = 50.0, -90.0
σ(x) = 1/(1+exp(-x)); softplus(x) = x > 30 ? x : log1p(exp(x))
Vh_of(u) = -80 + 80*σ(u);  u_of_Vh(Vh) = log((Vh+80)/(-Vh))
k_of(u)  =   2 + 28*σ(u);  u_of_k(k)   = log((k-2)/(30-k))
safediv(a, b) = abs(b) < 1e-9 ? a/(1e-9*sign(b+1e-30)) : a/b

# ----------------------------- base model -------------------------------------
# p = (gNa, gKd, gL, EL, VT, C)
@inline function rates(V, VT)
    x = V - VT
    am = -0.32*safediv(x-13, exp(-(x-13)/4)-1); bm = 0.28*safediv(x-40, exp((x-40)/5)-1)
    ah = 0.128*exp(-(x-17)/18);                  bh = 4/(1+exp(-(x-40)/5))
    an = -0.032*safediv(x-15, exp(-(x-15)/5)-1); bn = 0.5*exp(-(x-10)/40)
    return am, bm, ah, bh, an, bn
end
function gates_at(V, VT)   # steady-state gates for initial condition
    am,bm,ah,bh,an,bn = rates(V, VT); (am/(am+bm), ah/(ah+bh), an/(an+bn))
end

# ----------------------------- closures ---------------------------------------
nparams_mlp() = 3NH + 1
@inline function mlp(q, off, x)
    s = zero(x)*zero(eltype(q))
    @inbounds for i in 1:NH; s += q[off+2NH+i-1]*tanh(q[off+i-1]*x + q[off+NH+i-1]); end
    s + q[off+3NH]
end
const OFF_G = 5
# :struct q = [log g, uVh, uk, log τ, E]      :nn q = [uVh, uk, log τ, E, net(13)]
@inline function Ic(kind, q, V, z)
    kind === :none   && return zero(V)
    kind === :struct && return exp(q[1])*z*(V-q[5])
    return softplus(mlp(q, OFF_G, z))*(V-q[4])
end
@inline function zkin(kind, q, V, z)
    kind === :none && return zero(V)
    uVh, uk, τ = kind === :struct ? (q[2], q[3], exp(q[4])) : (q[1], q[2], exp(q[3]))
    (σ((V-Vh_of(uVh))/k_of(uk)) - z)/τ
end
@inline function rhs(x, I, p, kind, q, K, Vm)
    V, m, h, n, z = x
    gNa, gKd, gL, EL, VT, C = p
    am,bm,ah,bh,an,bn = rates(V, VT)
    Iion = gNa*m^3*h*(V-ENa) + gKd*n^4*(V-EK) + gL*(V-EL) + Ic(kind, q, V, z)
    ( (I-Iion)/C + K*(Vm-V), am*(1-m)-bm*m, ah*(1-h)-bh*h, an*(1-n)-bn*n, zkin(kind, q, V, z) )
end
function rk4g(x0, I, dt, p, kind=:none, q=Float64[]; K=0.0, Vmeas=nothing)
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
    X
end
spikes(V; thr=-10.0) = [i for i in 1:length(V)-1 if V[i+1] > thr && V[i] <= thr]

# ------------------------------ load data -------------------------------------
h = h5open(FILE, "r")
area = read_attribute(h, "area_cm2"); vrest = read_attribute(h, "vrest_mV")
@printf("cell %s  cre = %s  area = %.0f μm²  vrest = %.1f mV\n", read_attribute(h, "specimen_id"),
        read_attribute(h, "cre_line"), area*1e8, vrest)
sw_all = []
for name in keys(h)
    g = h[name]; nm = read_attribute(g, "stimulus_name"); fs = read_attribute(g, "fs_Hz")
    t = Float64.(read(g["t_ms"])); V = Float64.(read(g["V_mV"])); Ipa = Float64.(read(g["I_pA"]))
    dec = max(1, round(Int, DT_TARGET*fs/1000)); dt = dec*1000/fs
    t, V, Ipa = t[1:dec:end], V[1:dec:end], Ipa[1:dec:end]
    I = Ipa .* 1e-6 ./ area                                   # pA → μA/cm²
    on = findfirst(abs.(Ipa) .> 5.0); off = findlast(abs.(Ipa) .> 5.0)
    (on === nothing || off === nothing) && continue
    if startswith(nm, "Long Square")                           # crop: 100 ms pre, 300 ms post (AHP)
        a = max(1, on - round(Int, 100/dt)); b = min(length(V), off + round(Int, 300/dt))
        t, V, I = t[a:b] .- t[a], V[a:b], I[a:b]
    end
    push!(sw_all, (name=name, stim=nm, amp=read_attribute(g, "amplitude_pA"), dt=dt, t=t, V=V, I=I,
                   nspk=length(spikes(V))))
end
close(h)
ls = sort([s for s in sw_all if startswith(s.stim, "Long Square")], by=s -> s.amp)
noise = [s for s in sw_all if startswith(s.stim, "Noise")]
println("long-square sweeps (amp pA → spikes): ", join(["$(round(Int,s.amp))→$(s.nspk)" for s in ls], "  "))
println("noise sweeps: ", length(noise))
# training: the strongest subthreshold sweep + two spiking sweeps spread in amplitude; held-out: next spiking sweep
sub  = filter(s -> s.nspk == 0, ls); spk = filter(s -> s.nspk >= 3, ls)
length(spk) < 3 && error("need ≥3 spiking long squares")
train = vcat(isempty(sub) ? [] : [sub[end]], [spk[1], spk[end ÷ 2 + 1]])[1:min(end, N_TRAIN)]
test  = spk[end]
println("train: ", join(["$(round(Int,s.amp)) pA ($(s.nspk) spk)" for s in train], ", "), " | test: $(round(Int,test.amp)) pA ($(test.nspk) spk)")
V0 = mean(train[1].V[1:round(Int, 50/train[1].dt)])

# ------------------------------ fitting ---------------------------------------
function x0_for(p, kind)
    m,h_,n = gates_at(V0, p[5]); (V0, m, h_, n, 0.0)
end
splitθ(kind, θ) = kind === :none ? (θ, eltype(θ)[]) : (view(θ, 1:6), view(θ, 7:length(θ)))
function sweep_loss(sw, p, kind, q, K)
    X = rk4g(x0_for(p, kind), sw.I, sw.dt, p, kind, q; K=K, Vmeas=sw.V)
    mean((X[:,1] .- sw.V).^2)
end
train_loss(kind) = (θ, K) -> begin p, q = splitθ(kind, θ); mean(sweep_loss(sw, p, kind, q, K) for sw in train) end
function free_run(sw, kind, θ)  p, q = splitθ(kind, θ); rk4g(x0_for(p, kind), sw.I, sw.dt, p, kind, q)[:,1] end
function score(sw, kind, θ; tol=2.0)
    V = free_run(sw, kind, θ); s1 = spikes(sw.V); s2 = spikes(V)
    match = isempty(s2) || isempty(s1) ? 0.0 : mean(any(abs.(s2 .- i) .* sw.dt .< tol) for i in s1)
    sub = [i for i in 1:length(V) if sw.V[i] < -45]                      # subthreshold-only MSE
    (nspk=length(s2), ntrue=length(s1), match=match, submse=mean((V[sub] .- sw.V[sub]).^2))
end
function lbfgs(f, θ0; maxit=200)
    g!(G, θ) = (G .= ForwardDiff.gradient(f, θ))
    optimize(f, g!, θ0, LBFGS(), Optim.Options(iterations=maxit, g_tol=1e-6))
end
function anneal_grad(lossK, θ0; maxit=200, Ks=(3.0, 1.0, K_FINAL))
    θ = θ0; r = nothing
    for K in Ks; r = lbfgs(θ -> lossK(θ, K), θ; maxit=maxit); θ = Optim.minimizer(r); end
    r
end

# -------------------------- base model fit ------------------------------------
println("\nbase-model fit (Pospischil, no I_M)")
θb0 = [50.0, 5.0, 0.1, V0, -56.0, 1.0]
rb = anneal_grad(train_loss(:none), θb0; maxit=300); θb = Optim.minimizer(rb)
@printf("  gNa = %.1f  gKd = %.1f  gL = %.3f  EL = %.1f  VT = %.1f  C = %.2f   nudged train MSE (K=%.1f) = %.2f\n",
        θb..., K_FINAL, Optim.minimum(rb))
for sw in vcat(train, [test]); sc = score(sw, :none, θb)
    @printf("  free-run %4d pA: spikes %3d (true %3d)  2-ms match %.2f  subthreshold MSE %.2f\n", round(Int,sw.amp), sc.nspk, sc.ntrue, sc.match, sc.submse)
end

# --------------------- STA innovation diagnostic ------------------------------
function sta(kind, θ; K=3.0, pre=10.0, post=400.0)
    p, q = splitθ(kind, θ); acc = nothing; c = 0; τs = nothing
    for sw in train
        X = rk4g(x0_for(p, kind), sw.I, sw.dt, p, kind, q; K=K, Vmeas=sw.V); e = K .* (sw.V .- X[:,1])
        w0, w1 = round(Int, pre/sw.dt), round(Int, post/sw.dt)
        acc === nothing && (acc = zeros(w0+w1); τs = (-w0:w1-1) .* sw.dt)
        for i in spikes(sw.V)
            if i-w0 >= 1 && i+w1-1 <= length(e); acc .+= view(e, i-w0:i+w1-1); c += 1; end
        end
    end
    τs, acc ./ max(c,1), c
end
τs, S, nsta = sta(:none, θb)
win(a,b) = mean(S[(τs .>= a) .& (τs .< b)])
@printf("\nSTA of innovation (K=3, %d spikes): 5–20 ms %+.2f   20–100 ms %+.2f   100–300 ms %+.2f  μA/cm²\n",
        nsta, win(5,20), win(20,100), win(100,300))
println("  (sustained positive tail = missing outward current, e.g. I_M / I_AHP; negative = missing inward)")

# ----------------------------- closures ---------------------------------------
KIN_GRID = [(Vh, 5.0, τ) for Vh in (-45.0, -25.0) for τ in (100.0, 500.0)]
kin_idx(kind) = kind === :struct ? (8, 9, 10) : (7, 8, 9)
function multistart(kind, θ0; maxit=200)
    best = nothing; iv, ik, iτ = kin_idx(kind)
    starts = Any[("warm", θ0)]
    for (Vh,k,τ) in KIN_GRID
        θs = copy(θ0); θs[iv] = u_of_Vh(Vh); θs[ik] = u_of_k(k); θs[iτ] = log(τ); push!(starts, (@sprintf("Vh=%.0f τ=%.0f", Vh, τ), θs))
    end
    for (label, θs) in starts
        r = try anneal_grad(train_loss(kind), θs; maxit=maxit) catch e; @warn "start $label failed: $e"; nothing end
        r === nothing && continue
        @printf("    %-6s start %-14s: train = %9.3f%s\n", kind, label, Optim.minimum(r),
                abs(Optim.minimum(r)-Optim.minimum(rb)) < 1e-2 ? "  [dead latent]" : "")
        best = (best === nothing || Optim.minimum(r) < Optim.minimum(best)) ? r : best
    end
    best
end
println("\nclosure fits")
rng = MersenneTwister(1)
θ0 = vcat(θb, log(0.05), u_of_Vh(-35.0), u_of_k(10.0), log(300.0), -85.0)
rs = multistart(:struct, θ0); θs = Optim.minimizer(rs)
@printf("  struct : g = %.3f  Vh = %.1f  k = %.1f  τ = %.0f ms  E = %.1f   train = %.3f\n",
        exp(θs[7]), Vh_of(θs[8]), k_of(θs[9]), exp(θs[10]), θs[11], Optim.minimum(rs))
netg0 = vcat(10.0 .* ones(NH) .+ randn(rng, NH), randn(rng, NH), 0.1 .* ones(NH), -3.0)
θ0 = vcat(θb, u_of_Vh(-35.0), u_of_k(10.0), log(300.0), -85.0, netg0)
rn = multistart(:nn, θ0; maxit=600); θn = Optim.minimizer(rn)
@printf("  nn     : Vh = %.1f  k = %.1f  τ = %.0f ms  E = %.1f   train = %.3f\n",
        Vh_of(θn[7]), k_of(θn[8]), exp(θn[9]), θn[10], Optim.minimum(rn))
for z in (0.0, 0.1, 0.25, 0.5, 1.0); @printf("    g_nn(z=%.2f) = %.4f mS/cm²\n", z, softplus(mlp(θn[7:end], OFF_G, z))); end
println("\nfree-run scores (held-out and training sweeps)")
for (kind, θ, lab) in ((:none, θb, "base"), (:struct, θs, "struct"), (:nn, θn, "nn"))
    for sw in vcat([test], train); sc = score(sw, kind, θ)
        @printf("  %-6s %4d pA%s: spikes %3d (true %3d)  2-ms match %.2f  subthreshold MSE %.2f\n",
                lab, round(Int,sw.amp), sw === test ? " [test]" : "", sc.nspk, sc.ntrue, sc.match, sc.submse)
    end
end
for sw in noise[1:min(1,end)]; for (kind, θ, lab) in ((:none, θb, "base"), (:struct, θs, "struct"), (:nn, θn, "nn"))
    sc = score(sw, kind, θ); @printf("  %-6s noise sweep: spikes %3d (true %3d)  2-ms match %.2f  subthreshold MSE %.2f\n", lab, sc.nspk, sc.ntrue, sc.match, sc.submse)
end; end
# base-fit parameters after closure (they shift when the closure takes over adaptation)
@printf("\nbase params with nn closure: gNa = %.1f gKd = %.1f gL = %.3f EL = %.1f VT = %.1f C = %.2f\n", θn[1:6]...)

# ------------------------------- figures --------------------------------------
okabe = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000"]
set_theme!(fonts=(; regular="CMU Serif", bold="CMU Serif Bold"), fontsize=14)
outdir = joinpath(@__DIR__, "RUN_stage3_" * splitext(basename(FILE))[1]); mkpath(outdir)
f1 = Figure(size=(600,380)); ax = Axis(f1[1,1], xlabel="time from spike [ms]", ylabel="⟨K(Vmeas − V̂)⟩ [μA/cm²]", title="STA of innovation, base model")
lines!(ax, τs, S, color=okabe[6]); hlines!(ax, [0.0], color=:black, linestyle=:dash); save(joinpath(outdir,"fig_sta.pdf"), f1)
f2 = Figure(size=(760,760))
for (row, sw) in enumerate((test, train[end]))
    ax = Axis(f2[row,1], xlabel="t [ms]", ylabel="V [mV]", title=(sw === test ? "held-out " : "training ") * "$(round(Int,sw.amp)) pA")
    lines!(ax, sw.t, sw.V, color=(:gray,0.5), label="measured")
    lines!(ax, sw.t, free_run(sw, :none, θb), color=okabe[6], label="base")
    lines!(ax, sw.t, free_run(sw, :nn, θn), color=okabe[3], label="nn closure")
    row == 1 && axislegend(ax, position=:rt)
end
save(joinpath(outdir,"fig_freerun.pdf"), f2)
println("\nfigures written to ", outdir)