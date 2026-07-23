# reference_collision_model.jl
# Corrected reference implementation for the trapezoidal-electrode MEMS harvester
# collision model. Single C^1 vector field: no regime branches, no callbacks needed
# for the physics. Physics is identical to the Python-verified twin
# (verify_composite.py); this revision changes only the presentation layer:
# runs return a HarvesterRun whose REPL display is a ~12-line summary card
# (never the raw ODESolution dump), figures are saved to PNG rather than
# displayed, and zoom windows around contact episodes are sampled alias-free.
#
# Fixes relative to collision_model_v2.jl / test*.jl (all numerically verified):
#   F1  m2: removed stray Lf factor (was x2222 too light; 3.05 MHz vs 64.7 kHz).
#   F2  Fc: coupling spring -ke*(x1-x2) ALWAYS ON (v2's translational sign, the
#       correct one). Branch-switched contact force removed; boundary mechanics
#       handled by F4. Eliminates the +-9e-6 N force jump at |x2|=gp and both
#       flipped-branch instabilities.
#   F3  Fe: sign corrected to +(q^2/2C^2) dC/dx (attractive at constant charge).
#       Paper Eq. (26) and all five code files had it repulsive.
#   F4  Contact: unilateral Hunt-Crossley tip compliance on the electrode.
#       Required: with attractive Fe nothing else arrests pull-in.
#   F5  Gap floor hmin ~ 50 nm (mollified h~ = sqrt(h^2+hmin^2)) in C and film.
#       Removes the finite-time pull-in singularity at |x2|=gp; continuum
#       formulas lose validity below lambda ~ 70 nm anyway.
#   F6  Damping: corrected wedge closed form b(h) (== direct Reynolds to machine
#       precision). Paper Eq. (30)'s second term is spurious (x1794); v2's
#       transcription dimensionally broken; torsion-mirror rotational formula
#       (negative over the operating range) dropped. cdamp=1.0 physical; the
#       legacy 0.015 fudge is exposed as a parameter (it alone produces the
#       coherent ~65-68 kHz post-collision ringing; verified A/B).
#   F7  Electrical: q the only electrical state; Vout an observable.
#   F8  Solver: Rodas5P, state-scaled abstol, NO force_dtmin.
#
# Open modeling items:
#   T1  Post-contact film/capacitance from the tilted-wedge geometry with the
#       corrected conductance a^3 y^3 + 6 sigma_p lambda a^2 y^2 (journal rev.).
#   T2  cdamp: justify film suppression (venting/perforation/leakage) or accept
#       cushioned contacts.

module ReferenceCollisionModel

using DifferentialEquations   # matches your existing env; `using OrdinaryDiffEq`
using Printf                  # is sufficient if you add that package explicitly
using Plots

export build_params, rhs!, vout, run_pulse_demo, run_full_span,
       energy_ledger, HarvesterRun, save_figures, export_csv,
       sweep_gamma, mean_load_power

# ---------------------------------------------------------------- parameters --
function build_params(; g0=14e-6, Tp=120e-9, Tf=25e-6, wt=9e-6, wb=30e-6,
                        ws=14.7e-6, wss=14e-6, Leff=400e-6, Lf=450e-6,
                        Lsp=1400e-6, Lss=1000e-6, gss=14e-6,
                        m1=2.0933e-6, rho=2330.0, E=170e9,
                        e0=8.85e-12, ep=3.2, eta=1.849e-5,
                        N=160, cp=5e-12, Vbias=3.0, Rload=0.42e6,
                        cdamp=1.0, hmin=50e-9, kw=1e6, pw=1.5, cw=50.0,
                        aext = t -> 0.0)
    gp  = g0 - 2*Tp
    a   = (wb - wt)/Leff
    crl = e0*ep*Leff*Tf/Tp
    num = E*Tf*wt^2*(wb - wt)^3
    den = 6*Lf^3*((wb - 3*wt)*(wb - wt) + 2*wt^2*(log(Lf*wb) - log(Lf*wt)))
    ke  = num/den                                     # verified == direct integral
    modalcoeff = 0.236 + 0.045*(1.0 - wt/wb)
    m2 = modalcoeff * rho*Tf*Lf*(wb + wt)/2           # F1
    k1  = (2.0/3.0)*E*Tf*ws^3/Lsp^3
    k3  = (18.0/25.0)*E*Tf*ws/Lsp^3
    kss = E*Tf*wss^3/Lss^3
    Ne  = N/2
    events = Float64[]                                # guard-crossing times
    return (; g0, Tp, Tf, wt, wb, ws, wss, Leff, Lf, Lsp, Lss, gss,
              m1, rho, E, e0, ep, eta, N, cp, Vbias, Rload,
              cdamp, hmin, kw, pw, cw, aext,
              gp, a, crl, ke, m2, k1, k3, kss, Ne, events)
end

# ------------------------------------------------------------- constitutive --
@inline function cap(x2, p)                           # F5
    csum = 0.0
    dsum = 0.0
    for s in (1.0, -1.0)
        h0    = p.gp - s*x2
        ht    = sqrt(h0*h0 + p.hmin*p.hmin)
        dht   = -s*h0/ht
        Cair  = (p.e0*p.Tf/p.a)*log((ht + p.a*p.Leff)/ht)
        dCair = (p.e0*p.Tf/p.a)*(1/(ht + p.a*p.Leff) - 1/ht)*dht
        V     = 1/(2/p.crl + 1/Cair)
        csum += V
        dsum += (V/Cair)^2 * dCair
    end
    return p.Ne*csum, p.Ne*dsum
end

@inline bwedge(h, p) =                                # F6
    -(12*p.eta*p.Tf/p.a^3)*(2*p.a*p.Leff/(2h + p.a*p.Leff) + log(h/(h + p.a*p.Leff)))

@inline function forces(z, t, p)
    x1, v1, x2, v2, q = z
    ae = p.aext(t)
    Fs = -p.k1*x1 - 0.25*p.k3*x1^3
    Fss = abs(x1) <= p.gss ? 0.0 : -p.kss*(abs(x1) - p.gss)*sign(x1)
    Fc = -p.ke*(x1 - x2)                              # F2
    d  = abs(x2) - p.gp
    s  = sign(x2)
    Fw = d > 0 ? -s*p.kw*d^p.pw*max(1 + p.cw*v2*s, 0.0) : 0.0   # F4
    hr = sqrt((p.gp - x2)^2 + p.hmin^2)
    hl = sqrt((p.gp + x2)^2 + p.hmin^2)
    Fd = -p.cdamp*(bwedge(hr, p) + bwedge(hl, p))*v2  # F6
    Cv, dC = cap(x2, p)
    Ct = Cv + p.cp
    Fe = (q*q/(2*Ct*Ct))*dC/p.Ne                      # F3
    return ae, Fs + Fss, Fc, Fw, Fd, Fe, Ct
end

function rhs!(dz, z, p, t)
    ae, Fsp, Fc, Fw, Fd, Fe, Ct = forces(z, t, p)
    x1, v1, x2, v2, q = z
    dz[1] = v1
    dz[2] = (Fsp + p.Ne*Fc)/p.m1 - ae
    dz[3] = v2
    dz[4] = (-Fc + Fd + Fe + Fw)/p.m2 - ae
    dz[5] = (p.Vbias - q/Ct)/p.Rload                  # F7
    return nothing
end

vout(z, p) = p.Vbias - z[5]/(cap(z[3], p)[1] + p.cp)

function guard_callback(p)                            # diagnostic only
    cond(u, t, integ) = abs(u[3]) - p.gp
    aff!(integ) = (push!(p.events, integ.t); nothing)
    return ContinuousCallback(cond, aff!, aff!; save_positions = (false, false))
end

# ------------------------------------------------------------------- ledger --
function energy_ledger(sol, p)
    t = sol.t
    n = length(t)
    Etot = zeros(n); Pnet = zeros(n)
    for i in 1:n
        z = sol.u[i]
        x1, v1, x2, v2, q = z
        ae, _, _, _, _, _, Ct = forces(z, t[i], p)
        d  = max(abs(x2) - p.gp, 0.0)
        KE = 0.5*p.m1*v1^2 + p.Ne*0.5*p.m2*v2^2
        Vm = 0.5*p.k1*x1^2 + (0.25/4)*p.k3*x1^4 +
             (abs(x1) > p.gss ? 0.5*p.kss*(abs(x1) - p.gss)^2 : 0.0) +
             p.Ne*0.5*p.ke*(x1 - x2)^2 + p.Ne*(p.kw/(p.pw + 1))*d^(p.pw + 1)
        Etot[i] = KE + Vm + q^2/(2*Ct)
        qd  = (p.Vbias - q/Ct)/p.Rload
        hr  = sqrt((p.gp - x2)^2 + p.hmin^2)
        hl  = sqrt((p.gp + x2)^2 + p.hmin^2)
        s   = sign(x2)
        gate = (d > 0.0) && (1 + p.cw*v2*s > 0)
        Pw  = gate ? p.Ne*p.kw*d^p.pw*p.cw*v2^2 : 0.0
        Pnet[i] = -ae*(p.m1*v1 + p.Ne*p.m2*v2) + p.Vbias*qd - p.Rload*qd^2 -
                  p.Ne*p.cdamp*(bwedge(hr, p) + bwedge(hl, p))*v2^2 - Pw
    end
    work  = 0.0
    resid = 0.0
    for i in 2:n
        work += 0.5*(Pnet[i] + Pnet[i-1])*(t[i] - t[i-1])
        resid = max(resid, abs((Etot[i] - Etot[1]) - work))
    end
    thru = maximum(abs.(Etot)) + maximum(abs, Pnet)*(t[end] - t[1])
    return resid/thru
end

# ----------------------------------------------------- summaries and display --
struct RunSummary
    retcode::Any
    nsteps::Int
    t0::Float64
    t1::Float64
    nguard::Int
    episodes::Vector{Tuple{Float64,Float64}}
    duty::Float64
    maxpen::Float64
    ledger::Float64
    vpeak::Float64
    ctmin::Float64
    ctmax::Float64
end

struct HarvesterRun
    sol::Any
    p::Any
    summary::RunSummary
end

function contact_episodes(ev::Vector{Float64}; gap = 1e-3)
    isempty(ev) && return Tuple{Float64,Float64}[]
    evs = sort(ev)
    eps = Tuple{Float64,Float64}[]
    a = evs[1]; b = evs[1]
    for t in evs[2:end]
        if t - b > gap
            push!(eps, (a, b))
            a = t
        end
        b = t
    end
    push!(eps, (a, b))
    return eps
end

function summarize(sol, p; n_coarse = 40001, fine_dt = 2e-7)
    t0, t1 = sol.t[1], sol.t[end]
    tg = collect(range(t0, t1; length = n_coarse))
    U  = Array(sol(tg))
    Ct = [cap(x, p)[1] + p.cp for x in U[3, :]]
    V  = p.Vbias .- U[5, :] ./ Ct
    maxpen = maximum(abs.(U[3, :])) - p.gp
    vpeak  = maximum(abs.(V))
    # duty from the coarse grid over the FULL span (episode-count independent;
    # the earlier per-episode accumulation was silently capped at 8 episodes)
    duty = count(x -> abs(x) > p.gp, U[3, :]) / n_coarse
    eps = contact_episodes(p.events)
    for (i, (ta, tb)) in enumerate(eps)    # alias-free sharpening of peaks only
        i > 8 && break
        tz = collect(range(max(t0, ta - 5e-6), min(t1, tb + 5e-6); step = fine_dt))
        Uz = Array(sol(tz))
        pen = abs.(Uz[3, :]) .- p.gp
        maxpen = max(maxpen, maximum(pen))
        Ctz = [cap(x, p)[1] + p.cp for x in Uz[3, :]]
        vpeak = max(vpeak, maximum(abs.(p.Vbias .- Uz[5, :] ./ Ctz)))
    end
    return RunSummary(sol.retcode, length(sol.t), t0, t1, length(p.events),
                      eps, duty, maxpen,
                      energy_ledger(sol, p), vpeak,
                      minimum(Ct), maximum(Ct))
end

function Base.show(io::IO, ::MIME"text/plain", r::HarvesterRun)
    s = r.summary; p = r.p
    f2 = sqrt(p.ke/p.m2)/(2*pi)
    @printf(io, "HarvesterRun ──────────────────────────────────────────\n")
    @printf(io, "  retcode      %-12s steps        %d\n", string(s.retcode), s.nsteps)
    @printf(io, "  span         %.1f – %.1f μs   cdamp  %.3g   hmin  %.0f nm\n",
            s.t0*1e6, s.t1*1e6, p.cdamp, p.hmin*1e9)
    @printf(io, "  guard events %-4d  contact episodes %d   duty %.1f%%\n",
            s.nguard, length(s.episodes), 100*s.duty)
    if !isempty(s.episodes)
        (ta, tb) = s.episodes[1]
        @printf(io, "  first episode %.2f – %.2f μs   max penetration %.1f nm\n",
                ta*1e6, tb*1e6, s.maxpen*1e9)
    else
        @printf(io, "  no contact    (max |x2|-gp = %.1f nm)\n", s.maxpen*1e9)
    end
    @printf(io, "  ledger resid %.2e     peak |Vout| %.3f V\n", s.ledger, s.vpeak)
    @printf(io, "  Ct range     %.2f – %.2f pF    f2 %.1f kHz\n",
            s.ctmin*1e12, s.ctmax*1e12, f2/1e3)
    @printf(io, "  next: save_figures(run) → PNGs;  export_csv(run, \"run.csv\")\n")
end
Base.show(io::IO, r::HarvesterRun) = show(io, MIME"text/plain"(), r)

# ------------------------------------------------------------------ figures --
# Overview traces are decimated (fine for the slow states); episode zooms are
# sampled at zoom_dt = 0.2 μs, alias-free for the 64.7 kHz electrode mode.
function save_figures(run::HarvesterRun; dir = "figs", n_overview = 4000,
                      zoom_dt = 2e-7, max_zooms = 4)
    mkpath(dir)
    sol = run.sol; p = run.p
    t0, t1 = sol.t[1], sol.t[end]
    tg = collect(range(t0, t1; length = n_overview))
    U  = Array(sol(tg))
    Ct = [cap(x, p)[1] + p.cp for x in U[3, :]]
    V  = p.Vbias .- U[5, :] ./ Ct
    paths = String[]

    plt = plot(tg .* 1e3, U[1, :] .* 1e6; label = "x1", lw = 1,
               xlabel = "t (ms)", ylabel = "displacement (μm)",
               title = "shuttle and electrode (decimated overview)")
    plot!(plt, tg .* 1e3, U[3, :] .* 1e6; label = "x2", lw = 1)
    hline!(plt, [p.gp*1e6, -p.gp*1e6]; ls = :dash, lc = :gray, label = "±gp")
    push!(paths, joinpath(dir, "overview_displacements.png")); savefig(plt, paths[end])

    plt = plot(tg .* 1e3, V; label = "", lw = 1,
               xlabel = "t (ms)", ylabel = "Vout (V)",
               title = "output voltage (decimated overview)")
    push!(paths, joinpath(dir, "overview_vout.png")); savefig(plt, paths[end])

    for (i, (ta, tb)) in enumerate(run.summary.episodes)
        i > max_zooms && break
        tz = collect(range(max(t0, ta - 2e-5), min(t1, tb + 1.2e-4); step = zoom_dt))
        Uz = Array(sol(tz))
        rel = (Uz[1, :] .- Uz[3, :]) .* 1e9
        pen = (abs.(Uz[3, :]) .- p.gp) .* 1e9
        Ctz = [cap(x, p)[1] + p.cp for x in Uz[3, :]]
        Vz  = p.Vbias .- Uz[5, :] ./ Ctz
        p1 = plot(tz .* 1e6, rel; label = "", lw = 1,
                  ylabel = "x1 − x2 (nm)", title = @sprintf("episode %d (alias-free)", i))
        p2 = plot(tz .* 1e6, pen; label = "", lw = 1, ylabel = "|x2| − gp (nm)")
        hline!(p2, [0.0]; ls = :dash, lc = :gray, label = "")
        p3 = plot(tz .* 1e6, Vz; label = "", lw = 1,
                  xlabel = "t (μs)", ylabel = "Vout (V)")
        plt = plot(p1, p2, p3; layout = (3, 1), size = (900, 800))
        push!(paths, joinpath(dir, @sprintf("episode_%02d.png", i))); savefig(plt, paths[end])
    end
    return paths
end

# --------------------------------------------------------------------- csv --
function export_csv(run::HarvesterRun, path::AbstractString; n = 2000)
    tg = collect(range(run.sol.t[1], run.sol.t[end]; length = n))
    U  = Array(run.sol(tg))
    open(path, "w") do io
        println(io, "t,x1,v1,x2,v2,q,Vout")
        for (i, t) in enumerate(tg)
            Ct = cap(U[3, i], run.p)[1] + run.p.cp
            @printf(io, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                    t, U[1, i], U[2, i], U[3, i], U[4, i], U[5, i],
                    run.p.Vbias - U[5, i]/Ct)
        end
    end
    return path
end

# -------------------------------------------------------------- run helpers --
# Pulse acceptance test. Expected (Python-verified, cdamp=1.0): 2 guard events,
# one episode ~2–85 μs, max penetration ~22 nm, ledger ~1e-7.
function run_pulse_demo(; cdamp = 1.0, A = 2.1*9.81, Tpul = 150e-6, T = 600e-6)
    p  = build_params(cdamp = cdamp,
                      aext = t -> t < Tpul ? -A*sin(pi*t/Tpul) : 0.0)
    x2s = p.gp - 30e-9
    Cv0 = cap(x2s, p)[1]
    z0 = [x2s + 0.15e-6, 4e-3, x2s, 4e-3, p.Vbias*(Cv0 + p.cp)]
    prob = ODEProblem(rhs!, z0, (0.0, T), p)
    sol = solve(prob, Rodas5P();
                abstol = [1e-13, 1e-10, 1e-13, 1e-10, 1e-16],
                reltol = 1e-7, callback = guard_callback(p), maxiters = 10^7)
    return HarvesterRun(sol, p, summarize(sol, p))
end

# Full-span drive (v2-compatible sine ramp).
function run_full_span(; gamma = 2.1, f = 20.0, t_ramp = 0.2, T = 0.5,
                         cdamp = 1.0)
    A = gamma*9.81
    p = build_params(cdamp = cdamp,
                     aext = t -> A*min(t/t_ramp, 1.0)*sin(2*pi*f*t))
    Cv0 = cap(0.0, p)[1]
    z0 = [0.0, 0.0, 0.0, 0.0, p.Vbias*(Cv0 + p.cp)]
    prob = ODEProblem(rhs!, z0, (0.0, T), p)
    sol = solve(prob, Rodas5P();
                abstol = [1e-13, 1e-10, 1e-13, 1e-10, 1e-16],
                reltol = 1e-7, callback = guard_callback(p), maxiters = 10^8)
    return HarvesterRun(sol, p, summarize(sol, p))
end

# ------------------------------------------------------------------- sweep --
# Python-twin acceptance grid (verify_sweep.py: Radau, rtol 1e-6, T = 0.5 s,
# cdamp = 1.0, hmin = 50 nm). onset/duty/maxpen/pkV/Pload are the robust
# comparands; event counts can shift by a few for grazing taps (at gamma = 3.0
# they matched exactly across languages: 416 events, 25 episodes).
#  gamma  onset_ms  events  episodes  duty%  maxpen_nm  max_x1_um  pkV_mV  Pload_pW
#   1.00      --        0         0    0.00       0.0       6.10     0.01     0.000
#   1.50      --        0         0    0.00       0.0       8.74     0.02     0.001
#   2.00      --        0         0    0.00       0.0      11.09     0.05     0.003
#   2.10      --        0         0    0.00       0.0      11.54     0.06     0.004
#   2.50     213.9     48        12    0.46       4.3      13.78     2.45     0.264
#   3.00     185.4    416        25   13.07       7.6      13.80     4.03     0.408
#   3.50     160.3    372        27   24.19       8.7      13.81     4.86     0.521
#   4.00     136.0    390        29   31.15       9.8      13.82     5.58     0.585
function mean_load_power(run::HarvesterRun; tmin = 0.3, n = 20001)
    p = run.p
    lo = max(tmin, run.sol.t[1])
    tg = collect(range(lo, run.sol.t[end]; length = n))
    U  = Array(run.sol(tg))
    acc = 0.0
    for i in 1:n
        Ct = cap(U[3, i], p)[1] + p.cp
        qd = (p.Vbias - U[5, i]/Ct)/p.Rload
        acc += p.Rload*qd^2
    end
    return acc/n
end

function sweep_gamma(gammas = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
                     T = 0.5, cdamp = 1.0, csv = nothing)
    runs = HarvesterRun[]
    rows = String[]
    @printf(" gamma  onset_ms  events  episodes   duty%%  maxpen_nm  max_x1_um   pkV_mV  Pload_pW    steps\n")
    for g in gammas
        r = run_full_span(gamma = g, T = T, cdamp = cdamp)
        push!(runs, r)
        s = r.summary
        onset = isempty(r.p.events) ? "      -- " : @sprintf("%9.1f", minimum(r.p.events)*1e3)
        maxx1 = maximum(abs(u[1]) for u in r.sol.u)
        Pl = mean_load_power(r; tmin = 0.6*T)
        line = @sprintf("%6.2f %s %7d %9d %7.2f %10.1f %10.2f %8.2f %9.3f %8d",
                        g, onset, s.nguard, length(s.episodes), 100*s.duty,
                        max(s.maxpen, 0.0)*1e9, maxx1*1e6, s.vpeak*1e3,
                        Pl*1e12, s.nsteps)
        println(line)
        push!(rows, line)
    end
    if csv !== nothing
        open(csv, "w") do io
            println(io, "gamma,onset_ms,events,episodes,duty_pct,maxpen_nm,max_x1_um,pkV_mV,Pload_pW,steps")
            for line in rows
                println(io, join(split(line), ","))
            end
        end
    end
    return runs
end

end # module

# ------------------------------------------------------------------ driver --
# This footer executes whenever the file is run or include()d, so the file is
# self-driving: no `using` step needed. Qualified names (RCM.f) are deliberate:
# after editing and re-including you get "WARNING: replacing module ..." and any
# bare names previously imported with `using .ReferenceCollisionModel` keep
# pointing at the OLD module -- qualified calls always hit the fresh one.
# (For bare names, restart the REPL or use Revise.jl.)

RCM = ReferenceCollisionModel

res = RCM.run_pulse_demo()      # acceptance test: expect 2 guard events,
display(res)                    # one episode ~2-85 us, ~22 nm, ledger ~1e-7

# uncomment as needed:
RCM.save_figures(res; dir = "figs")            # PNG overviews + episode zooms
RCM.export_csv(res, "pulse_run.csv")           # decimated table
res_full = RCM.run_full_span(gamma = 3.0); display(res_full)
runs = RCM.sweep_gamma(csv = "sweep_julia.csv")  # acceptance grid; compare to
