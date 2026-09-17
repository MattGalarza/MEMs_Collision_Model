# =============================================================================
# mems_simulate.jl -- audited ODE runs (probe, driven) and the convergence study
# =============================================================================
# Included into module CorrectedMEMS. The solver packages (SciMLBase,
# OrdinaryDiffEqRosenbrock, ADTypes) are imported on the first call to
# `simulate`, so `--verify` keeps working with a bare Julia installation.
# Because of that late import, `_simulate` must be entered through
# `Base.invokelatest`; call `simulate`, never `_simulate`, from user code.
#
# Every run is accepted only if (i) the solver reports success and reaches the
# final time, (ii) the energy-ledger residual is below 1e-5 of the power
# throughput, and (iii) an independent Gauss quadrature of Vout^2/R agrees
# with the integrated load-energy state to 5e-4.
# =============================================================================

function simulate(; kwargs...)
    @eval import SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes
    return Base.invokelatest(_simulate; kwargs...)
end

"Four-cycle cosine amplitude ramp, then unity."
amplitude_ramp(t, freq) = 0.5 * (1 - cos(pi * min(t * freq / 4, 1.0)))

"""
    _simulate(; kind = :probe, ...)

* `kind = :probe` : unforced 600 microsecond contact transient from a fixed
  initial state just before tip contact.
* `kind = :drive` : `cycles` periods of ramped sinusoidal base acceleration.

Optional `closure` (and its parameters `theta`) are forwarded to `HarvesterRHS`;
the run then carries the 13th state so that the ledger still closes.
`autodiff = nothing` selects the validated finite-difference Jacobian; pass
e.g. `ADTypes.AutoForwardDiff()` to use the type-generic core instead.
Returns `(; metrics, sol, scale, model)`; `sol` is in scaled variables,
physical state = `sol(t) .* scale`.
"""
function _simulate(; p = Params(), panels = 512, kind = :probe, cycles = 8, freq = 20.0,
        acceleration = 4.95 * 9.80665, reltol = 1e-7, abstol = 1e-10,
        outdir = joinpath(@__DIR__, "results"), tag = string(kind),
        closure = nothing, theta = nothing, autodiff = nothing)
    kind in (:probe, :drive) || throw(ArgumentError("kind must be :probe or :drive"))
    (cycles > 4 && freq > 0) || throw(ArgumentError("need cycles > 4 and freq > 0"))
    m = Model(p; panels)
    mkpath(outdir)

    # ---- problem definition ------------------------------------------------
    nstates = closure === nothing ? N_LEDGER : N_LEDGER_CLOSURE
    scale = state_scales(m, nstates)
    u0 = zeros(nstates)
    if kind == :probe
        u0[1:5] = [m.gc + 0.15e-6, m.gc - 30e-9, 0.004, 0.004, 0.0]
    end
    accel = kind == :probe ? (t -> 0.0) :
            (t -> acceleration * amplitude_ramp(t, freq) * sin(2 * pi * freq * t))
    tend = kind == :probe ? 600e-6 : cycles / freq
    rhs = ScaledRHS(HarvesterRHS(m, accel; closure), scale)
    prob = theta === nothing ? SciMLBase.ODEProblem(rhs, u0 ./ scale, (0.0, tend)) :
           SciMLBase.ODEProblem(rhs, u0 ./ scale, (0.0, tend), theta)

    # ---- integration ---------------------------------------------------------
    jacobian = autodiff === nothing ? ADTypes.AutoFiniteDiff() : autodiff
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P(autodiff = jacobian);
        reltol, abstol, dtmax = kind == :probe ? 1e-6 : 2e-5, maxiters = 10^7,
        save_everystep = true, dense = true)
    SciMLBase.successful_retcode(sol) || error("Solver failed: $(sol.retcode)")
    isapprox(sol.t[end], tend; rtol = 1e-12) || error("Incomplete integration")
    physical(t) = sol(t) .* scale

    # ---- acceptance checks -----------------------------------------------------
    E0 = energy(m, u0)
    final = physical(tend)
    Es = scale[6]
    maxres = maximum(abs(ledger(m, z .* scale, E0)) for z in sol.u)
    throughput = final[12]
    relledger = maxres / max(throughput, eps() * Es)
    ERquad = resistor_energy_by_quadrature(sol, scale, p.Rload)
    erquad = abs(ERquad - final[8]) / max(abs(final[8]), eps() * Es)
    relledger < 1e-5 || error("Energy balance failed: residual/throughput = $relledger")
    erquad < 5e-4 || error("Resistor quadrature failed: relative error = $erquad")

    # ---- outputs ---------------------------------------------------------------
    times = range(0.0, tend; length = kind == :probe ? 2401 : 8001)
    rows = []
    for t in times
        u = physical(t)
        push!(rows, (t, t * 1e6, u[1] * 1e6, u[2] * 1e6, u[5] * 1e3, u[8],
            ledger(m, u, E0), max(u[2] - m.gc, -u[2] - m.gc) * 1e9))
    end
    writecsv(joinpath(outdir, tag * "_timeseries.csv"),
        ["t_s", "t_us", "x1_um", "x2_um", "Vout_mV", "ER_J", "energy_residual_J", "tip_overlap_nm"],
        rows)

    metrics = Dict{String,Any}(
        "ER_J" => final[8],
        "ER_quadrature_J" => ERquad,
        "ER_quadrature_relative_error" => erquad,
        "max_energy_residual_J" => maxres,
        "energy_residual_over_throughput" => relledger,
        "saved_steps" => length(sol.t),
        "duration_s" => tend,
        "film_panels" => panels,
        "reltol" => reltol,
        "scaled_abstol" => abstol)
    if kind == :drive
        add_drive_metrics!(metrics, m, sol, scale, freq, tend, outdir, tag)
    end

    open(joinpath(outdir, tag * "_summary.txt"), "w") do io
        println(io, "Julia ", VERSION, "; solver retcode = ", sol.retcode)
        for k in sort(collect(keys(metrics)))
            println(io, k, " = ", metrics[k])
        end
        println(io, "Finite transient; neither a periodic-attractor proof nor experimental validation.")
    end
    println(tag, ": ", sol.retcode, "; energy residual / throughput = ", relledger)
    return (; metrics, sol, scale, model = m)
end

"Independent 8-point Gauss integration of Vout^2/R over every accepted step, using the dense output."
function resistor_energy_by_quadrature(sol, scale, Rload)
    gx, gw = gauss_legendre(8)
    total = 0.0
    for j in 2:length(sol.t)
        a, b = sol.t[j-1], sol.t[j]
        mid = (a + b) / 2
        half = (b - a) / 2
        total += half * sum(gw[k] * (sol(mid + half * gx[k])[5] * scale[5])^2 / Rload
                            for k in eachindex(gx))
    end
    return total
end

"""
Driven-run diagnostics. The recurrence error compares the last two stroboscopic
states only; it is not a period classifier. A finely sampled window is written
around the last accepted-step sign change of the tip-contact condition, if any.
"""
function add_drive_metrics!(metrics, m, sol, scale, freq, tend, outdir, tag)
    physical(t) = sol(t) .* scale
    period = 1 / freq
    u = physical(tend)
    previous = physical(tend - period)
    metrics["last_cycle_scaled_recurrence_error"] = norm((u[1:5] - previous[1:5]) ./ scale[1:5])
    metrics["last_four_cycles_mean_power_W"] = (u[8] - physical(tend - 4 * period)[8]) / (4 * period)

    overlap(j) = abs(sol.u[j][2] * scale[2]) - m.gc
    crossings = [j for j in 2:length(sol.t) if overlap(j) * overlap(j - 1) < 0]
    if !isempty(crossings)
        tcontact = sol.t[last(crossings)]
        window = range(max(0, tcontact - 100e-6), min(tend, tcontact + 300e-6); length = 2001)
        writecsv(joinpath(outdir, tag * "_contact_window.csv"),
            ["t_s", "x1_um", "x2_um", "Vout_mV"],
            [(t, physical(t)[1] * 1e6, physical(t)[2] * 1e6, physical(t)[5] * 1e3) for t in window])
    end
    return metrics
end

"Format a number as LaTeX mantissa x 10^{exponent}."
function latexnum(v)
    mantissa, exponent = split(@sprintf("%.3e", v), "e")
    return mantissa * "\\times10^{" * string(parse(Int, exponent)) * "}"
end

"""
    convergence()

Three probe runs: reference, tighter solver tolerances, and tighter tolerances
with doubled film panels. Requires the integrated load energy to change by less
than 1 % at each refinement, and writes the summary table read by the TeX file.
"""
function convergence(; outdir = joinpath(@__DIR__, "results"))
    a = simulate(; outdir, tag = "probe")
    b = simulate(; outdir, tag = "probe_refined", reltol = 2e-8, abstol = 2e-11)
    c = simulate(; outdir, tag = "probe_grid_refined", reltol = 2e-8, abstol = 2e-11, panels = 1024)
    d1 = abs(a.metrics["ER_J"] / b.metrics["ER_J"] - 1)
    d2 = abs(b.metrics["ER_J"] / c.metrics["ER_J"] - 1)
    max(d1, d2) < 0.01 || error("Probe load energy is not converged to 1 percent")

    open(joinpath(outdir, "probe_convergence.txt"), "w") do io
        println(io, "ER relative change, tighter tolerances = ", d1)
        println(io, "ER relative change, doubled film panels = ", d2)
    end
    open(joinpath(outdir, "probe_summary.tex"), "w") do io
        println(io, "\\begin{tabular}{lr}\\toprule Quantity & Value ", "\\\\", "\\midrule")
        table = [("Integrated load energy (J)", a.metrics["ER_J"]),
            ("Energy residual / power throughput", a.metrics["energy_residual_over_throughput"]),
            ("Independent resistor quadrature error", a.metrics["ER_quadrature_relative_error"]),
            ("Load-energy change: solver refinement", d1),
            ("Load-energy change: film refinement", d2)]
        for (label, v) in table
            println(io, label, " & \$", latexnum(v), "\$ ", "\\\\")
        end
        println(io, "\\bottomrule\\end{tabular}")
    end
    return (; tolerance_change = d1, grid_change = d2)
end
