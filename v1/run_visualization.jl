# =============================================================================
# run_visualization.jl -- reference driven run and dense diagnostic exports
# =============================================================================
# Runs 8 cycles at 20 Hz, 4.95 g, Vbias from Params, then writes
#   results/overview.csv      whole run, 8001 samples
#   results/last_cycle.csv    final drive period, dense + every accepted step
#   results/contact.csv       -80 ... +320 microseconds around the first
#                             positive-wall contact entry of the final period
#   results/visualization_metrics.txt
# and finally the short contact probe. Column names are read by the figure
# TeX files; do not rename them.
# =============================================================================
using Printf, LinearAlgebra
using SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes
include(joinpath(@__DIR__, "collision_model_corrected.jl"))
using .CorrectedMEMS

const OUT = joinpath(@__DIR__, "results")
const FREQ = 20.0          # drive frequency [Hz]
const ACCEL_G = 4.95       # drive amplitude [g]
const CYCLES = 8

const HEADER = ["t_s", "t_ms", "tau_us", "x1_um", "x2_um", "tip_overlap_nm", "bend_nm",
    "v1_mm_s", "v2_mm_s", "Vout_mV", "C_pF", "PR_pW", "Qe2_uN", "Qwall2_uN",
    "Qbend2_uN", "Qfilm2_uN", "ER_pJ", "Wbase_nJ", "Wbias_nJ", "Dfilm_nJ",
    "Dwall_nJ", "Dstruct_nJ", "deltaE_nJ", "residual_aJ", "a_g"]

"One CSV row at time `t`; `tref` is the origin of the local time column tau."
function sample_row(m, physical, E0, t, tref)
    p = m.p
    u = physical(t)
    x1, x2, v1, v2, Vout = u[1:5]
    a_g = ACCEL_G * CorrectedMEMS.amplitude_ramp(t, FREQ) * sin(2 * pi * FREQ * t)
    c = constitutive(m, x1, x2)
    Q = generalized_forces(m, c, x1, x2, v1, v2, p.Vbias - Vout, a_g * 9.80665)
    return (t, t * 1e3, (t - tref) * 1e6, x1 * 1e6, x2 * 1e6, (abs(x2) - m.gc) * 1e9,
        (x2 - x1) * 1e9, v1 * 1e3, v2 * 1e3, Vout * 1e3, c.C * 1e12, Vout^2 / p.Rload * 1e12,
        Q.electrostatic[2] * 1e6, Q.wall[2] * 1e6, Q.spring[2] * 1e6, Q.film[2] * 1e6,
        u[8] * 1e12, u[6] * 1e9, u[7] * 1e9, u[9] * 1e9, u[11] * 1e9, u[10] * 1e9,
        (energy(m, u) - E0) * 1e9, ledger(m, u, E0) * 1e18, a_g)
end

"Time of the first upward crossing x2 = +gc in the final drive period (bisection on the dense output)."
function first_contact_entry(sol, scale, m, physical, tf, period)
    x2(j) = sol.u[j][2] * scale[2]
    entries = [j for j in 2:length(sol.t) if sol.t[j] >= tf - period && x2(j - 1) < m.gc && x2(j) >= m.gc]
    isempty(entries) && error("No positive contact entry in final cycle")
    j = first(entries)
    lo, hi = sol.t[j-1], sol.t[j]
    for _ in 1:50
        mid = (lo + hi) / 2
        if physical(mid)[2] > m.gc
            hi = mid
        else
            lo = mid
        end
    end
    return (lo + hi) / 2
end

function main()
    mkpath(OUT)
    println("Running reference model: $CYCLES cycles, $FREQ Hz, $ACCEL_G g.")
    flush(stdout)
    result = simulate(; kind = :drive, cycles = CYCLES, freq = FREQ,
        acceleration = ACCEL_G * 9.80665, outdir = OUT, tag = "visualization")
    sol, scale, m = result.sol, result.scale, result.model
    physical(t) = sol(t) .* scale
    period = 1 / FREQ
    tf = sol.t[end]
    E0 = energy(m, physical(0.0))

    tc = first_contact_entry(sol, scale, m, physical, tf, period)
    println("Contact entry at t = ", tc, " s. Exporting dense diagnostics.")
    flush(stdout)

    last_cycle_times = sort(unique(vcat(collect(range(tf - period, tf; length = 12001)),
        [t for t in sol.t if t >= tf - period])))
    exports = [("overview", collect(range(0.0, tf; length = 8001)), 0.0),
        ("last_cycle", last_cycle_times, tf - period),
        ("contact", collect(range(tc - 80e-6, tc + 320e-6; length = 4001)), tc)]
    for (name, times, tref) in exports
        CorrectedMEMS.writecsv(joinpath(OUT, name * ".csv"), HEADER,
            (sample_row(m, physical, E0, t, tref) for t in times))
    end

    # Extrema below are over accepted solver states, not optimized maxima.
    us = [z .* scale for z in sol.u]
    last_us = us[findall(>=(tf - period), sol.t)]
    mean_load_power = (us[end][8] - physical(tf - 4 * period)[8]) / (4 * period)
    open(joinpath(OUT, "visualization_metrics.txt"), "w") do io
        println(io, "Julia ", VERSION, "; successful native simulation; unchanged corrected model.")
        report = ["contact_entry_time_s" => tc,
            "last_four_cycles_mean_load_power_W" => mean_load_power,
            "last_cycle_sampled_max_abs_voltage_V" => maximum(abs(u[5]) for u in last_us),
            "last_cycle_sampled_max_tip_overlap_m" => maximum(abs(u[2]) - m.gc for u in last_us),
            "last_cycle_sampled_max_abs_beam_deflection_m" => maximum(abs(u[2] - u[1]) for u in last_us),
            "last_cycle_sampled_max_abs_shuttle_m" => maximum(abs(u[1]) for u in last_us),
            "total_film_loss_J" => us[end][9],
            "total_contact_loss_J" => us[end][11],
            "total_resistor_energy_J" => us[end][8],
            "energy_residual_over_throughput" => result.metrics["energy_residual_over_throughput"]]
        for (key, val) in report
            println(io, key, " = ", val)
        end
        println(io, "Peak quantities are sampled maxima, not independently optimized extrema.")
        println(io, "Finite reference run, not experimental validation or a stability proof.")
    end

    println("Running the short contact probe.")
    flush(stdout)
    simulate(; kind = :probe, outdir = OUT, tag = "probe")
    println("Data complete. Compile figure TeX files with pdflatex from this directory.")
    return nothing
end

main()
