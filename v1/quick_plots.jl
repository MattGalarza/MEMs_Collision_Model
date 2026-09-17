# =============================================================================
# quick_plots.jl -- run the model and look at it, nothing else
# =============================================================================
# In the Julia REPL, from the folder that holds these files:
#
#     include("quick_plots.jl")
#     quickplot()                              # 3 V drive, 8 cycles, original film closure
#     quickplot(film_model = :modal)           # same run with the 2-D film closure
#     quickplot(kind = :probe)                 # the 600 microsecond contact probe
#     compare()                                # :lengthwise vs :modal, last drive cycle
#
# A window opens and a PNG is saved in results/. Needs Plots (once: ] add Plots)
# plus the solver packages already in this folder's project environment.
# =============================================================================
using Plots
include(joinpath(@__DIR__, "collision_model_corrected.jl"))
using .CorrectedMEMS

default(fontfamily = "Computer Modern", palette = :okabe_ito, linewidth = 1.2,
    framestyle = :box, grid = true, legend = :topright)

"Run one case and return physical time histories (rows = accepted solver steps)."
function run_case(; film_model = :lengthwise, kind = :drive, cycles = 8, Vbias = 3.0, freq = 20.0)
    p = Params(; film_model, Vbias)
    tag = "quick_$(kind)_$(film_model)_$(Vbias)V"
    r = simulate(; p, kind, cycles, freq, tag)
    U = Matrix(reduce(hcat, [u .* r.scale for u in r.sol.u])')
    return (; t = r.sol.t, U, gc = r.model.gc, period = 1 / freq, tag, kind, film_model, metrics = r.metrics)
end

"At most `nmax` evenly spaced indices out of `idx` (keeps plotting fast)."
thin(idx, nmax = 40_000) = idx[1:max(1, length(idx) ÷ nmax):end]

"Indices and time axis of the detail window: last drive period, or the whole probe."
function detail_window(c)
    if c.kind == :probe
        return eachindex(c.t), c.t .* 1e6, "time (us)"
    end
    idx = findall(>=(c.t[end] - c.period), c.t)
    return idx, (c.t[idx] .- c.t[idx[1]]) .* 1e3, "time in last cycle (ms)"
end

function panels(c)
    t, U, gc = c.t, c.U, c.gc
    all_idx = thin(collect(eachindex(t)))
    tt = c.kind == :probe ? t[all_idx] .* 1e6 : t[all_idx] .* 1e3
    tlabel = c.kind == :probe ? "time (us)" : "time (ms)"
    idx, td, dlabel = detail_window(c)
    idx_thin = thin(collect(1:length(idx)), 60_000)
    didx, td = idx[idx_thin], td[idx_thin]

    p1 = plot(tt, U[all_idx, 1] .* 1e6; label = "shuttle x1", xlabel = tlabel, ylabel = "position (um)")
    plot!(p1, tt, U[all_idx, 2] .* 1e6; label = "tip x2", linestyle = :dash)
    hline!(p1, [gc, -gc] .* 1e6; color = :gray, linestyle = :dot, label = "contact")

    p2 = plot(tt, U[all_idx, 5] .* 1e3; label = "", xlabel = tlabel, ylabel = "Vout (mV)")

    p3 = plot(td, (abs.(U[didx, 2]) .- gc) .* 1e9; label = "", xlabel = dlabel,
        ylabel = "tip overlap (nm)", ylims = (-400, 60))
    hline!(p3, [0.0]; color = :gray, linestyle = :dot, label = "")

    p4 = plot(td, U[didx, 5] .* 1e3; label = "", xlabel = dlabel, ylabel = "Vout (mV)")

    p5 = plot(tt, U[all_idx, 6] .* 1e9; label = "shaker work", xlabel = tlabel, ylabel = "energy (nJ)",
        legend = :topleft)
    plot!(p5, tt, U[all_idx, 9] .* 1e9; label = "film loss")
    plot!(p5, tt, U[all_idx, 11] .* 1e9; label = "contact loss")

    p6 = plot(tt, U[all_idx, 8] .* 1e12; label = "", xlabel = tlabel, ylabel = "load energy (pJ)")

    title = "$(c.kind), film = $(c.film_model)"
    if haskey(c.metrics, "last_four_cycles_mean_power_W")
        power_pW = round(c.metrics["last_four_cycles_mean_power_W"] * 1e12; sigdigits = 3)
        title = string(title, ", mean load power = ", power_pW, " pW")
    end
    return plot(p1, p2, p3, p4, p5, p6; layout = (3, 2), size = (1300, 1000), plot_title = title,
        left_margin = 6Plots.mm, bottom_margin = 4Plots.mm)
end

"Run, show, and save one case. Keywords: film_model, kind, cycles, Vbias, freq."
function quickplot(; kwargs...)
    c = run_case(; kwargs...)
    fig = panels(c)
    display(fig)
    path = joinpath(@__DIR__, "results", c.tag * ".png")
    savefig(fig, path)
    println("Saved ", path)
    return fig
end

"Overlay the last drive cycle of the two film closures."
function compare(; cycles = 8, Vbias = 3.0)
    cases = [run_case(; film_model, cycles, Vbias) for film_model in (:lengthwise, :modal)]
    pv = plot(; xlabel = "time in last cycle (ms)", ylabel = "Vout (mV)")
    px = plot(; xlabel = "time in last cycle (ms)", ylabel = "tip overlap (nm)", ylims = (-400, 60))
    for c in cases
        idx, td, _ = detail_window(c)
        keep = thin(collect(1:length(idx)), 60_000)
        plot!(pv, td[keep], c.U[idx[keep], 5] .* 1e3; label = string(c.film_model))
        plot!(px, td[keep], (abs.(c.U[idx[keep], 2]) .- c.gc) .* 1e9; label = string(c.film_model))
    end
    fig = plot(pv, px; layout = (2, 1), size = (1200, 800), left_margin = 6Plots.mm)
    display(fig)
    path = joinpath(@__DIR__, "results", "quick_compare_$(Vbias)V.png")
    savefig(fig, path)
    println("Saved ", path)
    return fig
end
