# =============================================================================
# make_figures.jl -- thin the dense CSVs and render the PGFPlots figures
# =============================================================================
# Run after run_visualization.jl. Requires pdflatex; pdftoppm is optional and
# only adds PNG previews.   julia make_figures.jl [--prepare-only]
# =============================================================================
using Printf

const ROOT = @__DIR__
const FIGURES = ["01_driven_overview", "02_contact_closeup", "03_last_cycle", "04_energy_and_probe"]

"""
    thin_table(input, output, columns; bins = 600, keepcols = nothing)

Reduce a dense time series to a plottable size without losing peaks: the time
axis is cut into `bins` intervals, and in each one the first and last samples
are kept together with the minimum and maximum (plus one neighbor on each
side) of every column named in `columns`. `keepcols` selects the columns that
are written; all of them by default.
"""
function thin_table(input, output, columns; bins = 600, keepcols = nothing)
    lines = readlines(input)
    header = split(first(lines), ",")
    data = permutedims(reduce(hcat, [parse.(Float64, split(s, ",")) for s in lines[2:end]]))
    t = data[:, 1]
    tracked = [findfirst(==(c), header) for c in columns]
    all(!isnothing, tracked) || error("A tracked column is missing from $input")

    selected = Set([1, length(t)])
    for b in 1:bins
        lo = searchsortedfirst(t, t[1] + (t[end] - t[1]) * (b - 1) / bins)
        hi = min(length(t), searchsortedlast(t, t[1] + (t[end] - t[1]) * b / bins))
        lo > hi && continue
        push!(selected, lo, hi)
        for j in tracked
            vals = @view data[lo:hi, j]
            for k in (lo + argmin(vals) - 1, lo + argmax(vals) - 1)
                for n in max(lo, k - 1):min(hi, k + 1)
                    push!(selected, n)
                end
            end
        end
    end

    rows = sort!(collect(selected))
    out = isnothing(keepcols) ? collect(eachindex(header)) :
          [findfirst(==(c), header) for c in keepcols]
    open(output, "w") do io
        println(io, join(header[out], ","))
        for i in rows
            println(io, join(data[i, out], ","))
        end
    end
    println(basename(output), ": ", length(t), " dense samples -> ", length(rows), " display samples")
    return nothing
end

function render(name)
    println("Rendering ", name)
    flush(stdout)
    open(joinpath("figures", name * "_build.log"), "w") do io
        run(pipeline(`pdflatex -interaction=nonstopmode -halt-on-error -output-directory=figures $(name * ".tex")`;
            stdout = io, stderr = io))
    end
    if Sys.which("pdftoppm") !== nothing
        run(`pdftoppm -singlefile -r 300 -png $(joinpath("figures", name * ".pdf")) $(joinpath("figures", name))`)
    else
        @warn "PDF created. Install Poppler/pdftoppm if PNG previews are also wanted."
    end
    return nothing
end

function main()
    cd(ROOT)
    mkpath("figures")
    Sys.which("pdflatex") === nothing &&
        error("Install a LaTeX distribution with pdfLaTeX/PGFPlots, or use the supplied figures.")
    thin_table("results/last_cycle.csv", "results/last_cycle_plot.csv",
        ["x1_um", "v1_mm_s", "Vout_mV", "PR_pW", "C_pF"];
        keepcols = ["t_s", "tau_us", "x1_um", "v1_mm_s", "Vout_mV", "PR_pW", "C_pF"])
    thin_table("results/overview.csv", "results/energy_plot.csv",
        ["deltaE_nJ", "Wbase_nJ", "Wbias_nJ", "Dfilm_nJ", "Dwall_nJ", "ER_pJ", "residual_aJ"];
        keepcols = ["t_s", "t_ms", "deltaE_nJ", "Wbase_nJ", "Wbias_nJ", "Dfilm_nJ", "Dwall_nJ",
            "ER_pJ", "residual_aJ"])
    "--prepare-only" in ARGS && return nothing
    foreach(render, FIGURES)
    return nothing
end

main()
