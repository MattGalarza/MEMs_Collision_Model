#!/usr/bin/env julia
# Rebuild the two CSV views used by the dissertation chapter.
# Julia Base only; reads existing data and performs no physical simulation.

function chapter_transform_csv(source, destination, output_names, transform)
    isfile(source) || error("Missing source CSV: $source")
    written = 0
    open(source, "r") do input
        names = String.(split(chomp(readline(input)), ','))
        open(destination, "w") do output
            println(output, join(output_names, ','))
            for line in eachline(input)
                isempty(strip(line)) && continue
                values = parse.(Float64, split(chomp(line), ','))
                length(values) == length(names) || error("Invalid row in $source")
                row = Dict(zip(names, values)) 
                result = transform(row)
                result === nothing && continue
                length(result) == length(output_names) || error("Invalid output row")
                all(isfinite, result) || error("Nonfinite plotting value")
                println(output, join(result, ','))
                written += 1
            end
        end
    end
    written > 0 || error("No matching samples in $source")
    return written
end

function make_chapter_plot_data(;
    resultsdir=joinpath(@__DIR__, "results"),
    gc_um=13.71,
    late_start=0.3,
)
    nlate = chapter_transform_csv(
        joinpath(resultsdir, "drive_smoke_timeseries.csv"),
        joinpath(resultsdir, "chapter_drive_late.csv"),
        ["t_ms", "x1_um", "x2_um", "bending_nm", "Vout_mV", "ER_pJ"],
        row -> row["t_s"] < late_start ? nothing : (
            row["t_s"] * 1e3,
            row["x1_um"],
            row["x2_um"],
            (row["x2_um"] - row["x1_um"]) * 1e3,
            row["Vout_mV"],
            row["ER_J"] * 1e12,
        ),
    )

    first_time = Ref{Union{Nothing, Float64}}(nothing)
    nwindow = chapter_transform_csv(
        joinpath(resultsdir, "drive_smoke_contact_window.csv"),
        joinpath(resultsdir, "chapter_contact_window.csv"),
        ["tau_us", "clearance_nm", "Vout_mV"],
        row -> begin
            first_time[] === nothing && (first_time[] = row["t_s"])
            ((row["t_s"] - first_time[]) * 1e6,
             (gc_um - row["x2_um"]) * 1e3,
             row["Vout_mV"])
        end,
    )
    println("Wrote $nlate late-drive samples and $nwindow contact-window samples.")
    return (late_samples=nlate, window_samples=nwindow)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    make_chapter_plot_data()
end
