"""
Standalone MEMS contact diagnostic suite. In Julia:

    include("test.jl")
    result = test.main()

The original production model is preserved in baseline_reference.jl.
"""
module test
using LinearAlgebra, Printf, Dates, Statistics
using SHA, Serialization
import Pkg
include("baseline_reference.jl")
const BM = CorrectedMEMS
include("test_physics.jl")
include("test_static.jl")

function load_dependencies(;instantiate=true)
    Pkg.activate(@__DIR__)
    instantiate && Pkg.instantiate()
    @eval import SciMLBase, OrdinaryDiffEqRosenbrock, ADTypes, XLSX
    @eval using CairoMakie
    nothing
end

# Report and suite orchestration are included separately for readability.
include("test_reports.jl")

function main(;instantiate=true,kwargs...)
    load_dependencies(;instantiate)
    Base.invokelatest(run_suite;kwargs...)
end
end # module test

if abspath(PROGRAM_FILE)==@__FILE__
    test.main()
end
