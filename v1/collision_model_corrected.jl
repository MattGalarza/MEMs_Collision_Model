#!/usr/bin/env julia
# =============================================================================
# Corrected trapezoidal-electrode MEMS collision model (SI units throughout)
# Companion to Model_derivations_corrected.tex
# =============================================================================
#
# Layout (all files sit next to this one and form ONE module, CorrectedMEMS):
#
#   collision_model_corrected.jl   this file: module shell, exports, command line
#   mems_core.jl                   the model itself; LinearAlgebra only, AD-generic
#   mems_verify.jl                 assertions, legacy audits, audit CSVs (stdlib only)
#   mems_simulate.jl               audited ODE runs; loads the solver on first use
#
# Command line
#   julia collision_model_corrected.jl --verify
#   julia --project=. collision_model_corrected.jl --probe | --convergence | --drive
#
# From other code (figures, UDE, HDI)
#   include("collision_model_corrected.jl"); using .CorrectedMEMS 
#   m   = Model(Params(c1 = 1e-4))                     # rebuild after ANY parameter change
#   rhs = HarvesterRHS(m, t -> 0.0)                    # f(du, u, theta, t); 5, 12, or 13 states
#   c   = constitutive(m, x1, x2)                      # (; C, grad, D)
#   Q   = generalized_forces(m, c, x1, x2, v1, v2, Vc, a)   # every force term by name
#
# Adding a learned term (sketch; `net`, `st` from Lux, `theta` a ComponentArray):
#   s = reference_scales(m); zs = state_scales(m, 5)
#   closure = (u, t, theta, m) -> begin
#       out = first(net(u[1:5] ./ zs, theta, st))      # O(1) in, O(1) out
#       (s.F * out[1], s.F * out[2])                   # generalized force on (x1, x2) [N]
#   end
#   rhs  = ScaledRHS(HarvesterRHS(m, accel; closure), zs)
#   prob = ODEProblem(rhs, u0 ./ zs, tspan, theta)
# Integrate 13 states instead of 5 to audit the energy the learned term injects.
# =============================================================================
module CorrectedMEMS

using LinearAlgebra, Printf, Test

export Params, Model, constitutive, springs, wall, generalized_forces,
    HarvesterRHS, ScaledRHS, rhs!, energy, ledger,
    reference_scales, state_scales, observe_voltage, capacitor_charge,
    N_PHYSICAL, N_LEDGER, N_LEDGER_CLOSURE, STATE_NAMES,
    verify, simulate, convergence

include("mems_core.jl")
include("mems_verify.jl")
include("mems_simulate.jl")

function main(args = ARGS)
    mode = isempty(args) ? "--verify" : args[1]
    if mode == "--verify"
        verify()
    elseif mode == "--probe"
        simulate()
    elseif mode == "--convergence"
        convergence()
    elseif mode == "--drive"
        simulate(; kind = :drive, cycles = 30)
    elseif mode == "--help"
        println("--verify | --probe | --convergence | --drive\n",
            "See README.md for installation and parameter changes.")
    else
        error("Unknown mode: $mode")
    end
    return nothing
end

end # module CorrectedMEMS

if abspath(PROGRAM_FILE) == @__FILE__
    CorrectedMEMS.main()
end
