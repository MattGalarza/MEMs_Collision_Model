# =============================================================================
# regression_vs_original.jl
# =============================================================================
# Proves that the restructured model is the SAME model as the first corrected
# file. Standard libraries only.
#
#   1. Keep a copy of the first file under another name, e.g.
#          collision_model_corrected_original.jl
#   2. julia regression_vs_original.jl [path/to/original.jl]
#
# Precomputed quantities must agree to 1e-14 (the arithmetic is unchanged);
# constitutive outputs to 1e-12; accelerations to 1e-12 of the summed force
# magnitudes (the force terms are now added in a different order, so net
# forces near equilibrium can differ by roundoff of the large terms).
# =============================================================================
using LinearAlgebra, Random, Test

const ORIGINAL_PATH = isempty(ARGS) ?
                      joinpath(@__DIR__, "collision_model_corrected_original.jl") : ARGS[1]
isfile(ORIGINAL_PATH) || error("Original model not found at $ORIGINAL_PATH")

module Original
include(Main.ORIGINAL_PATH)
end

module Refactored
include(joinpath(@__DIR__, "collision_model_corrected.jl"))
end

const Old = Original.CorrectedMEMS
const New = Refactored.CorrectedMEMS

"Random state covering free flight, approach, and contact on both sides."
function random_state(rng, gc)
    side = rand(rng, (-1.0, 1.0))
    x1 = side * (gc + (rand(rng) - 0.7) * 1e-6)              # gc - 0.7 um ... gc + 0.3 um
    x2 = side * min(abs(x1) + (rand(rng) - 0.9) * 0.2e-6, gc + 20e-9)
    if rand(rng) < 0.25                                       # some free-flight states
        x1 *= rand(rng)
        x2 = x1 * (1 + 1e-3 * randn(rng))
    end
    return [x1, x2, 0.02 * randn(rng), 0.02 * randn(rng), 0.5 * randn(rng), zeros(7)...]
end

@testset "restructured model == original model" begin
    parameter_sets = [
        (;),
        (; c1 = 1e-4, ce = 1e-7, seal_at_contact = false),
        (; ghs = 14.3e-6, Vbias = 5.0, gap_slope = 0.0525, film_scale = 0.5),
    ]
    rng = MersenneTwister(20260917)
    accel = t -> 4.95 * 9.80665 * sin(2 * pi * 20.0 * t)

    for kw in parameter_sets
        mo = Old.Model(Old.Params(; kw...); panels = 256)
        mn = New.Model(New.Params(; kw...); panels = 256)
        label = isempty(kw) ? "reference parameters" : string(kw)

        @testset "precomputed model, $label" begin
            for name in (:ke, :k1, :k3, :kss, :gc, :alpha, :hd, :kp)
                @test isapprox(getfield(mn, name), getfield(mo, name); rtol = 1e-14)
            end
            for name in (:M, :Minv, :beta, :Cstruct, :y, :weights, :B1, :B2)
                @test isapprox(getfield(mn, name), getfield(mo, name); rtol = 1e-14)
            end
        end

        @testset "pointwise physics, $label" begin
            for _ in 1:40
                u = random_state(rng, mn.gc)
                t = rand(rng) * 0.05
                x1, x2, v1, v2, Vout = u[1:5]

                co = Old.constitutive(mo, x1, x2)
                cn = New.constitutive(mn, x1, x2)
                @test isapprox(cn.C, co.C; rtol = 1e-12)
                @test isapprox(cn.grad, co.grad; rtol = 1e-12)
                @test isapprox(cn.D, co.D; rtol = 1e-12)

                so = Old.springs(mo, x1, x2)
                sn = New.springs(mn, x1, x2)
                @test isapprox(sn.U, so.U; rtol = 1e-12)
                @test isapprox(sn.F, so.F; rtol = 1e-12)

                wo = Old.wall(mo, x2, v2)
                wn = New.wall(mn, x2, v2)
                @test isapprox(wn.U, wo.U; rtol = 1e-12, atol = 1e-40)
                @test isapprox(wn.F, wo.F; rtol = 1e-12, atol = 1e-30)
                @test isapprox(wn.loss, wo.loss; rtol = 1e-12, atol = 1e-30)

                duo = zero(u)
                dun = zero(u)
                Old.rhs!(duo, u, (mo, accel), t)
                New.rhs!(dun, u, (mn, accel), t)
                Q = New.generalized_forces(mn, cn, x1, x2, v1, v2, mn.p.Vbias - Vout, accel(t))
                force_scale = sum(norm, (Q.spring, Q.wall, Q.electrostatic, Q.film, Q.structural, Q.base))
                @test dun[1:2] == duo[1:2]
                @test norm(mn.M * (dun[3:4] - duo[3:4])) <= 1e-12 * force_scale
                @test isapprox(dun[5], duo[5]; rtol = 1e-12, atol = 1e-12 * abs(Vout) / (mn.p.Rload * cn.C))
                @test isapprox(dun[6:12], duo[6:12]; rtol = 1e-12)

                @test isapprox(New.energy(mn, u), Old.energy(mo, u); rtol = 1e-13)
            end
        end
    end
end
println("Restructured model reproduces the original at every sampled state.")
