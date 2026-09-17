# =============================================================================
# mems_verify.jl -- constitutive/energy verification and legacy audits
# =============================================================================
# Included into module CorrectedMEMS. Standard libraries only.
#
#   verify()  runs (a) the 48 assertions of the first corrected file, unchanged
#             in content and tolerance, (b) 7 interface checks for the
#             restructured right-hand side, and (c) 33 checks of the selectable
#             squeeze-film closures. It then writes the audit CSVs that the
#             TeX document plots.
#
# Numerical consistency is not experimental validation.
# =============================================================================

"Write a header and an iterable of row tuples as a comma-separated file."
function writecsv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
    return path
end

# -----------------------------------------------------------------------------
# Independent audits of the EARLIER formulation. Both wedge functions use the
# old slope (wb - wt)/Leff on purpose, so that each audit isolates one defect.
# -----------------------------------------------------------------------------

"""
    wedge_reference(p, h)

Open-open translational film coefficient of a rigid wedge with tip gap `h`,
b_tr = 12 eta Tf (I2 - I1^2/I0) (TeX Eq. `btr`), by Gauss-Legendre quadrature
in a logarithmic variable. One gap, one beam [N s/m].
"""
function wedge_reference(p, h)
    a = (p.wb - p.wt) / p.Leff
    kp = 6 * p.slip_coefficient * p.mean_free_path
    function moment(k)
        function integrand(z)
            y = h * expm1(z) / a
            gap = h * exp(z)
            return y^k / (a * gap * (gap + kp))
        end
        return integrate(integrand, 0.0, log1p(a * p.Leff / h))
    end
    return 12 * p.eta * p.Tf * (moment(2) - moment(1)^2 / moment(0))
end

"""
    wedge_legacy(p, h; n = 8001)

Reproduces the earlier uniform-grid pressure integration (rectangle-rule
boundary moments mixed with trapezoidal pressure integration). Kept only to
document the ~75 % discrepancy at a 50 nm gap.
"""
function wedge_legacy(p, h; n = 8001)
    a = (p.wb - p.wt) / p.Leff
    kp = 6 * p.slip_coefficient * p.mean_free_path
    y = range(0.0, p.Leff; length = n)
    dy = step(y)
    G = (h .+ a .* y) .^ 2 .* (h .+ a .* y .+ kp)
    cc = 12 * p.eta * sum(y ./ G) / sum(1 ./ G)
    dp = (cc .- 12 * p.eta .* y) ./ G
    pressure = 0.0
    result = 0.0
    for j in 2:n
        pressure += (dp[j-1] + dp[j]) * dy / 2
        result += pressure * dy
    end
    return abs(p.Tf * result)
end

"""
    restitution(vin, cw)

Coefficient of restitution of an isolated, unsmoothed Hunt-Crossley impact
(TeX Eq. `restitution`): solve a - log(1 + a) = -b - log(1 - b), a = cw vin,
and return e = b / a. Independent of stiffness and exponent.
"""
function restitution(vin, cw)
    cw == 0 && return 1.0
    a = cw * vin
    target = a - log1p(a)
    lo = 0.0
    hi = 1 - eps()
    for _ in 1:100
        b = (lo + hi) / 2
        if -b - log1p(-b) > target
            hi = b
        else
            lo = b
        end
    end
    return (lo + hi) / (2 * a)
end

# -----------------------------------------------------------------------------
# Verification
# -----------------------------------------------------------------------------

"States (x1, x2) used for pointwise checks: rest, approach, both contact sides, overlap."
function verification_states(m::Model)
    return [(0.0, 0.0),
        (m.gc - 1e-6, m.gc - 1e-6),
        (m.gc + 0.1e-6, m.gc - 2e-9),
        (m.gc + 0.2e-6, m.gc + 5e-9),
        (-m.gc - 0.1e-6, -m.gc - 2e-9)]
end

"Central-difference gradient of C_t, relative error against the analytic gradient."
function capacitance_gradient_error(m::Model, x1, x2; dx = 1e-12)
    C(a, b) = constitutive(m, a, b; film = false).C
    fd = [(C(x1 + dx, x2) - C(x1 - dx, x2)) / (2 * dx),
        (C(x1, x2 + dx) - C(x1, x2 - dx)) / (2 * dx)]
    grad = constitutive(m, x1, x2; film = false).grad
    return norm(fd - grad) / max(norm(grad), 1e-18)
end

"""
Pointwise power balance: grad E(u) . du/dt from finite differences of `energy`
against the ledger rates, normalized by the throughput rate.
"""
function point_energy_error(m::Model, x1, x2)
    u = [x1, x2, 0.003, -0.005, 0.2, zeros(7)...]
    du = zero(u)
    rhs!(du, u, (m, t -> 2.0), 0.0)
    scales = [m.gc, m.gc, 0.02, 0.02, 3.0]
    g = zeros(5)
    for j in 1:5
        step_j = scales[j] * 1e-6
        up = copy(u)
        um = copy(u)
        up[j] += step_j
        um[j] -= step_j
        g[j] = (energy(m, up) - energy(m, um)) / (2 * step_j)
    end
    exact = du[6] + du[7] - sum(du[8:11])
    return abs(dot(g, du[1:5]) - exact) / max(du[12], 1e-25)
end

function verify(; outdir = joinpath(@__DIR__, "results"))
    mkpath(outdir)
    m = Model()
    p = m.p
    refined = Model(p; panels = 1024)
    beam_mass = p.rho * p.Tf * p.Lf * (p.wb + p.wt) / 2
    shape_mass_fraction = m.M[2, 2] / p.n_beams / beam_mass
    caperr = 0.0
    filmerr = 0.0
    energyerr = 0.0

    @testset "Corrected MEMS model" begin
        # ---- the 48 assertions of the first corrected file ------------------
        @testset "structural benchmarks" begin
            @test isapprox(m.ke, 22.5969569824; rtol = 1e-9)
            @test isapprox(m.k1, 3.2799375; rtol = 1e-12)
            @test isapprox(m.k3, 3.035714285714e8; rtol = 1e-11)
            @test isapprox(m.kss, 5.831; rtol = 1e-12)
            @test isapprox(m.gc, 13.71e-6; rtol = 1e-12)
            @test isapprox(shape_mass_fraction, 0.369820966146; rtol = 1e-9)
            @test isposdef(Symmetric(m.M))
            @test isapprox(sum(m.M), p.m_shuttle + p.n_beams * beam_mass; rtol = 1e-12)
            @test isapprox(beam_shape(p, m.ke, p.Lf), 1.0; atol = 1e-12)
            @test beam_shape(p, m.ke, 0.0) == 0.0
        end
        @testset "regularization and quadrature" begin
            @test softpos(-1.0, 1e-12) > 0
            @test dsoftpos(-1.0, 1e-12) > 0
            @test isapprox(sum(m.weights), p.Leff; rtol = 1e-13)
            @test maximum(abs.(cumulative_simpson(m, ones(length(m.y))) - m.y)) < 1e-15
        end
        @testset "pointwise constitutive and energy checks" begin
            for (x1, x2) in verification_states(m)
                c = constitutive(m, x1, x2)
                cr = constitutive(refined, x1, x2)
                @test c.C > p.cp
                @test eigmin(Symmetric(c.D)) >= -1e-14 * norm(c.D)
                de = norm(c.D - cr.D) / max(norm(cr.D), 1e-30)
                filmerr = max(filmerr, de)
                @test de < 1e-5
                er = capacitance_gradient_error(m, x1, x2)
                caperr = max(caperr, er)
                @test er < 2e-5
                er = point_energy_error(m, x1, x2)
                energyerr = max(energyerr, er)
                @test er < 5e-5
            end
        end
        @testset "contact loss, legacy audits, symmetry" begin
            @test wall(m, m.gc + 10e-9, -0.03).loss > 0      # unloading cutoff branch
            @test wall(m, m.gc + 10e-9, 0.03).loss > 0
            @test isapprox(wedge_reference(p, 50e-9), 1.15476443946e-4; rtol = 1e-8)
            @test isapprox(wedge_legacy(p, 50e-9), 2.84470582048e-5; rtol = 1e-8)
            @test isapprox(restitution(0.005, 50.0), 0.8568510175; rtol = 1e-8)
            a = constitutive(m, 1e-6, 2e-6)
            b = constitutive(m, -1e-6, -2e-6)
            @test isapprox(a.C, b.C; rtol = 1e-14)
            @test norm(a.grad + b.grad) < 1e-20
            @test norm(a.D - b.D) < 1e-15
            @test norm(constitutive(m, 0.0, 0.0).grad) < 1e-20
        end

        # ---- interface checks added with the restructuring -------------------
        @testset "right-hand-side interface" begin
            x1, x2 = m.gc + 0.1e-6, m.gc - 2e-9
            u5 = [x1, x2, 0.003, -0.005, 0.2]
            accel = t -> 2.0

            # (1) 5-, 12-, and 13-state evaluations share the physical part.
            du5 = zeros(5)
            du12 = zeros(12)
            du13 = zeros(13)
            rhs = HarvesterRHS(m, accel)
            rhs(du5, u5, nothing, 0.0)
            rhs(du12, [u5; zeros(7)], nothing, 0.0)
            rhs(du13, [u5; zeros(8)], nothing, 0.0)
            @test du5 == du12[1:5] == du13[1:5]
            @test du13[13] == 0.0

            # (2) A closure force enters the dynamics and its work enters the ledger.
            Qc = (3e-7, -2e-7)
            with_closure = HarvesterRHS(m, accel; closure = (u, t, theta, model) -> Qc)
            dq = zeros(13)
            with_closure(dq, [u5; zeros(8)], nothing, 0.0)
            @test isapprox(dq[3:4] - du13[3:4], m.Minv * collect(Qc); rtol = 1e-9)
            @test isapprox(dq[13], Qc[1] * u5[3] + Qc[2] * u5[4]; rtol = 1e-14)

            # (3) The force decomposition sums to the force actually integrated.
            c = constitutive(m, x1, x2)
            Q = generalized_forces(m, c, x1, x2, u5[3], u5[4], p.Vbias - u5[5], 2.0)
            @test isapprox(m.Minv * Q.total, du5[3:4]; rtol = 1e-12)

            # (4) State-dependent code is generic in the number type (what AD needs).
            ub = big.(u5)
            dub = zeros(BigFloat, 5)
            rhs(dub, ub, nothing, big(0.0))
            @test isapprox(Float64.(dub), du5; rtol = 1e-9)

            # (5) The scaled wrapper is the same vector field in scaled variables.
            scale = state_scales(m, 5)
            dz = zeros(5)
            ScaledRHS(rhs, scale)(dz, u5 ./ scale, nothing, 0.0)
            @test isapprox(dz .* scale, du5; rtol = 1e-13)
        end
    end

    verify_film_closures(m)

    open(joinpath(outdir, "verification_summary.txt"), "w") do io
        println(io, "Julia ", VERSION, "; all native verification assertions passed.")
        entries = ("ke_N_per_m" => m.ke, "k1_N_per_m" => m.k1, "k3_N_per_m3" => m.k3,
            "kss_N_per_m" => m.kss, "contact_travel_m" => m.gc,
            "shape_mass_fraction" => shape_mass_fraction,
            "max_cap_gradient_relative_error" => caperr,
            "max_film_refinement_relative_error" => filmerr,
            "max_point_energy_relative_error" => energyerr)
        for (name, value) in entries
            println(io, name, " = ", value)
        end
        println(io, "Numerical consistency is not experimental validation.")
    end

    write_constitutive_audit(m, joinpath(outdir, "constitutive_audit.csv"))
    write_pinned_tip_path(m, joinpath(outdir, "pinned_tip_path.csv"))
    println("Verification outputs written to ", outdir)
    return (; caperr, filmerr, energyerr)
end

"""
    verify_film_closures(reference)

Checks of the selectable squeeze-film closures (33 assertions). Analytic
limits: the narrow-strip formula for `:thickness`, and the classical
finite-rectangle factor 1 - (192/pi^5)(W/L) sum_{n odd} tanh(n pi L / 2W)/n^5
for `:modal` with a uniform gap. Cross-language anchors come from the verified
Python twin (same discretization, 512 panels, 8 modes). Physical ordering:
opening an additional vent path can only reduce dissipation, so
D_modal <= D_thickness and D_modal <= D_lengthwise in the semidefinite order.
"""
function verify_film_closures(reference::Model)
    p = reference.p
    rebuilt(; kw...) = Model(Params(; (name => getfield(p, name) for name in fieldnames(Params))..., kw...);
        panels = reference.panels)
    thick = rebuilt(film_model = :thickness)
    modal = rebuilt(film_model = :modal)
    modal_refined = Model(modal.p; panels = 2 * reference.panels)
    total(D) = sum(D)                               # rigid translation, 1' D 1

    @testset "squeeze-film closures" begin
        @testset "analytic limits (uniform gap)" begin
            gap = p.g0 - 2 * p.Tp
            G = gap^2 * (gap + reference.kp)
            strip = 2 * p.n_beams * p.eta * p.Tf^3 * p.Leff / G
            flat_thick = rebuilt(film_model = :thickness, gap_slope = 0.0)
            flat_modal = rebuilt(film_model = :modal, gap_slope = 0.0)
            @test isapprox(total(constitutive(flat_thick, 0.0, 0.0).D), strip; rtol = 1e-6)
            aspect = p.Tf / p.Leff
            plate = 1 - 192 / pi^5 * aspect * sum(tanh(n * pi / (2 * aspect)) / n^5 for n in 1:2:199)
            @test isapprox(total(constitutive(flat_modal, 0.0, 0.0).D) / strip, plate; rtol = 1e-4)
        end
        @testset "cross-language anchors (Python twin)" begin
            entries(D) = [D[1, 1], D[1, 2], D[2, 2]]
            approach = (reference.gc - 1e-6, reference.gc - 1e-6)
            contact = (reference.gc + 0.1e-6, reference.gc - 2e-9)
            @test all(isapprox.(entries(constitutive(thick, approach...).D),
                [2.288636691328e-06, 8.286935353846e-06, 1.585295558390e-04]; rtol = 1e-8))
            @test all(isapprox.(entries(constitutive(modal, approach...).D),
                [2.232329670748e-06, 7.575106203150e-06, 9.241694091461e-05]; rtol = 1e-8))
            @test all(isapprox.(entries(constitutive(thick, contact...).D),
                [5.252861111351e-06, 9.286072624765e-05, 1.590764808095e-02]; rtol = 1e-8))
            @test all(isapprox.(entries(constitutive(modal, contact...).D),
                [4.645929955768e-06, 4.633850155302e-05, 1.922584610404e-03]; rtol = 1e-8))
        end
        @testset "passivity, refinement, energy identity, vent ordering" begin
            for (x1, x2) in verification_states(reference)
                Dm = constitutive(modal, x1, x2).D
                @test eigmin(Symmetric(Dm)) >= 0
                Dr = constitutive(modal_refined, x1, x2).D
                @test norm(Dm - Dr) / norm(Dr) < 2e-4
                @test point_energy_error(modal, x1, x2) < 5e-5
                @test eigmin(Symmetric(constitutive(thick, x1, x2).D - Dm)) >= 0
                @test eigmin(Symmetric(constitutive(reference, x1, x2).D - Dm)) >= 0
            end
            a = constitutive(modal, 1e-6, 2e-6)
            b = constitutive(modal, -1e-6, -2e-6)
            @test norm(a.D - b.D) < 1e-12 * norm(a.D)
        end
        @testset "number-type genericity of the modal solve" begin
            x1, x2 = reference.gc + 0.1e-6, reference.gc - 2e-9
            Db = constitutive(modal, big(x1), big(x2)).D
            @test isapprox(Float64.(Db), constitutive(modal, x1, x2).D; rtol = 1e-9)
        end
    end
    return nothing
end

"""
Audit table over tip gaps: legacy vs consistent film coefficient, and the ratio
of the earlier lumped-dielectric force to the local-stack force. Old slope
(wb - wt)/Leff throughout, so the two effects are not mixed with the slope change.
"""
function write_constitutive_audit(m::Model, path)
    p = m.p
    a = (p.wb - p.wt) / p.Leff
    K = p.eps0 * p.Tf
    C_coating = p.eps0 * p.epsr * p.Leff * p.Tf / p.Tp
    rows = []
    for h in exp.(range(log(1e-9), log(20e-6); length = 70))
        C_air = K / a * log1p(a * p.Leff / h)
        C_lumped = 1 / (2 / C_coating + 1 / C_air)
        force_old = (C_lumped / C_air)^2 * K / a * (1 / h - 1 / (h + a * p.Leff))
        force_new = K / a * (1 / (h + m.hd) - 1 / (h + m.hd + a * p.Leff))
        push!(rows, (h * 1e9, wedge_legacy(p, h), wedge_reference(p, h), force_old / force_new))
    end
    return writecsv(path, ["gap_nm", "legacy_film", "corrected_film", "force_ratio"], rows)
end

"""
Prescribed constitutive path (not an ODE trajectory): x1 = x2 up to nominal
contact, then x2 = gc while x1 continues.
"""
function write_pinned_tip_path(m::Model, path)
    p = m.p
    rows = []
    for z in range(-1e-6, 1e-6; length = 201)
        x1 = m.gc + z
        x2 = z <= 0 ? x1 : m.gc
        v = z <= 0 ? [1.0, 1.0] : [1.0, 0.0]
        c = constitutive(m, x1, x2)
        push!(rows, (z * 1e6, c.C * 1e12, 0.5 * p.Vbias^2 * c.grad[1] * 1e6,
            0.5 * p.Vbias^2 * c.grad[2] * 1e6, dot(v, c.D * v)))
    end
    return writecsv(path, ["travel_um", "C_pF", "Q1_uN", "Q2_uN", "D_path"], rows)
end
