#!/usr/bin/env julia
#=
verify_squeeze_films.jl

Numerical and visual verification of
    "Squeeze films between doubly tapered plates closing to contact".

Part 1  Derivations. Every identity used in §2–§6 is evaluated at random parameters in 256-bit
        arithmetic: derivatives by finite differences, integrals by adaptive Gauss–Kronrod.
        A correct identity leaves a residual near 1e-60; a wrong one leaves O(1).
Part 2  Cases. The film equation  ∇·(G∇p) = −12μ  (unit closing speed, G = h²(h + k_p)) is solved
        by finite volumes (a line-by-line port of the solver behind the paper's tables) and by
        Galerkin P1 finite elements (rigorous lower bounds), and compared with the paper's closed
        forms and tables: classes A–D, corners, strip law, scaling, monotonicity, venting.
Part 3  Figures, written to ./figures as PDF (Okabe–Ito palette, Computer Modern).

Run:    julia verify_squeeze_films.jl          (installs Plots, LaTeXStrings, QuadGK,
                                                SpecialFunctions on first use)
        Set QUICK = true below for a ~2 minute pass on coarser grids.
Exit code 0 if every check passes.
=#

const QUICK = false

import Pkg
for pkg in ("Plots", "LaTeXStrings", "QuadGK", "SpecialFunctions")
    Base.find_package(pkg) === nothing && Pkg.add(pkg)
end
using LinearAlgebra, SparseArrays, Printf, Random
using QuadGK, SpecialFunctions, Plots, LaTeXStrings

include(joinpath(@__DIR__, "reference_tables.jl"))

# ───────────────────────────── constants and bookkeeping ─────────────────────────────
const MU    = 1.849e-5          # air, Pa s
const W     = 25e-6             # device-layer thickness (film width), m
const ALPHA = 0.0467            # lengthwise wedge slope
const C0    = -1.109275         # constant of the screened logarithm
const DSTAR = MU * W / ALPHA^3  # damping scale D★ (N s/m)
const LFILM = 64W               # computational overlap length (paper, Appendix A)
const OKABE = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
const FIGDIR = joinpath(@__DIR__, "figures")
mkpath(FIGDIR)
default(fontfamily = "Computer Modern", framestyle = :box, grid = false, linewidth = 1.6,
        legendfontsize = 8, guidefontsize = 10, tickfontsize = 8, titlefontsize = 10,
        size = (520, 360), dpi = 300)

const RESULTS = NamedTuple{(:part, :name, :ok, :detail), Tuple{String, String, Bool, String}}[]
const RESIDUALS = Pair{String, Float64}[]          # for the derivation-residual figure
function check!(part, name, ok::Bool, detail = "")
    push!(RESULTS, (part = part, name = name, ok = ok, detail = detail))
    @printf("  [%s] %-60s %s\n", ok ? "PASS" : "FAIL", name, detail)
    flush(stdout)
    return ok
end
relerr(a, b) = abs(a / b - 1)
section(title) = (println("\n", "─"^100, "\n", title, "\n", "─"^100); flush(stdout))

# ═════════════════════════════════ PART 1: DERIVATIONS ═════════════════════════════════
section("PART 1  Derivations, checked at random parameters in 256-bit arithmetic")
setprecision(BigFloat, 256)
const RNG = MersenneTwister(20261009)
rb(lo, hi) = BigFloat(lo) + (BigFloat(hi) - BigFloat(lo)) * BigFloat(rand(RNG))
const QTOL = big"1e-50"
bquad(f, a, b; rtol = QTOL) = quadgk(f, BigFloat(a), BigFloat(b); rtol = rtol)[1]
# central differences in BigFloat (truncation O(δ²), rounding ~1e-77/δ)
d1(f, x; δ = big"1e-25") = (f(x + δ) - f(x - δ)) / (2δ)

function identity_check!(name, residual, scale; tol = big"1e-40")
    r = Float64(abs(residual) / max(abs(scale), big"1e-300"))
    push!(RESIDUALS, name => max(r, 1e-80))
    check!("derivation", name, r < Float64(tol), @sprintf("relative residual %.1e", r))
end

# 1.1  Profile (2.14): ODE, slip conditions; flux law (2.15); dissipation (2.16)
let worst_bc = big(0.0), worst_ode = big(0.0), worst_q = big(0.0), worst_d = big(0.0)
    for _ in 1:12
        P = rb(-5, 5); μ = rb(0.5, 2); b = rb(0.01, 0.5); hb = rb(-1, 1); h = rb(0.2, 2); ha = hb + h
        Ua = rb(-1, 1); Ub = rb(-1, 1)
        u(z)  = P / (2μ) * ((z - hb) * (z - ha) - b * h) + Ub + (Ua - Ub) * (z - hb + b) / (h + 2b)
        du(z) = d1(u, z)
        worst_bc  = max(worst_bc, abs((u(hb) - Ub) - b * du(hb)), abs((u(ha) - Ua) + b * du(ha)))
        z0 = hb + rb(0, 1) * h
        worst_ode = max(worst_ode, abs(μ * d1(du, z0; δ = big"1e-20") - P) / abs(P))
        q  = bquad(u, hb, ha); qf = -h^2 * (h + 6b) / (12μ) * P + h * (Ua + Ub) / 2
        worst_q = max(worst_q, abs(q - qf) / (abs(qf) + 1))
        up(z) = P / (2μ) * ((z - hb) * (z - ha) - b * h)      # Poiseuille part, both walls fixed
        diss = bquad(z -> μ * (P / (2μ) * (2z - ha - hb))^2, hb, ha) + μ / b * (up(ha)^2 + up(hb)^2)
        worst_d = max(worst_d, abs(diss - h^2 * (h + 6b) * P^2 / (12μ)) / (h^2 * (h + 6b) * P^2 / (12μ)))
    end
    identity_check!("profile (2.14): slip conditions at both walls", worst_bc, 1)
    identity_check!("profile (2.14): μ ∂²u/∂z² = ∂p/∂x", worst_ode, 1; tol = big"1e-25")
    identity_check!("flux law (2.15), slip + sliding walls", worst_q, 1)
    identity_check!("local dissipation (2.16), viscous + slip", worst_d, 1)
end

# 1.2  Proposition 2.1: the Reynolds equation is the compatibility of the two kinematic conditions
let worst = big(0.0)
    μ = big"1.3"
    ha(x, t) = 1 + big"0.3" * sin(x) + big"0.1" * t * x
    hb(x, t) = big"0.2" * cos(x) - big"0.05" * t
    Ua(x, t) = big"0.4" * cos(2x); Ub(x, t) = big"0.1" * x
    pp(x, t) = sin(big"1.3" * x) + t * x^2
    bs(x) = big"0.1" * (1 + sin(x) / 2)                       # slip length varying along x
    function u(x, z, t)
        h = ha(x, t) - hb(x, t); b = bs(x); px = d1(ξ -> pp(ξ, t), x)
        px / (2μ) * ((z - hb(x, t)) * (z - ha(x, t)) - b * h) + Ub(x, t) +
            (Ua(x, t) - Ub(x, t)) * (z - hb(x, t) + b) / (h + 2b)
    end
    for _ in 1:4
        x = rb(-1, 1); t = rb(0, 1)
        wb = d1(τ -> hb(x, τ), t) + u(x, hb(x, t), t) * d1(ξ -> hb(ξ, t), x)        # lower kinematic condition
        dudx_int = bquad(z -> d1(ξ -> u(ξ, z, t), x; δ = big"1e-18"), hb(x, t), ha(x, t); rtol = big"1e-25")
        wa = wb - dudx_int                                                          # continuity across the gap
        R_up = wa - d1(τ -> ha(x, τ), t) - u(x, ha(x, t), t) * d1(ξ -> ha(ξ, t), x)
        G(ξ) = (ha(ξ, t) - hb(ξ, t))^2 * ((ha(ξ, t) - hb(ξ, t)) + 6bs(ξ))
        flux(ξ) = G(ξ) / (12μ) * d1(s -> pp(s, t), ξ; δ = big"1e-18")
        hU(ξ) = (ha(ξ, t) - hb(ξ, t)) * (Ua(ξ, t) + Ub(ξ, t)) / 2
        R_rey = d1(flux, x; δ = big"1e-12") - d1(τ -> ha(x, τ) - hb(x, τ), t) - d1(hU, x)
        worst = max(worst, abs(R_up - R_rey) / (abs(R_rey) + 1))
    end
    identity_check!("Prop 2.1: upper kinematic residual = Reynolds residual", worst, 1; tol = big"1e-15")
end

# 1.3  Womersley mobility → h²(h + 6b) as λ → 0
let h = big"1.7", b = big"0.23", λ = big"1e-12"
    G = 3h^3 / λ^2 * (1 - tanh(λ) / (λ * (1 + (2b / h) * λ * tanh(λ))))
    identity_check!("Womersley mobility at λ = 1e-12 vs h²(h+6b)", G - h^2 * (h + 6b), h^3; tol = big"1e-20")
end

# 1.4  Corner ODE (Appendix B): L[f] = (cos³ψ f')' − 2cos³ψ f
let worst = big(0.0)
    L(f, ψ) = d1(s -> cos(s)^3 * d1(f, s; δ = big"1e-22"), ψ; δ = big"1e-18") - 2cos(ψ)^3 * f(ψ)
    for _ in 1:6
        ψ = rb(-1.2, 1.2)
        worst = max(worst, abs(L(s -> 1 / cos(s), ψ) + 1), abs(L(s -> sin(s) / cos(s)^2, ψ)), abs(L(s -> 1 / cos(s)^2, ψ)))
    end
    identity_check!("corner ODE: L[1/cos] = −1, L[sin/cos²] = L[1/cos²] = 0", worst, 1; tol = big"1e-20")
end

# 1.5  Proposition 5.1: the closed-form corner pressures solve ∇·(h³∇p) = −12μ with their edge conditions
function corner_closed(α, β, x, y; edge0 = :sealed, μ = MU)
    R = hypot(α, β); A = 12μ / R^3; h = α * x + β * y
    if edge0 == :sealed
        c1 = A * α * β * (α - β) / ((α + β) * (R^2 + α * β)); c2 = -(A * α * β + c1 * (R^2 + α^2)) / (2α * R)
    else
        c1 = A * α * (2α - β) / (R^2 + α^2 + 2α * β); c2 = (c1 * β - A * α) / R
    end
    return 12μ / (R^2 * h) + c1 * R * (α * y - β * x) / h^2 + c2 * R^2 * hypot(x, y) / h^2
end
let worst_pde = big(0.0), worst_tip = big(0.0), worst_face = big(0.0), worst_open = big(0.0)
    for edge0 in (:sealed, :open), _ in 1:3
        α = rb(0.02, 0.1); β = rb(0.01, 0.1); μ = big"1.0"
        p(x, y) = corner_closed(α, β, x, y; edge0 = edge0, μ = μ)
        x = rb(0.1, 1); y = rb(0.1, 1)
        Fx(ξ, η) = (α * ξ + β * η)^3 * d1(s -> p(s, η), ξ; δ = big"1e-22")
        Fy(ξ, η) = (α * ξ + β * η)^3 * d1(s -> p(ξ, s), η; δ = big"1e-22")
        div = d1(s -> Fx(s, y), x; δ = big"1e-18") + d1(s -> Fy(x, s), y; δ = big"1e-18")
        worst_pde = max(worst_pde, abs(div + 12μ) / 12)
        worst_tip = max(worst_tip, abs(Fx(big(0), y)) / abs(Fy(x, y)))           # sealed tip x = 0
        if edge0 == :sealed
            worst_face = max(worst_face, abs(Fy(x, big(0))) / abs(Fx(x, y)))     # blocked face y = 0
        else
            worst_open = max(worst_open, abs(p(x, big(0))) / abs(p(x, y)))       # open face y = 0
        end
    end
    identity_check!("Prop 5.1: corner pressures solve ∇·(h³∇p) = −12μ", worst_pde, 1; tol = big"1e-15")
    identity_check!("Prop 5.1: no flux through the sealed tip", worst_tip, 1; tol = big"1e-15")
    identity_check!("Prop 5.1: no flux through the blocked face (sealed corner)", worst_face, 1; tol = big"1e-15")
    identity_check!("Prop 5.1: p = 0 on the open face (open corner)", worst_open, 1; tol = big"1e-30")
end

# 1.6  Proposition 5.2: slip-corner solvability ∫₀^{π/2} 2 A_s k_p s(θ)² dθ = 12μ π/2
let α = big"0.0467", β = 2tan(big(pi) / 180), kp = big"4.565e-7", μ = big(MU)
    I = bquad(θ -> (α * cos(θ) + β * sin(θ))^2, 0, big(pi) / 2)
    As = 12 * big(pi) * μ / (kp * (big(pi) * (α^2 + β^2) + 4α * β))
    identity_check!("Prop 5.2: ∫ s² dθ = π(α²+β²)/4 + αβ", I - (big(pi) * (α^2 + β^2) / 4 + α * β), I)
    identity_check!("Prop 5.2: solvability gives A_s = 12πμ/(k_p(πR²+4αβ))", 2As * kp * I - 12μ * big(pi) / 2, 12μ)
end

# 1.7  Proposition 4.1: strip law, ∫₀^W p = 12μ ∫₀^W y²/G for one vent (order exchange)
let h0 = big"1e-6", β = big"0.035", kp = big"4.565e-7", Wb = big(W), μ = big(MU)
    G(y) = (h0 + β * y)^2 * (h0 + β * y + kp)
    p(y) = 12μ * bquad(s -> s / G(s), y, Wb)
    lhs = quadgk(p, big(0), Wb; rtol = big"1e-30")[1]
    rhs = 12μ * bquad(y -> y^2 / G(y), 0, Wb)
    identity_check!("Prop 4.1: ∫₀^W p dy = 12μ ∫₀^W y²/G dy", lhs - rhs, rhs; tol = big"1e-25")
end

# 1.8  Small-taper expansions (4.2)
let δ = big"1e-3"      # the remainder must fall as δ⁴: halving δ divides it by 16
    Iδ(d, k) = bquad(η -> η^k / (1 + d * (η - big"0.5"))^3, 0, 1)
    rem1(d) = 3Iδ(d, 2) - (1 - 3d / 4 + 3d^2 / 5 - 3d^3 / 8)
    rem2(d) = 12 * (Iδ(d, 2) - Iδ(d, 1)^2 / Iδ(d, 0)) - (1 + 3d^2 / 20)
    for (name, remf) in (("(4.2) one vent: 1 − 3δ/4 + 3δ²/5 − 3δ³/8 + O(δ⁴)", rem1), ("(4.2) two vents: 1 + 3δ²/20 + O(δ⁴)", rem2))
        r1 = remf(δ); r2 = remf(δ / 2)
        order = Float64(log2(abs(r1 / r2)))
        check!("derivation", name, abs(order - 4) < 0.05, @sprintf("remainder ∝ δ^%.3f, δ⁴ coefficient %.4f", order, Float64(r1 / δ^4)))
    end
end

# 1.9  Linear-taper closed form (4.3) and the finiteness bracket (Appendix C.3)
let h0 = big"0.3", β = big"0.7", Wb = big"2.1"
    hW = h0 + β * Wb
    num = bquad(y -> y^2 / (h0 + β * y)^3, 0, Wb)
    cf = (log(hW / h0) - big"1.5" + 2h0 / hW - h0^2 / (2hW^2)) / β^3
    identity_check!("(4.3) linear-taper strip integral, closed form", num - cf, cf)
    tmax = maximum(t -> -1.5 + 2t - t^2 / 2, range(0, 1; length = 10001))
    check!("derivation", "finiteness bracket −3/2 + 2t − t²/2 ≤ 1/2 on (0,1]", tmax <= 0.5, @sprintf("max %.6f", tmax))
end

# 1.10 Strip limits (4.4): J_B = 1/2, J_C = 1 − ln 2
let
    J0(ξ) = (1 / ξ^2 - 1 / (ξ + 1)^2) / 2
    J1(ξ) = (1 / ξ - 1 / (ξ + 1)) - ξ * J0(ξ)
    J2(ξ) = log((ξ + 1) / ξ) - 2ξ * (1 / ξ - 1 / (ξ + 1)) + ξ^2 * J0(ξ)
    JB = bquad(J2, 0, Inf)
    JC = quadgk(ξ -> J2(ξ) - J1(ξ)^2 / J0(ξ), big(0), big(Inf); rtol = big"1e-30")[1]
    identity_check!("(4.4) J_B = 1/2", JB - big"0.5", 1; tol = big"1e-30")
    identity_check!("(4.4) J_C = 1 − ln 2", JC - (1 - log(big(2))), 1; tol = big"1e-25")
end

# 1.11 Slice function F(a,k) (6.7): definition by a limit, its two limits, and F = Ψ − 1 − ln k
let a = big"0.37", k = big"1.9", U = big"1e30"
    F(a, k) = a * (2k + a) / k^2 * log(a) + a / k - (a + k)^2 / k^2 * log(a + k)
    Ψ(ϱ) = 1 + ϱ - log(1 + ϱ) - ϱ * (2 + ϱ) * log(1 + 1 / ϱ)
    integral = quadgk(u -> u^2 / ((u + a)^2 * (u + a + k)), big(0), big(1), big(10)^3, big(10)^6, big(10)^10, big(10)^15,
                      big(10)^20, big(10)^25, U; rtol = big"1e-35")[1]
    identity_check!("(6.7) F(a,k) = lim [∫₀^U u²/((u+a)²(u+a+k)) du − ln U]", integral - log(U) - F(a, k), 1; tol = big"1e-25")
    # F cancels terms of size a²/k² ≈ 1e60 here, so evaluate at 1024 bits; compare with the two-term expansion
    lim_k0 = setprecision(BigFloat, 1024) do
        kk = big"1e-30"; aa = BigFloat(a)
        F(aa, kk) - (-log(aa) - big"1.5" - kk / (3aa) + kk^2 / (12aa^2))
    end
    identity_check!("(6.7) F(a,k→0) = −ln a − 3/2 − k/(3a) + k²/(12a²) + O(k³)", lim_k0, 1)
    identity_check!("(6.7) F(a→0,k) → −ln k", F(big"1e-40", k) + log(k), 1; tol = big"1e-30")
    identity_check!("(6.7) F(a,k) = Ψ(a/k) − 1 − ln k", F(a, k) - (Ψ(a / k) - 1 - log(k)), 1)
end

# 1.12 Modal sum rules (3.5) and the screened-logarithm constants
let
    s2 = (1 - big(2)^-2) * zeta(big(2)); s4 = (1 - big(2)^-4) * zeta(big(4))
    identity_check!("(3.5) Σ_odd q⁻² = π²/8", s2 - big(pi)^2 / 8, 1)
    identity_check!("(3.5) Σ_odd q⁻⁴ = π⁴/96", s4 - big(pi)^4 / 96, 1)
    partial = sum(1 / big(2n + 1)^2 for n in 0:200000)
    check!("derivation", "(3.5) partial sum of q⁻² approaches π²/8 like 1/(4N)", abs(Float64(big(pi)^2 / 8 - partial) * 4 * 200001 - 1) < 1e-3, "")
    glaisher = big"1.28242712910062263687534256886979172776768892732500119206374"
    c0 = big"0.5" + log(big(pi)) + log(big(2)) / 3 - 12log(glaisher)
    C1 = log(2exp(big(1)) / big(pi)) + c0; C2 = log(exp(big(1)) / big(pi)) + c0
    check!("derivation", "c₀ = ½ + ln π + (ln 2)/3 − 12 ln A = −1.109275", abs(Float64(c0) - C0) < 5e-7, @sprintf("%.7f", Float64(c0)))
    check!("derivation", "C_X one vent = ln(2e/π) + c₀ = −0.5609 (as printed)", round(Float64(C1); digits = 4) == -0.5609, @sprintf("%.6f", Float64(C1)))
    check!("derivation", "C_X two vents = ln(e/π) + c₀ = −1.2540 (as printed)", round(Float64(C2); digits = 4) == -1.2540, @sprintf("%.6f", Float64(C2)))
end

# 1.13 Geometric means (6.8): linear βW/e, bow g_m/e², kink g_m/e, scallop g_m/4
let Wb = big"2.5", β = big"0.4", gm = big"0.8"
    mean_log(g) = quadgk(y -> log(g(y)), big(0), Wb / 2, Wb; rtol = big"1e-30")[1] / Wb
    identity_check!("(6.8) linear taper: exp⟨ln g⟩ = βW/e", mean_log(y -> β * y) - log(β * Wb / exp(big(1))), 1; tol = big"1e-25")
    identity_check!("(6.8) bow g_m(2y/W−1)²: g_m/e²", mean_log(y -> gm * (2y / Wb - 1)^2) - (log(gm) - 2), 1; tol = big"1e-25")
    identity_check!("(6.8) kink g_m|2y/W−1|: g_m/e", mean_log(y -> gm * abs(2y / Wb - 1)) - (log(gm) - 1), 1; tol = big"1e-25")
    scallop = quadgk(t -> log(big(2)) + 2log(abs(sin(t / 2 + big(pi) / 4))), big(0), 3big(pi) / 2, 2big(pi); rtol = big"1e-30")[1] / (2big(pi))
    identity_check!("(6.8) scallop g_m(1+sin)/2: g_m/4", (scallop - log(big(2))) - log(big"0.25"), 1; tol = big"1e-25")
end

# 1.14 Scaling coefficients (2.7)–(2.12) and homogeneity (Prop 6.5)
let ρ = rb(0.5, 2), μ = rb(1e-5, 2e-5), ω = rb(1e3, 1e5), V = rb(1e-3, 1), ld = rb(1e-5, 1e-4), hs = rb(1e-7, 1e-6), pa = rb(9e4, 1.1e5)
    U = V * ld / hs; Ps = μ * V * ld^2 / hs^3; scale = μ * U / hs^2
    res = maximum(abs, [ρ * ω * U / scale - ρ * ω * hs^2 / μ, ρ * U^2 / ld / scale - ρ * V * hs / μ, Ps / ld / scale - 1,
                        μ * U / ld^2 / scale - (hs / ld)^2, (μ * V / hs^2) / (Ps / hs) - (hs / ld)^2,
                        (ω * ρ * Ps / pa) / (ρ * V / hs) - μ * ω * ld^2 / (pa * hs^2)] ./
                       [ρ * ω * hs^2 / μ, ρ * V * hs / μ, 1, (hs / ld)^2, (hs / ld)^2, μ * ω * ld^2 / (pa * hs^2)])
    identity_check!("(2.7)–(2.12) every scaling coefficient", res, 1)
    c = rb(0.5, 3); a = rb(0.01, 0.1); b = rb(0.01, 0.1); kp = rb(1e-7, 1e-6); x = rb(0, 1e-5); y = rb(0, 1e-5)
    G(a, b, k) = (a * x + b * y)^2 * ((a * x + b * y) + k)
    identity_check!("Prop 6.5: G homogeneous of degree 3 in (α, β, k_p)", G(c * a, c * b, c * kp) - c^3 * G(a, b, kp), c^3 * G(a, b, kp))
end

runtime_failure!(name, err) = check!("cases", name * ": runtime error", false, first(sprint(showerror, err), 200))

# ═════════════════════════════════ PART 2: CASES ═════════════════════════════════
section("PART 2  Cases: film equation by finite volumes (paper's solver, ported) and P1 finite elements")

"Face coordinates on [0,L] with geometric grading from spacing d0_lo at 0 (and d0_hi at L)."
function graded_both(L, d0_lo; d0_hi = nothing, r = 1.08, dmax = L / 40)
    function one_side(d0, Lh)
        f = [0.0]; d = d0
        while f[end] < Lh
            push!(f, min(Lh, f[end] + d)); d = min(d * r, dmax)
        end
        return f
    end
    f = d0_hi === nothing ? one_side(d0_lo, L) :
        vcat(one_side(d0_lo, L / 2), L .- reverse(one_side(d0_hi, L / 2))[2:end])
    return unique(sort(f))
end

"""Cell-centred finite volumes for ∇·(G∇p) = 12μ ḣ with ḣ = −1 (closing). Faces xf, yf.
Boundaries: x = 0 tip (:sealed/:open), x = L root (:open), y = 0 and y = W faces (:blocked/:open).
Returns the damping D = ∫ p dA (N s/m), and optionally the field."""
function fv_solve(hfun, xf, yf; kp = 0.0, y0 = :blocked, yW = :open, x0 = :sealed, xL = :open, field = false)
    G(h) = h^2 * (h + kp)
    xc = (xf[1:end-1] .+ xf[2:end]) ./ 2; yc = (yf[1:end-1] .+ yf[2:end]) ./ 2
    nx, ny = length(xc), length(yc); dx = diff(xf); dy = diff(yf); N = nx * ny
    id(i, j) = (i - 1) * ny + j
    I = Int[]; J = Int[]; V = Float64[]; dg = zeros(nx, ny)
    sizehint!(I, 5N); sizehint!(J, 5N); sizehint!(V, 5N)
    for i in 1:nx-1, j in 1:ny
        c = G(hfun(xf[i+1], yc[j])) * dy[j] / (xc[i+1] - xc[i])
        dg[i, j] -= c; dg[i+1, j] -= c
        push!(I, id(i, j), id(i + 1, j)); push!(J, id(i + 1, j), id(i, j)); push!(V, c, c)
    end
    for i in 1:nx, j in 1:ny-1
        c = G(hfun(xc[i], yf[j+1])) * dx[i] / (yc[j+1] - yc[j])
        dg[i, j] -= c; dg[i, j+1] -= c
        push!(I, id(i, j), id(i, j + 1)); push!(J, id(i, j + 1), id(i, j)); push!(V, c, c)
    end
    if xL == :open; for j in 1:ny; dg[nx, j] -= G(hfun(xf[end], yc[j])) * dy[j] / (xf[end] - xc[nx]); end; end
    if x0 == :open; for j in 1:ny; dg[1, j]  -= G(hfun(xf[1], yc[j])) * dy[j] / (xc[1] - xf[1]); end; end
    if yW == :open; for i in 1:nx; dg[i, ny] -= G(hfun(xc[i], yf[end])) * dx[i] / (yf[end] - yc[ny]); end; end
    if y0 == :open; for i in 1:nx; dg[i, 1]  -= G(hfun(xc[i], yf[1])) * dx[i] / (yc[1] - yf[1]); end; end
    for i in 1:nx, j in 1:ny
        push!(I, id(i, j)); push!(J, id(i, j)); push!(V, dg[i, j])
    end
    A = sparse(I, J, V, N, N)
    rhs = [-12MU * dx[i] * dy[j] for i in 1:nx for j in 1:ny]
    P = permutedims(reshape((-A) \ (-rhs), ny, nx))     # −A is symmetric positive definite; P[i, j] at (xc[i], yc[j])
    D = sum(P .* (dx * dy'))
    return field ? (D, P, xc, yc) : D
end

"Bilinear interpolation of a cell-centred field."
function interp(P, xc, yc, x, y)
    i = clamp(searchsortedlast(xc, x), 1, length(xc) - 1); j = clamp(searchsortedlast(yc, y), 1, length(yc) - 1)
    tx = (x - xc[i]) / (xc[i+1] - xc[i]); ty = (y - yc[j]) / (yc[j+1] - yc[j])
    return (1 - tx) * (1 - ty) * P[i, j] + tx * (1 - ty) * P[i+1, j] + (1 - tx) * ty * P[i, j+1] + tx * ty * P[i+1, j+1]
end

# closure classes: (y = 0 face, y = W face, narrow side, drainage length)
const CLASSES = Dict('A' => (:blocked, :open, :none, W), 'B' => (:blocked, :open, :y0, W),
                     'C' => (:open, :open, :y0, W / 2), 'D' => (:blocked, :open, :yW, W))
const RG = QUICK ? 1.15 : 1.08                         # grading ratio (paper: 1.08)

function class_grids(cls, beta; r = RG, fine = false)
    narrow = CLASSES[cls][3]
    r = fine ? 1.05 : r; dmx = fine ? W / 16 : W / 8; dmy = fine ? W / 80 : W / 40
    xf = graded_both(LFILM, 1e-11; r = r, dmax = dmx)
    yf = (beta == 0 || narrow == :none) ? collect(range(0, W; length = fine ? 161 : 81)) :
         narrow == :y0 ? graded_both(W, 1e-11; d0_hi = W / 400, r = r, dmax = dmy) :
                         graded_both(W, W / 400; d0_hi = 1e-11, r = r, dmax = dmy)
    return xf, yf
end
function gapfun(cls, beta, hmin; alpha = ALPHA)
    CLASSES[cls][3] == :yW ? (x, y) -> hmin + alpha * x + beta * (W - y) : (x, y) -> hmin + alpha * x + beta * y
end

"Φ_X = D/D★ for class X at taper ratio β/α and slip parameter κ (contact: h_min = 1e-13 m)."
function Phi(cls, ba, kappa; hmin = 1e-13, tip = :sealed, alpha = ALPHA, r = RG, fine = false, xfyf = nothing)
    y0, yW, narrow, ld = CLASSES[cls]; beta = ba * alpha; kp = 2alpha * ld * kappa / π
    xf, yf = xfyf === nothing ? class_grids(cls, beta; r = r, fine = fine) : xfyf
    D = fv_solve(gapfun(cls, beta, hmin; alpha = alpha), xf, yf; kp = kp, y0 = y0, yW = yW, x0 = tip)
    return D / (MU * W / alpha^3)
end

"""Galerkin P1 finite elements on the triangulated tensor grid (nodes xn × yn). With the cubic
mobility integrated exactly (7-point degree-5 rule), the result is a rigorous lower bound on D."""
function fe_solve(hfun, xn, yn; kp = 0.0, y0 = :blocked, yW = :open, x0 = :sealed, xL = :open)
    G(h) = h^2 * (h + kp)
    nxn, nyn = length(xn), length(yn); node(i, j) = (i - 1) * nyn + j; N = nxn * nyn
    # Dunavant degree-5 rule on the reference triangle (barycentric coordinates, weights sum to 1)
    a1, b1 = 0.059715871789770, 0.470142064105115; a2, b2 = 0.797426985353087, 0.101286507323456
    bary = [(1/3, 1/3, 1/3), (a1, b1, b1), (b1, a1, b1), (b1, b1, a1), (a2, b2, b2), (b2, a2, b2), (b2, b2, a2)]
    wq = [0.225, fill(0.132394152788506, 3)..., fill(0.125939180544827, 3)...]
    I = Int[]; J = Int[]; V = Float64[]; F = zeros(N); M = zeros(N)
    for i in 1:nxn-1, j in 1:nyn-1
        for tri in (((i, j), (i + 1, j), (i + 1, j + 1)), ((i, j), (i + 1, j + 1), (i, j + 1)))
            v = [node(t...) for t in tri]; px = [xn[t[1]] for t in tri]; py = [yn[t[2]] for t in tri]
            bb = (py[2] - py[3], py[3] - py[1], py[1] - py[2]); cc = (px[3] - px[2], px[1] - px[3], px[2] - px[1])
            area = (bb[1] * cc[2] - bb[2] * cc[1]) / 2
            Gint = area * sum(wq[q] * G(hfun(sum(bary[q] .* px), sum(bary[q] .* py))) for q in 1:7)
            for k in 1:3, l in 1:3
                push!(I, v[k]); push!(J, v[l]); push!(V, Gint * (bb[k] * bb[l] + cc[k] * cc[l]) / (4area^2))
            end
            for k in 1:3
                F[v[k]] += 12MU * area / 3; M[v[k]] += area / 3
            end
        end
    end
    K = sparse(I, J, V, N, N)
    fixed = falses(N)
    for i in 1:nxn, j in 1:nyn
        ((xL == :open && i == nxn) || (x0 == :open && i == 1) || (yW == :open && j == nyn) || (y0 == :open && j == 1)) &&
            (fixed[node(i, j)] = true)
    end
    free = findall(!, fixed); p = zeros(N)
    p[free] = K[free, free] \ F[free]
    return dot(M, p)
end

# 2.1  Class A screened logarithm (Prop 6.2), sealed and open tip, with FE lower bounds
println("\n2.1  Class A (parallel sidewalls), no slip: screened logarithm 12[ln(1/Λ) + c₀ (− ½ if the tip is open)]")
screened = Tuple{Float64, Float64, Float64, Float64}[]
for hmin in (QUICK ? (1e-9, 1e-10) : (1e-8, 1e-9, 1e-10, 1e-11))
    Λ = π * hmin / (2ALPHA * W); exact = 12 * (log(1 / Λ) + C0)
    xf = graded_both(LFILM, hmin / ALPHA / 4; r = QUICK ? 1.1 : 1.05, dmax = W / 16); yf = collect(range(0, W; length = 41))
    fv = fv_solve((x, y) -> hmin + ALPHA * x, xf, yf) / DSTAR
    fe = fe_solve((x, y) -> hmin + ALPHA * x, xf, yf) / DSTAR
    push!(screened, (Λ, fv, fe, exact))
    if Λ < 2e-3      # at larger Λ the asymptote itself is off by O(Λ ln Λ): plotted, not asserted
        check!("cases", @sprintf("class A sealed tip, h_min = %.0e m: FE ≤ exact ≤ FV", hmin), fe <= exact <= fv && relerr(fv, exact) < 1e-2,
               @sprintf("FE %.3f  exact %.3f  FV %.3f", fe, exact, fv))
    else
        @printf("  [info] h_min = %.0e m (Λ = %.1e): FE %.3f, asymptote %.3f, FV %.3f; asymptote not yet accurate here\n", hmin, Λ, fe, exact, fv)
    end
end
let hmin = 1e-10
    Λ = π * hmin / (2ALPHA * W); exact = 12 * (log(1 / Λ) + C0 - 0.5)
    xf = graded_both(LFILM, hmin / ALPHA / 4; r = QUICK ? 1.1 : 1.05, dmax = W / 16); yf = collect(range(0, W; length = 41))
    fv = fv_solve((x, y) -> hmin + ALPHA * x, xf, yf; x0 = :open) / DSTAR
    fe = fe_solve((x, y) -> hmin + ALPHA * x, xf, yf; x0 = :open) / DSTAR
    check!("cases", "class A open tip, h_min = 1e-10 m: FE ≤ exact ≤ FV", fe <= exact <= fv && relerr(fv, exact) < 1e-2,
           @sprintf("FE %.3f  exact %.3f  FV %.3f", fe, exact, fv))
end

# 2.2  Contact maps (Appendix E): the port must reproduce the tables
println("\n2.2  Contact maps Φ_X(β/α, κ): ported solver vs the paper's tables")
map_points = [('B', 0.05, 0.0), ('B', 0.5, 0.0), ('B', 1.0, 0.3), ('B', 5.0, 3.0), ('C', 0.1, 0.0), ('C', 1.0, 1.0),
              ('D', 0.2, 0.0), ('D', 1.0, 0.1), ('D', 5.0, 3.0)]
for (cls, ba, k) in map_points
    v = Phi(cls, ba, k); ref = REF[(cls, ba, k)]
    check!("cases", @sprintf("Φ_%c(β/α = %g, κ = %g)", cls, ba, k), relerr(v, ref) < (QUICK ? 2e-2 : 1e-4), @sprintf("%.5f vs table %.5f", v, ref))
end

# 2.3  Class A slip saturation at contact
println("\n2.3  Class A with slip: finite contact damping")
for (k, ref) in zip(REF_A_KAPPA[1:(QUICK ? 2 : 4)], REF_A_SAT)
    v = Phi('A', 0.0, k)
    check!("cases", @sprintf("Φ_A(κ = %g) at contact", k), relerr(v, ref) < (QUICK ? 2e-2 : 5e-3), @sprintf("%.4f vs %.4f", v, ref))
end

# 2.4  Exact scaling (Prop 6.5): doubling α, β, k_p (and h_min) leaves Φ unchanged
println("\n2.4  Exact scaling and monotonicity")
let grids = class_grids('B', 0.7473 * ALPHA)
    p1 = Phi('B', 0.7473, 0.6142; xfyf = grids)
    p2 = Phi('B', 0.7473, 0.6142; alpha = 2ALPHA, hmin = 2e-13, xfyf = grids)
    check!("cases", "Prop 6.5: Φ_B unchanged when α, β, k_p are doubled", relerr(p1, p2) < 1e-9, @sprintf("%.10f vs %.10f", p1, p2))
end
# 2.5  Monotonicity (Prop 3.2)
let vr = [Phi('B', ba, 0.0) for ba in (0.1, 0.5, 1.0)], vk = [Phi('B', 0.5, k) for k in (0.0, 0.3, 1.0)]
    check!("cases", "Prop 3.2: Φ_B decreases with taper ratio", all(diff(vr) .< 0), string(round.(vr; digits = 3)))
    check!("cases", "Prop 3.2: Φ_B decreases with slip", all(diff(vk) .< 0), string(round.(vk; digits = 3)))
    sealed = vr[2]; opened = Phi('B', 0.5, 0.0; tip = :open)
    check!("cases", "Prop 3.2: venting the tip lowers the damping", opened < sealed, @sprintf("%.4f → %.4f", sealed, opened))
end

# 2.6  Tip venting table
println("\n2.6  Tip venting (open tip / sealed tip)")
for ((cls, ba), (s_ref, o_ref)) in REF_TIPVENT
    s = Phi(cls, ba, 0.0); o = Phi(cls, ba, 0.0; tip = :open)
    check!("cases", @sprintf("class %c, β/α = %g: open/sealed ratio", cls, ba), relerr(o / s, o_ref / s_ref) < (QUICK ? 2e-2 : 1e-3),
           @sprintf("%.4f vs %.4f", o / s, o_ref / s_ref))
end

# 2.7  Corners (Props 5.1, 5.2): field near the contact point vs closed forms
println("\n2.7  Corner fields at r = 10–100 nm vs closed forms")
corner_fields = Dict{Symbol, Any}()
let β = 2tan(deg2rad(1.0)), hmin = 1e-13, L = 400e-6
    xf = graded_both(L, 2e-11; r = QUICK ? 1.12 : 1.06); yf = graded_both(W, 2e-11; d0_hi = W / 400, r = QUICK ? 1.12 : 1.06, dmax = W / 60)
    for (edge0, y0) in ((:sealed, :blocked), (:open, :open))
        D, P, xc, yc = fv_solve((x, y) -> hmin + ALPHA * x + β * y, xf, yf; y0 = y0, field = true)
        corner_fields[edge0] = (P, xc, yc)
        worst = 0.0
        for ang in (10, 45, 80), r in (3e-8,)
            x, y = r * cosd(ang), r * sind(ang)
            worst = max(worst, relerr(interp(P, xc, yc, x, y), corner_closed(ALPHA, β, x, y; edge0 = edge0)))
        end
        t = deg2rad(45)
        e = log(interp(P, xc, yc, 1e-7cos(t), 1e-7sin(t)) / interp(P, xc, yc, 1e-8cos(t), 1e-8sin(t))) / log(10)
        check!("cases", "Prop 5.1 $(edge0) corner: amplitude at r = 30 nm", worst < (QUICK ? 3e-2 : 5e-3), @sprintf("max rel. error %.1e", worst))
        check!("cases", "Prop 5.1 $(edge0) corner: exponent d ln p/d ln r = −1", abs(e + 1) < (QUICK ? 2e-2 : 5e-3), @sprintf("%.4f", e))
    end
    kp = 456.5e-9; R = hypot(ALPHA, β); As = 12π * MU / (kp * (π * R^2 + 4ALPHA * β))
    D, P, xc, yc = fv_solve((x, y) -> hmin + ALPHA * x + β * y, xf, yf; kp = kp, field = true)
    t = deg2rad(45)
    s = (interp(P, xc, yc, 1e-9cos(t), 1e-9sin(t)) - interp(P, xc, yc, 1e-8cos(t), 1e-8sin(t))) / log(10)
    check!("cases", "Prop 5.2 slip corner: −dp/d ln r = A_s", relerr(s, As) < (QUICK ? 3e-2 : 5e-3), @sprintf("%.4e vs %.4e", s, As))
end

# 2.8  Strip law (Prop 4.1) at mid-length of a long film whose gap does not vary along x
println("\n2.8  Strip law, one and two vents")
strip_profiles = Dict{Symbol, Any}()
let β = 2tan(deg2rad(1.0)), hm = 1e-6, kp = 456.5e-9, L = 40W
    g(y) = (hm + β * (y - W / 2))^2 * (hm + β * (y - W / 2) + kp)
    xf = collect(range(0, L; length = QUICK ? 201 : 801)); yf = collect(range(0, W; length = QUICK ? 81 : 161))
    for (y0, label) in ((:blocked, :one), (:open, :two))
        ys = y0 == :blocked ? 0.0 : quadgk(y -> y / g(y), 0, W)[1] / quadgk(y -> 1 / g(y), 0, W)[1]
        strip = 12MU * quadgk(y -> (y - ys)^2 / g(y), 0, W)[1]
        D, P, xc, yc = fv_solve((x, y) -> hm + β * (y - W / 2), xf, yf; kp = kp, y0 = y0, x0 = :open, field = true)
        im = searchsortedfirst(xc, L / 2); line = P[im, :]
        mid = sum(line .* diff(yf))
        pexact = [12MU * quadgk(s -> (s - ys) / g(s), y, W)[1] for y in yc]
        strip_profiles[label] = (yc, line, pexact)
        check!("cases", "Prop 4.1 strip law, $(label) vent(s)", relerr(mid, strip) < (QUICK ? 5e-3 : 1e-3), @sprintf("%.6f vs %.6f N s/m²", mid, strip))
    end
end

# 2.9  Small-taper law (6.9): Φ_B − 12 ln(α/β) → C = −0.5609
println("\n2.9  Small taper: Φ_B − 12 ln(α/β) → ln(2e/π) + c₀")
const CPRED = log(2 * exp(1.0) / π) + C0            # −0.56086 (Julia spells Euler's number ℯ, not e)
bracket(ph, ba) = (ph - 12log(1 / ba)) / 12
function extrapolate_C(rs, bs)                       # fit C + a r ln(1/r) + b r through three points
    M = hcat(ones(3), rs .* log.(1 ./ rs), rs)
    return (M \ bs)[1]
end
small_rs = QUICK ? [0.002, 0.005, 0.01] : [0.001, 0.002, 0.005]
small_fe = Float64[]; small_fv = Float64[]; small_fvf = Float64[]
try
    for ba in small_rs
        g = class_grids('B', ba * ALPHA)
        push!(small_fv, bracket(Phi('B', ba, 0.0; xfyf = g), ba))
        push!(small_fe, bracket(fe_solve(gapfun('B', ba * ALPHA, 1e-13), g[1], g[2]) / DSTAR, ba))
        QUICK || push!(small_fvf, bracket(Phi('B', ba, 0.0; fine = true), ba))
        @printf("  β/α = %-6g bracket: FE %.5f  FV %.5f  %s\n", ba, small_fe[end], small_fv[end],
                QUICK ? "" : @sprintf("FV fine %.5f", small_fvf[end]))
    end
    let Cfe = extrapolate_C(small_rs, small_fe), Cfv = extrapolate_C(small_rs, small_fv)
        check!("cases", "(6.9) C_X = ln(2e/π) + c₀ lies between FE and FV extrapolations", Cfe <= CPRED <= Cfv,
               @sprintf("FE %.4f ≤ %.4f ≤ FV %.4f", Cfe, CPRED, Cfv))
        if !QUICK
            Cfvf = extrapolate_C(small_rs, small_fvf)
            check!("cases", "(6.9) refined FV extrapolation within 0.004 of the prediction", abs(Cfvf - CPRED) < 4e-3,
                   @sprintf("%.4f vs %.4f (coarse FV %.4f)", Cfvf, CPRED, Cfv))
        end
    end
catch err
    runtime_failure!("(6.9) small-taper constant", err)
end

# 2.10 Galerkin lower bound on a corner class (FE ≤ FV, within 0.5%)
println("\n2.10 Galerkin lower bound vs finite volumes on a corner class")
try
    let g = class_grids('B', 0.5ALPHA)
        xf, yf = g
        fv = Phi('B', 0.5, 0.0; xfyf = g)
        fe = fe_solve(gapfun('B', 0.5ALPHA, 1e-13), xf, yf) / DSTAR
        check!("cases", "class B, β/α = 0.5: FE (lower bound) ≤ FV, gap < 0.5%", fe <= fv && relerr(fe, fv) < 5e-3, @sprintf("FE %.4f  FV %.4f", fe, fv))
    end
catch err
    runtime_failure!("Galerkin lower bound, class B", err)
end

# ═════════════════════════════════ PART 3: FIGURES ═════════════════════════════════
section("PART 3  Figures → $(FIGDIR)")
const ASCII_MAP = ("μ" => "mu", "∂" => "d", "∇" => "grad ", "≤" => "<=", "∫" => "int ", "Σ" => "sum", "⟨" => "<", "⟩" => ">",
                   "²" => "^2", "³" => "^3", "⁴" => "^4", "⁻" => "^-", "₀" => "0", "π" => "pi", "α" => "alpha", "β" => "beta",
                   "δ" => "delta", "λ" => "lambda", "θ" => "theta", "ψ" => "psi", "Ψ" => "Psi", "−" => "-", "→" => "->",
                   "·" => "*", "½" => "1/2", "ₚ" => "p", "ℯ" => "e")
ascii_label(s) = replace(s, ASCII_MAP...)

# Fig 1: residuals of the derivation checks
try
    let names = ascii_label.(first.(RESIDUALS)), vals = log10.(last.(RESIDUALS))
        fig = bar(1:length(vals), vals; orientation = :h, yticks = (1:length(vals), names), legend = false,
                  color = OKABE[6], xlabel = L"\log_{10}(\mathrm{relative\ residual})", size = (760, 640),
                  title = "Derivations: residuals in 256-bit arithmetic", left_margin = 6Plots.mm)
        vline!(fig, [-15]; color = OKABE[7], linestyle = :dash)
        savefig(fig, joinpath(FIGDIR, "fig1_derivation_residuals.pdf"))
    end
catch err
    @warn "Fig 1 could not be drawn" exception = err
end

# Fig 2: screened logarithm with FE / FV bracket
try
    let Λ = [s[1] for s in screened]
        fig = plot(Λ, [s[4] for s in screened]; xscale = :log10, label = L"12[\ln(1/\Lambda)+c_0]", color = OKABE[1],
                   xlabel = L"\Lambda = \pi h_\mathrm{min}/(2\alpha W)", ylabel = L"D/D_\star", legend = :topright)
        scatter!(fig, Λ, [s[2] for s in screened]; label = "finite volumes", marker = :circle, color = OKABE[2])
        scatter!(fig, Λ, [s[3] for s in screened]; label = "finite elements (lower bound)", marker = :utriangle, color = OKABE[6])
        savefig(fig, joinpath(FIGDIR, "fig2_screened_log.pdf"))
    end
catch err
    @warn "Fig 2 could not be drawn" exception = err
end

# Fig 3: corner fields collapse onto the closed forms, r p vs angle
try
    let β = 2tan(deg2rad(1.0)), θs = range(1, 89; length = 89)
        fig = plot(xlabel = L"\theta\ (\mathrm{deg})", ylabel = L"r\,p\ \ (\mathrm{Pa\,m})", legend = :topleft)
        for (k, edge0) in enumerate((:sealed, :open))
            P, xc, yc = corner_fields[edge0]
            plot!(fig, θs, [1e-8 * corner_closed(ALPHA, β, 1e-8cosd(θ), 1e-8sind(θ); edge0 = edge0) for θ in θs];
                  color = OKABE[1], linestyle = k == 1 ? :solid : :dash, label = "closed form, $(edge0) corner")
            for (m, r) in enumerate((1e-8, 3e-8, 1e-7))
                scatter!(fig, θs[1:8:end], [r * interp(P, xc, yc, r * cosd(θ), r * sind(θ)) for θ in θs[1:8:end]];
                         color = OKABE[1+m], marker = k == 1 ? :circle : :diamond, markersize = 3,
                         label = k == 1 ? @sprintf("FV, r = %g nm", r * 1e9) : "")
            end
        end
        savefig(fig, joinpath(FIGDIR, "fig3_corner_collapse.pdf"))
    end
catch err
    @warn "Fig 3 could not be drawn" exception = err
end

# Fig 4: strip-law pressure profiles
try
    let fig = plot(xlabel = L"y/W", ylabel = L"p\ \ (\mathrm{Pa})", legend = :topright)
        for (k, label) in enumerate((:one, :two))
            yc, line, pexact = strip_profiles[label]
            plot!(fig, yc ./ W, pexact; color = OKABE[1], linestyle = k == 1 ? :solid : :dash, label = "strip law, $(label) vent(s)")
            scatter!(fig, yc[1:6:end] ./ W, line[1:6:end]; color = OKABE[1+k], markersize = 3, label = "FV mid-length, $(label)")
        end
        savefig(fig, joinpath(FIGDIR, "fig4_strip_law.pdf"))
    end
catch err
    @warn "Fig 4 could not be drawn" exception = err
end

# Fig 5: contact maps with the paper's table values
try
    let bas = QUICK ? [0.1, 0.5, 2.0] : [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
        fig = plot(xscale = :log10, yscale = :log10, xlabel = L"\beta/\alpha", ylabel = L"\Phi_X = D_c/D_\star", legend = :bottomleft)
        for (k, cls) in enumerate(('B', 'C', 'D'))
            for (m, κ) in enumerate((0.0, 1.0))
                plot!(fig, bas, [Phi(cls, ba, κ) for ba in bas]; color = OKABE[1+k], linestyle = m == 1 ? :solid : :dash,
                      label = @sprintf("class %c, κ = %g", cls, κ))
                scatter!(fig, REF_BA, [REF[(cls, Float64(ba), κ)] for ba in REF_BA]; color = OKABE[1+k], markersize = 2.5, label = "")
            end
        end
        savefig(fig, joinpath(FIGDIR, "fig5_contact_maps.pdf"))
    end
catch err
    @warn "Fig 5 could not be drawn" exception = err
end

# Fig 6: small-taper constant
try
    let fig = plot(xscale = :log10, xlabel = L"\beta/\alpha", ylabel = L"(\Phi_B - 12\ln(\alpha/\beta))/12", legend = :bottomright)
        scatter!(fig, small_rs, small_fv; label = "finite volumes", color = OKABE[2], marker = :circle)
        QUICK || scatter!(fig, small_rs, small_fvf; label = "finite volumes, refined", color = OKABE[4], marker = :diamond)
        scatter!(fig, small_rs, small_fe; label = "finite elements (lower bound)", color = OKABE[6], marker = :utriangle)
        hline!(fig, [CPRED]; label = L"\ln(2e/\pi)+c_0 = -0.5609", color = OKABE[1], linestyle = :dash)
        savefig(fig, joinpath(FIGDIR, "fig6_small_taper_constant.pdf"))
    end
catch err
    @warn "Fig 6 could not be drawn" exception = err
end

# ═════════════════════════════════ SUMMARY ═════════════════════════════════
section("SUMMARY")
for part in ("derivation", "cases")
    sel = filter(r -> r.part == part, RESULTS)
    @printf("  %-12s %3d / %3d passed\n", part, count(r -> r.ok, sel), length(sel))
end
failed = filter(r -> !r.ok, RESULTS)
isempty(failed) || (println("\n  Failed checks:"); foreach(r -> println("    ", r.name, "  ", r.detail), failed))
println("\n  Figures written to ", FIGDIR)
exit(isempty(failed) ? 0 : 1)
