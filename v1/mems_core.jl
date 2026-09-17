# =============================================================================
# mems_core.jl -- physics core of the corrected trapezoidal-electrode MEMS model
# =============================================================================
#
# This file is `include`d by collision_model_corrected.jl into the module
# CorrectedMEMS. It contains ONLY the model: parameters, reduced-order
# construction, constitutive laws, equations of motion, and energy accounting.
# It needs nothing beyond LinearAlgebra. Verification lives in mems_verify.jl
# and the ODE harness in mems_simulate.jl.
#
# Conventions
#   * SI units everywhere.
#   * Names mirror the symbols of Model_derivations_corrected.tex so that code
#     and derivation can be read side by side:
#         x1  shuttle coordinate            x2  electrode-tip coordinate
#         gc  tip travel to contact  g_c    alpha  gap slope  alpha
#         hd  dielectric air-equivalent thickness 2 Tp / eps_r
#         kp  slip length 6 sigma_p lambda  B1, B2  shape weights on the y grid
#         r   gap side: +1 closes for x > 0 ("right"), -1 opens ("left")
#   * Every function of the state is generic in its number type. Parameters
#     and the precomputed Model stay Float64; x1, x2, v, Vout, and the closure
#     parameters may be ForwardDiff duals, BigFloat, etc. No Float64 buffers
#     are allocated on a state-dependent path.
#
# State vector (see STATE_NAMES)
#     1:5    physical   x1, x2, v1, v2, Vout
#     6:12   ledger     Wbase, Wbias, ER, Dfilm, Dstruct, Dwall, throughput
#     13     optional   Wclosure (work done by a learned/added closure force)
# The right-hand side fills only as many entries as `du` provides, so the same
# function serves a 5-state UDE/HDI problem and the 12/13-state audited run.
# =============================================================================

const N_PHYSICAL = 5
const N_LEDGER = 12
const N_LEDGER_CLOSURE = 13
const STATE_NAMES = (:x1, :x2, :v1, :v2, :Vout,
                     :Wbase, :Wbias, :ER, :Dfilm, :Dstruct, :Dwall, :throughput,
                     :Wclosure)

"Gap sides r. The r = -1 gap opens and the r = +1 gap closes for positive displacement."
const GAP_SIDES = (-1.0, 1.0)

# -----------------------------------------------------------------------------
# 1. Parameters
# -----------------------------------------------------------------------------

"""
    Params(; kwargs...)

Immutable reference parameters (SI). Build a fresh `Model` after any change;
derived quantities (mass matrix, film grid, stiffnesses) are never updated in
place. "Status" remarks follow the parameter table of the TeX document.
"""
Base.@kwdef struct Params
    # --- electrode and gap geometry -----------------------------------------
    g0::Float64 = 14e-6               # bare (uncoated) tip gap [m]; confirm metrology datum
    Tp::Float64 = 120e-9              # dielectric coating thickness per face [m]
    Tf::Float64 = 25e-6               # device-layer thickness (out of plane) [m]
    wt::Float64 = 9e-6                # beam width at the narrow clamped root [m]
    wb::Float64 = 30e-6               # beam width at the wide free tip [m]
    Lf::Float64 = 450e-6              # full electrode length, root to tip [m]
    Leff::Float64 = 400e-6            # overlapped (active) length measured from the tip [m]
    gap_slope::Float64 = NaN          # gap slope alpha [-]; NaN -> (wb - wt)/Lf; override with metrology
    h_eff::Float64 = 50e-9            # residual mean air gap at asperity contact [m]; candidate
    n_beams::Int = 80                 # mobile beams (each has two gaps: 160 gap branches); confirm layout

    # --- suspension and stops -----------------------------------------------
    ws::Float64 = 14.7e-6             # suspension span width [m]
    Lsp::Float64 = 1400e-6            # suspension span length [m]
    n_parallel::Int = 4               # suspension branches in parallel
    n_series::Int = 6                 # spans in series per branch
    gamma3::Float64 = 1.0             # cubic-stiffness geometry correction [-]; not calibrated
    wss::Float64 = 14e-6              # soft-stop fiber width [m]
    Lss::Float64 = 1000e-6            # soft-stop fiber length [m]
    n_stop_beams::Int = 2             # soft-stop fibers per side
    gss::Float64 = 14e-6              # shuttle travel to soft stop [m]
    ghs::Float64 = Inf                # shuttle travel to hard stop [m]; Inf disables (unmeasured)
    khs::Float64 = 1e9                # hard-stop coefficient [N/m^phs]
    phs::Float64 = 1.5                # hard-stop exponent [-]

    # --- material and inertia -----------------------------------------------
    E::Float64 = 170e9                # Young's modulus [Pa]
    rho::Float64 = 2330.0             # density [kg/m^3]
    m_shuttle::Float64 = 2.0933e-6    # shuttle mass [kg]; MUST exclude the explicit mobile-beam mass

    # --- electrostatics and circuit -----------------------------------------
    eps0::Float64 = 8.85e-12          # vacuum permittivity [F/m]
    epsr::Float64 = 3.2               # coating relative permittivity [-]
    cp::Float64 = 5e-12               # parasitic capacitance [F]
    Vbias::Float64 = 3.0              # bias source [V]
    Rload::Float64 = 0.42e6           # load resistance [ohm]

    # --- gas film -------------------------------------------------------------
    eta::Float64 = 1.849e-5           # air viscosity [Pa s]
    mean_free_path::Float64 = 70e-9   # lambda [m]
    slip_coefficient::Float64 = 1.016 # sigma_p [-]
    film_scale::Float64 = 1.0         # calibration multiplier c_f on the film matrix [-]
    seal_at_contact::Bool = true      # candidate open-to-sealed tip boundary closure
    seal_width::Float64 = 25e-9       # sealing transition half-width l_s [m]; sensitivity parameter

    # --- structural damping (uncalibrated reference, NOT a measured zero) ----
    c1::Float64 = 0.0                 # shuttle damping [N s/m]
    ce::Float64 = 0.0                 # damping of relative beam motion, per beam [N s/m]

    # --- compliant tip contact (regularized Hunt-Crossley), per beam ----------
    kw::Float64 = 1e6                 # contact coefficient [N/m^pw]
    pw::Float64 = 1.5                 # contact exponent [-]
    cw::Float64 = 50.0                # Hunt-Crossley damping coefficient [s/m]

    # --- regularization widths ----------------------------------------------
    eps_gap::Float64 = 2e-9           # gap-floor smoothing [m]
    eps_wall::Float64 = 0.5e-9        # tip-contact smoothing [m]
    eps_stop::Float64 = 1e-9          # soft/hard-stop smoothing [m]
end

# -----------------------------------------------------------------------------
# 2. Small numerical building blocks
# -----------------------------------------------------------------------------

"""
    softpos(z, e)

Smooth positive part sigma_e(z) = (z + sqrt(z^2 + e^2))/2 (TeX Eq. `gap`).
The z < 0 branch is the algebraically identical form e^2 / (2 (sqrt(z^2+e^2) - z)),
which avoids cancellation. Strictly positive for every finite z.
"""
function softpos(z, e)
    r = hypot(z, e)
    return z >= 0 ? (z + r) / 2 : e^2 / (2 * (r - z))
end

"Exact derivative of `softpos` with respect to `z` (TeX Eq. `hgrad`); lies in (0, 1)."
function dsoftpos(z, e)
    r = hypot(z, e)
    return z >= 0 ? (1 + z / r) / 2 : e^2 / (2 * r * (r - z))
end

"C^2 smootherstep S(z) on [0, 1] (TeX Eq. `seal`)."
smootherstep(z) = z <= 0 ? zero(z) : z >= 1 ? one(z) : z^3 * (10 - 15 * z + 6 * z^2)

"Nodes and weights of n-point Gauss-Legendre quadrature on [-1, 1] (Golub-Welsch)."
function gauss_legendre(n)
    offdiagonal = [j / sqrt(4 * j * j - 1) for j in 1:n-1]
    F = eigen(SymTridiagonal(zeros(n), offdiagonal))
    return F.values, 2 .* F.vectors[1, :] .^ 2
end

const GL_NODES, GL_WEIGHTS = gauss_legendre(96)

"96-point Gauss-Legendre integral of `f` over [a, b]. Used for all beam (s-coordinate) integrals."
function integrate(f, a, b)
    mid = (a + b) / 2
    half = (b - a) / 2
    return half * sum(GL_WEIGHTS[j] * f(mid + half * GL_NODES[j]) for j in eachindex(GL_NODES))
end

# -----------------------------------------------------------------------------
# 3. Beam mechanics (s runs from the clamped root, s = 0, to the tip, s = Lf)
# -----------------------------------------------------------------------------

"In-plane beam width b(s) (TeX Eq. `width`)."
beam_width(p, s) = p.wt + (p.wb - p.wt) * s / p.Lf

"Bending rigidity E I(s) with I = Tf b^3 / 12."
bending_rigidity(p, s) = p.E * p.Tf * beam_width(p, s)^3 / 12

"Tip stiffness k_e of one tapered electrode from Castigliano's theorem (TeX Eq. `ke`)."
electrode_tip_stiffness(p) = 1 / integrate(s -> (p.Lf - s)^2 / bending_rigidity(p, s), 0.0, p.Lf)

"""
    beam_shape(p, ke, s)

Normalized static tip-load deflection phi(s), phi(0) = 0, phi(Lf) = 1 (TeX Eq. `shape`).
The single assumed shape behind the mass matrix, local gaps, electrostatic
projection, and film projection: w(s) = x1 + phi(s) (x2 - x1).
"""
function beam_shape(p, ke, s)
    s == 0 && return 0.0
    return ke * integrate(z -> (s - z) * (p.Lf - z) / bending_rigidity(p, z), 0.0, s)
end

"Suspension linear and cubic coefficients (k1, k3) (TeX Eq. `cubic`); 0.72 = 18/25."
function suspension_stiffness(p)
    k1 = p.n_parallel / p.n_series * p.E * p.Tf * p.ws^3 / p.Lsp^3
    k3 = p.gamma3 * p.n_parallel / p.n_series^3 * 0.72 * p.E * p.Tf * p.ws / p.Lsp^3
    return k1, k3
end

"Soft-stop stiffness: n_stop tip-loaded cantilever fibers in parallel, 3 E I / L^3 each."
soft_stop_stiffness(p) = p.n_stop_beams * p.E * p.Tf * p.wss^3 / (4 * p.Lss^3)

"Consistent 2x2 mass matrix of shuttle plus n_beams shaped electrodes (TeX Eq. `mass`)."
function consistent_mass_matrix(p, ke)
    function entry(i, j)
        function integrand(s)
            phi = beam_shape(p, ke, s)
            B = (1 - phi, phi)                 # (B1, B2) at this station
            return p.rho * p.Tf * beam_width(p, s) * B[i] * B[j]
        end
        return p.n_beams * integrate(integrand, 0.0, p.Lf)
    end
    M11 = p.m_shuttle + entry(1, 1)
    M12 = entry(1, 2)
    M22 = entry(2, 2)
    return [M11 M12; M12 M22]
end

# -----------------------------------------------------------------------------
# 4. Film/capacitance grid (y runs from the tip, y = 0, along the overlap)
# -----------------------------------------------------------------------------

"""
    graded_simpson_grid(p, alpha, panels)

Composite-Simpson nodes `y` (length 2 panels + 1) and positive weights on
[0, Leff]. Panel edges are geometrically graded toward the tip with grading
length h_eff / alpha, where gaps, conductance, and capacitance density vary
fastest; panel midpoints are arithmetic in y.
"""
function graded_simpson_grid(p, alpha, panels)
    grading_length = p.h_eff / max(alpha, 0.001)
    z = range(0.0, log1p(p.Leff / grading_length); length = panels + 1)
    edges = grading_length .* expm1.(z)
    edges[end] = p.Leff
    y = zeros(2 * panels + 1)
    weights = zero(y)
    for j in 1:panels
        i = 2 * j - 1
        a, b = edges[j], edges[j+1]
        d = b - a
        y[i] = a
        y[i+1] = (a + b) / 2
        y[i+2] = b
        weights[i] += d / 6
        weights[i+1] += 2 * d / 3
        weights[i+2] += d / 6
    end
    return y, weights
end

# -----------------------------------------------------------------------------
# 5. Model: parameters plus everything that is precomputed from them
# -----------------------------------------------------------------------------

"""
    Model(p = Params(); panels = 512)

Reduced-order model. Fields (TeX symbol in brackets):

| field     | meaning                                                         |
|:----------|:----------------------------------------------------------------|
| `p`       | the `Params` this model was built from                          |
| `ke`      | electrode tip stiffness per beam [k_e], N/m                     |
| `k1`,`k3` | suspension linear / cubic force coefficients, N/m, N/m^3        |
| `kss`     | soft-stop stiffness, N/m                                        |
| `gc`      | tip travel to contact [g_c = g0 - 2 Tp - h_eff], m              |
| `alpha`   | unperturbed gap slope [alpha]                                   |
| `hd`      | air-equivalent coating thickness [h_d = 2 Tp / eps_r], m        |
| `kp`      | slip length in G(h) = h^2 (h + kp) [k_p = 6 sigma_p lambda], m  |
| `M`,`Minv`| consistent mass matrix and inverse                              |
| `beta`    | base-excitation vector [beta = M 1]                             |
| `Cstruct` | structural damping matrix [C_s]                                 |
| `y`,`weights` | graded Simpson grid on the overlap                          |
| `B1`,`B2` | shape weights at the grid nodes, B1 + B2 = 1                    |
| `panels`  | number of Simpson panels                                        |
"""
struct Model
    p::Params
    ke::Float64
    k1::Float64
    k3::Float64
    kss::Float64
    gc::Float64
    alpha::Float64
    hd::Float64
    kp::Float64
    M::Matrix{Float64}
    Minv::Matrix{Float64}
    beta::Vector{Float64}
    Cstruct::Matrix{Float64}
    y::Vector{Float64}
    weights::Vector{Float64}
    B1::Vector{Float64}
    B2::Vector{Float64}
    panels::Int
end

require(condition, message) = condition || throw(ArgumentError(message))

function validate(p::Params, panels)
    require(p.n_beams > 0 && p.n_parallel > 0 && p.n_series > 0 && p.n_stop_beams > 0,
        "beam and span counts must be positive")
    require(0 < p.Leff <= p.Lf, "need 0 < Leff <= Lf")
    require(p.g0 > 2 * p.Tp + p.h_eff, "g0 must exceed 2 Tp + h_eff (positive travel to contact)")
    require(panels >= 16, "use at least 16 film panels")
    positive = (p.Tf, p.wt, p.wb, p.ws, p.wss, p.Lsp, p.Lss, p.m_shuttle, p.rho, p.E,
        p.eps0, p.epsr, p.eta, p.h_eff, p.eps_gap, p.eps_wall, p.eps_stop,
        p.seal_width, p.Rload, p.kw, p.pw, p.phs, p.gss)
    require(all(>(0), positive), "a strictly positive parameter is <= 0")
    nonnegative = (p.c1, p.ce, p.cw, p.film_scale, p.gamma3, p.cp, p.Tp,
        p.mean_free_path, p.slip_coefficient, p.khs)
    require(all(>=(0), nonnegative), "a nonnegative parameter is < 0")
    require(p.ghs >= p.gss, "hard stop must not engage before the soft stop")
    return nothing
end

function Model(p::Params = Params(); panels::Int = 512)
    validate(p, panels)
    alpha = isnan(p.gap_slope) ? (p.wb - p.wt) / p.Lf : p.gap_slope
    require(isfinite(alpha) && alpha >= 0, "gap slope must be finite and nonnegative")

    ke = electrode_tip_stiffness(p)
    k1, k3 = suspension_stiffness(p)
    kss = soft_stop_stiffness(p)
    M = consistent_mass_matrix(p, ke)
    require(isposdef(Symmetric(M)), "mass matrix is not positive definite")

    y, weights = graded_simpson_grid(p, alpha, panels)
    B2 = [beam_shape(p, ke, p.Lf - v) for v in y]   # y = Lf - s
    B1 = 1 .- B2

    Cstruct = [p.c1 0.0; 0.0 0.0] + p.n_beams * p.ce * [1.0 -1.0; -1.0 1.0]
    gc = p.g0 - 2 * p.Tp - p.h_eff
    hd = 2 * p.Tp / p.epsr
    kp = 6 * p.slip_coefficient * p.mean_free_path
    beta = M * ones(2)

    return Model(p, ke, k1, k3, kss, gc, alpha, hd, kp, M, inv(M), beta, Cstruct,
        y, weights, B1, B2, panels)
end

# -----------------------------------------------------------------------------
# 6. Constitutive laws
# -----------------------------------------------------------------------------

"""
    cumulative_simpson(m, f)

H(y) = integral of f from 0 to y at every grid node, fourth-order accurate at
panel midpoints and endpoints (TeX, "Spatial quadrature").
"""
function cumulative_simpson(m::Model, f)
    H = zeros(eltype(f), length(f))
    for j in 1:m.panels
        i = 2 * j - 1
        d = m.y[i+2] - m.y[i]
        H[i+1] = H[i] + d * (5 * f[i] + 8 * f[i+1] - f[i+2]) / 24
        H[i+2] = H[i] + d * (f[i] + 4 * f[i+1] + f[i+2]) / 6
    end
    return H
end

"""
    film_matrix_side(m, h, h1, h2, chi)

Generalized squeeze-film matrix of ONE gap of ONE beam, divided by 12 eta Tf,
in the centered Gram form (TeX Eq. `gram`):

    integral of (H - Hbar)(H - Hbar)' / G  +  chi I0 Hbar Hbar',   G = h^2 (h + kp).

`chi` in [0, 1] is the tip boundary closure (0 open, 1 sealed). The form is
positive semidefinite term by term for positive weights, so passivity does not
depend on cancellation. Returns the three independent entries (d11, d12, d22).
"""
function film_matrix_side(m::Model, h, h1, h2, chi)
    H1 = cumulative_simpson(m, h1)
    H2 = cumulative_simpson(m, h2)
    W = m.weights ./ (h .^ 2 .* (h .+ m.kp))   # quadrature weight / G(h)
    I0 = sum(W)
    mu1 = dot(W, H1) / I0                      # Hbar_1
    mu2 = dot(W, H2) / I0                      # Hbar_2
    Z1 = H1 .- mu1
    Z2 = H2 .- mu2
    d11 = dot(W, Z1 .^ 2) + chi * I0 * mu1^2
    d12 = dot(W, Z1 .* Z2) + chi * I0 * mu1 * mu2
    d22 = dot(W, Z2 .^ 2) + chi * I0 * mu2^2
    return d11, d12, d22
end

"Tip sealing fraction chi_r for gap side r (TeX Eq. `seal`)."
function seal_fraction(m::Model, r, x2)
    p = m.p
    p.seal_at_contact || return zero(x2)
    return smootherstep((r * x2 - m.gc + p.seal_width) / (2 * p.seal_width))
end

"""
    constitutive(m, x1, x2; film = true) -> (; C, grad, D)

Total capacitance `C` [F] (TeX Eq. `cap`), its exact coordinate gradient
`grad` = (dC/dx1, dC/dx2) [F/m], and the passive generalized film matrix `D`
[N s/m] summed over both gaps of all beams (TeX Eq. `gram`). With
`film = false` the film matrix is skipped and returned as zeros.
"""
function constitutive(m::Model, x1, x2; film::Bool = true)
    p = m.p
    T = promote_type(typeof(x1), typeof(x2), Float64)
    C = zero(T) + p.cp
    c1 = zero(T)
    c2 = zero(T)
    d11 = zero(T)
    d12 = zero(T)
    d22 = zero(T)
    cap_factor = p.n_beams * p.eps0 * p.Tf
    film_factor = 12 * p.eta * p.Tf * p.n_beams * p.film_scale

    for r in GAP_SIDES
        # Local gap h_r(y) = h_eff + softpos(d_r) and its coordinate derivatives.
        d = m.gc .+ m.alpha .* m.y .- r .* (m.B1 .* x1 .+ m.B2 .* x2)
        h = p.h_eff .+ softpos.(d, p.eps_gap)
        dh = dsoftpos.(d, p.eps_gap)
        h1 = -r .* dh .* m.B1
        h2 = -r .* dh .* m.B2

        # Local dielectric stack: air gap in series with two coatings, strips in parallel.
        stack = h .+ m.hd
        C += cap_factor * dot(m.weights, 1 ./ stack)
        c1 -= cap_factor * dot(m.weights, h1 ./ stack .^ 2)
        c2 -= cap_factor * dot(m.weights, h2 ./ stack .^ 2)

        if film
            s11, s12, s22 = film_matrix_side(m, h, h1, h2, seal_fraction(m, r, x2))
            d11 += film_factor * s11
            d12 += film_factor * s12
            d22 += film_factor * s22
        end
    end
    return (; C, grad = [c1, c2], D = [d11 d12; d12 d22])
end

"""
    springs(m, x1, x2) -> (; U, F)

Conservative structural potential `U` [J] (TeX Eq. `U0`) and generalized force
`F` = -grad U [N]: suspension (linear + cubic), electrode bending, two-sided
soft stops, and the optional hard stop.
"""
function springs(m::Model, x1, x2)
    p = m.p
    rel = x2 - x1
    bending = p.n_beams * m.ke
    U = m.k1 * x1^2 / 2 + m.k3 * x1^4 / 4 + bending * rel^2 / 2
    F1 = -m.k1 * x1 - m.k3 * x1^3 + bending * rel
    F2 = -bending * rel
    for r in GAP_SIDES
        z = r * x1 - p.gss
        s = softpos(z, p.eps_stop)
        U += m.kss * s^2 / 2
        F1 -= r * m.kss * s * dsoftpos(z, p.eps_stop)
        if isfinite(p.ghs)
            z = r * x1 - p.ghs
            s = softpos(z, p.eps_stop)
            U += p.khs * s^(p.phs + 1) / (p.phs + 1)
            F1 -= r * p.khs * s^p.phs * dsoftpos(z, p.eps_stop)
        end
    end
    return (; U, F = [F1, F2])
end

"""
    wall(m, x2, v2) -> (; U, F, loss)

Regularized Hunt-Crossley tip contact on coordinate x2, all beams, both sides
(TeX Eqs. `wallA`-`wallloss`): stored energy `U` [J], force `F` [N], and the
exact dissipated power `loss` >= 0 [W]. On the unloading branch where the gate
max(1 + cw v_r, 0) clips the force to zero the loss is -A v_r > 0, not zero.
"""
function wall(m::Model, x2, v2)
    p = m.p
    T = promote_type(typeof(x2), typeof(v2), Float64)
    U = zero(T)
    F = zero(T)
    loss = zero(T)
    for r in GAP_SIDES
        z = r * x2 - m.gc                       # penetration coordinate delta_r
        vr = r * v2                             # approach velocity
        s = softpos(z, p.eps_wall)
        A = p.kw * s^p.pw * dsoftpos(z, p.eps_wall)   # elastic force magnitude, = dU/d(delta)
        gate = max(1 + p.cw * vr, zero(vr))
        U += p.n_beams * p.kw * s^(p.pw + 1) / (p.pw + 1)
        F -= p.n_beams * r * A * gate
        loss += p.n_beams * A * vr * (gate - 1)
    end
    return (; U, F, loss)
end

# -----------------------------------------------------------------------------
# 7. Equations of motion
# -----------------------------------------------------------------------------

"""
    generalized_forces(m, c, x1, x2, v1, v2, Vc, a) -> NamedTuple

Every generalized mechanical force [N] on (x1, x2), term by term, given the
constitutive result `c = constitutive(m, x1, x2)`, capacitor voltage `Vc`, and
base acceleration `a`. Fields: `spring`, `wall` (length-2, acts on x2 only),
`electrostatic`, `film`, `structural`, `base`, their sum `total`, and the
contact dissipation `wall_loss` [W].

This is the single place where the force balance of TeX Eqs. `eom1`-`eom2` is
assembled. Diagnostics, figures, and discrepancy (UDE) targets should read the
individual terms from here instead of recomputing them.
"""
function generalized_forces(m::Model, c, x1, x2, v1, v2, Vc, a)
    v = [v1, v2]
    s = springs(m, x1, x2)
    w = wall(m, x2, v2)
    spring = s.F
    wall_force = [zero(w.F), w.F]
    electrostatic = 0.5 * Vc^2 * c.grad
    film = -(c.D * v)
    structural = -(m.Cstruct * v)
    base = -m.beta * a
    total = spring + wall_force + electrostatic + film + structural + base
    return (; spring, wall = wall_force, electrostatic, film, structural, base, total,
        wall_loss = w.loss)
end

"""
    HarvesterRHS(model, acceleration; closure = nothing)

Callable right-hand side `f(du, u, theta, t)` in the SciML in-place convention.

* `acceleration` : `t -> a(t)`, base acceleration [m/s^2].
* `closure`      : `nothing`, or `(u, t, theta, model) -> (Q1, Q2)`, an extra
  generalized force [N] added to the two mechanical equations. `theta` is the
  ODE parameter object and is passed through untouched, so it can hold neural
  network weights (e.g. a ComponentArray); the known physics never reads it.

How many states are integrated is decided by `length(du)`:
`N_PHYSICAL` (5) physical only; `N_LEDGER` (12) adds the energy ledger;
`N_LEDGER_CLOSURE` (13) also integrates the closure work, which is required
for the ledger residual to close whenever `closure !== nothing`.
"""
struct HarvesterRHS{A,C}
    model::Model
    acceleration::A
    closure::C
end
HarvesterRHS(model::Model, acceleration; closure = nothing) =
    HarvesterRHS(model, acceleration, closure)

function (f::HarvesterRHS)(du, u, theta, t)
    m = f.model
    p = m.p
    x1, x2, v1, v2, Vout = u[1], u[2], u[3], u[4], u[5]
    Vc = p.Vbias - Vout                       # capacitor voltage
    a = f.acceleration(t)

    c = constitutive(m, x1, x2)
    Q = generalized_forces(m, c, x1, x2, v1, v2, Vc, a)
    Qc1, Qc2 = f.closure === nothing ? (zero(x1), zero(x1)) : f.closure(u, t, theta, m)

    dv = m.Minv * [Q.total[1] + Qc1, Q.total[2] + Qc2]
    du[1] = v1
    du[2] = v2
    du[3] = dv[1]
    du[4] = dv[2]
    # Circuit, TeX Eq. `circuit`: q = C Vc and dq/dt = Vout / R.
    du[5] = -Vout / (p.Rload * c.C) + Vc / c.C * (c.grad[1] * v1 + c.grad[2] * v2)

    if length(du) >= N_LEDGER
        v = [v1, v2]
        P_base = -a * dot(m.beta, v)          # shaker work rate
        P_bias = p.Vbias * Vout / p.Rload     # bias-source work rate (either sign)
        P_load = Vout^2 / p.Rload             # resistor dissipation
        P_film = dot(v, c.D * v)
        P_struct = dot(v, m.Cstruct * v)
        P_closure = Qc1 * v1 + Qc2 * v2       # work rate of the closure force
        du[6] = P_base
        du[7] = P_bias
        du[8] = P_load
        du[9] = P_film
        du[10] = P_struct
        du[11] = Q.wall_loss
        du[12] = abs(P_base) + abs(P_bias) + P_load + P_film + P_struct + Q.wall_loss +
                 abs(P_closure)
        if length(du) >= N_LEDGER_CLOSURE
            du[13] = P_closure
        end
    end
    return nothing
end

"""
    rhs!(du, u, (model, acceleration), t)

Original calling convention, kept so existing scripts run unchanged.
Equivalent to `HarvesterRHS(model, acceleration)(du, u, nothing, t)`.
"""
function rhs!(du, u, context, t)
    m, acceleration = context
    return HarvesterRHS(m, acceleration)(du, u, nothing, t)
end

# -----------------------------------------------------------------------------
# 8. Energy accounting
# -----------------------------------------------------------------------------

"Stored relative-frame energy E(u) [J] (TeX Eq. `E`): kinetic + structural + contact + capacitor."
function energy(m::Model, u)
    c = constitutive(m, u[1], u[2]; film = false)
    v = [u[3], u[4]]
    return dot(v, m.M * v) / 2 + springs(m, u[1], u[2]).U + wall(m, u[2], u[4]).U +
           c.C * (m.p.Vbias - u[5])^2 / 2
end

"""
    ledger(m, u, E0)

Energy-balance residual r_E [J] (TeX Eq. `ledger`); zero for an exact solution.
If `u` carries the 13th state, the closure work is included.
"""
function ledger(m::Model, u, E0)
    residual = energy(m, u) - E0 - u[6] - u[7] + sum(u[8:11])
    return length(u) >= N_LEDGER_CLOSURE ? residual - u[13] : residual
end

# -----------------------------------------------------------------------------
# 9. Scaling and observation helpers (solver conditioning, UDE/HDI interfaces)
# -----------------------------------------------------------------------------

"""
    reference_scales(m) -> (; x, v, V, E, F, t)

Physical scales built from the model alone (no data): contact travel, the
corresponding suspension velocity, bias voltage (at least 1 V), suspension
energy, suspension force at contact, and the suspension time 1/omega_0.
Use these, not data-derived ranges, to nondimensionalize network inputs and
outputs; they stay O(1)-meaningful when a trajectory barely moves.
"""
function reference_scales(m::Model)
    x = m.gc
    omega0 = sqrt(m.k1 / sum(m.M))
    return (; x, v = x * omega0, V = max(abs(m.p.Vbias), 1.0), E = m.k1 * x^2,
        F = m.k1 * x, t = 1 / omega0)
end

"State-scale vector of length `n` (5, 12, or 13) matching `STATE_NAMES`."
function state_scales(m::Model, n::Integer = N_LEDGER)
    s = reference_scales(m)
    return [s.x, s.x, s.v, s.v, s.V, fill(s.E, n - N_PHYSICAL)...]
end

"""
    ScaledRHS(rhs, scale)

Wraps any in-place `rhs(du, u, theta, t)` so that the solver sees z = u ./ scale.
Conditioning matters here because nanometer-scale constitutive transitions
coexist with micrometer travel.
"""
struct ScaledRHS{F}
    rhs::F
    scale::Vector{Float64}
end

function (f::ScaledRHS)(dz, z, theta, t)
    f.rhs(dz, z .* f.scale, theta, t)
    dz ./= f.scale
    return nothing
end

"Measured output: the load-resistor voltage, which is state 5 itself."
observe_voltage(u) = u[5]

"Capacitor charge q = C_t(x) (Vbias - Vout) [C]; dq/dt = Vout/R is directly measurable."
capacitor_charge(m::Model, u) = constitutive(m, u[1], u[2]; film = false).C * (m.p.Vbias - u[5])

# -----------------------------------------------------------------------------
# 10. Names used by the first version of this file
# -----------------------------------------------------------------------------
const gausslegendre = gauss_legendre
const quad = integrate
const width = beam_width
const EI = bending_rigidity
const shape = beam_shape
