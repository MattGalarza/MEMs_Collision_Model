# ------------------------------------------ Libraries --------------------------------------
# collision_model_v3.jl -- the v2 script structure carrying the v3 (physics-complete)
# model. Same enclosed flow as collision_model_v2.jl: AnalyticalModel module with
# Params/create_params, individual force functions, CoupledSystem!, forcing-or-IC
# toggle, solve, then state and force plots.
#
# Physics relative to v2 (all numerically verified against the Python twin
# verify_full_model.py; see corrections memo for provenance):
#   F1  m2: stray Lf removed (was x2222 too light; electrode is 64.7 kHz, not 3 MHz).
#   F2  Coupling spring -ke*(x1-x2) is ALWAYS ON; no branch-switched contact force.
#       Boundary mechanics = unilateral Hunt-Crossley tip compliance (wall force).
#   F3  Electrostatic force sign corrected: Fe = +(q^2/2C^2) dC/dx2 (attractive).
#   F4  Capacitance: approach = translational wedge with root gap (gp-|x2|)+h_eff;
#       contact = wedge ROTATION about the sealed tip, u = |x2|-gp, au = a - u/Leff.
#       Exactly value-continuous at the seam; slope drops ~89x at contact
#       (translation -> rotation kinematics). Replaces the empirical Cmin/Cmax ramp.
#       Face-wrap toward Cmax needs u ~ a*Leff = 21 um; corrected mechanics caps u
#       at tens of nm, so that regime is dynamically unreachable below gamma~O(100).
#   F5  Film: slip-corrected conductance G = h^2(h + 6*sigmap*lambda). Approach =
#       translational squeeze (b = 1.155e-4 N s/m at 50 nm, zeta ~ 1.0); contact =
#       sealed-pivot rotational Reynolds (b_rot = 1.72e-5 N s/m, bounded pivot
#       pressure). Both tabulated at build from direct Reynolds quadrature and
#       C1-blended over the sealing layer |u| < h_eff/2. Replaces the Moy
#       transcription and the torsion-mirror rotational formula.
#   F6  c (damping scale) defaults to 1.0 = PHYSICAL film. The UDE discrepancy
#       model calibrates this later; insertion point is damping()/bside().
#   F7  Vout is an observable (Vbias - Q/Ct), not a state: 5 states, no stiff
#       relaxation equation. The plotting section computes it.
#   F9  (v3.3) Wall activation is C-infinity: d -> softpos(d, epsw), epsw = 0.5 nm.
#       Diagnosed via sol.stats: nreject = 963k vs naccept = 491k with njacs =
#       naccept -- the start-point Jacobian (d < 0 side, dFw/dx = 0) goes stale
#       across every d = 0 crossing, degrading the Rosenbrock step to low order
#       and pinning dt at ~ns around each of ~840 crossings. The softplus keeps
#       the Jacobian valid through the crossing. Physical cost: ~4 nN pre-contact
#       preload (1% of Fe at the floor); energy-ledger impact ~1e-17 J.
#   F8  h_eff (default 50 nm) is the single contact-gap parameter (residual air
#       gap at dielectric contact; measurable via the bias-sweep latch threshold).
#
# Operational notes (hard-won):
#   - NEVER attach a ContinuousCallback on |x2|-gp: at engagement the trajectory
#     skates the seam (thousands of micro-crossings); event root-finding + restart
#     multiplied cost ~50x. Count crossings in post-processing if needed.
#   - Engaged runs (alpha >~ 2.45) resolve ~200 kHz contact motion: expect minutes.
#     Sub-threshold runs are seconds. If the AD Jacobian complains about the film
#     tables, switch to Rodas5P(autodiff = false).
#   - Acceptance (Python-verified): alpha = 1.0 / 2.0 / 2.1 -> max|x1| = 6.098 /
#     11.094 / 11.54 um, no contact. Collision threshold alpha* = 2.45 +- 0.05.

using DifferentialEquations, Plots, Printf

# --------------------------------------- Analytical Model ----------------------------------

module AnalyticalModel
using DifferentialEquations
using Parameters
using LinearAlgebra
export Params, p, create_params, spring, collision, damping, electrostatic, CoupledSystem!

@with_kw mutable struct Params{T<:Real}
    # Fundamental geometric parameters
    g0::T = 14e-6        # Initial gap
    Tp::T = 120e-9       # Parylene-C thickness
    Tf::T = 25e-6        # Electrode thickness
    wt::T = 9e-6         # Electrode width, top
    wb::T = 30e-6        # Electrode width, bottom
    ws::T = 14.7e-6      # Suspension spring width
    wss::T = 14e-6       # Soft-stopper width
    Leff::T = 400e-6     # Effective electrode length
    Lf::T = 450e-6       # Full electrode length
    Lsp::T = 1400e-6     # Suspension spring length
    Lss::T = 1000e-6     # Soft-stopper length
    gss::T = 14e-6       # Soft-stopper position

    # Mass and material properties
    m1::T = 2.0933e-6    # Shuttle mass
    rho::T = 2330.0      # Density of silicon
    E::T = 170e9         # Young's modulus
    e::T = 8.85e-12      # Permittivity of free space
    ep::T = 3.2          # Relative permittivity of Parylene-C
    eta::T = 1.849e-5    # Viscosity of air
    c::T = 1.0           # F6: film scale -- PHYSICAL default (UDE calibrates later)
    c1::T = 0.0          # Optional shuttle damping (Eq. 1 omission; 0 = off)
    lambda::T = 70e-9    # Mean free path of air molecules (m)
    sigmap::T = 1.016    # Slip coefficient for rarefaction

    # Contact / boundary parameters (v3)
    h_eff::T = 50e-9     # F8: residual contact air gap (single parameter)
    epsw::T = 0.5e-9     # wall engagement width: C-infinity softplus activation
                         # (sub-roughness scale; pre-contact preload ~4 nN << Fe)
    kw::T = 1e6          # Hunt-Crossley wall stiffness (N/m^1.5)
    pw::T = 1.5          # Hunt-Crossley exponent
    cw::T = 50.0         # Hunt-Crossley dissipation (s/m)

    # Electrical parameters
    N::Int = 160         # Number of electrodes
    cp::T = 5e-12        # Parasitic capacitance
    Vbias::T = 3.0       # Bias voltage
    Rload::T = 0.42e6    # Load resistance

    # Derived parameters - calculated by create_params()
    gp::T = g0 - 2*Tp        # Initial electrode air gap (travel to contact)
    a::T = (wb - wt)/Leff    # Taper ratio
    crl::T = 0.0             # Dielectric layer capacitance
    kp::T = 0.0              # Slip conductance length 6*sigmap*lambda
    W::T = 0.0               # Sealing-blend half width h_eff/2
    a_min::T = 0.0           # Wedge-slope floor
    m2::T = 0.0              # Modal mass of electrode
    ke::T = 0.0              # Electrode spring constant
    k1::T = 0.0              # Linear spring constant
    k3::T = 0.0              # Cubic spring constant
    kss::T = 0.0             # Soft-stopper spring constant
    # Film tables (built by create_params from direct Reynolds quadrature)
    LH::Vector{T} = T[]      # log gap grid
    LB::Vector{T} = T[]      # log b_translational
    TRA::Vector{T} = T[]     # wedge-slope grid
    TRB::Vector{T} = T[]     # b_rotational (sealed pivot)
    # v3.2: analytic C-infinity film fits (replace piecewise-linear table lookup;
    # the lerp curvature jumps at 240 nodes -- 0.05 nm spacing at 1 nm gap --
    # were forcing sub-ns Jacobian re-resolution during seam skating)
    FTmu::T = 0.0; FTsd::T = 1.0; FTc::Vector{T} = T[]   # log b_trans vs log h
    FRmu::T = 0.0; FRsd::T = 1.0; FRc::Vector{T} = T[]   # log b_rot  vs log au
end

function create_params(p::Params{T}; verbose = true) where T<:Real
    # F1: corrected modal mass (trapezoid volume: rho*Tf*Lf*(wb+wt)/2; no Lf^2)
    modalcoeff = 0.236 + 0.045*(1.0 - p.wt/p.wb)
    m2Physical = 0.5 * p.Lf * p.rho * p.Tf * (p.wb + p.wt)
    p.m2 = modalcoeff * m2Physical

    # Electrode spring constant, ke (tapered Castigliano; verified == direct integral)
    num  = p.E * p.Tf * p.wt^2 * (p.wb - p.wt)^3
    dem  = 6 * p.Lf^3 * ((p.wb - 3*p.wt)*(p.wb - p.wt) +
                         2*p.wt^2*(log(p.Lf*p.wb) - log(p.Lf*p.wt)))
    p.ke = num/dem

    # Suspension / stopper spring constants
    p.k1  = (4.0/6.0) * p.E * p.Tf * p.ws^3 / p.Lsp^3
    p.k3  = (18.0/25.0) * p.E * p.Tf * p.ws / p.Lsp^3
    p.kss = p.E * p.Tf * p.wss^3 / p.Lss^3

    # Electrical / contact derived
    p.crl   = p.e * p.ep * p.Leff * p.Tf / p.Tp
    p.kp    = 6 * p.sigmap * p.lambda
    p.W     = 0.5 * p.h_eff
    p.a_min = p.a * 1e-3

    # F5: film tables from direct Reynolds quadrature (ground truth, built once)
    ny = 8001
    yy = range(0.0, p.Leff; length = ny); dy = step(yy)
    function b_trans_direct(h0)
        I1 = 0.0; I0 = 0.0
        for y in yy
            h = h0 + p.a*y; G = h*h*(h + p.kp)
            I1 += y/G; I0 += 1/G
        end
        I1 *= dy; I0 *= dy
        cc = 12*p.eta*I1/I0
        pr = cc/((h0)^2*(h0 + p.kp))
        pint = 0.0; F = 0.0
        for i in 2:ny
            y = yy[i]; h = h0 + p.a*y; G = h*h*(h + p.kp)
            cur = (-12*p.eta*y + cc)/G
            pint += 0.5*(pr + cur)*dy
            F += pint*dy
            pr = cur
        end
        return abs(p.Tf*F)
    end
    function b_rot_direct(au)
        ps = zeros(ny); pr = 0.0
        for i in 2:ny
            y = yy[i]; h = p.h_eff + au*y; G = h*h*(h + p.kp)
            cur = (-6*p.eta*y*y/p.Leff)/G
            ps[i] = ps[i-1] + 0.5*(pr + cur)*dy
            pr = cur
        end
        pL = ps[ny]; F = 0.0
        for i in 1:ny
            F += (ps[i] - pL)*dy
        end
        return abs(p.Tf*F)
    end
    TBH  = exp.(range(log(1e-9), log(4e-5); length = 240))
    p.LH = log.(TBH)
    p.LB = log.([b_trans_direct(h) for h in TBH])
    p.TRA = collect(range(p.a_min, p.a; length = 64))
    p.TRB = [b_rot_direct(x) for x in p.TRA]

    # v3.2: degree-10 / degree-6 log-log least-squares fits (C-infinity films)
    function _polyfit(xn, y, deg)
        V = [xn[i]^(deg - j) for i in eachindex(xn), j in 0:deg]
        return V \ y
    end
    p.FTmu = sum(p.LH)/length(p.LH)
    p.FTsd = sqrt(sum((p.LH .- p.FTmu).^2)/length(p.LH))
    p.FTc  = _polyfit((p.LH .- p.FTmu)./p.FTsd, p.LB, 10)
    lA     = log.(p.TRA)
    p.FRmu = sum(lA)/length(lA)
    p.FRsd = sqrt(sum((lA .- p.FRmu).^2)/length(lA))
    p.FRc  = _polyfit((lA .- p.FRmu)./p.FRsd, log.(p.TRB), 6)
    _h(x, c) = (r = c[1]; for k in 2:length(c); r = r*x + c[k]; end; r)
    ftres = maximum(abs.(expm1.([_h((p.LH[i]-p.FTmu)/p.FTsd, p.FTc) - p.LB[i] for i in eachindex(p.LH)])))
    frres = maximum(abs.(expm1.([_h((lA[i]-p.FRmu)/p.FRsd, p.FRc) - log(p.TRB[i]) for i in eachindex(lA)])))

    if verbose
        println("\n--- Modal Mass (F1 corrected) ---")
        println("Modal coefficient:  ", modalcoeff)
        println("Physical mass:      ", m2Physical)
        println("Modal mass, m2:     ", p.m2, "  (f2 = ",
                round(sqrt(p.ke/p.m2)/(2*pi)/1e3; digits = 1), " kHz)")
        println("\n--- Springs ---")
        println("ke = ", p.ke, "   k1 = ", p.k1, "   k3 = ", p.k3, "   kss = ", p.kss)
        println("\n--- Film tables (direct Reynolds, slip-corrected) ---")
        println("b_trans(gp)    = ", exp(p.LB[searchsortedlast(p.LH, log(p.gp))]),
                "  (expect ~2.30e-6)")
        println("b_trans(h_eff) = ", exp(p.LB[searchsortedlast(p.LH, log(p.h_eff))]),
                "  (expect ~1.16e-4)")
        println("b_rot(a)       = ", p.TRB[end], "  (expect ~1.72e-5)")
        println("film fit residuals: b_trans ", round(ftres*100; digits = 3),
                " %  b_rot ", round(frres*100; digits = 3), " %  (v3.2 analytic films)")
    end
    return p
end

# ------------------------------ constitutive helpers (v3) ------------------------------
@inline function lerp(x, xs, ys)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i])/(xs[i+1] - xs[i])
    return ys[i] + t*(ys[i+1] - ys[i])
end
@inline function _horner(x, c)
    r = c[1]
    @inbounds for k in 2:length(c)
        r = r*x + c[k]
    end
    return r
end
@inline btrans(h, p) = exp(_horner((log(clamp(h, 1e-9, 4e-5)) - p.FTmu)/p.FTsd, p.FTc))
@inline brot(au, p)  = exp(_horner((log(clamp(au, p.TRA[1], p.TRA[end])) - p.FRmu)/p.FRsd, p.FRc))

@inline function ctrans(u, p)      # u = penetration toward a wall (neg = separated)
    h0 = (u < 0 ? -u : 0.0) + p.h_eff
    C  = (p.e*p.Tf/p.a)*log((h0 + p.a*p.Leff)/h0)
    dC = u < 0 ? (p.e*p.Tf/p.a)*(1/h0 - 1/(h0 + p.a*p.Leff)) : 0.0
    return C, dC
end
@inline function crot(u, p)        # F4: wedge rotation about the sealed tip
    au = p.a - (u > 0 ? u : 0.0)/p.Leff
    au < p.a_min && (au = p.a_min)
    C  = (p.e*p.Tf/au)*log((p.h_eff + au*p.Leff)/p.h_eff)
    dC = u > 0 ? (p.e*p.Tf/p.Leff)*(log((p.h_eff + au*p.Leff)/p.h_eff)/au^2 -
                 p.Leff/(au*(p.h_eff + au*p.Leff))) : 0.0
    return C, dC
end
@inline function cap_side(u, p)    # C1 sealing blend over |u| < W
    u <= -p.W && return ctrans(u, p)
    u >=  p.W && return crot(u, p)
    x  = (u + p.W)/(2*p.W)
    s  = x*x*x*(10 - 15*x + 6*x*x)          # quintic smootherstep: C2 at edges
    ds = 30*x*x*(1 - x)*(1 - x)/(2*p.W)
    Ct_, dCt_ = ctrans(u, p)
    Cr_, dCr_ = crot(u, p)
    return (1 - s)*Ct_ + s*Cr_, (1 - s)*dCt_ + s*dCr_ + ds*(Cr_ - Ct_)
end
@inline function bside(u, p)       # blended film; UDE discrepancy multiplies here
    u <= -p.W && return btrans(-u + p.h_eff, p)
    u >=  p.W && return brot(p.a - u/p.Leff, p)
    x = (u + p.W)/(2*p.W); s = x*x*x*(10 - 15*x + 6*x*x)
    bt = btrans((u < 0 ? -u : 0.0) + p.h_eff, p)
    br = brot(p.a - (u > 0 ? u : 0.0)/p.Leff, p)
    return (1 - s)*bt + s*br
end

# --------------------------------- force functions ---------------------------------

# Suspension spring force, Fsp (+ soft stopper) -- unchanged from v2
function spring(x1, k1, k3, gss, kss)
    Fsp = -k1*x1 - 0.25*k3*x1^3
    Fss = abs(x1) <= gss ? 0.0 : -kss*(abs(x1) - gss)*sign(x1)
    return Fsp + Fss
end

# Electrode coupling + wall, Fc / Fw  (F2: spring always on; wall replaces the
# branch-switched contact force). collision_state kept for diagnostics/flow.
@inline softpos(d, e_) = 0.5*(d + sqrt(d*d + e_*e_))   # C-infinity positive part

function collision(x1, x2, x2dot, p)
    Fc = -p.ke*(x1 - x2)
    ur = x2 - p.gp
    ul = -x2 - p.gp
    dr = softpos(ur, p.epsw)      # ~ur when ur >> epsw; -> 0 like epsw^2/4|ur| far away
    dl = softpos(ul, p.epsw)
    gr = 1 + p.cw*x2dot           # HC gate: inert at |x2dot| << 20 mm/s; its g = 0
    gl = 1 - p.cw*x2dot           # kink is known-benign at operating speeds
    Fw = 0.0
    gr > 0 && (Fw += -p.kw*dr^p.pw*gr)
    gl > 0 && (Fw +=  p.kw*dl^p.pw*gl)
    collision_state = (ur > 0 || ul > 0) ? "contact" : "translational"
    return Fc, Fw, collision_state
end

# Viscous damping, Fd (F5: blended slip-corrected films, both walls).
# Signature keeps the v2 shape; only x2, x2dot enter the film physics.
function damping(x1, x1dot, x2, x2dot, p)
    ur = x2 - p.gp
    ul = -x2 - p.gp
    return -p.c*(bside(ur, p) + bside(ul, p))*x2dot
end

# Electrostatic coupling, Fe (F3: attractive; F4 capacitance; analytic dC)
function electrostatic(x1, x2, Qvar, p)
    Cr, dCr = cap_side(x2 - p.gp, p)
    Cl, dCl = cap_side(-x2 - p.gp, p)
    Vr = 1/(2/p.crl + 1/Cr)
    Vl = 1/(2/p.crl + 1/Cl)
    Cvar = (p.N/2)*(Vr + Vl)
    dC   = (p.N/2)*((Vr/Cr)^2*dCr - (Vl/Cl)^2*dCl)
    Ctotal = Cvar + p.cp
    Fe = +((Qvar^2/(2*Ctotal^2))*dC)/(p.N/2)     # attractive at constant charge
    return Ctotal, Fe
end

# 5 states: x1, x1dot, x2, x2dot, Qvar. Vout = Vbias - Q/Ct is an observable (F7).
function CoupledSystem!(dz, z, p, t, current_acceleration)
    z1, z2, z3, z4, z5 = z
    Fs = spring(z1, p.k1, p.k3, p.gss, p.kss) - p.c1*z2
    Fc, Fw, _ = collision(z1, z3, z4, p)
    Fd = damping(z1, z2, z3, z4, p)
    Ctotal, Fe = electrostatic(z1, z3, z5, p)
    Fext = current_acceleration
    dz[1] = z2
    dz[2] = (Fs + (p.N/2)*Fc)/p.m1 - Fext
    dz[3] = z4
    dz[4] = (-Fc + Fd + Fe + Fw)/p.m2 - Fext
    dz[5] = (p.Vbias - z5/Ctotal)/p.Rload
    return nothing
end

# Initialize a default Params instance and calculate dependent parameters
p = Params{Float64}()
p = create_params(p; verbose = false)

end # module AnalyticalModel

import .AnalyticalModel

# --------------------------------------- External Force ------------------------------------

# Sine Wave External Force
f = 20.0        # Frequency (Hz)
alpha = 2.1     # Applied acceleration constant; contact threshold alpha* = 2.45 +- 0.05
g = 9.81        # Gravitational constant (m/s^2)
A = alpha*g
t_ramp = 0.2    # Ramp-up duration (s)
ramp(t) = t < t_ramp ? t/t_ramp : 1.0
Fext_sine = t -> A*ramp(t)*sin(2*pi*f*t)

# ------------------------------------- Set Input Force ------------------------------------

# Set to `true` to use sine forcing, `false` for a near-contact displaced IC
# (free evolution: one contact episode probe, no external force)
use_sine = true
Fext_input = use_sine ? Fext_sine : (t -> 0.0)

# ------------------------------------ Initialize Parameters --------------------------------

p_new = deepcopy(AnalyticalModel.p)

# Initial conditions
if use_sine
    x10, x10dot, x20, x20dot = 0.0, 0.0, 0.0, 0.0
else
    # Near-contact probe IC (episode-scale evaluation without forcing)
    x20    = p_new.gp - 30e-9
    x20dot = 4e-3
    x10    = x20 + 0.15e-6
    x10dot = 4e-3
end

Ctotal0, _ = AnalyticalModel.electrostatic(x10, x20, 0.0, p_new)
Q0 = p_new.Vbias * Ctotal0
z0 = [x10, x10dot, x20, x20dot, Q0]
tspan  = use_sine ? (0.0, 0.5) : (0.0, 600e-6)
abstol = [1e-13, 1e-10, 1e-13, 1e-10, 1e-16]   # state-scaled (F8 op-notes)
reltol = use_sine ? 1e-5 : 1e-6

# ---------------------------------- Solve Analytical Model ---------------------------------

function CoupledSystem_wrapper!(dz, z, p, t)
    AnalyticalModel.CoupledSystem!(dz, z, p, t, Fext_input(t))
end

eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)

# Rodas5P; NO seam callbacks (see header). If AD complains about the film tables:
# sol = solve(eqn, Rodas5P(autodiff = false); abstol, reltol, maxiters = Int(1e9))
sol = solve(eqn, Rodas5P(); abstol = abstol, reltol = reltol, maxiters = Int(1e9))

println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)
println("Solver stats:  ", sol.stats)   # nreject >> naccept near taps => rejection storm

# --- step-floor diagnosis A/B (run one at a time; compare near-seam min/med dt) ---
# The tap-boundary dt collapse to ~5e-13 s persists after the C-infinity film fix,
# so it is a localized rejection storm, not distributed table roughness. Suspects:
# (a) AD-Jacobian branch mismatch when a step straddles d=0 / the blend edges;
# (b) the estimator at the wall's divergent curvature (d^{3/2} at d=0+).
# Experiment 1:  sol = solve(eqn, Rodas5P(autodiff = false); abstol, reltol, maxiters = Int(1e9))
# Experiment 2:  sol = solve(eqn, FBDF(); abstol = abstol, reltol = reltol, maxiters = Int(1e9))
# If (1) removes the floor -> AD-branch interaction; if only (2) does -> estimator/order.

# ----------------------------------------- Plotting -----------------------------------------

x1    = [u[1] for u in sol.u]
x1dot = [u[2] for u in sol.u]
x2    = [u[3] for u in sol.u]
x2dot = [u[4] for u in sol.u]
Qvar  = [u[5] for u in sol.u]

# Vout and Ctotal as observables (F7)
Ctot = [AnalyticalModel.electrostatic(x1[i], x2[i], Qvar[i], p_new)[1] for i in eachindex(x1)]
V    = [p_new.Vbias - Qvar[i]/Ctot[i] for i in eachindex(Qvar)]

p3 = plot(sol.t, x1,    xlabel = "Time (s)", ylabel = "x1 (m)",    title = "Shuttle Mass Displacement (x1)");    display(p3)
p4 = plot(sol.t, x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)");   display(p4)
p5 = plot(sol.t, x2,    xlabel = "Time (s)", ylabel = "x2 (m)",    title = "Mobile Electrode Displacement (x2)"); display(p5)
p6 = plot(sol.t, x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Velocity (x2dot)"); display(p6)
p7 = plot(sol.t, Qvar,  xlabel = "Time (s)", ylabel = "Qvar (C)",  title = "Charge (Qvar)");                     display(p7)
p8 = plot(sol.t, V,     xlabel = "Time (s)", ylabel = "Vout (V)",  title = "Output Voltage (observable)");       display(p8)

# v3 diagnostics: penetration and total capacitance
pen = [abs(x2[i]) - p_new.gp for i in eachindex(x2)]
p8b = plot(sol.t, pen .* 1e9, xlabel = "Time (s)", ylabel = "|x2|-gp (nm)",
           title = "Penetration (contact when > 0)"); hline!(p8b, [0.0]; ls = :dash, lc = :gray, label = "")
display(p8b)
p8c = plot(sol.t, Ctot .* 1e12, xlabel = "Time (s)", ylabel = "Ctotal (pF)",
           title = "Total Capacitance"); display(p8c)

# Force reconstruction over the solution
Fs_array = Float64[]; Fc_array = Float64[]; Fw_array = Float64[]
Fd_array = Float64[]; Fe_array = Float64[]; CollisionState_array = String[]
for i in eachindex(sol.t)
    z1, z2, z3, z4, z5 = sol.u[i]
    push!(Fs_array, AnalyticalModel.spring(z1, p_new.k1, p_new.k3, p_new.gss, p_new.kss))
    Fc, Fw, st = AnalyticalModel.collision(z1, z3, z4, p_new)
    push!(Fc_array, Fc); push!(Fw_array, Fw); push!(CollisionState_array, st)
    push!(Fd_array, AnalyticalModel.damping(z1, z2, z3, z4, p_new))
    push!(Fe_array, AnalyticalModel.electrostatic(z1, z3, z5, p_new)[2])
end

p9  = plot(sol.t, Fs_array, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension + Soft-stopper Spring Force"); display(p9)
p10 = plot(sol.t, Fc_array, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Electrode Coupling Force (always-on spring)"); display(p10)
p10b = plot(sol.t, Fw_array, xlabel = "Time (s)", ylabel = "Fw (N)", title = "Tip Contact (Hunt-Crossley wall) Force"); display(p10b)
p11 = plot(sol.t, Fd_array, xlabel = "Time (s)", ylabel = "Fd (N)", title = "Viscous Film Force (slip-corrected, blended)"); display(p11)
p12 = plot(sol.t, Fe_array, xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force (attractive)"); display(p12)
p13 = plot(sol.t, Fext_input, xlabel = "Time (s)", ylabel = "Fext (m/s^2)", title = "Applied Base Acceleration"); display(p13)

# =========================================================================================
# ============================== EVALUATION SUITE (v3.2) ==================================
# =========================================================================================
# (1) force-function characterization sweeps      (2) two-cycle steady-state zooms
# (3) composite multi-panels (states / forces)    (4) collision/discontinuity metrics
# All sections are read-only with respect to the solve above.

AM = AnalyticalModel

# ------------------------- shared sampling helpers (interpolant-based) -------------------
# Sample states + observables + reconstructed forces on a uniform grid from sol's
# dense output. Used by the zooms, composites, and metrics so every figure draws
# from the same machinery.
function sample_window(sol, p, t0, t1; dt = 2e-6, Fext = Fext_input)
    tg = collect(max(t0, sol.t[1]):dt:min(t1, sol.t[end]))
    U  = Array(sol(tg))
    n  = length(tg)
    Ct = zeros(n); Fe = zeros(n); Fs = zeros(n); Fc = zeros(n); Fw = zeros(n); Fd = zeros(n)
    for i in 1:n
        z1, z2, z3, z4, z5 = U[1,i], U[2,i], U[3,i], U[4,i], U[5,i]
        Fs[i] = AM.spring(z1, p.k1, p.k3, p.gss, p.kss)
        fc, fw, _ = AM.collision(z1, z3, z4, p)
        Fc[i] = fc; Fw[i] = fw
        Fd[i] = AM.damping(z1, z2, z3, z4, p)
        ct, fe = AM.electrostatic(z1, z3, z5, p)
        Ct[i] = ct; Fe[i] = fe
    end
    V   = [p.Vbias - U[5,i]/Ct[i] for i in 1:n]
    pen = [abs(U[3,i]) - p.gp for i in 1:n]
    ae  = [Fext(t) for t in tg]
    return (; t = tg, x1 = U[1,:], x1dot = U[2,:], x2 = U[3,:], x2dot = U[4,:],
              Q = U[5,:], V, Ct, pen, Fs, Fc, Fw, Fd, Fe, ae)
end

# =========================== (1) FORCE CHARACTERIZATION ==================================
# Each internal force swept over the state parameters it depends on.
#   Fs(x1)            : 1-D (shows cubic hardening + stopper kink at gss)
#   Fc(x1 - x2)       : 1-D in relative displacement (pure linear, slope -ke)
#   Fe(x2)            : 1-D under charge equilibration q = Vb*Ct(x2), for which
#                       Fe = (Vb^2/2)*dC/(N/2) exactly (tau_RC << drive period)
#   Fd(x2, x2dot)     : SEPARABLE: Fd = -c*B(x2)*x2dot -> characterize B(x2) once;
#                       overlay Fd at representative velocities (linearity check)
#   Fw(pen, x2dot)    : affine in x2dot up to the Hunt-Crossley clip -> 3 slices

# dense composite grid for x2: coarse across the span + fine near both walls
x2grid = sort(vcat(range(-p_new.gp + 1e-7, p_new.gp - 1e-7; length = 1200),
                   p_new.gp .- exp10.(range(-9.3, -6.3; length = 400)),
                  -p_new.gp .+ exp10.(range(-9.3, -6.3; length = 400)),
                   p_new.gp .+ range(1e-9, 2e-7; length = 200),
                  -p_new.gp .- range(1e-9, 2e-7; length = 200)))

Fe_ch = zeros(length(x2grid)); Ct_ch = zeros(length(x2grid)); B_ch = zeros(length(x2grid))
for (i, xv) in enumerate(x2grid)
    ct, _ = AM.electrostatic(0.0, xv, 0.0, p_new)
    qeq   = p_new.Vbias*ct
    _, fe = AM.electrostatic(0.0, xv, qeq, p_new)
    Fe_ch[i] = fe; Ct_ch[i] = ct
    B_ch[i]  = AM.bside(xv - p_new.gp, p_new) + AM.bside(-xv - p_new.gp, p_new)
end
x1grid = range(-16e-6, 16e-6; length = 1601)
Fs_ch  = [AM.spring(xv, p_new.k1, p_new.k3, p_new.gss, p_new.kss) for xv in x1grid]
rgrid  = range(-0.5e-6, 0.5e-6; length = 401)
Fc_ch  = [-p_new.ke*r for r in rgrid]
pengrid = range(0.0, 60e-9; length = 301)
Fw_m10 = [AM.collision(0.0,  p_new.gp + d, -10e-3, p_new)[2] for d in pengrid]
Fw_0   = [AM.collision(0.0,  p_new.gp + d,   0.0,  p_new)[2] for d in pengrid]
Fw_p10 = [AM.collision(0.0,  p_new.gp + d, +10e-3, p_new)[2] for d in pengrid]
vgrid  = [-10e-3, -3e-3, 3e-3, 10e-3]

c1a = plot(x1grid .* 1e6, Fs_ch .* 1e6; lw = 1.2, legend = false,
           xlabel = "x1 (um)", ylabel = "Fs (uN)",
           title = "Suspension + stopper: cubic hardening, kink at gss")
vline!(c1a, [p_new.gss*1e6, -p_new.gss*1e6]; ls = :dash, lc = :gray)
c1b = plot(rgrid .* 1e9, Fc_ch .* 1e6; lw = 1.2, legend = false,
           xlabel = "x1 - x2 (nm)", ylabel = "Fc (uN)",
           title = "Coupling spring (slope -ke)")
c1c = plot(x2grid .* 1e6, Fe_ch .* 1e6; lw = 1.0, legend = false,
           xlabel = "x2 (um)", ylabel = "Fe (uN)",
           title = "Electrostatic, q equilibrated (full span)")
useam = (x2grid .- p_new.gp)
msk = (useam .> -5e-7) .& (useam .< 2e-7)
c1d = plot(useam[msk] .* 1e9, Fe_ch[msk] .* 1e6; lw = 1.2, legend = false,
           xlabel = "u = x2 - gp (nm)", ylabel = "Fe (uN)",
           title = "Fe across the seam (89x slope drop)")
vline!(c1d, [0.0]; ls = :dash, lc = :gray)
c1e = plot(x2grid .* 1e6, B_ch; yscale = :log10, lw = 1.0, legend = false,
           xlabel = "x2 (um)", ylabel = "B(x2) (N s/m)",
           title = "Film coefficient (Fd = -c*B(x2)*x2dot)")
c1f = plot(pengrid .* 1e9, Fw_m10 .* 1e6; lw = 1.2, label = "x2dot = -10 mm/s",
           xlabel = "penetration (nm)", ylabel = "Fw (uN)",
           title = "Wall: affine in x2dot up to HC clip")
plot!(c1f, pengrid .* 1e9, Fw_0 .* 1e6; lw = 1.2, label = "0")
plot!(c1f, pengrid .* 1e9, Fw_p10 .* 1e6; lw = 1.2, label = "+10 mm/s")
char_fig = plot(c1a, c1b, c1c, c1d, c1e, c1f; layout = (3, 2), size = (1100, 1100),
                titlefontsize = 9, guidefontsize = 8, tickfontsize = 7)
display(char_fig)

# damping linearity check at three gaps (separability made explicit)
vv = range(-20e-3, 20e-3; length = 41)
cl = plot(; xlabel = "x2dot (mm/s)", ylabel = "Fd (uN)",
          title = "Fd linear in velocity at fixed x2 (separability)",
          titlefontsize = 10)
for x2v in (0.0, p_new.gp - 1e-6, p_new.gp - 1e-7)
    plot!(cl, vv .* 1e3, [AM.damping(0.0, 0.0, x2v, v, p_new)*1e6 for v in vv];
          lw = 1.2, label = string("gap = ", round((p_new.gp - x2v)*1e9; digits = 0), " nm"))
end
display(cl)

# ---- capacitance review: monotonicity at the seam + deep-wrap asymptotes ----
# Characterization-only helper: the paper/v2 branch capacitance (Cmin/Cmax ramp)
function C_paper(x2, p)
    Cmin = 7.105299639935359e-14
    Cmax = 9.95200974248769e-11
    kk   = 2*p.Tp/(p.a*p.Leff)
    if abs(x2) < p.gp
        Cair_r = (p.e*p.Tf/p.a)*log((p.gp - x2 + p.a*p.Leff)/(p.gp - x2))
        Cair_l = (p.e*p.Tf/p.a)*log((p.gp + x2 + p.a*p.Leff)/(p.gp + x2))
        return (p.N/2)*(1/(2/p.crl + 1/Cair_r) + 1/(2/p.crl + 1/Cair_l)) + p.cp
    else
        u   = abs(x2) - p.gp
        Cc  = Cmin + (Cmax - Cmin)*log(1 + kk*u)/log(1 + kk*p.a*p.Leff)
        Cnc = (p.e*p.Tf/p.a)*log((p.gp + abs(x2) + p.a*p.Leff)/(p.gp + abs(x2)))
        return (p.N/2)*(1/(2/p.crl + 1/Cc) + 1/(2/p.crl + 1/Cnc)) + p.cp
    end
end

# (i) seam window: is C monotone? (answers the "up then down" reading directly)
useam2 = range(-3e-7, 3e-7; length = 2001)
Cv3_s  = [AM.electrostatic(0.0, p_new.gp + u, 0.0, p_new)[1] for u in useam2]
Cpp_s  = [C_paper(p_new.gp + u, p_new) for u in useam2]
cc1 = plot(useam2 .* 1e9, Cv3_s .* 1e12; lw = 1.4, label = "v3 (wedge + blend)",
           xlabel = "u = x2 - gp (nm)", ylabel = "C_total (pF)",
           title = "C at the seam: both models MONOTONE (the up/down is Fe, i.e. dC)")
plot!(cc1, useam2 .* 1e9, Cpp_s .* 1e12; lw = 1.2, ls = :dash, label = "paper Cmin/Cmax branches")
vspan!(cc1, [-p_new.W*1e9, p_new.W*1e9]; alpha = 0.12, color = :orange, label = "sealing blend")
display(cc1)

# (ii) deep-wrap comparison: asymptote set by the PRESSED-contact residual gap
uw = range(1e-9, p_new.a*p_new.Leff*0.999; length = 1200)
Cv3_w = [AM.electrostatic(0.0, p_new.gp + u, 0.0, p_new)[1] for u in uw]
Cpp_w = [C_paper(p_new.gp + u, p_new) for u in uw]
Ccrl_limit = (p_new.N/2)*(p_new.crl/2) + p_new.cp     # air gap -> 0 idealization
cc2 = plot(uw .* 1e6, Cv3_w .* 1e12; lw = 1.4, label = "v3 (pressed gap = h_eff)",
           xlabel = "root advance u (um)", ylabel = "C_total (pF)",
           title = "Deep wrap: monotone rise; asymptote = pressed-contact air gap")
plot!(cc2, uw .* 1e6, Cpp_w .* 1e12; lw = 1.2, ls = :dash, label = "paper Cmin/Cmax ramp")
hline!(cc2, [Ccrl_limit*1e12]; ls = :dot, lc = :red,
       label = "dielectric-only limit (pressed gap -> 0)")
vspan!(cc2, [0.0, 0.06]; alpha = 0.18, color = :green,
       label = "dynamically reachable (u <= ~60 nm)")
display(cc2)

# (iii) capacitance through the displacement range up to collision and past
xs_a = range(-p_new.gp + 1e-7, p_new.gp + 1e-7; length = 3000)
cc3a = plot(xs_a .* 1e6, [AM.electrostatic(0.0, x, 0.0, p_new)[1]*1e12 for x in xs_a];
            lw = 1.3, legend = false, xlabel = "x2 (um)", ylabel = "C_total (pF)",
            title = "v3 capacitance, full travel (rest 5.62 -> contact 7.18 pF)")
vline!(cc3a, [p_new.gp*1e6, -p_new.gp*1e6]; ls = :dash, lc = :gray)
xs_b = range(p_new.gp - 2e-6, p_new.gp + 1e-7; length = 3000)
cc3b = plot((xs_b .- p_new.gp) .* 1e9,
            [AM.electrostatic(0.0, x, 0.0, p_new)[1]*1e12 for x in xs_b];
            lw = 1.4, legend = false, xlabel = "u = x2 - gp (nm)", ylabel = "C_total (pF)",
            title = "approach + past contact: max slope AT contact, then 89x kink")
vspan!(cc3b, [-p_new.W*1e9, p_new.W*1e9]; alpha = 0.15, color = :orange, label = "")
vline!(cc3b, [0.0]; ls = :dash, lc = :gray)
xs_c = range(p_new.gp - 1.5e-7, p_new.gp + 1e-7; length = 2500)
cc3c = plot((xs_c .- p_new.gp) .* 1e9,
            [AM.electrostatic(0.0, x, 0.0, p_new)[1]*1e12 for x in xs_c];
            lw = 1.5, label = "v3", xlabel = "u (nm)", ylabel = "C_total (pF)",
            title = "micro-window: monotone through the seam")
plot!(cc3c, (xs_c .- p_new.gp) .* 1e9, [C_paper(x, p_new)*1e12 for x in xs_c];
      lw = 1.1, ls = :dash, label = "paper branches")
vspan!(cc3c, [0.0, 60.0]; alpha = 0.12, color = :green, label = "reached (u <= ~7 nm)")
display(plot(cc3a, cc3b, cc3c; layout = (3, 1), size = (900, 1050),
             titlefontsize = 9, guidefontsize = 8))

# ==================== (2)+(3) TWO-CYCLE ZOOMS AND COMPOSITE PANELS =======================
Tdrive = 1/f
Wfull = sample_window(sol, p_new, sol.t[1], sol.t[end]; dt = 2e-6)
Wzoom = sample_window(sol, p_new, sol.t[end] - 2*Tdrive, sol.t[end]; dt = 2e-6)
# NOTE: 2 us sampling resolves episode structure but ALIASES the ~200 kHz contact
# micro-ring; the alias-free view is the metrics-section episode window (0.1 us).

function state_panels(Wd, tag)
    ts = Wd.t .* 1e3
    ps = [plot(ts, Wd.x1 .* 1e6;    ylabel = "x1 (um)",    legend = false),
          plot(ts, Wd.x1dot .* 1e3; ylabel = "x1dot (mm/s)", legend = false),
          plot(ts, Wd.x2 .* 1e6;    ylabel = "x2 (um)",    legend = false),
          plot(ts, Wd.x2dot .* 1e3; ylabel = "x2dot (mm/s)", legend = false),
          plot(ts, Wd.Q .* 1e12;    ylabel = "Q (pC)",     legend = false),
          plot(ts, Wd.V .* 1e3;     ylabel = "Vout (mV)",  legend = false),
          plot(ts, Wd.pen .* 1e9;   ylabel = "|x2|-gp (nm)", legend = false),
          plot(ts, Wd.Ct .* 1e12;   ylabel = "Ct (pF)", xlabel = "t (ms)", legend = false)]
    hline!(ps[7], [0.0]; ls = :dash, lc = :gray)
    hline!(ps[3], [p_new.gp*1e6, -p_new.gp*1e6]; ls = :dash, lc = :gray)
    fig = plot(ps...; layout = (4, 2), size = (1100, 1300),
               plot_title = string("States - ", tag),
               titlefontsize = 9, guidefontsize = 8, tickfontsize = 7)
    return fig
end
function force_panels(Wd, tag)
    ts = Wd.t .* 1e3
    ps = [plot(ts, Wd.Fs .* 1e6; ylabel = "Fs (uN)", legend = false),
          plot(ts, Wd.Fc .* 1e6; ylabel = "Fc (uN)", legend = false),
          plot(ts, Wd.Fw .* 1e6; ylabel = "Fw (uN)", legend = false),
          plot(ts, Wd.Fd .* 1e6; ylabel = "Fd (uN)", legend = false),
          plot(ts, Wd.Fe .* 1e6; ylabel = "Fe (uN)", legend = false),
          plot(ts, Wd.ae;        ylabel = "a_ext (m/s^2)", xlabel = "t (ms)", legend = false)]
    fig = plot(ps...; layout = (3, 2), size = (1100, 1100),
               plot_title = string("Forces - ", tag),
               titlefontsize = 9, guidefontsize = 8, tickfontsize = 7)
    return fig
end

display(state_panels(Wfull, "full span"))
display(force_panels(Wfull, "full span"))
display(state_panels(Wzoom, "last 2 drive cycles (steady state)"))
display(force_panels(Wzoom, "last 2 drive cycles (steady state)"))

# individual two-cycle versions of each original plot (p3z..p13z)
tzs = Wzoom.t .* 1e3
p3z  = plot(tzs, Wzoom.x1;         xlabel = "t (ms)", ylabel = "x1 (m)",    title = "x1 - 2 cycles");    display(p3z)
p4z  = plot(tzs, Wzoom.x1dot;      xlabel = "t (ms)", ylabel = "x1dot",     title = "x1dot - 2 cycles"); display(p4z)
p5z  = plot(tzs, Wzoom.x2;         xlabel = "t (ms)", ylabel = "x2 (m)",    title = "x2 - 2 cycles");    display(p5z)
p6z  = plot(tzs, Wzoom.x2dot;      xlabel = "t (ms)", ylabel = "x2dot",     title = "x2dot - 2 cycles"); display(p6z)
p7z  = plot(tzs, Wzoom.Q;          xlabel = "t (ms)", ylabel = "Q (C)",     title = "Q - 2 cycles");     display(p7z)
p8z  = plot(tzs, Wzoom.V;          xlabel = "t (ms)", ylabel = "Vout (V)",  title = "Vout - 2 cycles");  display(p8z)
p9z  = plot(tzs, Wzoom.Fs;         xlabel = "t (ms)", ylabel = "Fs (N)",    title = "Fs - 2 cycles");    display(p9z)
p10z = plot(tzs, Wzoom.Fc;         xlabel = "t (ms)", ylabel = "Fc (N)",    title = "Fc - 2 cycles");    display(p10z)
p10bz= plot(tzs, Wzoom.Fw;         xlabel = "t (ms)", ylabel = "Fw (N)",    title = "Fw - 2 cycles");    display(p10bz)
p11z = plot(tzs, Wzoom.Fd;         xlabel = "t (ms)", ylabel = "Fd (N)",    title = "Fd - 2 cycles");    display(p11z)
p12z = plot(tzs, Wzoom.Fe;         xlabel = "t (ms)", ylabel = "Fe (N)",    title = "Fe - 2 cycles");    display(p12z)
p13z = plot(tzs, Wzoom.ae;         xlabel = "t (ms)", ylabel = "a_ext",     title = "Forcing - 2 cycles"); display(p13z)

# ======================= (4) COLLISION / DISCONTINUITY METRICS ===========================
println("\n================ COLLISION / DISCONTINUITY METRICS ================")
pen_c = Wfull.pen
engaged = any(pen_c .> 0)
if !engaged
    println("No contact this run (closest approach ",
            round(-maximum(pen_c)*1e9; digits = 1), " nm). Metrics skipped.")
else
    # --- episode structure from the coarse grid ---
    inc = pen_c .> 0
    onset = Wfull.t[findfirst(inc)]
    # cluster engaged samples into episodes (gap > 1 ms)
    tclist = Wfull.t[inc]
    epb = [1]
    for i in 2:length(tclist)
        (tclist[i] - tclist[i-1]) > 1e-3 && push!(epb, i)
    end
    n_ep = length(epb)
    duty = count(inc)/length(inc)

    # --- alias-free episode window: 0.1 us sampling around the FIRST episode ---
    tw0 = max(onset - 2e-5, sol.t[1])
    # extend window to episode end: first 0.2 ms of sustained separation, cap 6 ms
    iend = findfirst(i -> Wfull.t[i] > onset &&
                          all(.!inc[i:min(i + 100, length(inc))]),
                     1:length(inc))
    tw1 = min(iend === nothing ? onset + 6e-3 : Wfull.t[iend] + 1e-4,
              min(onset + 6e-3, sol.t[end]))
    Wep = sample_window(sol, p_new, tw0, tw1; dt = 1e-7)
    penz = Wep.pen
    crossings = count(!=(0), diff(sign.(penz)))
    # impact velocities at inward crossings
    vimp = Float64[]
    for i in 2:length(penz)
        if penz[i] > 0 && penz[i-1] <= 0
            push!(vimp, abs(Wep.x2dot[i]))
        end
    end
    # time-in-blend-layer fraction (seam skating occupancy)
    ublend = abs.(abs.(Wep.x2) .- p_new.gp) .< p_new.W
    fblend = count(ublend)/length(ublend)
    # contact-phase ring frequency: zero crossings of detrended relative coord
    mR = (Wep.t .> onset + 2e-6) .& (Wep.t .< onset + 6e-5)
    if count(mR) > 20
        r = (Wep.x1 .- Wep.x2)[mR]
        tt_ = Wep.t[mR]
        A_ = [ones(length(tt_)) tt_ tt_.^2]
        r = r .- A_*(A_\r)
        f_ring = count(!=(0), diff(sign.(r)))/(2*(tt_[end] - tt_[1]))
    else
        f_ring = NaN
    end
    # per-episode effective restitution (electrode KE out/in) within the window
    KEe = 0.5*p_new.m2 .* Wep.x2dot.^2
    iin  = findfirst(penz .> 0)
    iout = findlast(penz .> 0)
    e_eff = (iin === nothing || iout === nothing || iout >= length(penz)) ? NaN :
            sqrt(max(KEe[min(iout+1, end)], 0)/max(KEe[max(iin-1, 1)], 1e-300))

    # --- solver behavior across the seam (from accepted steps) ---
    dts = diff(sol.t)
    x2s = [u[3] for u in sol.u]
    nearseam = [abs(abs(0.5*(x2s[i] + x2s[i+1])) - p_new.gp) < 5*p_new.W for i in 1:length(dts)]
    dt_near = dts[nearseam]; dt_far = dts[.!nearseam]
    @printf("onset               : %.3f ms\n", onset*1e3)
    @printf("episodes / duty     : %d / %.2f %%   (coarse grid, 2 us)\n", n_ep, 100*duty)
    @printf("max penetration     : %.1f nm\n", maximum(pen_c)*1e9)
    @printf("first-episode window: %.3f-%.3f ms at 0.1 us (alias-free)\n", tw0*1e3, tw1*1e3)
    @printf("seam crossings (win): %d      blend-layer occupancy: %.1f %%\n", crossings, 100*fblend)
    if !isempty(vimp)
        @printf("impact |x2dot|      : n=%d  mean %.2f mm/s  max %.2f mm/s\n",
                length(vimp), 1e3*sum(vimp)/length(vimp), 1e3*maximum(vimp))
    end
    mc = penz .> 0
    tin = count(mc)*1e-7
    f_contact = tin > 2e-6 ? count(!=(0), diff(sign.(Wep.x2dot[mc])))/(2*tin) : NaN
    f_bounce  = count(!=(0), diff(sign.(penz .- 0.5*(maximum(penz) + minimum(penz)))))/
                (2*(Wep.t[end] - Wep.t[1]))
    @printf("ring, window (zc)   : %.1f kHz  | contact-conditioned: %.1f kHz (predicted 100-160 kHz at <=7 nm pen)\n", f_ring/1e3, f_contact/1e3)
    @printf("stopper bounce      : %.2f kHz  (rigid-pin hybrid limit 4.68 kHz; duty-cycled tether)\n", f_bounce/1e3)
    @printf("per-tap restitution : e_eff = %.3f  (electrode KE out/in; window = first tap,\n", e_eff)
    @printf("                      since 0.2 ms separation splits taps of one episode)\n")
    @printf("accepted steps      : %d total | near-seam min/med dt = %.1e / %.1e s | far = %.1e / %.1e s\n",
            length(sol.t), minimum(dt_near), sort(dt_near)[div(end,2)+1],
            minimum(dt_far), sort(dt_far)[div(end,2)+1])

    # --- energy ledger (global correctness of the blended discontinuity) ---
    function ledger(Wd, p; Fext = Fext_input)
        n = length(Wd.t)
        E = zeros(n); Pn = zeros(n)
        for i in 1:n
            x1_, v1_, x2_, v2_, q_ = Wd.x1[i], Wd.x1dot[i], Wd.x2[i], Wd.x2dot[i], Wd.Q[i]
            d = max(abs(x2_) - p.gp, 0.0)
            KE = 0.5*p.m1*v1_^2 + (p.N/2)*0.5*p.m2*v2_^2
            Vm = 0.5*p.k1*x1_^2 + (0.25/4)*p.k3*x1_^4 +
                 (abs(x1_) > p.gss ? 0.5*p.kss*(abs(x1_) - p.gss)^2 : 0.0) +
                 (p.N/2)*0.5*p.ke*(x1_ - x2_)^2 + (p.N/2)*(p.kw/(p.pw + 1))*d^(p.pw + 1)
            E[i] = KE + Vm + q_^2/(2*Wd.Ct[i])
            qd = (p.Vbias - q_/Wd.Ct[i])/p.Rload
            Bv = AM.bside(x2_ - p.gp, p) + AM.bside(-x2_ - p.gp, p)
            s_ = sign(x2_)
            gate = (d > 0) && (1 + p.cw*v2_*s_ > 0)
            Pw = gate ? (p.N/2)*p.kw*d^p.pw*p.cw*v2_^2 : 0.0
            Pn[i] = -Fext(Wd.t[i])*(p.m1*v1_ + (p.N/2)*p.m2*v2_) +
                    p.Vbias*qd - p.Rload*qd^2 - (p.N/2)*p.c*Bv*v2_^2 - p.c1*v1_^2 - Pw
        end
        work = 0.0; resid = 0.0
        for i in 2:n
            work += 0.5*(Pn[i] + Pn[i-1])*(Wd.t[i] - Wd.t[i-1])
            resid = max(resid, abs((E[i] - E[1]) - work))
        end
        thru = maximum(abs.(E)) + maximum(abs, Pn)*(Wd.t[end] - Wd.t[1])
        return resid/thru
    end
    @printf("energy ledger       : episode window %.2e | full span (2 us, aliased) %.2e\n",
            ledger(Wep, p_new), ledger(Wfull, p_new))

    # --- discontinuity-handling figures ---
    m4 = plot((Wep.t .- onset) .* 1e6, penz .* 1e9; lw = 0.7, legend = false,
              xlabel = "t - onset (us)", ylabel = "|x2|-gp (nm)",
              title = "Episode (alias-free): seam skate + contact ring")
    hline!(m4, [0.0]; ls = :dash, lc = :gray)
    hline!(m4, [p_new.W*1e9, -p_new.W*1e9]; ls = :dot, lc = :orange)
    display(m4)
    m5 = plot(penz .* 1e9, Wep.x2dot .* 1e3; lw = 0.5, legend = false,
              xlabel = "|x2|-gp (nm)", ylabel = "x2dot (mm/s)",
              title = "Phase portrait at the seam (crossing dynamics)")
    vline!(m5, [0.0]; ls = :dash, lc = :gray)
    display(m5)
    m6 = plot(sol.t[2:end] .* 1e3, log10.(dts); lw = 0.4, legend = false,
              xlabel = "t (ms)", ylabel = "log10 dt (s)",
              title = "Accepted step size (cost of the discontinuity)")
    display(m6)
    if !isempty(vimp)
        m7 = histogram(vimp .* 1e3; bins = 20, legend = false,
                       xlabel = "impact |x2dot| (mm/s)", ylabel = "count",
                       title = "Impact-velocity distribution (first episode)")
        display(m7)
    end
end