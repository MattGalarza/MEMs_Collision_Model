# ------------------------------------------ Libraries --------------------------------------

using DifferentialEquations, Plots, Printf

# Define high-quality theme for journal publication
function set_journal_theme()
    default(
        fontfamily="Computer Modern",  # LaTeX-like font
        linewidth=1.15,                 # Thicker lines
        foreground_color_legend=nothing, # Transparent legend background
        background_color_legend=nothing, # Transparent legend background
        legendfontsize=10,             # Legend font size
        guidefontsize=12,              # Axis label font size
        tickfontsize=10,               # Tick label font size
        titlefontsize=14,              # Title font size
        size=(800, 600),               # Figure size
        dpi=600,                       # High DPI for print quality
        grid=false,                    # No grid by default
        framestyle=:box,               # Box-style frame
        foreground_color_axis=:black,  # Black axes
        tick_direction=:out,           # Ticks pointing outward
        palette=:default               # Color palette
    )
end
set_journal_theme()

# --------------------------------------- Analytical Model ----------------------------------

module AnalyticalModel
using DifferentialEquations
using Parameters
using LinearAlgebra
export Params, p, create_params, spring, collision, damping, electrostatic, CoupledSystem!

@with_kw mutable struct Params{T<:Real}
    # Fundamental geometric parameters
    g0::T = 14e-6  # Initial gap
    Tp::T = 120e-9  # Parylene-C thickness
    Tf::T = 25e-6  # Electrode thickness
    wt::T = 9e-6  # Electrode width, top
    wb::T = 30e-6  # Electrode width, bottom
    ws::T = 14.7e-6  # Suspension spring width
    wss::T = 14e-6  # Soft-stopper width
    Leff::T = 400e-6  # Effective electrode length
    Lf::T = 450e-6  # Full electrode length
    Lsp::T = 1400e-6  # Suspension spring length
    Lss::T = 1000e-6  # Soft-stopper length
    gss::T = 14e-6  # Soft-stopper position

    # Mass and material properties
    m1::T = 2.0933e-6  # Shuttle mass
    rho::T = 2330.0  # Density of silicon
    E::T = 170e9  # Young's modulus
    e::T = 8.85e-12  # Permittivity of free space
    ep::T = 3.2  # Relative permittivity of Parylene-C
    eta::T = 1.849e-5  # Viscosity of air
    c::T = 1.0  # F6: film scale -- PHYSICAL default (UDE calibrates later)
    c1::T = 0.001  # Shuttle damping [N s/m]. CAUTION: 0.05 gives
    lambda::T = 70e-9  # Mean free path of air molecules (m)
    sigmap::T = 1.016  # Slip coefficient for rarefaction

    # Model-variant switches (v3.8) -----------------------------------------
    cap_model::Symbol = :v3   # :v3 = derived rotation branch; :ramp = paper
                              # Eq. 23 Cmin->Cmax log ramp, continuity-repaired
                              # (A/B experiment; prediction: :ramp latches at
                              # first contact since Fe_contact ~ 19 uN >> ke*Rc)
    Wfac::T = 0.5             # blend half-width W = Wfac*h_eff. 0.5 = legacy;
                              # 0.2-0.3 tightens the window for crisper
                              # stick-unstick (mild near-seam step-cost rise;
                              # keep epsc >= 1 nm when tightening)
    Cmax_ramp::T = 9.952e-11  # ramp asymptote (paper's Cmax, per pair side)
    # Film-hypothesis hooks (v3.9): both default to 1.0 = physics as derived.
    crot_scale::T = 1.0       # multiplies the in-contact (sealed-pivot) film.
                              # >>1 encodes trapped-gas compressibility at
                              # contact (sigma_sq >> 1 there): kills ring-in,
                              # gives the clean touch of the press hypothesis.
    vent_open::T = 1.0        # scales film damping on the OPENING stroke
                              # (venting / inlet-limit asymmetry). <1 frees the
                              # post-detach ring-out at f2 ~ 64.7 kHz.
    veps::T = 1e-5            # C-infinity velocity gate width for vent_open
                              # (softsign over 10 um/s; a hard sign switch on
                              # x2dot would be a derivative jump on a manifold
                              # the trajectory crosses constantly)

    # Contact / boundary parameters (v3) ------------------------------------
    # COLLISION-WINDOW TUNING GUIDE (v3.7). Ordering constraint:
    #   epsw < epsc < W < h_eff   (activation < sealing gate < blend < floor).
    # h_eff : residual pressed-contact air gap. Sets THREE things at once: the
    #         Fe adhesion floor (~1/h_eff), the film floor b_tr(h_eff), and the
    #         seam location of all contact physics. Raising it calms episodes
    #         (weaker adhesion), lowers pkV, delays onset slightly; it is the
    #         latch-bias-measurable parameter. W = h_eff/2 follows automatically.
    # epsc  : width over which the Fe slope makes its 89x translational->
    #         rotational transition. Too small (<1 nm): step cost returns near
    #         the seam. Too large (>~W/3): smears the kinematic kink and
    #         inflates in-contact coupling. 1-5 nm is the sane band.
    # epsw  : wall activation width (sub-asperity). Keep smallest; it only
    #         rounds the last fraction-of-nm of Fw onset. 0.2-1 nm.
    # kw,pw : wall stiffness/exponent. Penetration depth scales like
    #         (F_press/kw)^(1/pw); contact-ring frequency like
    #         sqrt(pw*kw*delta^(pw-1)/m2)/2pi. Changing pw requires re-tuning
    #         kw to match stiffness at the operating penetration (e.g. pw=2
    #         with kw~1.2e10 matches the 7 nm point).
    # cw    : Hunt-Crossley loss. Directly sets per-tap restitution e_eff and
    #         how fast the grazing bounce decays; 10-200 s/m explorable.
    # Wfac  : window-width factor, W = Wfac*h_eff. The direct 'tighter window'
    #         control, now decoupled from the physics floors: h_eff moves the
    #         adhesion well and film floor; Wfac moves ONLY the blend span.
    # cap_model : capacitance A/B. :v3 (derived) vs :ramp (paper Eq. 23 shape,
    #         value-anchored at the seam). :ramp raises the in-contact slope
    #         ~58x -- expect latch-on-first-contact; that is the test.
    # Derived, do not hand-set: W (= Wfac*h_eff), a_min, kp.
    h_eff::T = 50e-9     # F8: residual contact air gap (single parameter)
    epsc::T = 2e-9       # v3.5: sealing-transition width for the Fe derivative
                         # gates (removes the ~0.18 uN Fe JUMP at u = 0 that pinned
                         # dt at ~1e-13 s; diagnosed by the dt localizer)
    epsw::T = 0.5e-9     # wall engagement width: C-infinity softplus activation
                         # (sub-roughness scale; pre-contact preload ~4 nN << Fe)
    kw::T = 1e6          # Hunt-Crossley wall stiffness (N/m^1.5)
    pw::T = 1.5          # Hunt-Crossley exponent
    cw::T = 50.0         # Hunt-Crossley dissipation (s/m)

    # Electrical parameters
    N::Int = 160  # Number of electrodes
    cp::T = 5e-12  # Parasitic capacitance
    Vbias::T = 3.0  # Bias voltage
    Rload::T = 0.42e6  # Load resistance

    # Derived parameters - calculated by create_params()
    gp::T = g0 - 2*Tp  # Initial electrode air gap (travel to contact)
    a::T = (wb - wt)/Leff  # Taper ratio
    crl::T = 0.0  # Dielectric layer capacitance
    kp::T = 0.0  # Slip conductance length 6*sigmap*lambda
    W::T = 0.0  # Sealing-blend half width h_eff/2
    a_min::T = 0.0  # Wedge-slope floor
    m2::T = 0.0  # Modal mass of electrode
    ke::T = 0.0  # Electrode spring constant
    k1::T = 0.0  # Linear spring constant
    k3::T = 0.0  # Cubic spring constant
    kss::T = 0.0  # Soft-stopper spring constant
    # Film tables (built by create_params from direct Reynolds quadrature)
    LH::Vector{T} = T[]  # log gap grid
    LB::Vector{T} = T[]  # log b_translational
    TRA::Vector{T} = T[]  # wedge-slope grid
    TRB::Vector{T} = T[]  # b_rotational (sealed pivot)
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
    p.W     = p.Wfac * p.h_eff
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

@inline function ctrans(u, p)      # v3.5: SMOOTH derivative gate over epsc
    m  = -u
    h0 = softpos(m, p.epsc) + p.h_eff
    sm = 0.5*(1 + m/sqrt(m*m + p.epsc*p.epsc))
    C  = (p.e*p.Tf/p.a)*log((h0 + p.a*p.Leff)/h0)
    dC = (p.e*p.Tf/p.a)*(1/h0 - 1/(h0 + p.a*p.Leff))*sm
    return C, dC
end
@inline function crot(u, p)        # F4 + v3.5 smooth gate (fixes the Fe jump)
    up = softpos(u, p.epsc)
    sp = 0.5*(1 + u/sqrt(u*u + p.epsc*p.epsc))
    au = p.a - up/p.Leff
    au < p.a_min && (au = p.a_min)
    C  = (p.e*p.Tf/au)*log((p.h_eff + au*p.Leff)/p.h_eff)
    dC = (p.e*p.Tf/p.Leff)*(log((p.h_eff + au*p.Leff)/p.h_eff)/au^2 -
         p.Leff/(au*(p.h_eff + au*p.Leff)))*sp
    return C, dC
end
@inline function cramp(u, p)       # paper Eq. 23 log ramp, continuity-repaired:
    # anchored to the true seam air value (kills the 16x value jump), softpos-
    # gated (C1), but keeping the paper's slope structure -- which is 58x the
    # approach slope at u = 0+. Inherits the paper's dimensional oddity in kk*u.
    up = softpos(u, p.epsc)
    sp = 0.5*(1 + u/sqrt(u*u + p.epsc*p.epsc))
    Cs  = (p.e*p.Tf/p.a)*log((p.h_eff + p.a*p.Leff)/p.h_eff)   # seam air value
    kk  = 2*p.Tp/(p.a*p.Leff)
    den = log1p(kk*p.a*p.Leff)
    C  = Cs + (p.Cmax_ramp - Cs)*log1p(kk*up)/den
    dC = (p.Cmax_ramp - Cs)*kk/((1 + kk*up)*den)*sp
    return C, dC
end
@inline capR(u, p) = p.cap_model === :ramp ? cramp(u, p) : crot(u, p)

@inline function cap_side(u, p)    # C1 sealing blend over |u| < W
    u <= -p.W && return ctrans(u, p)
    u >=  p.W && return capR(u, p)
    x  = (u + p.W)/(2*p.W)
    s  = x*x*x*(10 - 15*x + 6*x*x)          # quintic smootherstep: C2 at edges
    ds = 30*x*x*(1 - x)*(1 - x)/(2*p.W)
    Ct_, dCt_ = ctrans(u, p)
    Cr_, dCr_ = capR(u, p)
    return (1 - s)*Ct_ + s*Cr_, (1 - s)*dCt_ + s*dCr_ + ds*(Cr_ - Ct_)
end
@inline function bside(u, p)       # blended film; UDE discrepancy multiplies here
    u <= -p.W && return btrans(-u + p.h_eff, p)
    u >=  p.W && return p.crot_scale*brot(p.a - u/p.Leff, p)
    x = (u + p.W)/(2*p.W); s = x*x*x*(10 - 15*x + 6*x*x)
    bt = btrans(softpos(-u, p.epsc) + p.h_eff, p)
    br = p.crot_scale*brot(p.a - softpos(u, p.epsc)/p.Leff, p)
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
    spr = 0.5*(1 + ur/sqrt(ur*ur + p.epsw*p.epsw))   # v3.6: sigma'(d) factors make
    spl = 0.5*(1 + ul/sqrt(ul*ul + p.epsw*p.epsw))   # Fw = -dV/dx exactly, with V =
    gr = 1 + p.cw*x2dot                              # (kw/(pw+1))*softpos(d)^(pw+1);
    gl = 1 - p.cw*x2dot                              # also zeroes the preload identically
    Fw = 0.0
    gr > 0 && (Fw += -p.kw*dr^p.pw*spr*gr)
    gl > 0 && (Fw +=  p.kw*dl^p.pw*spl*gl)
    collision_state = (ur > 0 || ul > 0) ? "contact" : "translational"
    return Fc, Fw, collision_state
end

# Viscous damping, Fd (F5: blended slip-corrected films, both walls).
# Signature keeps the v2 shape; only x2, x2dot enter the film physics.
function damping(x1, x1dot, x2, x2dot, p)
    ur = x2 - p.gp
    ul = -x2 - p.gp
    # smooth open/close gates: closing stroke -> 1, opening -> vent_open
    sv = x2dot/sqrt(x2dot*x2dot + p.veps*p.veps)
    gr = p.vent_open + (1 - p.vent_open)*0.5*(1 + sv)   # right wall closes when x2dot > 0
    gl = p.vent_open + (1 - p.vent_open)*0.5*(1 - sv)   # left wall closes when x2dot < 0
    return -p.c*(gr*bside(ur, p) + gl*bside(ul, p))*x2dot
end

# Electrostatic coupling, Fe (F3: attractive; F4 capacitance; analytic dC).
# v3.7: third argument is Vout (the new state). With Vc = Vbias - Vout the
# capacitor voltage, Fe = (Vc^2/2)*dCt/(N/2) exactly -- no charge, no
# large-number subtraction anywhere. q = Vc*Ct is an observable.
function electrostatic(x1, x2, Vout, p)
    Cr, dCr = cap_side(x2 - p.gp, p)
    Cl, dCl = cap_side(-x2 - p.gp, p)
    Vr = 1/(2/p.crl + 1/Cr)
    Vl = 1/(2/p.crl + 1/Cl)
    Cvar = (p.N/2)*(Vr + Vl)
    dC   = (p.N/2)*((Vr/Cr)^2*dCr - (Vl/Cl)^2*dCl)
    Ctotal = Cvar + p.cp
    Vc = p.Vbias - Vout
    Fe = +((Vc*Vc/2)*dC)/(p.N/2)                 # attractive at constant charge
    return Ctotal, Fe, dC
end

# 5 states: x1, x1dot, x2, x2dot, Vout (v3.7 cancellation-free formulation).
# Derivation: Vout = Vbias - q/Ct and qdot = Vout/R give
#   dVout/dt = -Vout/(R*Ct) + ((Vbias - Vout)/Ct)*(dCt/dx2)*x2dot,
# mathematically identical to the charge formulation; error control now acts
# on the mV-scale output directly. q = (Vbias - Vout)*Ct is the observable.
function CoupledSystem!(dz, z, p, t, current_acceleration)
    z1, z2, z3, z4, z5 = z
    Fs = spring(z1, p.k1, p.k3, p.gss, p.kss) - p.c1*z2
    Fc, Fw, _ = collision(z1, z3, z4, p)
    Fd = damping(z1, z2, z3, z4, p)
    Ctotal, Fe, dC = electrostatic(z1, z3, z5, p)
    Fext = current_acceleration
    dz[1] = z2
    dz[2] = (Fs + (p.N/2)*Fc)/p.m1 - Fext
    dz[3] = z4
    dz[4] = (-Fc + Fd + Fe + Fw)/p.m2 - Fext
    dz[5] = -z5/(p.Rload*Ctotal) + ((p.Vbias - z5)/Ctotal)*dC*z4
    return nothing
end

# Initialize a default Params instance and calculate dependent parameters
p = Params{Float64}()
p = create_params(p; verbose = false)

end # module AnalyticalModel

import .AnalyticalModel

# --------------------------------------- External Force ------------------------------------

# Sine Wave External Force
f = 20.0  # Frequency (Hz)
alpha = 2.8  # Applied acceleration constant; contact threshold alpha* = 2.45 +- 0.05
g = 9.81  # Gravitational constant (m/s^2)
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

# Equilibrated start: q = Vbias*Ct  <=>  Vout = 0 exactly.
z0 = [x10, x10dot, x20, x20dot, 0.0]
tspan  = use_sine ? (0.0, 0.5) : (0.0, 600e-6)
abstol = [1e-13, 1e-10, 1e-13, 1e-10, 1e-9]    # state-scaled; z5 is Vout [V]
reltol = use_sine ? 1e-5 : 1e-6

# ---------------------------------- Solve Analytical Model ---------------------------------

function CoupledSystem_wrapper!(dz, z, p, t)
    AnalyticalModel.CoupledSystem!(dz, z, p, t, Fext_input(t))
end

eqn = ODEProblem(CoupledSystem_wrapper!, z0, tspan, p_new)

# Rodas5P; NO seam callbacks (see header). If AD complains about the film tables:
# sol = solve(eqn, Rodas5P(autodiff = false); abstol, reltol, maxiters = Int(1e9))
sol = solve(eqn, Rodas5P(); abstol = abstol, reltol = reltol, maxiters = Int(1e9))

println(">>> collision_model version v3.8 (cap_model A/B switch; Wfac window control; experiment harness) <<<")
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
AM = AnalyticalModel

# Sample states + observables + reconstructed forces on a uniform grid from the
# dense interpolant. Used by ALL plots, zooms, and metrics: uniform grids avoid
# the rendering artifacts of plotting at raw adaptive steps (sparse in cruise,
# ultra-dense in taps).
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
        ct, fe, _ = AM.electrostatic(z1, z3, z5, p)
        Ct[i] = ct; Fe[i] = fe
    end
    V   = U[5,:]                                  # Vout is the state (v3.7)
    Q   = [(p.Vbias - V[i])*Ct[i] for i in 1:n]   # charge as observable
    pen = [abs(U[3,i]) - p.gp for i in 1:n]
    ae  = [Fext(t) for t in tg]
    return (; t = tg, x1 = U[1,:], x1dot = U[2,:], x2 = U[3,:], x2dot = U[4,:],
              Q, V, Ct, pen, Fs, Fc, Fw, Fd, Fe, ae)
end

# Uniform 10 us overview grid for the state and force plots
Wover = sample_window(sol, p_new, sol.t[1], sol.t[end]; dt = 1e-5)
to = Wover.t

p3  = plot(to, Wover.x1,    xlabel = "Time (s)", ylabel = "x1 (m)",     title = "Shuttle Mass Displacement (x1)");    display(p3)
p4  = plot(to, Wover.x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)");    display(p4)
p5  = plot(to, Wover.x2,    xlabel = "Time (s)", ylabel = "x2 (m)",     title = "Mobile Electrode Displacement (x2)")
hline!(p5, [p_new.gp, -p_new.gp]; ls = :dash, lc = :gray, label = ""); display(p5)
p6  = plot(to, Wover.x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Velocity (x2dot)"); display(p6)
p7  = plot(to, Wover.Q,     xlabel = "Time (s)", ylabel = "Q (C)",      title = "Charge (observable, v3.7)");          display(p7)
p8  = plot(to, Wover.V,     xlabel = "Time (s)", ylabel = "Vout (V)",   title = "Output Voltage (STATE, cancellation-free)"); display(p8)

# v3 diagnostics: penetration and total capacitance
p8b = plot(to, Wover.pen .* 1e9, xlabel = "Time (s)", ylabel = "|x2|-gp (nm)",
           title = "Penetration (contact when > 0)")
hline!(p8b, [0.0]; ls = :dash, lc = :gray, label = ""); display(p8b)
p8c = plot(to, Wover.Ct .* 1e12, xlabel = "Time (s)", ylabel = "Ctotal (pF)",
           title = "Total Capacitance"); display(p8c)

p9   = plot(to, Wover.Fs, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension + Soft-stopper Spring Force"); display(p9)
p10  = plot(to, Wover.Fc, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Electrode Coupling Force (always-on spring)"); display(p10)
p10b = plot(to, Wover.Fw, xlabel = "Time (s)", ylabel = "Fw (N)", title = "Tip Contact (Hunt-Crossley wall) Force"); display(p10b)
p11  = plot(to, Wover.Fd, xlabel = "Time (s)", ylabel = "Fd (N)", title = "Viscous Film Force (slip-corrected, blended)"); display(p11)
p12  = plot(to, Wover.Fe, xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force (attractive)"); display(p12)
p13  = plot(to, Wover.ae, xlabel = "Time (s)", ylabel = "a_ext (m/s^2)", title = "Applied Base Acceleration"); display(p13)

# ============================== EVALUATION SUITE (v3.7) ==================================
# =========================================================================================
# (1) force-function characterization sweeps      (2) two-cycle steady-state zooms
# (3) composite multi-panels (states / forces)    (4) collision/discontinuity metrics
# All sections are read-only with respect to the solve above.

# =========================== (1) FUNCTION VISUALIZER (v3.7) ==============================
# Every constitutive function plotted individually over its governing range,
# with the collision window marked, plus branch-decomposition views showing
# exactly how each law bridges the gap.
function visualize_functions(p)
    W_ = p.W; hc = p.h_eff
    markband!(pl) = vspan!(pl, [-W_*1e9, W_*1e9]; alpha = 0.10, color = :orange, label = "")
    # V1 suspension + stopper
    x1g = range(-16e-6, 16e-6; length = 1601)
    v = plot(x1g .* 1e6, [AM.spring(x, p.k1, p.k3, p.gss, p.kss)*1e6 for x in x1g];
             legend = false, xlabel = "x1 (um)", ylabel = "Fs (uN)",
             title = "V1  Suspension + stopper")
    vline!(v, [p.gss*1e6, -p.gss*1e6]; ls = :dash, lc = :gray); display(v)
    # V2 coupling spring
    rg = range(-0.5e-6, 0.5e-6; length = 301)
    display(plot(rg .* 1e9, [-p.ke*r*1e6 for r in rg]; legend = false,
                 xlabel = "x1 - x2 (nm)", ylabel = "Fc (uN)",
                 title = "V2  Coupling spring, slope -ke"))
    # V3 capacitance, full travel
    xg = range(-p.gp + 1e-7, p.gp + 1e-7; length = 2500)
    display(plot(xg .* 1e6, [AM.electrostatic(0.0, x, 0.0, p)[1]*1e12 for x in xg];
                 legend = false, xlabel = "x2 (um)", ylabel = "C_t (pF)",
                 title = "V3  Total capacitance, full travel"))
    # V4 capacitance through the collision window
    ug = range(-2e-6, 1e-7; length = 2500)
    v = plot(ug .* 1e9, [AM.electrostatic(0.0, p.gp + u, 0.0, p)[1]*1e12 for u in ug];
             legend = false, xlabel = "u = x2 - gp (nm)", ylabel = "C_t (pF)",
             title = "V4  Capacitance through the collision window")
    markband!(v); vline!(v, [0.0]; ls = :dash, lc = :gray); display(v)
    # V5 capacitance slope (89x drop)
    ug2 = range(-5e-7, 1e-7; length = 3000)
    v = plot(ug2 .* 1e9, [abs(AM.electrostatic(0.0, p.gp + u, 0.0, p)[3]) for u in ug2];
             yscale = :log10, legend = false, xlabel = "u (nm)",
             ylabel = "|dC_t/dx2| (F/m)", title = "V5  Capacitance slope: 89x kinematic drop")
    markband!(v); display(v)
    # V6 electrostatic force, full travel (equilibrated: Vout = 0)
    display(plot(xg .* 1e6, [AM.electrostatic(0.0, x, 0.0, p)[2]*1e6 for x in xg];
                 legend = false, xlabel = "x2 (um)", ylabel = "Fe (uN)",
                 title = "V6  Electrostatic force, q equilibrated"))
    # V7 Fe through the window
    v = plot(ug2 .* 1e9, [AM.electrostatic(0.0, p.gp + u, 0.0, p)[2]*1e6 for u in ug2];
             legend = false, xlabel = "u (nm)", ylabel = "Fe (uN)",
             title = "V7  Fe through the collision window")
    markband!(v); display(v)
    # V8 film coefficient, full travel
    display(plot(xg .* 1e6, [AM.bside(x - p.gp, p) + AM.bside(-x - p.gp, p) for x in xg];
                 yscale = :log10, legend = false, xlabel = "x2 (um)",
                 ylabel = "B (N s/m)", title = "V8  Film coefficient, full travel"))
    # V9 film branches and blend across the window
    ub = range(-3*W_, 3*W_; length = 800)
    v = plot(ub .* 1e9, [AM.btrans(AM.softpos(-u, p.epsc) + hc, p) for u in ub];
             label = "translational branch", xlabel = "u (nm)", ylabel = "b (N s/m)",
             title = "V9  Film branches and blend")
    plot!(v, ub .* 1e9, [AM.brot(p.a - AM.softpos(u, p.epsc)/p.Leff, p) for u in ub];
          label = "rotational branch")
    plot!(v, ub .* 1e9, [AM.bside(u, p) for u in ub]; lw = 2.0, label = "blended")
    markband!(v); display(v)
    # V10 capacitance branches and blend (per pair side)
    v = plot(ub .* 1e9, [AM.ctrans(u, p)[1]*1e15 for u in ub];
             label = "translational branch", xlabel = "u (nm)",
             ylabel = "C_air (fF, per pair side)",
             title = "V10  Capacitance branches and blend")
    plot!(v, ub .* 1e9, [AM.capR(u, p)[1]*1e15 for u in ub];
          label = p.cap_model === :ramp ? "paper ramp branch" : "rotational branch")
    plot!(v, ub .* 1e9, [AM.cap_side(u, p)[1]*1e15 for u in ub]; lw = 2.0, label = "blended")
    markband!(v); display(v)
    # V11 wall force with velocity slices
    dg = range(0.0, 60e-9; length = 400)
    v = plot(dg .* 1e9, [AM.collision(0.0, p.gp + d, 0.0, p)[2]*1e6 for d in dg];
             label = "x2dot = 0", xlabel = "penetration (nm)", ylabel = "Fw (uN)",
             title = "V11  Wall force (Hunt-Crossley)")
    plot!(v, dg .* 1e9, [AM.collision(0.0, p.gp + d, -10e-3, p)[2]*1e6 for d in dg]; label = "-10 mm/s")
    plot!(v, dg .* 1e9, [AM.collision(0.0, p.gp + d, +10e-3, p)[2]*1e6 for d in dg]; label = "+10 mm/s")
    display(v)
    # V12 wall activation zoom (zero preload)
    dz2 = range(-5e-9, 5e-9; length = 600)
    display(plot(dz2 .* 1e9, [AM.collision(0.0, p.gp + d, 0.0, p)[2]*1e9 for d in dz2];
                 legend = false, xlabel = "penetration (nm)", ylabel = "Fw (nN)",
                 title = "V12  Wall activation zoom (softpos, zero preload)"))
    # V13 damping separability
    vg = range(-20e-3, 20e-3; length = 41)
    v = plot(; xlabel = "x2dot (mm/s)", ylabel = "Fd (uN)",
             title = "V13  Fd linear in velocity (separability)")
    for x2v in (0.0, p.gp - 1e-6, p.gp - 1e-7)
        plot!(v, vg .* 1e3, [AM.damping(0.0, 0.0, x2v, w, p)*1e6 for w in vg];
              label = string("gap = ", round((p.gp - x2v)*1e9; digits = 0), " nm"))
    end
    display(v)
    # V14 the regularizers themselves
    xr = range(-4.0, 4.0; length = 400)
    v = plot(xr, [AM.softpos(x, 1.0) for x in xr]; label = "softpos(x, eps)",
             xlabel = "x/eps (and shifted blend axis)", ylabel = "",
             title = "V14  Regularizers")
    plot!(v, xr, [max(x, 0.0) for x in xr]; ls = :dash, label = "(x)+")
    xs2 = range(0.0, 1.0; length = 200)
    plot!(v, 4 .* xs2 .- 4, [x^3*(10 - 15*x + 6*x^2)*4 for x in xs2];
          label = "smootherstep (scaled)")
    display(v)
    return nothing
end
visualize_functions(p_new)


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
    # nm-resolution localization of the crawl (accepted steps binned by u)
    us = [abs(0.5*(x2s[i] + x2s[i+1])) - p_new.gp for i in 1:length(dts)]
    edges = [-1.0, -200e-9, -50e-9, -25e-9, -5e-9, -1e-9, 1e-9, 5e-9, 25e-9, 1.0]
    println("dt localization by u = |x2|-gp (accepted steps):")
    for k in 1:(length(edges) - 1)
        m = [(us[i] > edges[k]) && (us[i] <= edges[k+1]) for i in 1:length(us)]
        n = count(m); n == 0 && continue
        dd = sort(dts[m])
        @printf("  u in (%9.1f, %9.1f] nm : n = %8d  min = %.1e  med = %.1e\n",
                edges[k]*1e9, edges[k+1]*1e9, n, dd[1], dd[div(end, 2) + 1])
    end
    idx = findall(i -> abs(us[i]) < 6e-8, 1:length(us))
    idx = idx[1:max(1, div(length(idx), 20000)):end]
    m8 = scatter([us[i]*1e9 for i in idx], [log10(dts[i]) for i in idx];
                 ms = 1, msw = 0, legend = false, xlabel = "u (nm)",
                 ylabel = "log10 dt (s)",
                 title = "Accepted dt vs distance to seam (crawl localization)")
    display(m8)

    # --- energy ledger (global correctness of the blended discontinuity) ---
    function ledger(Wd, p; Fext = Fext_input)
        n = length(Wd.t)
        E = zeros(n); Pn = zeros(n)
        for i in 1:n
            x1_, v1_, x2_, v2_, q_ = Wd.x1[i], Wd.x1dot[i], Wd.x2[i], Wd.x2dot[i], Wd.Q[i]
            d = AM.softpos(abs(x2_) - p.gp, p.epsw)   # v3.6: matches wall potential
            KE = 0.5*p.m1*v1_^2 + (p.N/2)*0.5*p.m2*v2_^2
            Vm = 0.5*p.k1*x1_^2 + (0.25/4)*p.k3*x1_^4 +
                 (abs(x1_) > p.gss ? 0.5*p.kss*(abs(x1_) - p.gss)^2 : 0.0) +
                 (p.N/2)*0.5*p.ke*(x1_ - x2_)^2 + (p.N/2)*(p.kw/(p.pw + 1))*d^(p.pw + 1)
            E[i] = KE + Vm + q_^2/(2*Wd.Ct[i])
            qd = (p.Vbias - q_/Wd.Ct[i])/p.Rload
            sv_ = v2_/sqrt(v2_*v2_ + p.veps*p.veps)
            Bv = (p.vent_open + (1 - p.vent_open)*0.5*(1 + sv_))*AM.bside(x2_ - p.gp, p) +
                 (p.vent_open + (1 - p.vent_open)*0.5*(1 - sv_))*AM.bside(-x2_ - p.gp, p)
            s_ = sign(x2_)
            spd  = 0.5*(1 + (abs(x2_) - p.gp)/sqrt((abs(x2_) - p.gp)^2 + p.epsw^2))
            gate = (1 + p.cw*v2_*s_ > 0)
            Pw = gate ? (p.N/2)*p.kw*d^p.pw*spd*p.cw*v2_^2 : 0.0
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


# =========================== STEP-FLOOR PROBE (opt-in, v3.4) ===========================
# Six short A/B runs over (0, 216 ms) -- cruise + first tap -- each seconds-to-a-
# minute. Interpretation: floor gone under Rosenbrock23/FBDF but not noAD =>
# high-order estimator vs d^(3/2) temporal regularity at d ~ 0; gone under
# "quadratic wall" => the 3/2 power law itself (legitimately replaceable: line
# contact is not Hertz-3/2 anyway); gone only with kw = 0 => wall pathway confirmed
# but activation-independent; gone under noAD => AD after all, mechanism revised.
run_probe = false
if run_probe
    function probe_run(tag; alg = Rodas5P(), pmod! = identity)
        pp = deepcopy(AnalyticalModel.p); pmod!(pp)
        zz = [0.0, 0.0, 0.0, 0.0, 0.0]        # v3.7: z5 = Vout, equilibrated
        pr = ODEProblem(CoupledSystem_wrapper!, zz, (0.0, 0.2160), pp)
        ss = solve(pr, alg; abstol = abstol, reltol = reltol, maxiters = Int(1e9))
        dts_ = diff(ss.t); x2_ = [u[3] for u in ss.u]
        ns_ = [abs(abs(0.5*(x2_[i] + x2_[i+1])) - pp.gp) < 5*pp.W for i in 1:length(dts_)]
        dn_ = dts_[ns_]
        @printf("%-24s acc=%9d rej=%9d  near-seam min/med dt = %.1e / %.1e  maxpen = %.1f nm\n",
                tag, ss.stats.naccept, ss.stats.nreject,
                isempty(dn_) ? NaN : minimum(dn_),
                isempty(dn_) ? NaN : sort(dn_)[div(end, 2) + 1],
                (maximum(abs.(x2_)) - pp.gp)*1e9)
        return nothing
    end
    probe_run("baseline Rodas5P")
    probe_run("Rosenbrock23 (order 2/3)"; alg = Rosenbrock23())
    probe_run("FBDF (variable order)";    alg = FBDF())
    probe_run("Rodas5P autodiff=false";   alg = Rodas5P(autodiff = false))
    probe_run("wall off (kw = 0)";        pmod! = pp -> (pp.kw = 0.0; nothing))
    probe_run("quadratic wall (pw = 2)";  pmod! = pp -> (pp.kw = 1.2e10; pp.pw = 2.0; nothing))
end

# =========================================================================================
# ===================== STANDALONE EXPERIMENT HARNESS (opt-in, v3.8) ======================
# =========================================================================================
# Automates the three-pass protocol:
#   Pass 1  invariance: fresh default build must reproduce the main run above
#   Pass 2  capacitance A/B: cap_model = :ramp latch test (pre-registered
#           prediction: Fe_contact >> ke*Rc  =>  latch-on-first-contact)
#   Pass 3  window study: Wfac, h_eff, kw moved ONE at a time
# Each experiment is a fresh solve with rebuilt derived parameters (film tables
# depend on h_eff; W depends on Wfac). Uses the top-level forcing: set
# use_sine = true. Cost: 8 solves of t_exp seconds each (t_exp = 0.35 is the
# quick pass; 0.5 gives more settled steady-state cycles).
run_experiments = false
t_exp = 0.5

if run_experiments && !use_sine
    @warn "Experiment harness expects sine forcing; set use_sine = true."
end
if run_experiments && use_sine

# ---- compact per-run metrics (the five numbers we discuss, plus signatures) ----
function ledger_exp(Wd, p)
    n = length(Wd.t); E = zeros(n); Pn = zeros(n)
    for i in 1:n
        x1_, v1_, x2_, v2_, q_ = Wd.x1[i], Wd.x1dot[i], Wd.x2[i], Wd.x2dot[i], Wd.Q[i]
        d  = AM.softpos(abs(x2_) - p.gp, p.epsw)
        KE = 0.5*p.m1*v1_^2 + (p.N/2)*0.5*p.m2*v2_^2
        Vm = 0.5*p.k1*x1_^2 + (0.25/4)*p.k3*x1_^4 +
             (abs(x1_) > p.gss ? 0.5*p.kss*(abs(x1_) - p.gss)^2 : 0.0) +
             (p.N/2)*0.5*p.ke*(x1_ - x2_)^2 + (p.N/2)*(p.kw/(p.pw + 1))*d^(p.pw + 1)
        E[i] = KE + Vm + q_^2/(2*Wd.Ct[i])
        qd = (p.Vbias - q_/Wd.Ct[i])/p.Rload
        Bv = AM.bside(x2_ - p.gp, p) + AM.bside(-x2_ - p.gp, p)
        s_ = sign(x2_)
        spd  = 0.5*(1 + (abs(x2_) - p.gp)/sqrt((abs(x2_) - p.gp)^2 + p.epsw^2))
        gate = (1 + p.cw*v2_*s_ > 0)
        Pw = gate ? (p.N/2)*p.kw*d^p.pw*spd*p.cw*v2_^2 : 0.0
        Pn[i] = -Fext_input(Wd.t[i])*(p.m1*v1_ + (p.N/2)*p.m2*v2_) +
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

function experiment_metrics(sol_, p; label = "")
    Wf  = sample_window(sol_, p, sol_.t[1], sol_.t[end]; dt = 2e-6)
    inc = Wf.pen .> 0
    stats = sol_.stats
    rej = stats.nreject/max(stats.naccept, 1)
    half_period_ms = 0.5e3/f
    if !any(inc)
        return (; label, engaged = false, onset = NaN, duty = 0.0,
                  hold = 0.0, maxpen = -maximum(Wf.pen)*1e9, e_eff = NaN,
                  f_c = NaN, Vpk = maximum(abs.(Wf.V))*1e3, Vdc = NaN,
                  ledg = NaN, rej, latch = false)
    end
    onset = Wf.t[findfirst(inc)]
    duty  = 100*count(inc)/length(inc)
    maxpen = maximum(Wf.pen)*1e9
    # longest contiguous hold (latch discriminator; tap holds << half period)
    hold = 0; run_ = 0
    for b in inc
        run_ = b ? run_ + 1 : 0
        hold = max(hold, run_)
    end
    hold_ms = hold*2e-6*1e3
    latch = hold_ms > 0.25*half_period_ms
    # alias-free first-episode window: e_eff, contact-conditioned ring
    tw0 = max(onset - 2e-5, sol_.t[1]); tw1 = min(onset + 6e-3, sol_.t[end])
    Wep_ = sample_window(sol_, p, tw0, tw1; dt = 1e-7)
    penz = Wep_.pen
    KEe = 0.5*p.m2 .* Wep_.x2dot.^2
    iin = findfirst(penz .> 0); iout = findlast(penz .> 0)
    e_eff = (iin === nothing || iout === nothing || iout >= length(penz)) ? NaN :
            sqrt(max(KEe[min(iout+1, end)], 0)/max(KEe[max(iin-1, 1)], 1e-300))
    mc = penz .> 0; tin = count(mc)*1e-7
    f_c = tin > 2e-6 ? count(!=(0), diff(sign.(Wep_.x2dot[mc])))/(2*tin)/1e3 : NaN
    # output signatures (post-ramp), incl. the latch DC shift over last 2 cycles
    mpr = Wf.t .> t_ramp
    Vpk = maximum(abs.(Wf.V[mpr]))*1e3
    m2c = Wf.t .> (Wf.t[end] - 2/f)
    Vdc = (sum(Wf.V[m2c])/count(m2c))*1e3
    ledg = ledger_exp(Wep_, p)
    return (; label, engaged = true, onset = onset*1e3, duty, hold = hold_ms,
              maxpen, e_eff, f_c, Vpk, Vdc, ledg, rej, latch)
end

function print_exp_header()
    @printf("%-22s | eng | onset  | duty  | hold  | maxpen |  e_eff | f_c    | Vpk   | Vdc    | ledger  | rej   | latch\n",
            "experiment")
    @printf("%-22s |     | (ms)   | (%%)   | (ms)  | (nm)   |        | (kHz)  | (mV)  | (mV)   |         |       |\n", "")
    println(repeat("-", 128))
end
function print_exp_row(m)
    @printf("%-22s |  %s  | %6.2f | %5.2f | %5.2f | %6.1f | %6.3f | %6.1f | %5.1f | %6.2f | %.1e | %5.3f |  %s\n",
            m.label, m.engaged ? "Y" : "n", m.onset, m.duty, m.hold, m.maxpen,
            m.e_eff, m.f_c, m.Vpk, m.Vdc, m.ledg, m.rej, m.latch ? "YES" : "no")
end

# ---- driver: fresh params, single-knob modification, rebuild, solve, measure ----
function run_experiment(label; mod! = (pp -> nothing), keepsol = false, alg = Rodas5P())
    pp = Params{Float64}()
    mod!(pp)
    AnalyticalModel.create_params(pp; verbose = false)   # rebuilds W and film tables
    z0e = [0.0, 0.0, 0.0, 0.0, 0.0]
    pre = ODEProblem(CoupledSystem_wrapper!, z0e, (0.0, t_exp), pp)
    se  = solve(pre, alg; abstol = abstol, reltol = reltol, maxiters = Int(1e9))
    m   = experiment_metrics(se, pp; label)
    print_exp_row(m)
    return keepsol ? (m, se, pp) : (m, nothing, pp)
end

println("\n================== EXPERIMENT HARNESS (alpha = ", alpha, ", t_exp = ", t_exp, " s) ==================")
print_exp_header()

# Pass 1: invariance -- the main run above, measured with the same instrument,
# then a fresh default build. The two rows must agree to plotting precision.
m_main = experiment_metrics(sol, p_new; label = "main run (above)")
print_exp_row(m_main)
mB, sB, pB = run_experiment("baseline (fresh build)"; keepsol = true)

# Pass 2: capacitance A/B with the pre-registered latch prediction
mR, sR, pR = run_experiment("cap_model = :ramp";
                            mod! = pp -> (pp.cap_model = :ramp), keepsol = true)
FeC = AM.electrostatic(0.0, pR.gp + 5e-9, 0.0, pR)[2]
Rc  = (pR.k1*pR.gp + pR.m1*alpha*g)/(pR.k1 + (pR.N/2)*pR.ke)
@printf("\n[A/B verdict] ramp Fe(u=+5nm) = %.2f uN vs release capacity ke*Rc = %.2f uN\n",
        FeC*1e6, pR.ke*Rc*1e6)
println(mR.latch ?
    "[A/B verdict] LATCH observed (hold $(round(mR.hold; digits=2)) ms) -- prediction CONFIRMED; ramp branch dynamically self-refutes." :
    "[A/B verdict] NO latch (hold $(round(mR.hold; digits=2)) ms) -- prediction REFUTED; report this table back for investigation.")

# Pass 3: window study, one knob at a time (v3 capacitance throughout)
println()
mW3, _, _  = run_experiment("Wfac = 0.3";    mod! = pp -> (pp.Wfac = 0.3))
mW2, sW2, pW2 = run_experiment("Wfac = 0.2"; mod! = pp -> (pp.Wfac = 0.2), keepsol = true)
mH8, _, _  = run_experiment("h_eff = 80 nm"; mod! = pp -> (pp.h_eff = 80e-9))
mH10, _, _ = run_experiment("h_eff = 100 nm";mod! = pp -> (pp.h_eff = 100e-9))
mK3, _, _  = run_experiment("kw = 3e6";      mod! = pp -> (pp.kw = 3e6))
mK10, sK10, pK10 = run_experiment("kw = 1e7";mod! = pp -> (pp.kw = 1e7), keepsol = true)

println("\nReading guide: crisper stick-unstick = hold DOWN, e_eff UP, f_c UP at rej <~ 0.1;")
println("h_eff rows move the adhesion well (duty/hold down, Vpk down is the trade);")
println("latch discriminator: hold > ", round(0.25*0.5e3/f; digits = 1), " ms. Ledger must stay ~1e-8..1e-7.")

# ---- comparison figures ----
Td = 1/f
WzB = sample_window(sB, pB, sB.t[end] - 2*Td, sB.t[end]; dt = 2e-6)
WzR = sample_window(sR, pR, sR.t[end] - 2*Td, sR.t[end]; dt = 2e-6)
e1 = plot(WzB.t .* 1e3, WzB.pen .* 1e9; label = "baseline :v3",
          xlabel = "t (ms)", ylabel = "|x2|-gp (nm)",
          title = "E1  A/B penetration, last 2 cycles (latch test)")
plot!(e1, WzR.t .* 1e3, WzR.pen .* 1e9; label = "cap_model = :ramp")
hline!(e1, [0.0]; ls = :dash, lc = :gray, label = ""); display(e1)
e2 = plot(WzB.t .* 1e3, WzB.V .* 1e3; label = "baseline :v3",
          xlabel = "t (ms)", ylabel = "Vout (mV)",
          title = "E2  A/B output voltage, last 2 cycles")
plot!(e2, WzR.t .* 1e3, WzR.V .* 1e3; label = "cap_model = :ramp"); display(e2)
# episode-aligned crispness comparison
function first_episode(sol_, p, m)
    isnan(m.onset) && return nothing
    t0 = m.onset*1e-3
    return sample_window(sol_, p, max(t0 - 2e-5, 0.0), min(t0 + 3e-3, sol_.t[end]); dt = 1e-7)
end
EB = first_episode(sB, pB, mB); EW = first_episode(sW2, pW2, mW2); EK = first_episode(sK10, pK10, mK10)
e3 = plot(; xlabel = "t - onset (us)", ylabel = "|x2|-gp (nm)",
          title = "E3  First episode, aligned: window crispness")
EB !== nothing && plot!(e3, (EB.t .- mB.onset*1e-3) .* 1e6, EB.pen .* 1e9; label = "baseline")
EW !== nothing && plot!(e3, (EW.t .- mW2.onset*1e-3) .* 1e6, EW.pen .* 1e9; label = "Wfac = 0.2")
EK !== nothing && plot!(e3, (EK.t .- mK10.onset*1e-3) .* 1e6, EK.pen .* 1e9; label = "kw = 1e7")
hline!(e3, [0.0]; ls = :dash, lc = :gray, label = ""); display(e3)

end # run_experiments


# =========================== MODEL EXPERIMENTS (opt-in, v3.8.1) ==========================
# Harness rev v3.10: live-parameter echo per run (catches stale-module runs:
# three consecutive bit-identical tables under 'changed' c1 motivated this);
# Lambda prints n/a when no complete taps; film battery F1-F3 tests the
# press hypothesis (in-contact damping up, opening suction down).
# Prior rev v3.9: impact ratio Lambda = KE_in/V_wall(delta_max) classifies
# each tap PRESS vs BALLISTIC (restitution reported only when ballistic; a
# pressed contact has no meaningful e_eff); per-experiment 4-panel figures
# saved as exp_<tag>.png; run_voltage_sweep battery maps the pull-in boundary
# V_PI for both capacitance models -- the model discriminator.
# Prior rev v3.8.2: per-tap segmentation -- e_eff is now the MEDIAN PER-TAP
# restitution over complete contacts (the previous episode-level value measured
# e^n across n taps), f_ring/tap counts x2dot crossings in the first 20 us of
# each tap (the damped entry chirp; ~tau 15 us at baseline), and n_tap / t_tap
# report tap count and median single-contact duration.
# Prior rev v3.8.1: episode window sized to the first full episode; e_eff only
# when a genuine separation exists; regime column TAP / HOLD / LATCH (HOLD =
# longest contact > 0.15*T_drive, i.e. drive-enforced press-and-hold); ring
# split into f_hold (hold-phase grind) and f_entry (first 100 us electrode
# transient -- the frequency-up channel); P in watts (scientific) + Vrms.
# Automated three-pass protocol from the capacitance / collision-window review:
#   E0      invariance baseline (:v3, Wfac 0.5, h_eff 50 nm, kw 1e6) -- must
#           reproduce the main run above (same params, forcing, tolerances)
#   E1      cap_model = :ramp   paper Eq. 23 A/B; PREDICTION: latch-on-contact
#   E2-E3   Wfac 0.3, 0.2       tighter blend window (crisper stick-unstick?)
#   E4-E5   h_eff 80, 100 nm    shallower adhesion well ("air gap potential less")
#   E6-E7   kw 3e6, 1e7         stiffer wall (shorter taps, higher ring freq)
# Each run: FRESH Params + create_params (film tables rebuilt; struct-default
# edits propagate, runtime p_new edits do NOT), same forcing/solver/tolerances.
# LATCH criterion: longest continuous contact interval > half a drive period.
# Cost: 8 solves at main-run size; expect several minutes total.
run_experiments = true
make_plots = true          # per-experiment 4-panel figure, saved to exp_<tag>.png
if run_experiments

    function experiment_metrics(sol, p)
        Wc  = sample_window(sol, p, sol.t[1], sol.t[end]; dt = 2e-6)
        pen = Wc.pen
        inc = pen .> 0
        Td  = 1/f
        duty = count(inc)/length(inc)
        lh = 0; cur = 0
        for b in inc
            b ? (cur += 1; lh = max(lh, cur)) : (cur = 0)
        end
        hold  = lh*2e-6                     # longest continuous contact (s)
        latch = hold > 0.5*Td
        pkV   = maximum(abs.(Wc.V))
        mL    = Wc.t .>= (Wc.t[end] - 0.2)  # harvest over the last 0.2 s
        Pavg  = sum(Wc.V[mL].^2)/count(mL)/p.Rload
        Vrms  = sqrt(sum(Wc.V[mL].^2)/count(mL))
        regime = latch ? "LATCH" : (hold > 0.15*Td ? "HOLD" : "TAP")
        if !any(inc)
            return (; engaged = false, duty, hold, latch, onset = NaN,
                      maxpen = maximum(pen), vimp = NaN, eeff = NaN,
                      fring = NaN, lam = NaN, ballistic = false, ntap = 0,
                      ttap = NaN, occW = NaN, pkV, Pavg, Vrms, regime = "FREE",
                      ledg = NaN)
        end
        onset = Wc.t[findfirst(inc)]
        tw0 = max(onset - 2e-5, sol.t[1])
        dur = min(max(3*hold, 2e-3), 2.5e-2)   # cover the first FULL episode
        tw1 = min(onset + dur, sol.t[end])
        We  = sample_window(sol, p, tw0, tw1; dt = 1e-7)   # alias-free window
        pz  = We.pen
        vimp = Float64[]
        for i in 2:length(pz)
            (pz[i] > 0 && pz[i-1] <= 0) && push!(vimp, abs(We.x2dot[i]))
        end
        occW = count(abs.(abs.(We.x2) .- p.gp) .< p.W)/length(pz)
        mc = pz .> 0
        tin = count(mc)*1e-7
        fring = tin > 2e-6 ?
            count(!=(0), diff(sign.(We.x2dot[mc])))/(2*tin) : NaN
        KEe  = 0.5*p.m2 .* We.x2dot.^2
        # ---- per-tap segmentation (complete contacts only; v3.8.2) ----
        starts = Int[]; stops = Int[]
        for i in 2:length(mc)
            (mc[i] && !mc[i-1]) && push!(starts, i)
            (!mc[i] && mc[i-1]) && push!(stops, i - 1)
        end
        (!isempty(stops) && !isempty(starts) && stops[1] < starts[1]) && popfirst!(stops)
        etaps = Float64[]; ttaps = Float64[]; lams = Float64[]
        for k in 1:min(length(starts), length(stops))
            i1, i2 = starts[k], stops[k]
            (i1 > 1 && i2 + 1 <= length(pz)) || continue
            push!(etaps, sqrt(max(KEe[i2 + 1], 0)/max(KEe[i1 - 1], 1e-300)))
            push!(ttaps, (i2 - i1 + 1)*1e-7)
            # impact ratio: entry KE vs wall potential at this tap's peak
            # penetration. Lambda << 1 = quasi-static PRESS (drive pushes the
            # electrode in and out; no strike, no ring, restitution undefined);
            # Lambda ~ 1 = BALLISTIC impact (ring physics active).
            dk = maximum(pz[i1:i2])
            Vwk = (p.kw/(p.pw + 1))*dk^(p.pw + 1)
            push!(lams, KEe[i1 - 1]/max(Vwk, 1e-300))
        end
        _med(v) = isempty(v) ? NaN : sort(v)[div(length(v), 2) + 1]
        ntap = length(etaps)
        ttap = _med(ttaps)           # median single-tap contact time
        lam  = _med(lams)            # median impact ratio Lambda
        ballistic = !isnan(lam) && lam > 0.1
        eeff = ballistic ? _med(etaps) : NaN   # restitution only if ballistic
        # energy ledger on the fine window (global acceptance test, v3.6 form)
        n = length(We.t); E = zeros(n); Pn = zeros(n)
        for i in 1:n
            x1_, v1_, x2_, v2_, q_ = We.x1[i], We.x1dot[i], We.x2[i], We.x2dot[i], We.Q[i]
            d  = AM.softpos(abs(x2_) - p.gp, p.epsw)
            KE = 0.5*p.m1*v1_^2 + (p.N/2)*0.5*p.m2*v2_^2
            Vm = 0.5*p.k1*x1_^2 + (0.25/4)*p.k3*x1_^4 +
                 (abs(x1_) > p.gss ? 0.5*p.kss*(abs(x1_) - p.gss)^2 : 0.0) +
                 (p.N/2)*0.5*p.ke*(x1_ - x2_)^2 + (p.N/2)*(p.kw/(p.pw + 1))*d^(p.pw + 1)
            E[i] = KE + Vm + q_^2/(2*We.Ct[i])
            qd = (p.Vbias - q_/We.Ct[i])/p.Rload
            sv_ = v2_/sqrt(v2_*v2_ + p.veps*p.veps)
            Bv = (p.vent_open + (1 - p.vent_open)*0.5*(1 + sv_))*AM.bside(x2_ - p.gp, p) +
                 (p.vent_open + (1 - p.vent_open)*0.5*(1 - sv_))*AM.bside(-x2_ - p.gp, p)
            s_ = sign(x2_)
            spd  = 0.5*(1 + (abs(x2_) - p.gp)/sqrt((abs(x2_) - p.gp)^2 + p.epsw^2))
            gate = (1 + p.cw*v2_*s_ > 0)
            Pw = gate ? (p.N/2)*p.kw*d^p.pw*spd*p.cw*v2_^2 : 0.0
            Pn[i] = -Fext_input(We.t[i])*(p.m1*v1_ + (p.N/2)*p.m2*v2_) +
                    p.Vbias*qd - p.Rload*qd^2 - (p.N/2)*p.c*Bv*v2_^2 - p.c1*v1_^2 - Pw
        end
        work = 0.0; resid = 0.0
        for i in 2:n
            work += 0.5*(Pn[i] + Pn[i-1])*(We.t[i] - We.t[i-1])
            resid = max(resid, abs((E[i] - E[1]) - work))
        end
        ledg = resid/(maximum(abs.(E)) + maximum(abs, Pn)*(We.t[end] - We.t[1]))
        return (; engaged = true, duty, hold, latch, onset,
                  maxpen = maximum(pen),
                  vimp = isempty(vimp) ? NaN : maximum(vimp),
                  eeff, ntap, ttap, fring, lam, ballistic, occW, pkV, Pavg, Vrms,
                  regime, ledg)
    end

    results = NamedTuple[]
    function run_experiment(tag; kwargs...)
        pp = AnalyticalModel.Params{Float64}(; kwargs...)
        AnalyticalModel.create_params(pp; verbose = false)
        prob = ODEProblem(CoupledSystem_wrapper!, [0.0, 0.0, 0.0, 0.0, 0.0],
                          (0.0, tspan[2]), pp)
        tel = @elapsed ss = solve(prob, Rodas5P(); abstol = abstol,
                                  reltol = reltol, maxiters = Int(1e9))
        dts_ = diff(ss.t)
        x2_  = [u[3] for u in ss.u]
        ns_  = [abs(abs(0.5*(x2_[i] + x2_[i+1])) - pp.gp) < 5*pp.W for i in 1:length(dts_)]
        dmed = any(ns_) ? sort(dts_[ns_])[div(count(ns_), 2) + 1] : NaN
        M  = experiment_metrics(ss, pp)
        rr = ss.stats.nreject/max(ss.stats.naccept, 1)
        @printf("%-14s %-10s acc %8d  rej/acc %.3f  seam-med-dt %.1e  wall %.0f s\n",
                tag, string(ss.retcode), ss.stats.naccept, rr, dmed, tel)
        @printf("    params: c1 %.1e  h_eff %.0f nm  Wfac %.2f  kw %.1e  cw %.0f  Vb %.2f V  cap %s  crot_scale %.1f  vent %.2f\n",
                pp.c1, pp.h_eff*1e9, pp.Wfac, pp.kw, pp.cw, pp.Vbias,
                string(pp.cap_model), pp.crot_scale, pp.vent_open)
        if M.engaged
            @printf("    onset %.1f ms  duty %5.2f %%  hold %6.2f ms  [%s]  maxpen %5.1f nm  occW %4.1f %%\n",
                    M.onset*1e3, 100*M.duty, M.hold*1e3, M.regime,
                    M.maxpen*1e9, 100*M.occW)
            @printf("    vimp %5.2f mm/s  tap: n=%d  t_tap %5.1f us  Lambda %.1e [%s]  e_eff %s  f_hold %5.1f kHz\n",
                    1e3*M.vimp, M.ntap, 1e6*M.ttap, M.lam,
                    M.ntap == 0 ? "n/a" : (M.ballistic ? "BALLISTIC" : "PRESS"),
                    M.ballistic ? @sprintf("%5.3f", M.eeff) : "  n/a",
                    M.fring/1e3)
            @printf("    pkV %7.2f mV  Vrms %6.3f mV  P %.2e W  ledger %.1e\n",
                    1e3*M.pkV, 1e3*M.Vrms, M.Pavg, M.ledg)
        else
            @printf("    NO CONTACT (closest %.1f nm)  pkV %7.2f mV  Vrms %6.3f mV  P %.2e W\n",
                    -M.maxpen*1e9, 1e3*M.pkV, 1e3*M.Vrms, M.Pavg)
        end
        push!(results, (; tag, M.duty, M.maxpen, M.lam, M.ballistic, M.pkV, M.Vrms, M.regime))
        if make_plots && M.engaged
            Wp = sample_window(ss, pp, ss.t[1], ss.t[end]; dt = 1e-5)
            Wt = sample_window(ss, pp, M.onset - 2e-5, M.onset + 4e-4; dt = 1e-7)
            Wl = sample_window(ss, pp, ss.t[end] - 2/f, ss.t[end]; dt = 2e-6)
            pa = plot(Wp.t, Wp.pen .* 1e9; legend = false, xlabel = "t (s)",
                      ylabel = "|x2|-gp (nm)", title = "penetration, full span")
            hline!(pa, [0.0]; ls = :dash, lc = :gray)
            pb = plot((Wt.t .- M.onset) .* 1e6, Wt.pen .* 1e9; legend = false,
                      xlabel = "t - onset (us)", ylabel = "pen (nm)",
                      title = "first taps (alias-free)")
            hline!(pb, [0.0]; ls = :dash, lc = :gray)
            hline!(pb, [pp.W*1e9, -pp.W*1e9]; ls = :dot, lc = :orange)
            pc = plot((Wt.t .- M.onset) .* 1e6, Wt.x2dot .* 1e3; legend = false,
                      xlabel = "t - onset (us)", ylabel = "x2dot (mm/s)",
                      title = "entry velocity (press vs ballistic at a glance)")
            pd = plot(Wl.t .* 1e3, Wl.V .* 1e3; legend = false,
                      xlabel = "t (ms)", ylabel = "Vout (mV)",
                      title = "output, last 2 drive cycles")
            fig = plot(pa, pb, pc, pd; layout = (2, 2), size = (1100, 800),
                       plot_title = tag, titlefontsize = 9, guidefontsize = 8,
                       tickfontsize = 7)
            display(fig)
            savefig(fig, "exp_" * replace(tag, " " => "_") * ".png")
        end
        return nothing
    end

    println("\n================ MODEL EXPERIMENTS (alpha = $(alpha), t_end = $(tspan[2]) s) ================")
    run_experiment("E0 baseline")
    run_experiment("E1 :ramp";     cap_model = :ramp)
    run_experiment("E2 Wfac 0.3";  Wfac = 0.3)
    run_experiment("E3 Wfac 0.2";  Wfac = 0.2)
    run_experiment("E4 heff 80n";  h_eff = 80e-9)
    run_experiment("E5 heff 100n"; h_eff = 100e-9)
    run_experiment("E6 kw 3e6";    kw = 3e6)
    run_experiment("E7 kw 1e7";    kw = 1e7)

    # -------- PULL-IN / RELEASE VOLTAGE SWEEP (the real discriminator) --------
    # Release condition: Fe(h_eff; Vb) > ke*Rc  =>  V_PI = Vb*sqrt(ke*Rc/Fe(Vb)).
    # Pre-registered predictions at alpha = 2.8 (Rc ~ 57 nm):
    #   :v3   rotation branch  V_PI ~ 5.7 V  (taps at 3-5 V, latches by ~6 V)
    #   :ramp paper branch     V_PI ~ 0.8 V  (taps only below ~0.8 V)
    # A bench measurement of the tap->latch bias boundary calibrates h_eff.
    run_voltage_sweep = true
    if run_voltage_sweep
        println("\n================ PULL-IN VOLTAGE SWEEP (both cap models) ================")
        for cm in (:v3, :ramp), Vb in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0)
            run_experiment(string(cm, " V=", Vb); cap_model = cm, Vbias = Vb)
        end
        println("\n---- pull-in regime grid ----")
        for r in results
            (startswith(r.tag, "v3") || startswith(r.tag, "ramp")) &&
                @printf("%-12s %7s\n", r.tag, r.regime)
        end
    end

    # -------- FILM-HYPOTHESIS TESTS (ring-in vs ring-out) --------
    # Panels to watch: entry-velocity wiggles (ring-in) and post-episode tails.
    # Pre-registered predictions:
    #   F1 crot x10 : entry ring GONE (clean touch); t_tap, duty ~unchanged;
    #                 pkV unchanged (output is approach-generated)
    #   F2 vent 0.3 : post-detach ring-out appears at ~f2 = 64.7 kHz in x2dot
    #                 and as a high-frequency tail on the Vout spikes
    #   F3 both     : the full press-hypothesis phenomenology -- silent capture,
    #                 oscillatory release. Bench Vout spectrum discriminates:
    #                 ~155 kHz content = capture ring real; ~65 kHz = ring-out.
    run_film_tests = true
    if run_film_tests
        println("\n================ FILM-HYPOTHESIS TESTS ================")
        run_experiment("F1 crot x10"; crot_scale = 10.0)
        run_experiment("F2 vent 0.3"; vent_open = 0.3)
        run_experiment("F3 both";     crot_scale = 10.0, vent_open = 0.3)
    end

    println("\n---- comparison (decision metrics) ----")
    @printf("%-14s %7s %8s %10s %10s %8s %9s %7s\n",
            "tag", "duty%", "maxpen", "Lambda", "tap-type", "pkV", "Vrms", "regime")
    for r in results
        @printf("%-14s %7.2f %7.1fn %10.1e %10s %7.2fm %8.3fm %7s\n",
                r.tag, 100*r.duty, r.maxpen*1e9, r.lam,
                isnan(r.lam) ? "n/a" : (r.ballistic ? "BALLISTIC" : "press"),
                1e3*r.pkV, 1e3*r.Vrms, r.regime)
    end
end