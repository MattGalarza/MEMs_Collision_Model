# ===========================================================================================
# collision_model_v5 : the CORRECTED two-coordinate model (collision_model_corrected.jl)
# written in the syntax, style and outline of collision_model_v4.jl.
#
# Mathematics : collision_model_corrected.jl, unchanged (common displacement field, consistent
#               mass matrix, local-stack capacitance, projected Reynolds film, Hunt-Crossley
#               tip wall, Vout state), plus ONE addition:
#                 * thickness-vented squeeze film (chapter App. D, Eq. D.2), vent_faces = 1 | 2
# Defaults    : the closures behind comparison panels (d) and (e):
#                 vent_faces = 2, ce = 0.75e-4 N s/m per beam, c1 = 5.29e-5 N s/m (Q = 50
#                 placeholder -- measure the shuttle ring-down before trusting it)
#               alpha = 4.95 -> panel (d);  alpha = 2.7 -> panel (e)
# Regression  : vent_faces = 0, c1 = 0, ce = 0, panels = 512, alpha = 4.95 is the submitted
#               corrected model; its last-four-cycle load power is 9.8853e-13 W.
# Expected    : (Python twin, these defaults, 10 cycles, ledger on; agree to ~1e-3)
#               4.95 g: residual/throughput 5e-8, ER = 5.248e-12 J, last-four-cycle power
#                       1.307e-11 W, Vout -22.7/+32.1 mV, 12 contacts per half cycle with
#                       flights 1479, 978, 706, 525, 394 ... us, release at 1.54 g
#               2.7 g : residual/throughput 4e-9, ER = 2.696e-12 J, power 8.00e-12 W,
#                       Vout -15.8/+21.0 mV, 19 contacts, flights 2050, 1495, 1165, 922 ... us
# Status      : written WITHOUT a Julia runtime. This text was syntax-checked only; the
#               algorithm as written here was executed through a line-by-line Python mirror
#               (all selfcheck items pass; right-hand side equals the validated twin to 2e-13
#               with the lengthwise film and 2e-5 with the vented film; vented film vs the
#               classical rectangular plate 4e-5, vs a Bessel closed form on a wedge 1e-6..3e-4).
#               Run AnalyticalModel.selfcheck() FIRST (run_selfcheck below).
# ===========================================================================================

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
        palette=:okabe_ito             # Colour-blind-safe palette (v4 used :default)
    )
end
set_journal_theme()

# --------------------------------------- Analytical Model ----------------------------------

module AnalyticalModel
using DifferentialEquations
using Parameters
using LinearAlgebra
export Params, p, create_params, spring, collision, damping, electrostatic, CoupledSystem!,
       energy, ledger, selfcheck

@with_kw mutable struct Params{T<:Real}
    # Fundamental geometric parameters
    g0::T = 14e-6        # Initial (bare) tip gap
    Tp::T = 120e-9       # Parylene-C thickness, per face
    Tf::T = 25e-6        # Electrode (device-layer) thickness
    wt::T = 9e-6         # Electrode width at the clamped root (narrow)
    wb::T = 30e-6        # Electrode width at the free tip (wide)
    ws::T = 14.7e-6      # Suspension spring width
    wss::T = 14e-6       # Soft-stopper width
    Leff::T = 400e-6     # Effective (overlap) electrode length
    Lf::T = 450e-6       # Full electrode length
    Lsp::T = 1400e-6     # Suspension spring length
    Lss::T = 1000e-6     # Soft-stopper length
    gss::T = 14e-6       # Soft-stopper position
    gap_slope::T = NaN   # Facing-wall gap slope. NaN -> (wb - wt)/Lf (mobile and fixed
                         # faces both inclined). Override with metrology.

    # Array / suspension topology
    N::Int = 160         # Number of gap branches (N/2 mobile electrodes, two faces each)
    nsp_par::Int = 4     # Suspension: parallel chains
    nsp_ser::Int = 6     # Suspension: series spans per chain
    nss::Int = 2         # Soft-stopper cantilevers acting in parallel, per side
    gamma3::T = 1.0      # Cubic-stiffness geometry correction (not calibrated)

    # Mass and material properties
    m1::T = 2.0933e-6    # Shuttle-only mass. MUST exclude the explicit mobile beams:
                         # their distributed inertia enters through the mass matrix M
    rho::T = 2330.0      # Density of silicon
    E::T = 170e9         # Young's modulus
    e::T = 8.85e-12      # Permittivity of free space
    ep::T = 3.2          # Relative permittivity of Parylene-C
    eta::T = 1.849e-5    # Viscosity of air
    lambda::T = 70e-9    # Mean free path of air molecules (m)
    sigmap::T = 1.016    # Slip coefficient for rarefaction

    # Dissipation closures -------------------------------------------------------
    # These three lines are what separates panels (d)/(e) from the submitted model.
    # c1 : shuttle damping [N s/m]. 5.29e-5 = mtot*w0/Q with Q = 50 (PLACEHOLDER).
    #      It matters: the flight-time ratio of the chatter is 0.73-0.79 with it and
    #      0.76-0.85 without. Replace by a measured ring-down before fitting ce.
    # ce : damping of RELATIVE tip/root velocity, per beam [N s/m]. The only loss that
    #      acts while the tips rest on the wall and the shuttle bounces on the N/2
    #      beam springs (contact mode ~4.3 kHz). zeta_c = (N/2)*ce/(2*M11*w_c) = 0.052
    #      here; chosen so the simulated flight-time ratio matches the bench (~0.73).
    #      Stands in for an interface loss: it also damps the free 55 kHz beam mode
    #      (zeta ~ 0.57), which no measurement here constrains.
    # vent_faces : squeeze-film drainage. 0 = lengthwise only (submitted model, gas
    #      leaves through the two ends of the 400 um overlap); 1 = one device-layer
    #      face open as well; 2 = both faces open (top surface + cavity underneath).
    #      Lengthwise-only overdamps by ~210x at rest, ~40x at a 1 um gap, 6.5x at
    #      contact relative to 2 (57x / 14x / 3.7x relative to 1).
    c::T = 1.0           # Film scale -- physical default
    c1::T = 5.29e-5      # Shuttle damping [N s/m]   (submitted model: 0)
    ce::T = 0.75e-4      # Beam relative damping, per beam [N s/m]   (submitted model: 0)
    vent_faces::Int = 2  # Film drainage: 0 lengthwise | 1 one face open | 2 both faces open
    vent_modes::Int = 4  # Thickness modes solved exactly; higher modes use the strip limit
                         # (1, 2, 4, 8 modes -> 2.213, 2.083, 2.068, 2.066 mN s/m at contact)

    # Hard stop (disabled until the as-fabricated gap is measured) ----------------
    ghs::T = Inf         # Hard-stop engagement position [m]; must satisfy ghs >= gss
    khs::T = 1.0e9       # Hard-stop stiffness [N/m^phs]
    phs::T = 1.5         # Hard-stop exponent

    # Contact / boundary parameters -----------------------------------------------
    # h_eff : residual pressed-contact air gap; floor of every local gap. Sets the
    #         electrostatic force at contact (~1/(h_eff + hd)) and the film floor.
    # epsg  : smoothing of the positive part in the gap law h = h_eff + softpos(d).
    # epsw  : wall activation width. epss : stopper activation width.
    # ls    : tip-vent sealing half-width; chi goes 0 -> 1 over -ls < delta < ls.
    # kw,pw : wall stiffness / exponent, per beam.  cw : Hunt-Crossley loss.
    #         cw acts on the TIP velocity, which is nearly zero in the contact mode,
    #         so it does NOT set the shuttle-level restitution; ce does.
    # Tested insensitive (wall x100, epsw/10, epsg/4, ls/5): none of these widths
    # creates or removes the chatter; they are numerics, not physics.
    h_eff::T = 50e-9     # Residual contact air gap
    epsg::T = 2e-9       # Gap positive-part smoothing
    epsw::T = 0.5e-9     # Wall engagement smoothing
    epss::T = 1e-9       # Soft/hard stopper engagement smoothing
    ls::T = 25e-9        # Tip-vent sealing half-width
    seal::Bool = true    # Robin tip closure: vent seals as the tip lands (false: always open)
    kw::T = 1e6          # Hunt-Crossley wall stiffness, per beam (N/m^1.5)
    pw::T = 1.5          # Hunt-Crossley exponent
    cw::T = 50.0         # Hunt-Crossley dissipation (s/m)

    # Electrical parameters
    cp::T = 5e-12        # Parasitic capacitance (parallel to the variable capacitor)
    Vbias::T = 3.0       # Bias voltage
    Rload::T = 0.42e6    # Load resistance

    # Numerical resolution
    panels::Int = 256    # Graded Simpson panels along the overlap (2*panels + 1 nodes).
                         # vs 1024 panels the vented film differs by 1.3e-3 / 3.4e-4 / 7e-5
                         # at 128 / 256 / 512; dC by < 1e-6. Panels (d),(e) were run at 128.

    # Derived parameters - calculated by create_params()
    nb::Int = 0              # Number of mobile electrodes, N/2
    a::T = 0.0               # Gap slope actually used
    gc::T = 0.0              # Tip travel to nominal contact, g0 - 2*Tp - h_eff
    hd::T = 0.0              # Effective dielectric thickness, 2*Tp/ep
    kp::T = 0.0              # Slip conductance length 6*sigmap*lambda
    ke::T = 0.0              # Electrode tip stiffness, per beam
    k1::T = 0.0              # Linear spring constant
    k3::T = 0.0              # Cubic spring constant
    kss::T = 0.0             # Soft-stopper spring constant
    mtot::T = 0.0            # Total moving mass, m1 + (N/2)*m_beam
    M::Matrix{T} = zeros(T, 2, 2)      # Consistent mass matrix in (x1, x2)
    Minv::Matrix{T} = zeros(T, 2, 2)   # Its inverse
    beta::Vector{T} = zeros(T, 2)      # Base-excitation vector, M*[1, 1]
    # Quadrature grid along the overlap (built once by create_params)
    y::Vector{T} = T[]       # Graded nodes, y = 0 at the mobile tip
    wq::Vector{T} = T[]      # Simpson weights
    vol::Vector{T} = T[]     # Control-volume widths (vented film)
    B1::Vector{T} = T[]      # Weight of x1 in the face displacement, 1 - phi
    B2::Vector{T} = T[]      # Weight of x2 in the face displacement, phi
end

# ------------------------------ constitutive helpers ------------------------------
# Stable C-infinity positive part, its exact derivative, and the C2 quintic step
@inline softpos(z, e_)  = z >= 0 ? (z + hypot(z, e_))/2 : e_^2/(2*(hypot(z, e_) - z))
@inline dsoftpos(z, e_) = z >= 0 ? (1 + z/hypot(z, e_))/2 :
                                   e_^2/(2*hypot(z, e_)*(hypot(z, e_) - z))
@inline smootherstep(z) = z <= 0 ? 0.0 : z >= 1 ? 1.0 : z^3*(10 - 15*z + 6*z^2)

# 96-point Gauss-Legendre rule (Golub-Welsch) for the beam integrals
function gausslegendre(n)
    F = eigen(SymTridiagonal(zeros(n), [j/sqrt(4*j*j - 1) for j in 1:n-1]))
    return F.values, 2 .* F.vectors[1, :].^2
end
const GLX, GLW = gausslegendre(96)
glquad(f, lo, hi) = (hi - lo)/2*sum(GLW[j]*f((lo + hi)/2 + (hi - lo)/2*GLX[j]) for j in eachindex(GLX))

# Tapered section and the unit-tip-load shape phi(s): phi(0) = 0, phi(Lf) = 1
bwidth(s, p) = p.wt + (p.wb - p.wt)*s/p.Lf
bEI(s, p)    = p.E*p.Tf*bwidth(s, p)^3/12
bshape(s, ke, p) = s == 0 ? 0.0 : ke*glquad(z -> (s - z)*(p.Lf - z)/bEI(z, p), 0.0, s)

function create_params(p::Params{T}; verbose = true) where T<:Real
    @assert p.N > 0 && iseven(p.N) && p.nsp_par > 0 && p.nsp_ser > 0 && p.nss > 0
    @assert 0 < p.Leff <= p.Lf && p.g0 > 2*p.Tp + p.h_eff && p.panels >= 16
    @assert all(>(0), (p.Tf, p.wt, p.wb, p.ws, p.wss, p.Lsp, p.Lss, p.m1, p.rho, p.E, p.e,
                       p.ep, p.eta, p.h_eff, p.epsg, p.epsw, p.epss, p.ls, p.Rload, p.kw,
                       p.pw, p.phs, p.gss))
    @assert all(>=(0), (p.c, p.c1, p.ce, p.cw, p.gamma3, p.cp, p.Tp, p.lambda, p.sigmap, p.khs))
    @assert p.ghs >= p.gss && p.vent_faces in (0, 1, 2) && p.vent_modes >= 1

    p.nb = div(p.N, 2)
    p.a  = isnan(p.gap_slope) ? (p.wb - p.wt)/p.Lf : p.gap_slope
    @assert isfinite(p.a) && p.a >= 0

    # Electrode tip stiffness: unit tip load on the tapered clamped beam (Castigliano)
    p.ke = 1/glquad(s -> (p.Lf - s)^2/bEI(s, p), 0.0, p.Lf)

    # Suspension / stopper spring constants
    p.k1  = p.nsp_par/p.nsp_ser*p.E*p.Tf*p.ws^3/p.Lsp^3
    p.k3  = p.gamma3*p.nsp_par/p.nsp_ser^3*0.72*p.E*p.Tf*p.ws/p.Lsp^3
    p.kss = p.nss*p.E*p.Tf*p.wss^3/(4*p.Lss^3)

    # Consistent mass matrix from the common displacement field w = (1 - phi)*x1 + phi*x2
    mass(i, j) = p.nb*glquad(s -> begin
        phi = bshape(s, p.ke, p); B = (1 - phi, phi)
        p.rho*p.Tf*bwidth(s, p)*B[i]*B[j]
    end, 0.0, p.Lf)
    m11 = p.m1 + mass(1, 1); m12 = mass(1, 2); m22 = mass(2, 2)
    p.M = [m11 m12; m12 m22]
    @assert isposdef(Symmetric(p.M))
    p.Minv = inv(p.M)
    p.beta = p.M*ones(2)
    p.mtot = sum(p.M)

    # Electrical / contact derived
    p.gc = p.g0 - 2*p.Tp - p.h_eff
    p.hd = 2*p.Tp/p.ep
    p.kp = 6*p.sigmap*p.lambda

    # Graded panel endpoints (uniform in log(1 + y/lg)); Simpson midpoints are
    # arithmetic in physical y. lg = h_eff/a is the length over which the wedge opens by h_eff.
    lg = p.h_eff/max(p.a, 0.001)
    zz = range(0.0, log1p(p.Leff/lg); length = p.panels + 1)
    edges = lg .* expm1.(zz)
    edges[end] = p.Leff
    p.y = zeros(2*p.panels + 1); p.wq = zero(p.y)
    for j in 1:p.panels
        i = 2*j - 1; ya, yb = edges[j], edges[j+1]; d = yb - ya
        p.y[i] = ya; p.y[i+1] = (ya + yb)/2; p.y[i+2] = yb
        p.wq[i] += d/6; p.wq[i+1] += 2*d/3; p.wq[i+2] += d/6
    end
    n = length(p.y)
    p.vol = zeros(n)
    p.vol[1] = (p.y[2] - p.y[1])/2; p.vol[n] = (p.y[n] - p.y[n-1])/2
    for i in 2:n-1
        p.vol[i] = (p.y[i+1] - p.y[i-1])/2
    end
    p.B2 = [bshape(p.Lf - v, p.ke, p) for v in p.y]
    p.B1 = 1 .- p.B2

    if verbose
        nbke = p.nb*p.ke
        k11 = p.k1 + nbke; k12 = -nbke
        K  = [k11 k12; k12 nbke]
        fr = sqrt.(eigvals(Symmetric(K), Symmetric(p.M)))./(2*pi)
        d0 = film(0.0, 0.0, p); dc = film(p.gc, p.gc, p)
        Cc, Fe1c, Fe2c, _, _ = electrostatic(p.gc, p.gc, 0.0, p)
        Fsc = p.k1*p.gc + p.k3*p.gc^3
        println("\n--- Springs ---")
        println("ke = ", p.ke, "   k1 = ", p.k1, "   k3 = ", p.k3, "   kss = ", p.kss)
        println("(expect ke = 22.5969569824, k1 = 3.2799375, k3 = 3.0357e8, kss = 5.831)")
        println("\n--- Mass matrix (x1, x2) ---")
        println("M11 = ", p.M[1,1], "   M12 = ", p.M[1,2], "   M22 = ", p.M[2,2], "   mtot = ", p.mtot)
        println("\n--- Modes ---")
        println("shuttle f1 = ", round(fr[1]; digits = 1), " Hz   free beam f2 = ",
                round(fr[2]/1e3; digits = 2), " kHz   tips-pinned contact mode = ",
                round(sqrt((p.k1 + nbke)/p.M[1,1])/(2*pi); digits = 0), " Hz")
        println("(expect 197.3 Hz, 55.43 kHz, 4667 Hz)")
        println("zeta_c of ce in the contact mode = ",
                round(p.nb*p.ce/(2*sqrt((p.k1 + nbke)*p.M[1,1])); sigdigits = 3),
                "   shuttle Q of c1 = ", (p.c1 > 0 ? round(sqrt(p.k1*p.mtot)/p.c1; sigdigits = 3) : Inf))
        println("\n--- Film, rigid translation, whole array (vent_faces = ", p.vent_faces, ") ---")
        println("b(rest)    = ", d0[1] + 2*d0[2] + d0[3], "  (expect 1.957e-6 | 7.32e-6 | 4.187e-4 for 2 | 1 | 0)")
        println("b(contact) = ", dc[1] + 2*dc[2] + dc[3], "  (expect 2.068e-3 | 3.60e-3 | 1.347e-2 for 2 | 1 | 0)")
        println("\n--- Electrostatics at nominal contact, Vbias = ", p.Vbias, " V ---")
        println("Ct = ", Cc*1e12, " pF   Fe1 = ", Fe1c*1e6, " uN   Fe2 = ", Fe2c*1e6, " uN")
        println("(expect 7.0974 pF, 0.377 uN, 13.05 uN at 3 V)")
        println("suspension force at gc = ", Fsc*1e6, " uN   static hold voltage = ",
                round(p.Vbias*sqrt(Fsc/(Fe1c + Fe2c)); digits = 2), " V  (expect 45.75 uN, 5.54 V)")
    end
    return p
end

# Local gap and its coordinate sensitivities on wall r = +-1:
#   d = gc + a*y - r*(B1*x1 + B2*x2),  h = h_eff + softpos(d),  h_i = dh/dx_i
function gapfield(x1, x2, r, p)
    d  = p.gc .+ p.a .* p.y .- r .* (p.B1 .* x1 .+ p.B2 .* x2)
    h  = p.h_eff .+ softpos.(d, p.epsg)
    dh = dsoftpos.(d, p.epsg)
    return h, -r .* dh .* p.B1, -r .* dh .* p.B2
end

# Half-panel Simpson primitive H(y) = int_0^y f on the graded grid
function cumulative_simpson(f, p)
    H = zeros(length(f))
    for j in 1:p.panels
        i = 2*j - 1; d = p.y[i+2] - p.y[i]
        H[i+1] = H[i] + d*(5*f[i] + 8*f[i+1] - f[i+2])/24
        H[i+2] = H[i] + d*(f[i] + 4*f[i+1] + f[i+2])/6
    end
    return H
end

# Lengthwise-only film (submitted model; chapter Eq. 3.88), one wall, all beams.
# (G p')' = 12 eta hdot, p(Leff) = 0, Robin tip vent; centred-moment form is PSD.
function film_lengthwise(h, h1, h2, chi, p)
    H1 = cumulative_simpson(h1, p); H2 = cumulative_simpson(h2, p)
    W  = p.wq ./ (h.^2 .* (h .+ p.kp)); I0 = sum(W)
    mu1 = dot(W, H1)/I0; mu2 = dot(W, H2)/I0
    Z1 = H1 .- mu1; Z2 = H2 .- mu2
    fac = 12*p.eta*p.Tf*p.nb*p.c
    return fac*(dot(W, Z1.^2)   + chi*I0*mu1^2),
           fac*(dot(W, Z1.*Z2)  + chi*I0*mu1*mu2),
           fac*(dot(W, Z2.^2)   + chi*I0*mu2^2)
end

# Thickness-vented film (chapter App. D, Eq. D.2), one wall, all beams.
#   d/dy(G dp/dy) + G d2p/dzeta2 = 12 eta hdot,   p = 0 on the open device-layer faces.
# Cosine modes in zeta, k_n = (2n+1)pi/W. Mode n solves (G P')' - G k_n^2 P = h_j with
# P(Leff) = 0 and the same Robin tip vent as the lengthwise film. Modes >= vent_modes use
# their strip limit P = -h_j/(G k_n^2). Conservative finite volumes on the graded nodes;
# solve and projection share the control-volume weights, so the result is a sum of
# B'(-A_n)^{-1}B terms: symmetric positive semidefinite for ANY grid. One open face is
# the mirrored strip: width 2*Tf, half the force.
function film_vented(h, h1, h2, chi, p)
    y = p.y; n = length(y)
    Wv   = p.vent_faces == 2 ? p.Tf : 2*p.Tf
    frac = p.vent_faces == 2 ? 1.0 : 0.5
    G    = h.^2 .* (h .+ p.kp)
    I0   = dot(p.wq, 1 ./ G)
    gface = [2*G[i]*G[i+1]/(G[i] + G[i+1])/(y[i+1] - y[i]) for i in 1:n-1]  # face conductances
    open_tip = chi <= 1e-12
    kappa    = open_tip ? 0.0 : (1 - chi)/(chi*I0)
    H = hcat(h1, h2)
    B = H .* p.vol
    B[n, :] .= 0.0                        # p = 0 at the open far end
    open_tip && (B[1, :] .= 0.0)          # p = 0 at a fully open tip
    V = zeros(2, 2); csum = 0.0
    for k in 0:p.vent_modes-1
        q = 2*k + 1; kn = q*pi/Wv; w = Wv*(16/(q*pi)^2)/2; csum += 1/q^4
        dg = -(G .* p.vol) .* kn^2        # leak to the faces
        dg[1:n-1] .-= gface               # right-face conductance of node i
        dg[2:n]   .-= gface               # left-face conductance of node i
        dl = copy(gface); du = copy(gface)
        if open_tip
            dg[1] = 1.0; du[1] = 0.0
        else
            dg[1] -= kappa                # Robin vent: G p' = kappa p at the tip
        end
        dg[n] = 1.0; dl[n-1] = 0.0
        P = Tridiagonal(dl, dg, du) \ B
        V .-= w .* (B' * P)
    end
    tail = (pi^4/96 - csum)*8/pi^4*Wv^3   # strip limit of the unsolved modes
    S = (H .* p.vol)' * (H ./ G)
    fac = 12*p.eta*p.nb*p.c*frac
    return fac*(V[1,1] + tail*S[1,1]),
           fac*((V[1,2] + V[2,1])/2 + tail*S[1,2]),
           fac*(V[2,2] + tail*S[2,2])
end

# Generalized film matrix D(x) = [d11 d12; d12 d22], both walls, all beams
function film(x1, x2, p)
    d11 = 0.0; d12 = 0.0; d22 = 0.0
    for r in (-1.0, 1.0)
        h, h1, h2 = gapfield(x1, x2, r, p)
        chi = p.seal ? smootherstep((r*x2 - p.gc + p.ls)/(2*p.ls)) : 0.0
        a11, a12, a22 = p.vent_faces == 0 ? film_lengthwise(h, h1, h2, chi, p) :
                                            film_vented(h, h1, h2, chi, p)
        d11 += a11; d12 += a12; d22 += a22
    end
    return d11, d12, d22
end

# --------------------------------- force functions ---------------------------------
# ALL forces below are generalized forces of the WHOLE array (N/2 beams, both walls) on
# the coordinates x1 (shuttle) and x2 (actual tip displacement, same frame as x1).

# Suspension spring force, Fsp (+ soft stopper, + optional hard stop); acts on x1.
# Exact gradients of their potentials, so the ledger stays an exact acceptance test.
function spring(x1, p)
    Fsp = -p.k1*x1 - p.k3*x1^3
    Fss = 0.0; Fhs = 0.0
    for r in (-1.0, 1.0)
        zs = r*x1 - p.gss
        Fss -= r*p.kss*softpos(zs, p.epss)*dsoftpos(zs, p.epss)
        if isfinite(p.ghs)
            zh = r*x1 - p.ghs
            Fhs -= r*p.khs*softpos(zh, p.epss)^p.phs*dsoftpos(zh, p.epss)
        end
    end
    return Fsp + Fss + Fhs
end

# Electrode coupling + wall, Fc / Fw. Fc (beam bending, always on) acts on x1 and -Fc on
# x2; Fw (Hunt-Crossley tip wall with tensile cutoff) acts on x2. Pw >= 0 is the exact
# wall loss including the energy removed on the zero-force unloading branch.
function collision(x1, x2, x2dot, p)
    Fc = -p.nb*p.ke*(x1 - x2)
    Fw = 0.0; Pw = 0.0
    for r in (-1.0, 1.0)
        zw = r*x2 - p.gc; vr = r*x2dot
        s  = softpos(zw, p.epsw)
        A  = p.kw*s^p.pw*dsoftpos(zw, p.epsw)
        gate = max(1 + p.cw*vr, 0.0)
        Fw -= p.nb*r*A*gate
        Pw += p.nb*A*vr*(gate - 1)
    end
    collision_state = abs(x2) > p.gc ? "contact" : "translational"
    return Fc, Fw, collision_state, Pw
end

# Viscous film damping, Fd: projected Reynolds film on BOTH coordinates,
#   [Fd1, Fd2] = -D(x)*[x1dot, x2dot].  Also returns D = (d11, d12, d22) for the ledger.
function damping(x1, x1dot, x2, x2dot, p)
    d11, d12, d22 = film(x1, x2, p)
    Fd1 = -(d11*x1dot + d12*x2dot)
    Fd2 = -(d12*x1dot + d22*x2dot)
    return Fd1, Fd2, (d11, d12, d22)
end

# Electrostatic coupling, Fe: local dielectric stack (two coatings + air in series at each
# strip, strips / faces / beams in parallel), Fe_i = (Vc^2/2)*dCt/dx_i with Vc = Vbias - Vout.
# x1 changes the interior gaps through B1, so Fe1 is nonzero even at a fixed tip.
function electrostatic(x1, x2, Vout, p)
    Ctotal = p.cp; dC1 = 0.0; dC2 = 0.0
    for r in (-1.0, 1.0)
        h, h1, h2 = gapfield(x1, x2, r, p)
        f = h .+ p.hd
        Ctotal += p.nb*p.e*p.Tf*dot(p.wq, 1 ./ f)
        dC1    -= p.nb*p.e*p.Tf*dot(p.wq, h1 ./ f.^2)
        dC2    -= p.nb*p.e*p.Tf*dot(p.wq, h2 ./ f.^2)
    end
    Vc  = p.Vbias - Vout
    Fe1 = 0.5*Vc*Vc*dC1
    Fe2 = 0.5*Vc*Vc*dC2
    return Ctotal, Fe1, Fe2, dC1, dC2
end

# 5 states: x1, x1dot, x2, x2dot, Vout.
#   M*[x1ddot, x2ddot] = [f1, f2],   M = consistent mass matrix, base input -beta*a(t)
#   dVout/dt = -Vout/(R*Ct) + ((Vbias - Vout)/Ct)*(dC1*x1dot + dC2*x2dot)
# Optional states 6-12 (when length(z) >= 12) are work integrals, not physics:
#   Wbase, Wbias, ER, Dfilm, Dstruct, Dwall, throughput.
function CoupledSystem!(dz, z, p, t, current_acceleration)
    z1, z2, z3, z4, z5 = z
    Fs = spring(z1, p) - p.c1*z2
    Fc, Fw, _, Pw = collision(z1, z3, z4, p)
    Fb = -p.nb*p.ce*(z2 - z4)                 # beam relative damping: on x1, and -Fb on x2
    Fd1, Fd2, D = damping(z1, z2, z3, z4, p)
    Ctotal, Fe1, Fe2, dC1, dC2 = electrostatic(z1, z3, z5, p)
    Fext = current_acceleration
    f1 = Fs + Fc + Fb + Fd1 + Fe1 - p.beta[1]*Fext
    f2 =     -Fc - Fb + Fd2 + Fe2 + Fw - p.beta[2]*Fext
    dz[1] = z2
    dz[2] = p.Minv[1,1]*f1 + p.Minv[1,2]*f2
    dz[3] = z4
    dz[4] = p.Minv[2,1]*f1 + p.Minv[2,2]*f2
    dz[5] = -z5/(p.Rload*Ctotal) + ((p.Vbias - z5)/Ctotal)*(dC1*z2 + dC2*z4)
    if length(dz) >= 12
        Pb = -Fext*(p.beta[1]*z2 + p.beta[2]*z4)
        Pe = p.Vbias*z5/p.Rload
        PR = z5*z5/p.Rload
        Pf = D[1]*z2*z2 + 2*D[2]*z2*z4 + D[3]*z4*z4
        Ps = p.c1*z2*z2 + p.nb*p.ce*(z4 - z2)^2
        dz[6] = Pb; dz[7] = Pe; dz[8] = PR; dz[9] = Pf; dz[10] = Ps; dz[11] = Pw
        dz[12] = abs(Pb) + abs(Pe) + PR + Pf + Ps + Pw
    end
    return nothing
end

# Stored energy T + U0 + Uw + Ue in base-relative coordinates, and the exact residual
#   ledger = E - E0 - Wbase - Wbias + ER + Dfilm + Dstruct + Dwall  (= 0 analytically)
function energy(z, p)
    x1, v1, x2, v2, Vout = z[1], z[2], z[3], z[4], z[5]
    Tk = 0.5*(p.M[1,1]*v1*v1 + 2*p.M[1,2]*v1*v2 + p.M[2,2]*v2*v2)
    U  = p.k1*x1^2/2 + p.k3*x1^4/4 + p.nb*p.ke*(x2 - x1)^2/2
    for r in (-1.0, 1.0)
        U += p.kss*softpos(r*x1 - p.gss, p.epss)^2/2
        isfinite(p.ghs) && (U += p.khs*softpos(r*x1 - p.ghs, p.epss)^(p.phs + 1)/(p.phs + 1))
        U += p.nb*p.kw*softpos(r*x2 - p.gc, p.epsw)^(p.pw + 1)/(p.pw + 1)
    end
    Ctotal, _, _, _, _ = electrostatic(x1, x2, Vout, p)
    return Tk + U + 0.5*Ctotal*(p.Vbias - Vout)^2
end
ledger(z, E0, p) = energy(z, p) - E0 - z[6] - z[7] + z[8] + z[9] + z[10] + z[11]

# Acceptance tests against the corrected-model regression constants and the Python twin.
# Builds its own 512-panel parameter sets, so it is independent of the run settings.
function selfcheck(; verbose = true)
    allok = true
    function chk(name, val, ref, rtol)
        pass = isapprox(val, ref; rtol = rtol)
        verbose && println(rpad(name, 52), pass ? "PASS  " : "FAIL  ", val, "   (ref ", ref, ")")
        allok = allok && pass
        return pass
    end
    function chkmax(name, val, tol)
        pass = val <= tol
        verbose && println(rpad(name, 52), pass ? "PASS  " : "FAIL  ", val, "   (tol ", tol, ")")
        allok = allok && pass
        return pass
    end
    sumD(d) = d[1] + 2*d[2] + d[3]
    verbose && println("\n--- selfcheck ---")
    q = create_params(Params{Float64}(panels = 512); verbose = false)
    mphys = q.rho*q.Tf*q.Lf*(q.wb + q.wt)/2
    chk("ke [N/m]", q.ke, 22.5969569824, 1e-9)
    chk("k1 [N/m]", q.k1, 3.2799375, 1e-12)
    chk("k3 [N/m^3]", q.k3, 3.035714285714e8, 1e-11)
    chk("kss [N/m]", q.kss, 5.831, 1e-12)
    chk("gc [m]", q.gc, 13.71e-6, 1e-12)
    chk("tip mass fraction M22/(nb*m_beam)", q.M[2,2]/q.nb/mphys, 0.369820966146, 1e-9)
    chk("sum(M) = m1 + nb*m_beam", sum(q.M), q.m1 + q.nb*mphys, 1e-12)
    chk("shape at the tip", bshape(q.Lf, q.ke, q), 1.0, 1e-12)
    chk("sum of Simpson weights = Leff", sum(q.wq), q.Leff, 1e-13)
    # capacitance and its gradient at nominal contact
    Cc, _, _, g1, g2 = electrostatic(q.gc, q.gc, 0.0, q)
    chk("Ct at contact [F]", Cc, 7.0974456842e-12, 1e-8)
    chk("dCt/dx1 at contact [F/m]", g1, 8.3856419104e-08, 1e-7)
    chk("dCt/dx2 at contact [F/m]", g2, 2.9006020264e-06, 1e-7)
    # vented film, both faces / one face, rigid translation (x1 = x2)
    xs = (0.0, q.gc - 0.95e-6, q.gc)
    for (x, ref) in zip(xs, (1.9570803057e-06, 1.1783055953e-04, 2.0677254530e-03))
        chk("vented film (2 faces) sum(D), x = $(x)", sumD(film(x, x, q)), ref, 1e-6)
    end
    q1 = create_params(Params{Float64}(panels = 512, vent_faces = 1); verbose = false)
    for (x, ref) in zip(xs, (7.3237497784e-06, 3.3537785951e-04, 3.5991030680e-03))
        chk("vented film (1 face) sum(D), x = $(x)", sumD(film(x, x, q1)), ref, 1e-6)
    end
    # lengthwise film = submitted model
    q0 = create_params(Params{Float64}(panels = 512, vent_faces = 0); verbose = false)
    for (x, ref) in zip(xs, (4.1873376307e-04, 4.7557756074e-03, 1.3468200914e-02))
        chk("lengthwise film sum(D), x = $(x)", sumD(film(x, x, q0)), ref, 1e-8)
    end
    # classical rectangular plate: uniform gap, no slip, all four edges open
    qr = create_params(Params{Float64}(panels = 512, gap_slope = 0.0, sigmap = 0.0, seal = false);
                       verbose = false)
    hr = qr.gc + qr.h_eff
    bt = 1 - 192/pi^5*(qr.Tf/qr.Leff)*sum(tanh(k*pi*qr.Leff/(2*qr.Tf))/k^5 for k in 1:2:199)
    chk("vented film vs rectangular-plate solution", sumD(film(0.0, 0.0, qr)),
        qr.N*qr.eta*qr.Leff*qr.Tf^3/hr^3*bt, 2e-4)
    # passivity, monotonicity, gradient and point-wise power balance at five states
    states = [(0.0, 0.0), (q.gc - 1e-6, q.gc - 1e-6), (q.gc + 0.1e-6, q.gc - 2e-9),
              (q.gc + 0.2e-6, q.gc + 5e-9), (-q.gc - 0.1e-6, -q.gc - 2e-9)]
    worstpsd = 0.0; worstmono = 0.0; worstgrad = 0.0; worstpow = 0.0
    for (x1, x2) in states
        d = film(x1, x2, q); d0 = film(x1, x2, q0)
        worstpsd  = max(worstpsd, -d[1]/abs(d[1]), -(d[1]*d[3] - d[2]^2)/(d[1]*d[3]))
        e11 = d0[1] - d[1]; e12 = d0[2] - d[2]; e22 = d0[3] - d[3]
        worstmono = max(worstmono, -e11/abs(e11), -(e11*e22 - e12^2)/(d0[1]*d0[3]))
        dx = 1e-12
        _, _, _, c1, c2 = electrostatic(x1, x2, 0.0, q)
        f1 = (electrostatic(x1 + dx, x2, 0.0, q)[1] - electrostatic(x1 - dx, x2, 0.0, q)[1])/(2*dx)
        f2 = (electrostatic(x1, x2 + dx, 0.0, q)[1] - electrostatic(x1, x2 - dx, 0.0, q)[1])/(2*dx)
        worstgrad = max(worstgrad, hypot(f1 - c1, f2 - c2)/max(hypot(c1, c2), 1e-18))
        z = [x1, 0.003, x2, -0.005, 0.2, zeros(7)...]; dz = zero(z)
        CoupledSystem!(dz, z, q, 0.0, 2.0)
        sc = [q.gc, 0.02, q.gc, 0.02, 3.0]; g = zeros(5)
        for j in 1:5
            hj = sc[j]*1e-6; zp = copy(z); zm = copy(z); zp[j] += hj; zm[j] -= hj
            g[j] = (energy(zp, q) - energy(zm, q))/(2*hj)
        end
        exact = dz[6] + dz[7] - (dz[8] + dz[9] + dz[10] + dz[11])
        worstpow = max(worstpow, abs(dot(g, dz[1:5]) - exact)/max(dz[12], 1e-25))
    end
    chkmax("film matrix positive semidefinite (5 states)", worstpsd, 1e-12)
    chkmax("D_lengthwise - D_vented positive semidefinite", worstmono, 1e-12)
    chkmax("dCt vs central difference (rel.)", worstgrad, 2e-5)
    chkmax("dE/dt vs power ledger, point-wise (rel.)", worstpow, 5e-5)
    # symmetry of the two walls
    Ca = electrostatic(1e-6, 2e-6, 0.0, q)[1]; Cb = electrostatic(-1e-6, -2e-6, 0.0, q)[1]
    chk("wall symmetry of Ct", Ca, Cb, 1e-13)
    verbose && println(allok ? ">>> selfcheck: ALL PASS" : ">>> selfcheck: FAILURES ABOVE -- do not trust the run")
    return allok
end

# Initialize a default Params instance and calculate dependent parameters
p = Params{Float64}()
p = create_params(p; verbose = false)

end # module AnalyticalModel

import .AnalyticalModel

# --------------------------------------- External Force ------------------------------------

# Sine Wave External Force
f = 20.0        # Frequency (Hz)
alpha = 4.95    # Applied acceleration constant (g). 4.95 -> panel (d); 2.7 -> panel (e).
                # Quasi-static contact threshold at 3 V is between 2.0 and 2.1
g = 9.80665     # Gravitational constant (m/s^2)
A = alpha*g
n_ramp = 4      # Ramp-up duration in drive cycles (C1 cosine ramp, zero end slopes)
ramp(t) = 0.5*(1 - cos(pi*min(t*f/n_ramp, 1.0)))
Fext_sine = t -> A*ramp(t)*sin(2*pi*f*t)

# ------------------------------------- Set Input Force ------------------------------------

# Set to `true` to use sine forcing, `false` for a near-contact displaced IC
# (free evolution: one contact episode probe, no external force)
use_sine = true
Fext_input = use_sine ? Fext_sine : (t -> 0.0)

# ------------------------------------ Initialize Parameters --------------------------------

p_new = deepcopy(AnalyticalModel.p)

# To change parameters, set the fields and REBUILD the derived quantities, e.g.
#   p_new.Vbias = 5.0;  AnalyticalModel.create_params(p_new; verbose = false)
# Submitted corrected model (lengthwise film, no structural loss):
#   p_new.vent_faces = 0; p_new.c1 = 0.0; p_new.ce = 0.0; p_new.panels = 512
#   AnalyticalModel.create_params(p_new; verbose = false)
AnalyticalModel.create_params(p_new; verbose = true)

run_selfcheck = true     # acceptance tests (a few seconds); see the header
run_selfcheck && AnalyticalModel.selfcheck()

# Initial conditions
if use_sine
    x10, x10dot, x20, x20dot = 0.0, 0.0, 0.0, 0.0
else
    # Near-contact probe IC (episode-scale evaluation without forcing)
    x20    = p_new.gc - 30e-9
    x20dot = 4e-3
    x10    = p_new.gc + 0.15e-6
    x10dot = 4e-3
end

# Equilibrated start: q = Vbias*Ct  <=>  Vout = 0 exactly.
# use_ledger appends the seven work integrals (states 6-12): an exact energy acceptance
# test for every run, at roughly twice the Jacobian cost. Set false for speed.
use_ledger = true
z0 = [x10, x10dot, x20, x20dot, 0.0]
use_ledger && (z0 = vcat(z0, zeros(7)))
n_cycles = 10
tspan  = use_sine ? (0.0, n_cycles/f) : (0.0, 600e-6)

# State scaling (corrected-model numerics). The finite-difference Jacobian perturbs each
# state by ~1.5e-8*max(|z|, 1): in SI metres that is 15 nm, far wider than the contact
# physics, so the solver integrates O(1) scaled states and abstol is one scalar.
xs = p_new.gc                                  # displacement scale
vs = xs*sqrt(p_new.k1/p_new.mtot)              # velocity scale
Es = p_new.k1*xs^2                             # energy scale
zscale = [xs, vs, xs, vs, max(abs(p_new.Vbias), 1.0)]
use_ledger && (zscale = vcat(zscale, fill(Es, 7)))
abstol = 1e-10                                 # on the scaled states
reltol = 1e-7
dtmax  = use_sine ? 2e-5 : 1e-6                # brackets the ~100 us contact events

# ---------------------------------- Solve Analytical Model ---------------------------------

function CoupledSystem_wrapper!(dz, z, p, t)
    AnalyticalModel.CoupledSystem!(dz, z .* zscale, p, t, Fext_input(t))
    dz ./= zscale
    return nothing
end

eqn = ODEProblem(CoupledSystem_wrapper!, z0 ./ zscale, tspan, p_new)

# Rodas5P with a FINITE-DIFFERENCE Jacobian (the quadrature arrays are Float64, not dual
# numbers); NO seam callbacks: the model is one continuous vector field through contact.
fdjac = isdefined(@__MODULE__, :AutoFiniteDiff) ? AutoFiniteDiff() : false
sol = solve(eqn, Rodas5P(autodiff = fdjac); abstol = abstol, reltol = reltol, dtmax = dtmax,
            maxiters = Int(1e7))

println(">>> collision_model version v5 (corrected two-coordinate model, v4 outline; vent_faces = ",
        p_new.vent_faces, ", c1 = ", p_new.c1, ", ce = ", p_new.ce, ") <<<")
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)
println("Solver stats:  ", sol.stats)

# Energy acceptance test: max |ledger| / throughput should be < 1e-5 (submitted model: ~1e-9..1e-7)
if use_ledger
    E0     = AnalyticalModel.energy(z0, p_new)
    stride = max(1, div(length(sol.u), 20000))
    maxres = maximum(abs(AnalyticalModel.ledger(u .* zscale, E0, p_new)) for u in sol.u[1:stride:end])
    zend   = sol.u[end] .* zscale
    @printf("energy residual / throughput = %.3e   (accept < 1e-5)\n", maxres/max(zend[12], eps()*Es))
    @printf("load energy ER(t_end) = %.6e J\n", zend[8])
    if use_sine && n_cycles > 4
        ER4 = (sol(tspan[2] - 4/f) .* zscale)[8]
        @printf("last-four-cycle mean load power = %.6e W   (submitted model at 4.95 g: 9.8853e-13)\n",
                (zend[8] - ER4)/(4/f))
    end
end

# ----------------------------------------- Plotting -----------------------------------------
AM = AnalyticalModel

# Sample states + observables + reconstructed forces on a uniform grid from the
# dense interpolant. Uniform grids avoid the rendering artifacts of plotting at raw
# adaptive steps (sparse in cruise, ultra-dense in taps).
function sample_window(sol, p, t0, t1; dt = 2e-6, Fext = Fext_input)
    tg = collect(max(t0, sol.t[1]):dt:min(t1, sol.t[end]))
    U  = Array(sol(tg)) .* zscale                 # back to SI units
    n  = length(tg)
    Ct = zeros(n); Fs = zeros(n); Fc = zeros(n); Fw = zeros(n); Fb = zeros(n)
    Fe1 = zeros(n); Fe2 = zeros(n); Fd1 = zeros(n); Fd2 = zeros(n)
    for i in 1:n
        z1, z2, z3, z4, z5 = U[1,i], U[2,i], U[3,i], U[4,i], U[5,i]
        Fs[i] = AM.spring(z1, p) - p.c1*z2
        fc, fw, _, _ = AM.collision(z1, z3, z4, p)
        Fc[i] = fc; Fw[i] = fw
        Fb[i] = -p.nb*p.ce*(z2 - z4)
        fd1, fd2, _ = AM.damping(z1, z2, z3, z4, p)
        Fd1[i] = fd1; Fd2[i] = fd2
        ct, fe1, fe2, _, _ = AM.electrostatic(z1, z3, z5, p)
        Ct[i] = ct; Fe1[i] = fe1; Fe2[i] = fe2
    end
    V   = U[5,:]                                  # Vout is the state
    Q   = [(p.Vbias - V[i])*Ct[i] for i in 1:n]   # charge as observable
    pen = [abs(U[3,i]) - p.gc for i in 1:n]       # nominal tip overlap
    ae  = [Fext(t) for t in tg]
    return (; t = tg, x1 = U[1,:], x1dot = U[2,:], x2 = U[3,:], x2dot = U[4,:],
              Q, V, Ct, pen, Fs, Fc, Fb, Fw, Fd1, Fd2, Fe1, Fe2, ae)
end

# Uniform 10 us overview grid for the state and force plots
Wover = sample_window(sol, p_new, sol.t[1], sol.t[end]; dt = 1e-5)
to = Wover.t

p3  = plot(to, Wover.x1,    xlabel = "Time (s)", ylabel = "x1 (m)",     title = "Shuttle Mass Displacement (x1)", label = "");    display(p3)
p4  = plot(to, Wover.x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)", label = "");    display(p4)
p5  = plot(to, Wover.x2,    xlabel = "Time (s)", ylabel = "x2 (m)",     title = "Mobile Electrode Tip Displacement (x2)", label = "")
hline!(p5, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); display(p5)
p6  = plot(to, Wover.x2dot, xlabel = "Time (s)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Tip Velocity (x2dot)", label = ""); display(p6)
p7  = plot(to, Wover.Q,     xlabel = "Time (s)", ylabel = "Q (C)",      title = "Charge (observable)", label = "");               display(p7)
p8  = plot(to, Wover.V,     xlabel = "Time (s)", ylabel = "Vout (V)",   title = "Output Voltage (state)", label = "");            display(p8)

# Diagnostics: penetration and total capacitance
p8b = plot(to, Wover.pen .* 1e9, xlabel = "Time (s)", ylabel = "|x2|-gc (nm)",
           title = "Penetration (contact when > 0)", label = "")
hline!(p8b, [0.0]; ls = :dash, lc = :gray, label = ""); display(p8b)
p8c = plot(to, Wover.Ct .* 1e12, xlabel = "Time (s)", ylabel = "Ctotal (pF)",
           title = "Total Capacitance", label = ""); display(p8c)

p9   = plot(to, Wover.Fs, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension + Stopper Force on x1 (incl. -c1*x1dot)", label = ""); display(p9)
p10  = plot(to, Wover.Fc, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Electrode Coupling Force (on x1; -Fc on x2)", label = ""); display(p10)
p10a = plot(to, Wover.Fb, xlabel = "Time (s)", ylabel = "Fb (N)", title = "Beam Relative-Damping Force (on x1; -Fb on x2)", label = ""); display(p10a)
p10b = plot(to, Wover.Fw, xlabel = "Time (s)", ylabel = "Fw (N)", title = "Tip Contact (Hunt-Crossley wall) Force on x2", label = ""); display(p10b)
p11  = plot(to, [Wover.Fd1 Wover.Fd2], xlabel = "Time (s)", ylabel = "Fd (N)", title = "Squeeze-Film Force (projected Reynolds film)",
            label = ["on x1" "on x2"], legend = :topright); display(p11)
p12  = plot(to, [Wover.Fe1 Wover.Fe2], xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force (attractive)",
            label = ["on x1" "on x2"], legend = :topright); display(p12)
p13  = plot(to, Wover.ae, xlabel = "Time (s)", ylabel = "a_ext (m/s^2)", title = "Applied Base Acceleration", label = ""); display(p13)