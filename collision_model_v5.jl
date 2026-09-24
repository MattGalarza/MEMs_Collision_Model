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
       energy, energy_parts, ledger, selfcheck
 
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
    Vbias::T = 3       # Bias voltage
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
 
# Stored energy by reservoir, in base-relative coordinates:
#   Tk  kinetic (consistent mass matrix)      Usp suspension + soft/hard stoppers
#   Ube electrode bending                     Uw  tip-wall contact potential
#   Ue  electrical field energy, Ct*Vc^2/2
# and the exact residual
#   ledger = E - E0 - Wbase - Wbias + ER + Dfilm + Dstruct + Dwall  (= 0 analytically)
function energy_parts(z, p, Ctotal)
    x1, v1, x2, v2, Vout = z[1], z[2], z[3], z[4], z[5]
    Tk  = 0.5*(p.M[1,1]*v1*v1 + 2*p.M[1,2]*v1*v2 + p.M[2,2]*v2*v2)
    Usp = p.k1*x1^2/2 + p.k3*x1^4/4
    Ube = p.nb*p.ke*(x2 - x1)^2/2
    Uw  = 0.0
    for r in (-1.0, 1.0)
        Usp += p.kss*softpos(r*x1 - p.gss, p.epss)^2/2
        isfinite(p.ghs) && (Usp += p.khs*softpos(r*x1 - p.ghs, p.epss)^(p.phs + 1)/(p.phs + 1))
        Uw  += p.nb*p.kw*softpos(r*x2 - p.gc, p.epsw)^(p.pw + 1)/(p.pw + 1)
    end
    Ue = 0.5*Ctotal*(p.Vbias - Vout)^2
    return (; Tk, Usp, Ube, Uw, Ue, E = Tk + Usp + Ube + Uw + Ue)
end
energy_parts(z, p) = energy_parts(z, p, electrostatic(z[1], z[3], z[5], p)[1])
energy(z, p) = energy_parts(z, p).E
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
f = 200.0        # Frequency (Hz)
alpha = 2.0    # Applied acceleration constant (g). 4.95 -> panel (d); 2.7 -> panel (e).
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
 
println(">>> collision_model version v5.1 (corrected two-coordinate model, v4 outline; vent_faces = ",
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
 
# Sample states + observables + reconstructed forces, energies and powers on a uniform grid
# from the dense interpolant. Uniform grids avoid the rendering artifacts of plotting at raw
# adaptive steps (sparse in cruise, ultra-dense in taps).
function sample_window(sol, p, t0, t1; dt = 2e-6, Fext = Fext_input)
    tg = collect(max(t0, sol.t[1]):dt:min(t1, sol.t[end]))
    U  = Array(sol(tg)) .* zscale                 # back to SI units
    n  = length(tg)
    Ct = zeros(n); Fs = zeros(n); Fc = zeros(n); Fw = zeros(n); Fb = zeros(n)
    Fe1 = zeros(n); Fe2 = zeros(n); Fd1 = zeros(n); Fd2 = zeros(n)
    Tk = zeros(n); Usp = zeros(n); Ube = zeros(n); Uw = zeros(n); Ue = zeros(n); E = zeros(n)
    Pb = zeros(n); Pe = zeros(n); PR = zeros(n); Pf = zeros(n); Ps = zeros(n); Pw = zeros(n)
    for i in 1:n
        z1, z2, z3, z4, z5 = U[1,i], U[2,i], U[3,i], U[4,i], U[5,i]
        ai = Fext(tg[i])
        Fs[i] = AM.spring(z1, p) - p.c1*z2
        fc, fw, _, pwl = AM.collision(z1, z3, z4, p)
        Fc[i] = fc; Fw[i] = fw
        Fb[i] = -p.nb*p.ce*(z2 - z4)
        fd1, fd2, _ = AM.damping(z1, z2, z3, z4, p)
        Fd1[i] = fd1; Fd2[i] = fd2
        ct, fe1, fe2, _, _ = AM.electrostatic(z1, z3, z5, p)
        Ct[i] = ct; Fe1[i] = fe1; Fe2[i] = fe2
        en = AM.energy_parts(U[1:5, i], p, ct)
        Tk[i] = en.Tk; Usp[i] = en.Usp; Ube[i] = en.Ube; Uw[i] = en.Uw; Ue[i] = en.Ue; E[i] = en.E
        # instantaneous powers: dE/dt = Pb + Pe - PR - Pf - Ps - Pw
        Pb[i] = -ai*(p.beta[1]*z2 + p.beta[2]*z4)     # base excitation (signed)
        Pe[i] = p.Vbias*z5/p.Rload                    # bias source (signed)
        PR[i] = z5*z5/p.Rload                         # load resistor
        Pf[i] = -(fd1*z2 + fd2*z4)                    # squeeze film
        Ps[i] = p.c1*z2*z2 + p.nb*p.ce*(z4 - z2)^2    # structure: shuttle c1 + beam ce
        Pw[i] = pwl                                   # tip wall (Hunt-Crossley)
    end
    V    = U[5,:]                                  # Vout is the state
    Q    = [(p.Vbias - V[i])*Ct[i] for i in 1:n]   # charge as observable
    pen  = [abs(U[3,i]) - p.gc for i in 1:n]       # nominal tip overlap
    htip = [p.h_eff + AM.softpos(p.gc - abs(U[3,i]), p.epsg) for i in 1:n]   # near-wall tip air gap
    ae   = [Fext(t) for t in tg]
    led  = size(U, 1) >= 12                        # work-integral states present?
    col(k) = led ? U[k,:] : fill(NaN, n)
    return (; t = tg, x1 = U[1,:], x1dot = U[2,:], x2 = U[3,:], x2dot = U[4,:],
              Q, V, Ct, pen, htip, Fs, Fc, Fb, Fw, Fd1, Fd2, Fe1, Fe2, ae,
              Tk, Usp, Ube, Uw, Ue, E, Pb, Pe, PR, Pf, Ps, Pw,
              Wb = col(6), Wbias = col(7), ER = col(8), Df = col(9), Ds = col(10),
              Dw = col(11), thr = col(12), led)
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
 
# ==================== (2) LAST-TWO-CYCLE TWIN OF EVERY STATE / FORCE PLOT ====================
Tdrive = 1/f
Wzoom  = sample_window(sol, p_new, sol.t[end] - 2*Tdrive, sol.t[end]; dt = 2e-6)
# NOTE: 2 us sampling resolves the chatter structure and the ~100 us contacts; the alias-free
# view of a single impact is the 0.1 us close-up in section (4).
tzs = Wzoom.t .* 1e3
 
p3z  = plot(tzs, Wzoom.x1,    xlabel = "t (ms)", ylabel = "x1 (m)",      title = "Shuttle Mass Displacement (x1) - last 2 cycles", label = "")
hline!(p3z, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); display(p3z)
p4z  = plot(tzs, Wzoom.x1dot, xlabel = "t (ms)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot) - last 2 cycles", label = ""); display(p4z)
p5z  = plot(tzs, Wzoom.x2,    xlabel = "t (ms)", ylabel = "x2 (m)",      title = "Mobile Electrode Tip Displacement (x2) - last 2 cycles", label = "")
hline!(p5z, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); display(p5z)
p6z  = plot(tzs, Wzoom.x2dot, xlabel = "t (ms)", ylabel = "x2dot (m/s)", title = "Mobile Electrode Tip Velocity (x2dot) - last 2 cycles", label = ""); display(p6z)
p7z  = plot(tzs, Wzoom.Q,     xlabel = "t (ms)", ylabel = "Q (C)",       title = "Charge (observable) - last 2 cycles", label = ""); display(p7z)
p8z  = plot(tzs, Wzoom.V,     xlabel = "t (ms)", ylabel = "Vout (V)",    title = "Output Voltage (state) - last 2 cycles", label = ""); display(p8z)
p8bz = plot(tzs, Wzoom.pen .* 1e9, xlabel = "t (ms)", ylabel = "|x2|-gc (nm)", title = "Penetration (contact when > 0) - last 2 cycles", label = "")
hline!(p8bz, [0.0]; ls = :dash, lc = :gray, label = ""); display(p8bz)
p8cz = plot(tzs, Wzoom.Ct .* 1e12, xlabel = "t (ms)", ylabel = "Ctotal (pF)", title = "Total Capacitance - last 2 cycles", label = ""); display(p8cz)
p9z   = plot(tzs, Wzoom.Fs, xlabel = "t (ms)", ylabel = "Fs (N)", title = "Suspension + Stopper Force on x1 - last 2 cycles", label = ""); display(p9z)
p10z  = plot(tzs, Wzoom.Fc, xlabel = "t (ms)", ylabel = "Fc (N)", title = "Electrode Coupling Force - last 2 cycles", label = ""); display(p10z)
p10az = plot(tzs, Wzoom.Fb, xlabel = "t (ms)", ylabel = "Fb (N)", title = "Beam Relative-Damping Force - last 2 cycles", label = ""); display(p10az)
p10bz = plot(tzs, Wzoom.Fw, xlabel = "t (ms)", ylabel = "Fw (N)", title = "Tip Contact (Hunt-Crossley wall) Force - last 2 cycles", label = ""); display(p10bz)
p11z  = plot(tzs, hcat(Wzoom.Fd1, Wzoom.Fd2), xlabel = "t (ms)", ylabel = "Fd (N)", title = "Squeeze-Film Force - last 2 cycles",
             label = ["on x1" "on x2"], legend = :topright); display(p11z)
p12z  = plot(tzs, hcat(Wzoom.Fe1, Wzoom.Fe2), xlabel = "t (ms)", ylabel = "Fe (N)", title = "Electrostatic Force - last 2 cycles",
             label = ["on x1" "on x2"], legend = :topright); display(p12z)
p13z  = plot(tzs, Wzoom.ae, xlabel = "t (ms)", ylabel = "a_ext (m/s^2)", title = "Applied Base Acceleration - last 2 cycles", label = ""); display(p13z)
 
# ================================ (1) ENERGY RELATIONS ====================================
# Identity being displayed (chapter Eq. 4.9):  dE/dt = Pb + Pbias - PR - Pf - Ps - Pw, with
#   E = Tk + Usp + Ube + Uw + Ue. The cumulative terms are the ledger states 6-12, so this
# section needs use_ledger = true; the stored energies and powers are reconstructed from the
# states and do not. Inputs (base, bias) are signed; the four losses are nonnegative.
function energy_plots(Wd, tx, xl, tag)
    Ubw = Wd.Ube .+ Wd.Uw
    dUe = Wd.Ue .- Wd.Ue[1]
    e1 = plot(tx, hcat(Wd.Tk, Wd.Usp, Ubw, dUe) .* 1e12, xlabel = xl, ylabel = "Energy (pJ)",
              title = string("Stored Energy by Reservoir - ", tag),
              label = ["kinetic" "suspension + stoppers" "beam bending + wall" "electrical (change)"],
              legend = :topleft)
    e2 = plot(tx, hcat(Wd.Wb, Wd.Wbias, Wd.ER, Wd.Df, Wd.Ds, Wd.Dw) .* 1e12, xlabel = xl, ylabel = "Cumulative energy (pJ)",
              title = string("Work In and Losses Out (ledger) - ", tag),
              label = ["W base (in)" "W bias (in)" "E load" "D film" "D struct (c1, ce)" "D wall"],
              legend = :topleft)
    dE  = Wd.E .- Wd.E[1]
    net = (Wd.Wb .- Wd.Wb[1]) .+ (Wd.Wbias .- Wd.Wbias[1]) .- (Wd.ER .- Wd.ER[1]) .-
          (Wd.Df .- Wd.Df[1]) .- (Wd.Ds .- Wd.Ds[1]) .- (Wd.Dw .- Wd.Dw[1])
    thr = max(Wd.thr[end] - Wd.thr[1], 1e-300)
    e3a = plot(tx, hcat(dE, net) .* 1e12, ylabel = "Energy (pJ)", title = string("Conservation Check - ", tag),
               label = ["E(t) - E(start)" "work in - losses"], ls = [:solid :dash], legend = :topleft)
    e3b = plot(tx, max.(abs.(dE .- net) ./ thr, 1e-16), xlabel = xl, ylabel = "|residual| / throughput",
               yscale = :log10, label = "")
    e3  = plot(e3a, e3b; layout = (2, 1), size = (800, 800))
    dEdt = Wd.Pb .+ Wd.Pe .- Wd.PR .- Wd.Pf .- Wd.Ps .- Wd.Pw
    e4a = plot(tx, hcat(Wd.Pb, dEdt) .* 1e9, ylabel = "Power (nW)", title = string("Power Flows - ", tag),
               label = ["base input Pb (signed)" "dE/dt"], legend = :topleft)
    e4b = plot(tx, max.(hcat(Wd.PR, Wd.Pf, Wd.Ps, Wd.Pw) .* 1e12, 1e-3), xlabel = xl, ylabel = "Loss power (pW)",
               yscale = :log10, label = ["load PR" "film Pf" "structure Ps" "wall Pw"], legend = :topleft)
    e4  = plot(e4a, e4b; layout = (2, 1), size = (800, 800))
    return e1, e2, e3, e4
end
 
if use_ledger
    e1, e2, e3, e4 = energy_plots(Wover, to, "Time (s)", "full run")
    display(e1); display(e2); display(e3); display(e4)
    e1z, e2z, e3z, e4z = energy_plots(Wzoom, tzs, "t (ms)", "last 2 cycles")
    display(e1z); display(e2z); display(e3z); display(e4z)
 
    # Energy budget over the last two cycles: inputs = sinks (+ change of stored energy)
    dlt(v) = v[end] - v[1]
    budget = [dlt(Wzoom.Wb), dlt(Wzoom.Wbias), -dlt(Wzoom.E), dlt(Wzoom.ER), dlt(Wzoom.Df), dlt(Wzoom.Ds), dlt(Wzoom.Dw)]
    e5 = bar(["W base", "W bias", "-dE stored", "E load", "D film", "D struct", "D wall"], budget .* 1e12;
             legend = false, ylabel = "Energy over the last 2 cycles (pJ)", xrotation = 30,
             title = "Energy Budget: first three bars = last four")
    display(e5)
    # Cross-check of the load energy from the electrical plane: int Vout dQ = int Vout^2/R dt
    WQ = sum(0.5*(Wzoom.V[i] + Wzoom.V[i+1])*(Wzoom.Q[i+1] - Wzoom.Q[i]) for i in 1:length(Wzoom.t)-1)
    println("\n================ ENERGY BUDGET, LAST TWO CYCLES ================")
    @printf("inputs : W base %+.4e J   W bias %+.4e J   release of stored energy %+.4e J\n", budget[1], budget[2], budget[3])
    @printf("sinks  : load %.4e J   film %.4e J   structure %.4e J   wall %.4e J\n", budget[4], budget[5], budget[6], budget[7])
    @printf("balance: inputs - sinks = %+.3e J   (%.2e of the inputs)\n", sum(budget[1:3]) - sum(budget[4:7]),
            abs(sum(budget[1:3]) - sum(budget[4:7]))/max(abs(sum(budget[1:3])), 1e-300))
    @printf("load share of the base work  E_load/W_base = %.4f   (an efficiency only if the window is periodic)\n",
            budget[4]/max(budget[1], 1e-300))
    @printf("loop area of the (Q, Vout) plane = %.4e J   vs ledger E_load = %.4e J\n", WQ, budget[4])
else
    println("Energy plots skipped: set use_ledger = true to carry the work integrals.")
end
 
# ================================ (3) PHASE SPACE ========================================
# Drawn from the ACCEPTED solver steps (dense in the taps, sparse in cruise), so every point
# lies on the computed trajectory and contact loops are not aliased by a uniform grid.
Zs  = Array(sol) .* zscale                          # nstate x nsteps, SI units
izz = findall(>=(sol.t[end] - 2*Tdrive), sol.t)     # last two drive cycles
gcu = p_new.gc*1e6
 
ph1  = plot(Zs[1,:] .* 1e6, Zs[2,:] .* 1e3, xlabel = "x1 (um)", ylabel = "x1dot (mm/s)",
            title = "Shuttle Phase Plane - full trajectory", label = "", lw = 0.5)
vline!(ph1, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); display(ph1)
ph1z = plot(Zs[1,izz] .* 1e6, Zs[2,izz] .* 1e3, xlabel = "x1 (um)", ylabel = "x1dot (mm/s)",
            title = "Shuttle Phase Plane - last 2 cycles", label = "", lw = 0.8)
vline!(ph1z, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); display(ph1z)
 
ph2  = plot(Zs[3,:] .* 1e6, Zs[4,:] .* 1e3, xlabel = "x2 (um)", ylabel = "x2dot (mm/s)",
            title = "Tip Phase Plane - full trajectory", label = "", lw = 0.5)
vline!(ph2, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); display(ph2)
ph2z = plot(Zs[3,izz] .* 1e6, Zs[4,izz] .* 1e3, xlabel = "x2 (um)", ylabel = "x2dot (mm/s)",
            title = "Tip Phase Plane - last 2 cycles", label = "", lw = 0.8)
vline!(ph2z, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); display(ph2z)
 
# Electrode bending plane: loops appear only while the tips are on a wall (contact mode)
ph3  = plot((Zs[3,:] .- Zs[1,:]) .* 1e9, (Zs[4,:] .- Zs[2,:]) .* 1e3, xlabel = "x2 - x1 (nm)",
            ylabel = "x2dot - x1dot (mm/s)", title = "Electrode Bending Plane - full trajectory", label = "", lw = 0.5); display(ph3)
ph3z = plot((Zs[3,izz] .- Zs[1,izz]) .* 1e9, (Zs[4,izz] .- Zs[2,izz]) .* 1e3, xlabel = "x2 - x1 (nm)",
            ylabel = "x2dot - x1dot (mm/s)", title = "Electrode Bending Plane - last 2 cycles", label = "", lw = 0.8); display(ph3z)
 
# Near-wall portrait: nominal overlap vs wall-normal tip velocity, both walls folded together
nearw(idx) = [i for i in idx if abs(Zs[3,i]) - p_new.gc > -200e-9]
inear = nearw(1:size(Zs, 2)); inearz = nearw(izz)
ph4  = scatter((abs.(Zs[3,inear]) .- p_new.gc) .* 1e9, sign.(Zs[3,inear]) .* Zs[4,inear] .* 1e3,
               xlabel = "|x2| - gc (nm)", ylabel = "wall-normal tip velocity (mm/s)", ms = 1.2, msw = 0,
               title = "Phase Portrait at the Contact Boundary - full trajectory", label = "")
vline!(ph4, [0.0]; ls = :dash, lc = :gray, label = ""); display(ph4)
ph4z = scatter((abs.(Zs[3,inearz]) .- p_new.gc) .* 1e9, sign.(Zs[3,inearz]) .* Zs[4,inearz] .* 1e3,
               xlabel = "|x2| - gc (nm)", ylabel = "wall-normal tip velocity (mm/s)", ms = 1.5, msw = 0,
               title = "Phase Portrait at the Contact Boundary - last 2 cycles", label = "")
vline!(ph4z, [0.0]; ls = :dash, lc = :gray, label = ""); display(ph4z)
 
# Electrical plane: the area enclosed per cycle is the energy delivered to the load
ph5  = plot(Wover.Q .* 1e12, Wover.V .* 1e3, xlabel = "Q (pC)", ylabel = "Vout (mV)",
            title = "Electrical Plane (Q, Vout) - full trajectory", label = "", lw = 0.5); display(ph5)
ph5z = plot(Wzoom.Q .* 1e12, Wzoom.V .* 1e3, xlabel = "Q (pC)", ylabel = "Vout (mV)",
            title = "Electrical Plane (Q, Vout) - last 2 cycles", label = "", lw = 0.8); display(ph5z)
 
# ========================== (4) LABELLED COLLISION CLOSE-UPS =============================
# Same-side contact sequences from the accepted steps: entries closer together than `gap`
# belong to one sequence (one wall, one half cycle).
function contact_sequences(tt, x2, x1dot, gc; gap = 8e-3)
    dls  = abs.(x2) .- gc
    ient = [i for i in 2:length(dls) if dls[i-1] < 0 && dls[i] >= 0]
    iext = [i for i in 2:length(dls) if dls[i-1] >= 0 && dls[i] < 0]
    seqs = NamedTuple[]
    isempty(ient) && return seqs, ient, iext
    k0 = 1
    for k in 2:length(ient)+1
        if k > length(ient) || tt[ient[k]] - tt[ient[k-1]] > gap
            ks   = k0:k-1
            jx   = findfirst(>(ient[ks[end]]), iext)
            tend = jx === nothing ? tt[end] : tt[iext[jx]]
            vin  = [abs(x1dot[ient[j]]) for j in ks]
            push!(seqs, (t0 = tt[ient[k0]], t1 = tend, side = sign(x2[ient[k0]]), n = length(ks),
                         vmax = maximum(vin), ihard = ient[ks[argmax(vin)]]))
            k0 = k
        end
    end
    return seqs, ient, iext
end
 
# Linear-interpolated zero crossings of d(t): upward = entries, downward = exits
function crossings(t, d)
    tin = Float64[]; tout = Float64[]
    for i in 2:length(d)
        if d[i-1] < 0 && d[i] >= 0
            push!(tin, t[i-1] + (t[i] - t[i-1])*(-d[i-1])/(d[i] - d[i-1]))
        elseif d[i-1] >= 0 && d[i] < 0
            push!(tout, t[i-1] + (t[i] - t[i-1])*d[i-1]/(d[i-1] - d[i]))
        end
    end
    return tin, tout
end
 
# Eight labelled panels for one window. All quantities are projected on the wall normal
# r = +-1, so "positive" always means "toward / into the wall".
function collision_figure(Wd, p; tag = "", unit = :us)
    r    = sign(Wd.x2[argmax(abs.(Wd.x2))])
    dl   = r .* Wd.x2 .- p.gc                       # nominal overlap delta; contact when > 0
    tin, tout = crossings(Wd.t, dl)
    tref = isempty(tin) ? Wd.t[1] : tin[1]
    sc   = unit === :us ? 1e6 : 1e3
    ul   = unit === :us ? "us" : "ms"
    tu   = (Wd.t .- tref) .* sc
    xlab = string("time from first nominal contact (", ul, ")")
    xL   = tu[1] + 0.02*(tu[end] - tu[1]); xR = tu[end] - 0.02*(tu[end] - tu[1])
    nm(v) = round(Int, v*1e9)
    function marks!(pl)
        isempty(tin)  || vline!(pl, (tin .- tref) .* sc;  ls = :dot, lc = :green, lw = 0.8, label = "")
        isempty(tout) || vline!(pl, (tout .- tref) .* sc; ls = :dot, lc = :red,   lw = 0.8, label = "")
        return pl
    end
    function hlabel!(pl, yv, str, yspan; side = :left)
        hline!(pl, [yv]; ls = :dash, lc = :gray, label = "")
        annotate!(pl, side === :left ? xL : xR, yv + 0.05*yspan, text(str, 7, side))
        return pl
    end
    span(v) = max(maximum(v) - minimum(v), 1e-30)
 
    # (a) shuttle and tip against the contact boundary
    s1 = (r .* Wd.x1 .- p.gc) .* 1e9; s2 = dl .* 1e9
    pa = plot(tu, hcat(s1, s2); label = ["shuttle  r*x1 - gc" "tip  r*x2 - gc"], ylabel = "nm", legend = :bottomright,
              title = "Shuttle and tip relative to the contact boundary")
    hlabel!(pa, 0.0, "x = gc: nominal contact (tip air gap = h_eff)", span(s1)); marks!(pa)
 
    # (b) the physical tip air gap on a log axis, with the three lengths that matter there
    hn = Wd.htip .* 1e9
    pb = plot(tu, hn; yscale = :log10, label = "", ylabel = "tip air gap (nm)",
              title = "Gap closure: h(tip) = h_eff + softpos(gc - r*x2)")
    hline!(pb, [p.h_eff, p.h_eff + p.ls, p.hd] .* 1e9; ls = :dash, lc = :gray, label = "")
    annotate!(pb, xR, p.h_eff*1e9*1.12, text(string("h_eff = ", nm(p.h_eff), " nm: residual gap (floor)"), 7, :right))
    annotate!(pb, xL, (p.h_eff + p.ls)*1e9*1.12, text(string("h_eff + ls = ", nm(p.h_eff + p.ls), " nm: tip vent starts to seal"), 7, :left))
    annotate!(pb, xR, p.hd*1e9*1.30, text(string("hd = 2Tp/ep = ", nm(p.hd), " nm: coating-equivalent gap"), 7, :right))
    marks!(pb)
 
    # (c) nominal overlap with the sealing window
    top = 1.4*max(maximum(dl)*1e9, 1.0)
    pc = plot(tu, dl .* 1e9; label = "", ylabel = "delta = r*x2 - gc (nm)", ylims = (-3*p.ls*1e9, top),
              title = "Nominal overlap and the tip-vent sealing window")
    hlabel!(pc, 0.0, "delta = 0", top; side = :right)
    hlabel!(pc, p.ls*1e9, string("+ls = ", nm(p.ls), " nm: vent sealed (chi = 1)"), top)
    hlabel!(pc, -p.ls*1e9, string("-ls: vent open (chi = 0)"), top)
    marks!(pc)
 
    # (d) electrode bending
    bn = r .* (Wd.x2 .- Wd.x1) .* 1e9
    pd = plot(tu, bn; label = "", ylabel = "r*(x2 - x1) (nm)", title = "Electrode bending (negative: shuttle pushes past the tip)")
    hline!(pd, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pd)
 
    # (e) wall-normal velocities
    pe = plot(tu, hcat(r .* Wd.x1dot, r .* Wd.x2dot) .* 1e3; label = ["shuttle  r*x1dot" "tip  r*x2dot"], ylabel = "mm/s",
              legend = :topright, title = "Velocities normal to the wall")
    hline!(pe, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pe)
 
    # (f) forces on the tip coordinate x2, projected on the wall normal
    Ft = hcat((-r) .* Wd.Fc, r .* Wd.Fw, r .* Wd.Fe2, r .* Wd.Fd2, (-r) .* Wd.Fb) .* 1e6
    pf = plot(tu, Ft; label = ["bending" "wall contact" "electrostatic" "squeeze film" "beam damping"], ylabel = "uN",
              legend = :bottomright, title = "Forces on the tip coordinate (+ = into the wall)")
    hline!(pf, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pf)
 
    # (g) output voltage
    pg = plot(tu, Wd.V .* 1e3; label = "", ylabel = "Vout (mV)", xlabel = xlab, title = "Output voltage")
    marks!(pg)
 
    # (h) phase portrait at the contact boundary
    ph = plot(dl .* 1e9, r .* Wd.x2dot .* 1e3; label = "", xlabel = "delta = r*x2 - gc (nm)", ylabel = "r*x2dot (mm/s)",
              xlims = (-3*p.ls*1e9, top), title = "Phase portrait at the contact boundary")
    vline!(ph, [0.0]; ls = :dash, lc = :gray, label = "")
 
    # metrics of the FIRST contact in the window (shuttle-level restitution, contact time, peaks)
    met = (; side = r, n_contacts = length(tin), t_in = NaN, t_contact = NaN, v_in = NaN, v_out = NaN,
             e_shuttle = NaN, overlap_max = NaN, bending_max = NaN, Fw_max = NaN)
    jo = isempty(tin) ? nothing : findfirst(>(tin[1]), tout)
    if jo !== nothing
        i1 = findfirst(>=(tin[1]), Wd.t); i2 = findfirst(>=(tout[jo]), Wd.t)
        vin = r*Wd.x1dot[i1]; vout = r*Wd.x1dot[i2]
        met = (; side = r, n_contacts = length(tin), t_in = tin[1], t_contact = tout[jo] - tin[1], v_in = vin, v_out = vout,
                 e_shuttle = -vout/vin, overlap_max = maximum(dl[i1:i2]), bending_max = maximum(abs, bn[i1:i2])*1e-9,
                 Fw_max = maximum(abs, Wd.Fw[i1:i2]))
        annotate!(pe, xL, minimum(r .* Wd.x1dot)*1e3 + 0.08*span(r .* Wd.x1dot)*1e3,
                  text(string("shuttle in ", round(vin*1e3; digits = 2), ", out ", round(vout*1e3; digits = 2),
                              " mm/s:  e = ", round(-vout/vin; digits = 3)), 7, :left))
        annotate!(pc, xL, 0.88*top, text(string("first contact ", round((tout[jo] - tin[1])*1e6; digits = 1),
                              " us, peak overlap ", round(maximum(dl[i1:i2])*1e9; digits = 2), " nm"), 7, :left))
    end
    fig = plot(pa, pb, pc, pd, pe, pf, pg, ph; layout = (4, 2), size = (1200, 1500),
               plot_title = string("Collision close-up - ", tag, "   (green / red dotted: contact made / lost)"),
               titlefontsize = 9, guidefontsize = 8, tickfontsize = 7, legendfontsize = 7)
    return fig, met
end
 
seqs, ient, iext = contact_sequences(sol.t, Zs[3,:], Zs[2,:], p_new.gc)
if isempty(seqs)
    println("\nNo contact this run (closest approach ",
            round(-maximum(abs.(Zs[3,:]) .- p_new.gc)*1e9; digits = 1), " nm). Collision close-ups skipped.")
else
    println("\n================ CONTACT SEQUENCES (one line per wall visit) ================")
    println("   start (s)    wall   contacts   hardest shuttle entry (mm/s)   duration (ms)")
    for s in seqs
        @printf("   %9.5f     %+d     %5d        %8.3f                     %8.3f\n",
                s.t0, Int(s.side), s.n, s.vmax*1e3, (s.t1 - s.t0)*1e3)
    end
    @printf("wall +1: %d visits, %d contacts     wall -1: %d visits, %d contacts\n",
            count(s -> s.side > 0, seqs), sum([s.n for s in seqs if s.side > 0]; init = 0),
            count(s -> s.side < 0, seqs), sum([s.n for s in seqs if s.side < 0]; init = 0))
 
    recent = [s for s in seqs if s.t0 >= sol.t[end] - 2*Tdrive]
    isempty(recent) && (recent = seqs)
 
    # (4a) the hardest single impact of the last two cycles, 0.1 us sampling
    sA  = recent[argmax([s.vmax for s in recent])]
    jA  = findfirst(>(sA.ihard), iext)
    tA1 = jA === nothing ? sol.t[end] : sol.t[iext[jA]]
    Wimp = sample_window(sol, p_new, sol.t[sA.ihard] - 40e-6, tA1 + 80e-6; dt = 1e-7)
    figA, mA = collision_figure(Wimp, p_new; unit = :us,
                                tag = string("hardest impact, t = ", round(sol.t[sA.ihard]; digits = 5), " s, wall ", Int(sA.side)))
    display(figA)
 
    # (4b) the longest chatter sequence of the last two cycles: approach - chatter - dwell - release
    sB   = recent[argmax([s.n for s in recent])]
    Wep  = sample_window(sol, p_new, sB.t0 - 0.5e-3, sB.t1 + 0.5e-3; dt = max(1e-7, (sB.t1 - sB.t0 + 1e-3)/40000))
    figB, mB = collision_figure(Wep, p_new; unit = :ms,
                                tag = string("longest sequence, ", sB.n, " contacts from t = ", round(sB.t0; digits = 5), " s, wall ", Int(sB.side)))
    display(figB)
 
    println("\n================ CLOSE-UP METRICS (first contact of each window) ================")
    for (nm_, m) in (("hardest impact  ", mA), ("longest sequence", mB))
        @printf("%s: wall %+d, %d contacts in window; first contact %.1f us; shuttle %.3f -> %.3f mm/s (e = %.3f); peak overlap %.2f nm; peak bending %.1f nm; peak wall force %.1f uN\n",
                nm_, Int(m.side), m.n_contacts, m.t_contact*1e6, m.v_in*1e3, m.v_out*1e3, m.e_shuttle,
                m.overlap_max*1e9, m.bending_max*1e9, m.Fw_max*1e6)
    end
end




#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

using Plots, Printf

"""
    animate_electrode_contact(sol, p, zscale; frequency, outfile, kwargs...)

Animate the last two forcing cycles of the existing SCALED solution. No new
simulation is run. Root = state 1, tip = state 3; both are restored to SI units.
Beam shape: w(s,t) = (1-phi(s))*x1(t) + phi(s)*x2(t), using AnalyticalModel.bshape.

`window=(ta,tb)` selects an explicit interval in seconds, useful for slow motion.
`seconds=10, fps=30` makes 300 frames. The movie is an overview, not a count of
contact crossings; very fast rebounds need a narrow window / more frames.
`frame_times` optionally supplies increasing physical times for variable-speed
playback. Use `animate_electrode_cycle` below for both wall visits of one cycle.
The diagram compresses the vertical dimension; horizontal geometry is in µm.
Orange means nominal contact, r*x2 >= gc. The tiny residual film is reported
numerically, and the tip is NEVER clamped or snapped to a wall for display.
`outfile` accepts .gif or .mp4; a first-frame PNG is also saved beside the movie.
"""
function animate_electrode_contact(sol, p, zscale;
        frequency, cycles=2, window=nothing, seconds=10, fps=30,
        outfile="electrode_last_two_cycles.gif", trace_dt=2e-6,
        frame_times=nothing, visits=NamedTuple[])
    @assert length(zscale) == length(sol.u[1]) >= 5
    @assert frequency > 0 && cycles > 0 && seconds > 0 && fps >= 1 && trace_dt > 0
    ta, tb = window === nothing ? (max(sol.t[1], sol.t[end]-cycles/frequency), sol.t[end]) : window
    sol.t[1] <= ta < tb <= sol.t[end] || error("Animation window must lie inside sol.t.")
    ext = lowercase(splitext(outfile)[2])
    ext in (".gif", ".mp4") || error("Use a .gif or .mp4 output filename.")
    window === nothing && tb-ta < (cycles/frequency)*(1-1e-10) &&
        @warn "Solution is shorter than the requested cycle window; showing available time."
    dst = abspath(outfile); mkpath(dirname(dst)); png = splitext(dst)[1]*"_preview.png"
    si(t) = sol(t) .* zscale                     # IMPORTANT: recover physical states
    times = frame_times === nothing ? collect(range(ta,tb;length=max(2,round(Int,seconds*fps)))) : collect(frame_times)
    length(times)>=2 && all(isfinite,times) && all(diff(times).>0) ||
        error("frame_times must be finite and strictly increasing.")
    ta<=first(times)<last(times)<=tb || error("frame_times must lie inside the window.")
    nf=length(times); variable_speed=frame_times!==nothing
    nt = clamp(ceil(Int, (tb-ta)/trace_dt)+1, 1001, 100001)
    tt = range(ta, tb; length=nt); U = hcat(si.(tt)...)
    tx = (tt .- ta).*1e3; X1=U[1,:].*1e6; X2=U[3,:].*1e6
    bend = (U[3,:].-U[1,:]).*1e9
    # Full mobile beam, clamped root at s=0 and free tip at s=Lf.
    s = collect(range(0, p.Lf; length=101)); yy = 100 .* s ./ p.Lf
    phi = [AnalyticalModel.bshape(q, p.ke, p) for q in s]
    halfwidth = [AnalyticalModel.bwidth(q, p)*1e6/2 for q in s]
    gp = p.g0-2p.Tp; ybottom = 100*(1-p.Leff/p.Lf)
    # Fixed facing walls reproduce gp + a*y - r*w over the actual overlap.
    inner_top = (p.wb/2 + gp)*1e6
    inner_bottom = (AnalyticalModel.bwidth(p.Lf-p.Leff,p)/2 + gp + p.a*p.Leff)*1e6
    fixed_center = max(inner_bottom,inner_top) + 9.0
    shuttle_half = 2fixed_center
    extent = max(shuttle_half + maximum(abs, X1) + 8, 2fixed_center-inner_top+8)
    pad(v) = max(0.12*(maximum(v)-minimum(v)), 0.01)
    xlim = 1.12max(maximum(abs, X1), maximum(abs, X2), p.gc*1e6)
    blim = (minimum(bend)-pad(bend), maximum(bend)+pad(bend))
    teal="#187E89"; gray="#C7C9CC"; orange="#D66B22"; red="#C83B43"
    box(xa,xb,ya,yb) = Shape([xa,xb,xb,xa], [ya,ya,yb,yb])
    common = (; fontfamily="Computer Modern", titlefontsize=11, guidefontsize=10,
        tickfontsize=9, legendfontsize=9, background_color=:white, grid=false,
        linewidth=1.5, dpi=100, margin=4*Plots.mm)
    @printf("Animating %.6f to %.6f s; %d frames; physical step %.3f to %.3f µs/frame.\n",
        ta,tb,nf,1e6*minimum(diff(times)),1e6*maximum(diff(times)))
    mktempdir() do frame_dir
        anim = Animation(frame_dir)
        for (k,t) in enumerate(times)
            z=si(t); x1=z[1]*1e6; x2=z[3]*1e6; tm=(t-ta)*1e3
            w=(1 .- phi).*x1 .+ phi.*x2
            onleft=-z[3]>=p.gc; onright=z[3]>=p.gc
            side=z[3]<0 ? -1 : z[3]>0 ? 1 : (z[4]<0 ? -1 : 1)
            wall=side<0 ? "LEFT" : "RIGHT"
            j=isempty(visits) ? nothing : argmin([abs(t-(v.entry+v.exit)/2) for v in visits])
            visit=j===nothing ? nothing : visits[j]
            state=onleft ? "LEFT CONTACT" : onright ? "RIGHT CONTACT" :
                visit!==nothing && visit.entry<=t<=visit.exit ? "BETWEEN REBOUNDS - "*(visit.side<0 ? "LEFT" : "RIGHT") :
                side*z[4]>0 ? "APPROACH "*wall : "RELEASE / AWAY FROM "*wall
            hm=p.h_eff+AnalyticalModel.softpos(p.gc+z[3],p.epsg)
            hp=p.h_eff+AnalyticalModel.softpos(p.gc-z[3],p.epsg)
            a=plot(; common..., xlims=(-extent,extent), ylims=(-48,135),
                framestyle=:none, ticks=false, legend=false,
                title=@sprintf("2D electrode schematic  |  t = %.6f s",t))
            plot!(a,box(-extent+3,extent-3,118,130); c=gray,lc=:black,label="")
            annotate!(a,0,124,text("Fixed support",10,:black))
            for (r,active) in ((-1,onleft),(1,onright))
                ix=r.*[inner_bottom,inner_top,inner_top]
                ox=r.*[2fixed_center-inner_bottom,2fixed_center-inner_top,2fixed_center-inner_top]
                sy=[ybottom,100.,118.]
                plot!(a,Shape(vcat(ix,reverse(ox)),vcat(sy,reverse(sy)));
                    c=active ? orange : gray,lc=:black,label="")
            end
            plot!(a,box(x1-shuttle_half,x1+shuttle_half,-9,0);c=teal,lc=:black,label="")
            plot!(a,Shape(vcat(w.-halfwidth,reverse(w.+halfwidth)),vcat(yy,reverse(yy)));
                c=teal,lc=:black,label="")
            plot!(a,w,yy;c=:white,ls=:dash,lw=1,label="")
            scatter!(a,[x1,x2],[0,100];c=red,ms=4,msc=red,label="")
            annotate!(a,x1,-4.5,text("Shuttle  x1(t)",10,:white))
            annotate!(a,x2,105,text("x2(t)",10,:black))
            # Tip-gap brackets; the numerical readout resolves nanometre clearances.
            for (xa,xb) in ((-inner_top,x2-p.wb*1e6/2),(x2+p.wb*1e6/2,inner_top))
                plot!(a,[xa,xb],[110,110];c=teal,lw=1.5,label="")
                plot!(a,[xa,xa,NaN,xb,xb],[108,112,NaN,108,112];c=teal,lw=1,label="")
            end
            if onleft || onright
                r=onright ? 1 : -1
                scatter!(a,[x2+r*p.wb*1e6/2],[100];c=orange,ms=6,msc=:black,label="")
            end
            annotate!(a,0,-16,text(state,10,(onleft||onright) ? orange : teal))
            annotate!(a,0,-24,text(@sprintf("x1 = %.4f µm    x2 = %.4f µm",x1,x2),9,:black))
            annotate!(a,0,-31,text(@sprintf("Fluid h₋ = %.1f nm    h₊ = %.1f nm",hm*1e9,hp*1e9),9,:black))
            annotate!(a,0,-38,text(@sprintf("Nominal overlap |x2| − gc = %.2f nm",(abs(z[3])-p.gc)*1e9),9,:black))
            annotate!(a,0,-45,text("Horizontal scale retained; vertical dimension compressed",8,:gray35))
            b=plot(tx,hcat(X1,X2); common...,c=[teal red],label=["Shuttle x1" "Tip x2"],
                xlims=(0,(tb-ta)*1e3),ylims=(-xlim,xlim),ylabel="Displacement (µm)",
                xlabel="",title=variable_speed ? "One cycle; playback slower around wall visits" : "Motion and nominal contact limits",
                legend=:outertop,legend_column=2)
            for v in visits
                vspan!(b,([v.lo,v.hi].-ta).*1e3;c=orange,alpha=0.09,label="")
            end
            hline!(b,[-p.gc,p.gc].*1e6;c=:gray,ls=:dash,lw=1,label="")
            vline!(b,[tm];c=:black,lw=1,label="")
            scatter!(b,[tm,tm],[x1,x2];c=[teal,red],ms=4,msc=:white,label="")
            c=plot(tx,bend;common...,c=teal,label="",xlims=(0,(tb-ta)*1e3),ylims=blim,
                xlabel="Time from window start (ms)",ylabel="x2 − x1 (nm)",title="Relative electrode bending",legend=false)
            vline!(c,[tm];c=:black,lw=1,label="")
            scatter!(c,[tm],[(z[3]-z[1])*1e9];c=red,ms=4,msc=:white,label="")
            if visit!==nothing
                # Wall-relative coordinates resolve nanometre contact on the same
                # trajectory. This close-up switches between the two wall visits.
                v=visit; sel=findall(q->v.lo<=q<=v.hi,tt)
                td=(tt[sel].-v.entry).*1e3
                d1=(v.side.*U[1,sel].-p.gc).*1e9
                d2=(v.side.*U[3,sel].-p.gc).*1e9
                d=plot(td,hcat(d1,d2);common...,c=[teal red],label=["Shuttle" "Tip"],
                    legend=false,xlims=((v.lo-v.entry)*1e3,(v.hi-v.entry)*1e3),
                    xlabel="Time from first contact in this visit (ms)",ylabel="r*x - gc (nm)",
                    title="Visit $j/$(length(visits)): "*(v.side<0 ? "left" : "right")*" approach and release")
                hline!(d,[0.];c=:gray,ls=:dash,label="")
                vline!(d,([v.entry,v.exit].-v.entry).*1e3;c=orange,ls=:dot,label="")
                if v.lo<=t<=v.hi
                    vline!(d,[(t-v.entry)*1e3];c=:black,label="")
                    scatter!(d,fill((t-v.entry)*1e3,2),(v.side.*[z[1],z[3]].-p.gc).*1e9;
                        c=[teal,red],ms=4,msc=:white,label="")
                end
                right=plot(b,c,d;layout=(3,1))
            else
                right=plot(b,c;layout=(2,1))
            end
            fig=plot(a,right;layout=grid(1,2; widths=[0.50,0.50]),
                size=(1160,isempty(visits) ? 640 : 820),dpi=100,background_color=:white)
            k==1 && savefig(fig,png)
            frame(anim,fig)
            (k==1 || k%60==0 || k==nf) && println("  animation frame ",k,"/",nf)
        end
        ext==".gif" ? gif(anim,dst;fps=fps) : mp4(anim,dst;fps=fps)
    end
    println("Saved animation: ",dst)
    return (; path=dst, preview=png, tspan=(ta,tb), frames=nf,
        simulation_dt=variable_speed ? nothing : (tb-ta)/(nf-1),
        frame_times=times, simulation_dt_range=extrema(diff(times)))
end

"""
    animate_electrode_cycle(sol, p, zscale; frequency, outfile, kwargs...)

Replay the final complete forcing cycle: both wall visits, their approaches,
and releases. Default cycle boundaries are measured from `sol.t[1]`; set
`cycle_start` to choose a different start time. This does not assert steady state.

Playback allocates more frames to each visit, from `pre` seconds before first
contact to `post` seconds after final release. Physical time stays on screen.
One visit includes all rebounds on that side before the tip crosses the center.
Every detected contiguous contact episode contributes a frame at its sampled
maximum overlap. This helps show brief contacts; it is not an exact crossing
counter. No contact is invented when the solution does not reach a wall.
"""
function animate_electrode_cycle(sol,p,zscale; frequency, cycle_start=nothing,
        seconds=16, fps=30, pre=1e-3, post=2e-3, slowdown=8.0,
        trace_dt=2e-6, outfile="electrode_one_cycle_two_hits.gif")
    frequency>0 && seconds>0 && fps>=1 && pre>=0 && post>=0 && slowdown>=1 && trace_dt>0 ||
        error("Invalid frequency, playback, or sampling settings.")
    T=1/frequency; t0=first(sol.t); tend=last(sol.t)
    if cycle_start===nothing
        n=floor(Int,(tend-t0)/T+1e-10)
        n>=1 || error("A complete forcing cycle is needed. Use animate_electrode_contact for a short probe.")
        ta=t0+(n-1)*T; tb=min(tend,t0+n*T)
    else
        ta=Float64(cycle_start); tb=ta+T
    end
    t0<=ta<tb<=tend || error("Requested cycle must lie inside the solution.")
    # Accepted solver times retain short events even when the regular trace grid
    # is coarser. Dense interpolation supplies intermediate visual samples.
    tt=sort!(unique!(vcat(collect(range(ta,tb;length=max(1001,ceil(Int,T/trace_dt)+1))),
        [t for t in sol.t if ta<=t<=tb])))
    x2=[sol(t)[3]*zscale[3] for t in tt]; delta=abs.(x2).-p.gc
    hit=delta.>=0; visits=NamedTuple[]; peaks=Float64[]
    # Refine observed brackets in the dense solution; this is for movie timing,
    # not an independently converged physical event-detection calculation.
    function crossing(a,b,r)
        ga=r*sol(a)[3]*zscale[3]-p.gc
        for _ in 1:35
            m=(a+b)/2; gm=r*sol(m)[3]*zscale[3]-p.gc
            if (ga>=0)==(gm>=0); a=m; ga=gm; else; b=m; end
        end
        return (a+b)/2
    end
    i=1
    while i<=length(tt)
        r=x2[i]<0 ? -1 : 1; j=i
        while j<length(tt) && (x2[j+1]<0 ? -1 : 1)==r; j+=1; end
        ci=[k for k in i:j if hit[k]]
        if !isempty(ci)
            a=first(ci); b=last(ci)
            entry=a==1 ? ta : crossing(tt[a-1],tt[a],r)
            leave=b==length(tt) ? tb : crossing(tt[b],tt[b+1],r)
            push!(visits,(;side=r,entry,exit=leave,lo=max(ta,entry-pre),hi=min(tb,leave+post),
                complete_entry=a>1,complete_exit=b<length(tt)))
        end
        i=j+1
    end
    i=1
    while i<=length(tt)
        if !hit[i]; i+=1; continue; end
        j=i
        while j<length(tt) && hit[j+1]; j+=1; end
        k=i-1+argmax(@view delta[i:j]); push!(peaks,tt[k]); i=j+1
    end
    length(visits)==2 && sort([v.side for v in visits])==[-1,1] ||
        @warn "This computed cycle does not contain exactly one visit to each wall; showing the actual trajectory." visits=length(visits)
    any(v->!v.complete_entry || !v.complete_exit,visits) &&
        @warn "A wall visit crosses the cycle boundary. Set cycle_start to a free-travel instant to see both ends."
    # Uniform movie frames in weighted time => smaller physical increments near
    # contacts. Required peak samples are inserted without changing any state.
    weight(t)=any(v->v.lo<=t<=v.hi,visits) ? slowdown : 1.0
    clock=vcat(0.,cumsum(diff(tt).*[weight((tt[k]+tt[k+1])/2) for k in 1:length(tt)-1]))
    nbase=max(2,round(Int,seconds*fps)-length(peaks)); times=Float64[]
    for q in range(0.,last(clock);length=nbase)
        k=min(searchsortedlast(clock,q),length(clock)-1)
        push!(times,tt[k]+(q-clock[k])/(clock[k+1]-clock[k])*(tt[k+1]-tt[k]))
    end
    times[1]=ta; times[end]=tb; times=sort!(unique!(vcat(times,peaks)))
    println("One-cycle wall visits: ",length(visits),"; detected contact episodes shown: ",length(peaks))
    for (j,v) in enumerate(visits)
        @printf("  Visit %d (%s): first contact %.9f s, final release %.9f s\n",j,v.side<0 ? "left" : "right",v.entry,v.exit)
    end
    movie=animate_electrode_contact(sol,p,zscale;frequency,cycles=1,window=(ta,tb),
        seconds,fps,outfile,trace_dt,frame_times=times,visits)
    return (;movie...,visits,contact_peak_times=peaks)
end

animate_electrode_cycle(sol, p_new, zscale;
    frequency=f,
    seconds=16,
    fps=30,
    outfile="electrode_one_cycle_two_hits.gif")
