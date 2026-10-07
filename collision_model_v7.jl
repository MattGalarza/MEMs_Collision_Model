# collision_model_v7.jl -- v6 with the gas-film closure of the companion gas-film paper (Oct 2026).
# Every change is in the squeeze film; the other forces, the states and the energy ledger are untouched.
#   * one device-layer face vents by default (the device lies on its PCB): vent_faces = 1
#   * the film conductance carries the linearized-BGK kinetic flow factor (kinetic = true)
#   * the face drainage passes through the Stokes-cell exit factor f0(h) (exit_factor = true)
#   * the tip end vents through a pocket calibrated on the 3-D Stokes solution, sealing fraction 0.75 at
#     rest (chi_rest), in series with the contact sealing, which now starts at nominal contact (seal_shift)
#   * free path: the BGK equivalent length 74.87 nm (with sigmap = 1.016191)
# The v6 closure is recovered to round-off with vent_faces = 2, kinetic = false, exit_factor = false,
# seal_shift = false, chi_rest = 0, lambda = 70e-9, sigmap = 1.016.
# ce was fitted with the v6 film; check the printed chatter-interval ratio against the bench 0.73 (Params notes).

# ------------------------------------------ Libraries --------------------------------------
 
using DifferentialEquations, Plots, Printf, XLSX
 
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
using Printf
export Params, p, create_params, spring, collision, damping, electrostatic, CoupledSystem!,
       energy, energy_parts, ledger, report, energy_check, forces, capacitance, classtip, normq
 
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
    Lsp::T = 1600e-6     # Suspension spring length
    Lss::T = 600e-6     # Soft-stopper length
    gss::T = 14e-6       # Soft-stopper position
    gap_slope::T = NaN   # Facing-wall gap slope. NaN -> (wb - wt)/Lf
 
    # Array / suspension topology
    N::Int = 160         # Number of gap branches (N/2 compliant electrodes, two faces each)
    nsp_par::Int = 4     # Suspension: parallel chains
    nsp_ser::Int = 6     # Suspension: series spans per chain
    nss::Int = 2         # Soft-stopper cantilevers acting in parallel, per side
    gamma3::T = 1.0      # Cubic-stiffness geometry correction (not calibrated)
 
    # Mass and material properties
    m1::T = 2.0933e-6    # Shuttle-only mass. MUST exclude the explicit electrodes: the compliant
                         # ones enter through M; the stiff ones are added when orient = 1
    rho::T = 2330.0      # Density of silicon
    E::T = 170e9         # Young's modulus
    e::T = 8.85e-12      # Permittivity of free space
    ep::T = 3.2          # Relative permittivity of Parylene-C
    eta::T = 1.849e-5    # Viscosity of air
    lambda::T = 74.866e-9 # Equivalent free path ell = eta*sqrt(2/(rho_air*p_a)), 20 C and 1 atm: the length the
                          # BGK slip coefficient and the kinetic factor refer to (v6: 70e-9)
    sigmap::T = 1.016191  # BGK slip coefficient, diffuse walls (v6: 1.016)
 
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
    #      v7 CHECK: ce was fitted with the v6 film. The v7 film removes 12-16% of the shuttle speed per chatter
    #      flight (v6 film: 5-8%) and barely damps the contact mode itself (< 0.5%). With the contact keeping
    #      e ~ 0.86 (from ce) and a flight keeping ~0.85, a decoupled estimate puts the v7 chatter-interval
    #      ratio near the bench 0.73 with ce unchanged. The run prints the ratio; refit ce only if it misses.
    #      c1 must exclude the comb film, which is 21.0 uN s/m at rest in v7 (1.96 in v6).
    # vent_faces : squeeze-film drainage. 0 = lengthwise only (submitted model, gas
    #      leaves through the two ends of the 400 um overlap); 1 = one device-layer
    #      face open as well; 2 = both faces open (top surface + cavity underneath).
    #      Lengthwise-only overdamps by ~210x at rest, ~40x at a 1 um gap, 6.5x at
    #      contact relative to 2 (57x / 14x / 3.7x relative to 1) [ratios for the v6 closure].
    c::T = 1.0           # Film scale -- physical default
    c1::T = 5.29e-5      # Shuttle damping [N s/m]   (submitted model: 0)
    ce::T = 0.75e-4      # Beam relative damping, per beam [N s/m]   (submitted model: 0)
    vent_faces::Int = 1  # Film drainage: 0 lengthwise | 1 one face open (as built, on the PCB) | 2 both faces open
    vent_modes::Int = 4  # Thickness modes solved exactly; higher modes use the strip limit. v7, one face:
                         # 1, 2, 4, 8, 24 modes -> 3.563, 3.059, 2.989, 2.981, 2.980 mN s/m at contact (21.003 uN s/m
                         # at rest from 2 modes); use 8 for final figures (v6: 2.213, 2.083, 2.068, 2.066)

    # Gas-film closure (v7, gas-film paper). v6 closure: vent_faces = 2, kinetic = false, exit_factor = false,
    # seal_shift = false, chi_rest = 0, lambda = 70e-9, sigmap = 1.016 (reproduces v6 to round-off).
    kinetic::Bool = true     # Linearized-BGK flow factor on the film conductance (false: first-order slip only)
    exit_factor::Bool = true # Stokes-cell exit factor f0(h) on the face drainage (false: p = 0 exactly at the face)
    seal_shift::Bool = true  # Contact sealing over 0 < delta < 2*ls after nominal contact (false: centred on contact)
    chi_rest::T = 0.75       # Tip-pocket sealing fraction at rest, from the 3-D Stokes solution (0: open tip)

    # Electrode orientation (v5.2) ---------------------------------------------------
    # orient : 1 = compliant (tapered) electrodes anchored to the SUBSTRATE, stiff electrodes on
    #          the shuttle (fabricated devices, from micrographs); 0 = compliant electrodes on the
    #          shuttle (design intent, v5.1). In both, x1 is the shuttle and x2 the compliant tips.
    #          Changes M, beta, B1, the beam forces and the contact reaction; statics are identical.
    orient::Int = 1

    # Gap disorder, model-form correction (v5.2) --------------------------------------
    # Each compliant finger sits offset by a FIXED eps within its channel (fabrication), so the
    # fingers engage at staggered times. eps ~ N(mu_off, sig_off^2) is represented by n_cls
    # equal-population classes, each with its own tip coordinate (states appended after the work
    # integrals). n_cls = 1 and mu_off = 0 recover the synchronized model exactly. The measured
    # 1.85 kHz ring and the 150 Hz output suggest sig_off ~ 90 nm (to be measured by SEM).
    n_cls::Int = 1       # Offset classes (odd keeps a median class as the reference tip, states 3-4)
    mu_off::T = 0.0      # Mean finger offset [m] (systematic: selects the dominant wall)
    sig_off::T = 0.0     # Spread of finger offsets [m] (staggered engagement)
 
    # Hard stop (disabled until the as-fabricated gap is measured) ----------------
    ghs::T = Inf         # Hard-stop engagement position [m]; must satisfy ghs >= gss
    khs::T = 1.0e9       # Hard-stop stiffness [N/m^phs]
    phs::T = 1.5         # Hard-stop exponent
 
    # Contact / boundary parameters -----------------------------------------------
    # h_eff : residual pressed-contact air gap; floor of every local gap. Sets the
    #         electrostatic force at contact (~1/(h_eff + hd)) and the film floor.
    # epsg  : smoothing of the positive part in the gap law h = h_eff + softpos(d).
    # epsw  : wall activation width. epss : stopper activation width.
    # ls    : tip-vent sealing width; chi goes 0 -> 1 over 0 < delta < 2*ls (seal_shift) or -ls < delta < ls (v6).
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
    seal::Bool = true    # Robin tip closure: vent seals as the tip lands (false: no contact sealing; the tip
                         # then vents through the pocket alone, fully open if chi_rest = 0)
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
    kpocket::T = Inf         # Tip-pocket conductance [m^2] (Robin end condition G dp/dy = kappa p)
    ke::T = 0.0              # Electrode tip stiffness, per beam
    k1::T = 0.0              # Linear spring constant
    k3::T = 0.0              # Cubic spring constant
    kss::T = 0.0             # Soft-stopper spring constant
    mtot::T = 0.0            # Suspended mass, m1 + (N/2)*m_beam (both orientations)
    M::Matrix{T} = zeros(T, 2, 2)      # Consistent mass matrix in (x1, x2)
    Minv::Matrix{T} = zeros(T, 2, 2)   # Its inverse
    beta::Vector{T} = zeros(T, 2)      # Base-excitation vector (x2 entry split over the classes)
    # Quadrature grid along the overlap (built once by create_params)
    y::Vector{T} = T[]       # Graded nodes, y = 0 at the mobile tip
    wq::Vector{T} = T[]      # Simpson weights
    vol::Vector{T} = T[]     # Control-volume widths (vented film)
    B1::Vector{T} = T[]      # Weight of x1 in the relative face displacement: 1 - phi | -1 (orient 1)
    B2::Vector{T} = T[]      # Weight of x2 in the relative face displacement, phi
    ecls::Vector{T} = T[]    # Class offsets, median class first
    wcls::Vector{T} = T[]    # Class weights, 1/n_cls
end
 
# ------------------------------ constitutive helpers ------------------------------
# Stable C-infinity positive part, its exact derivative, and the C2 quintic step
@inline softpos(z, e_)  = z >= 0 ? (z + hypot(z, e_))/2 : e_^2/(2*(hypot(z, e_) - z))
@inline dsoftpos(z, e_) = z >= 0 ? (1 + z/hypot(z, e_))/2 :
                                   e_^2/(2*hypot(z, e_)*(hypot(z, e_) - z))
@inline smootherstep(z) = z <= 0 ? 0.0 : z >= 1 ? 1.0 : z^3*(10 - 15*z + 6*z^2)

# Standard-normal quantile (Acklam's rational approximation, |error| < 1.5e-9 vs a reference);
# builds the equal-population offset classes without a Distributions dependency
function normq(q)
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
          1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
          6.680131188771972e+01, -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00)
    if q < 0.02425
        u = sqrt(-2*log(q))
        return (((((c[1]*u + c[2])*u + c[3])*u + c[4])*u + c[5])*u + c[6])/((((d[1]*u + d[2])*u + d[3])*u + d[4])*u + 1)
    elseif q > 1 - 0.02425
        u = sqrt(-2*log(1 - q))
        return -(((((c[1]*u + c[2])*u + c[3])*u + c[4])*u + c[5])*u + c[6])/((((d[1]*u + d[2])*u + d[3])*u + d[4])*u + 1)
    end
    u = q - 0.5; v = u*u
    return (((((a[1]*v + a[2])*v + a[3])*v + a[4])*v + a[5])*v + a[6])*u/(((((b[1]*v + b[2])*v + b[3])*v + b[4])*v + b[5])*v + 1)
end
 
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

# ------------------------------ gas-film closure (v7) ------------------------------
# Linearized-BGK plane Poiseuille flow, diffuse walls (validated against Barichello et al. 2001 to seven digits):
# reduced flow rate G_P(delta), delta = h/lambda, tabulated for 1e-5 <= delta <= 60; ln(delta*G_P) is monotone and
# is interpolated by PCHIP in ln(delta). Table and interpolation copied from the paper's solver (squeeze_film_fd.jl).
const KIN_D = [9.999999999999999e-06, 1.1388973796922415e-05, 1.2970872414698539e-05, 1.4772492605422543e-05, 1.682435311983875e-05, 1.916121168320134e-05, 2.1822653777726372e-05, 2.4853763205383563e-05, 2.8305885790102787e-05, 3.2237499156215916e-05, 3.671520331684516e-05, 4.181484885242285e-05, 4.7622821790251516e-05, 5.423750695046804e-05, 6.177095454692779e-05, 7.03507782745846e-05, 8.012231703623429e-05, 9.125109692743828e-05, 0.0001039256351847022, 0.00011836063359470917, 0.00013480061545972778, 0.00015352406772798543, 0.00017484815845509683, 0.00019913410950852365, 0.00022679331552660548, 0.0002582943127849667, 0.00029417071602020684, 0.0003350302576576041, 0.0003815650825638619, 0.0004345634727140361, 0.0004949232003839766, 0.0005636667360662092, 0.0006419585687254839, 0.0007311249317924354, 0.0008326762690460736, 0.0009483328209484852, 0.0010800537648543815, 0.0012300704027193954, 0.0014009239584940997, 0.0015955086254770129, 0.00181712059283214, 0.002069513881761337, 0.002356963937174706, 0.002684340052077382, 0.0030571878515138653, 0.003481823233316095, 0.003965439356975268, 0.004516228492987621, 0.0051435207967550425, 0.005857942357816869, 0.006671595201705823, 0.007598262293590093, 0.008653641016384118, 0.00985560907835718, 0.011224527354612058, 0.012783584792451562, 0.014559191223196672, 0.01658142473453697, 0.01888454118172828, 0.021507554468560564, 0.024494897427831785, 0.02789717449638785, 0.03177201893475335, 0.036185069112322873, 0.04121108039600719, 0.046935191477298896, 0.05345436658884934, 0.06087903804114902, 0.06933497690324891, 0.07896542351613227, 0.08993351392881116, 0.10242504336003874, 0.11665161349761234, 0.13285421694930283, 0.15130731956462556, 0.1723235097804087, 0.19625879374827782, 0.2235186259414737, 0.25456477739715466, 0.28992315793955825, 0.3301927248894628, 0.37605562917005037, 0.42828877068028764, 0.4877769586793909, 0.5555279001142092, 0.632689269786006, 0.72056815151868, 0.8206531796543067, 0.9346397559443963, 1.0644587690012692, 1.2123093028059746, 1.3806958883422527, 1.5724709293848431, 1.7908830211186217, 2.0396319800873237, 2.3229315176579513, 2.6455806186651625, 3.013044834362333, 3.431548866750505, 3.908182012628031, 4.451018253542416, 5.069253025921794, 5.773358988219298, 6.57526342370561, 7.488550284044557, 8.528690296191936, 9.713303030539645, 11.062455369638311, 12.599001433443439, 14.348969719287528, 16.342004014579885, 18.611865551125124, 21.1970049073607, 24.141213346316686, 27.494364622711448, 31.31325982510913, 35.66258956443913, 40.6160298079796, 46.25748992180995, 52.68253406308962, 60.0]
const KIN_G = [6.853741530642733, 6.780417036395427, 6.70709864625514, 6.633787027380521, 6.560482917199392, 6.48718713044266, 6.41390056682413, 6.340624219433568, 6.267359183879465, 6.194106668252924, 6.120868003978412, 6.047644657601969, 5.974438243602629, 5.901250538281386, 5.8280834948168305, 5.7549392595535895, 5.6818201896005185, 5.608728871830426, 5.535668143342408, 5.462641113487963, 5.389651187520524, 5.316702091962343, 5.243797901751522, 5.170943069249882, 5.0981424551656565, 5.0254013614606245, 4.9527255662769125, 4.8801213609322325, 4.8075955889972075, 4.735155687468187, 4.662809730019578, 4.590566472306283, 4.518435399255351, 4.446426774265302, 4.3745516901962125, 4.302822122003805, 4.231250980833252, 4.159852169350721, 4.088640638050669, 4.017632442232379, 3.946844799300257, 3.8762961459934604, 3.806006195114114, 3.7359959912760687, 3.6662879651664415, 3.5969059857725214, 3.527875410009863, 3.459223129163532, 3.3909776115586183, 3.323168940877811, 3.2558288495786902, 3.1889907469017262, 3.122689741036616, 3.056962655101532, 2.9918480367196945, 2.927386161125368, 2.863619027925489, 2.8005903518623176, 2.7383455481908214, 2.6769317135839907, 2.6163976038294536, 2.556793609967119, 2.498171734955004, 2.440585573428212, 2.3840902976452587, 2.328742653289812, 2.2746009694234677, 2.221725187563767, 2.1701769156012727, 2.1200195130720405, 2.0713182151818263, 2.0241403039436197, 1.9785553358634986, 1.9346354368070882, 1.8924556760352085, 1.8520945329388177, 1.8136344717810189, 1.7771626418067477, 1.742771722478073, 1.7105609363879222, 1.6806372556853115, 1.6531168316761178, 1.628126681746306, 1.605806672958995, 1.5863118477019629, 1.5698151436621617, 1.556510568243731, 1.5466168963235072, 1.5403819699286554, 1.538087688900252, 1.5400557926824134, 1.5466545447246993, 1.5583064421830306, 1.5754970841082905, 1.5987853405392412, 1.628814972327681, 1.6663278568362379, 1.7121789780994436, 1.7673533428030044, 1.8329849880382287, 1.9103782575637671, 2.0010315465136164, 2.1066637580670333, 2.2292437879428237, 2.371023623157809, 2.534574541208286, 2.7228299396244884, 2.9391323933134523, 3.1872874929151367, 3.471626484771623, 3.7970783108401807, 4.169252301257675, 4.594533428667247, 5.080191711431437, 5.634507196738212, 6.266911469811236, 6.98814994844532, 7.810465596892156, 8.747807880198096, 9.816070350514547, 11.033360902484304]
function pchip_slopes1(x, y)
    hx = diff(x); s = diff(y) ./ hx
    d = zeros(length(x)); d[1] = s[1]; d[end] = s[end]
    for k in 2:length(x)-1
        w1 = 2*hx[k] + hx[k-1]; w2 = hx[k] + 2*hx[k-1]
        d[k] = s[k-1]*s[k] > 0 ? (w1 + w2)/(w1/s[k-1] + w2/s[k]) : 0.0
    end
    return d
end
const LKD = log.(KIN_D)
const LKG = log.(KIN_D .* KIN_G)
const DKG = pchip_slopes1(LKD, LKG)
function gp_bgk(dl)
    dl < KIN_D[1] && return 0.35887 + 0.56410*log(1/dl)
    dl > KIN_D[end] && return dl/6 + 1.016191 + 1.0650/dl - 2.1246/dl^2
    lw = log(dl)
    j  = clamp(searchsortedlast(LKD, lw), 1, length(LKD) - 1)
    dx = LKD[j+1] - LKD[j]; t = (lw - LKD[j])/dx
    v  = (2t^3 - 3t^2 + 1)*LKG[j] + (t^3 - 2t^2 + t)*dx*DKG[j] + (-2t^3 + 3t^2)*LKG[j+1] + (t^3 - t^2)*dx*DKG[j+1]
    return exp(v)/dl
end
# Film conductance [m^3]: first-order slip h^2 (h + kp), times the BGK factor Q(delta)/(1 + 6 sigmap/delta) >= 1
function filmG(h, p)
    G = h^2*(h + p.kp)
    p.kinetic || return G
    dl = h/p.lambda
    return G*6*gp_bgk(dl)/dl/(1 + 6*p.sigmap/dl)
end
# Exit factor of the face drainage, f0 = (1 + 2 (c0 + c1 X/(X + x0)) X)^3 with X = h/(2 l_d), l_d = Tf (one open face)
# or Tf/2 (two): the gas turning out of the slot, fitted to the paper's 2-D Stokes cells within 0.21%
const F0_ONE = (0.4207, 0.2449, 1.0015)
const F0_TWO = (0.4066, 0.1279, 0.2733)
function exitf(h, p)
    c0, c1, x0 = p.vent_faces == 2 ? F0_TWO : F0_ONE
    X = h/(2*(p.vent_faces == 2 ? p.Tf/2 : p.Tf))
    return (1 + 2*(c0 + c1*X/(X + x0))*X)^3
end
 
function create_params(p::Params{T}; verbose = true) where T<:Real
    @assert p.N > 0 && iseven(p.N) && p.nsp_par > 0 && p.nsp_ser > 0 && p.nss > 0
    @assert 0 < p.Leff <= p.Lf && p.g0 > 2*p.Tp + p.h_eff && p.panels >= 16
    @assert all(>(0), (p.Tf, p.wt, p.wb, p.ws, p.wss, p.Lsp, p.Lss, p.m1, p.rho, p.E, p.e,
                       p.ep, p.eta, p.h_eff, p.epsg, p.epsw, p.epss, p.ls, p.Rload, p.kw,
                       p.pw, p.phs, p.gss))
    @assert all(>=(0), (p.c, p.c1, p.ce, p.cw, p.gamma3, p.cp, p.Tp, p.lambda, p.sigmap, p.khs))
    @assert p.ghs >= p.gss && p.vent_faces in (0, 1, 2) && p.vent_modes >= 1
    @assert p.orient in (0, 1) && p.n_cls >= 1 && p.sig_off >= 0 && isfinite(p.mu_off)
    @assert 0 <= p.chi_rest <= 1 && (!p.kinetic || p.lambda > 0)
 
    p.nb = div(p.N, 2)
    p.a  = isnan(p.gap_slope) ? (p.wb - p.wt)/p.Lf : p.gap_slope
    @assert isfinite(p.a) && p.a >= 0
 
    # Electrode tip stiffness: unit tip load on the tapered clamped beam (Castigliano)
    p.ke = 1/glquad(s -> (p.Lf - s)^2/bEI(s, p), 0.0, p.Lf)
 
    # Suspension / stopper spring constants
    p.k1  = p.nsp_par/p.nsp_ser*p.E*p.Tf*p.ws^3/p.Lsp^3
    p.k3  = p.gamma3*p.nsp_par/p.nsp_ser^3*0.72*p.E*p.Tf*p.ws/p.Lsp^3
    p.kss = p.nss*p.E*p.Tf*p.wss^3/(4*p.Lss^3)
 
    # Consistent mass matrix. orient 0: the compliant electrodes ride on the shuttle,
    # w = (1 - phi)*x1 + phi*x2. orient 1: they are anchored to the substrate, w = phi*x2, and the
    # stiff electrodes (same planform, mass m_beam) ride on the shuttle.
    mass(i, j) = p.nb*glquad(s -> begin
        phi = bshape(s, p.ke, p); B = (1 - phi, phi)
        p.rho*p.Tf*bwidth(s, p)*B[i]*B[j]
    end, 0.0, p.Lf)
    mbeam = p.rho*p.Tf*p.Lf*(p.wt + p.wb)/2
    if p.orient == 0
        m11 = p.m1 + mass(1, 1); m12 = mass(1, 2); m22 = mass(2, 2)
        p.M = [m11 m12; m12 m22]
        p.beta = p.M*ones(2)
        p.mtot = sum(p.M)
    else
        m11 = p.m1 + p.nb*mbeam; m22 = mass(2, 2)
        p.M = [m11 0.0; 0.0 m22]
        p.beta = [m11, p.nb*glquad(s -> p.rho*p.Tf*bwidth(s, p)*bshape(s, p.ke, p), 0.0, p.Lf)]
        p.mtot = m11
    end
    @assert isposdef(Symmetric(p.M))
    p.Minv = inv(p.M)
 
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
    p.B1 = p.orient == 0 ? 1 .- p.B2 : -ones(length(p.B2))   # stiff face moves with x1 when orient = 1

    # Offset classes: centres of n_cls equal-population bands of N(mu_off, sig_off^2), median first
    nc = p.n_cls
    ec = [p.mu_off + p.sig_off*normq((j - 0.5)/nc) for j in 1:nc]
    p.ecls = ec[sortperm(collect(1:nc); by = j -> (abs(j - (nc + 1)/2), j))]
    p.wcls = fill(1/nc, nc)

    # Tip-pocket conductance, calibrated at rest: sealing fraction chi_rest measured against the film's own
    # lengthwise resistance I0 (gas-film paper: 3-D Stokes solution, one open face; doubled for two faces)
    if p.chi_rest <= 0
        p.kpocket = Inf
    elseif p.chi_rest >= 1
        p.kpocket = 0.0
    else
        hr, _, _ = gapfield(0.0, 0.0, 1.0, p)
        p.kpocket = (1 - p.chi_rest)/(p.chi_rest*dot(p.wq, 1 ./ filmG.(hr, Ref(p))))*(p.vent_faces == 2 ? 2.0 : 1.0)
    end
 
    verbose && report(p)
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
#   d/dy(G dp/dy) + (G/f0) d2p/dzeta2 = 12 eta hdot,   p = 0 on the open device-layer faces.
# v7: G carries the BGK factor (filmG), f0(h) is the exit factor of the face drainage (exitf), the strip tail is
# P = -h_j f0/(G k_n^2), and the tip end is the contact sealing in series with the calibrated pocket.
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
    G    = filmG.(h, Ref(p))                              # slip, times the BGK factor when p.kinetic
    fx   = p.exit_factor ? exitf.(h, Ref(p)) : ones(n)      # exit factor of the face drainage
    I0   = dot(p.wq, 1 ./ G)
    gface = [2*G[i]*G[i+1]/(G[i] + G[i+1])/(y[i+1] - y[i]) for i in 1:n-1]  # face conductances
    # tip end: contact sealing (chi, against the film's own resistance I0) in series with the pocket
    kseal = chi <= 1e-12 ? Inf : chi >= 1 ? 0.0 : (1 - chi)/(chi*I0)
    if isinf(kseal)
        kappa = p.kpocket
    elseif isinf(p.kpocket)
        kappa = kseal
    elseif kseal == 0 || p.kpocket == 0
        kappa = 0.0
    else
        kappa = 1/(1/kseal + 1/p.kpocket)
    end
    open_tip = isinf(kappa)
    H = hcat(h1, h2)
    B = H .* p.vol
    B[n, :] .= 0.0                        # p = 0 at the open far end
    open_tip && (B[1, :] .= 0.0)          # p = 0 at a fully open tip
    V = zeros(2, 2); csum = 0.0
    for k in 0:p.vent_modes-1
        q = 2*k + 1; kn = q*pi/Wv; w = Wv*(16/(q*pi)^2)/2; csum += 1/q^4
        dg = -(G ./ fx .* p.vol) .* kn^2  # leak to the faces, through the exit factor
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
    S = (H .* p.vol)' * (H .* fx ./ G)
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
        tip = p.orient == 0 ? x2 : x2 - x1        # compliant tip relative to the facing stiff face
        delta = r*tip - p.gc                      # nominal tip overlap (contact when > 0)
        chi = !p.seal ? 0.0 : p.seal_shift ? smootherstep(delta/(2*p.ls)) : smootherstep((delta + p.ls)/(2*p.ls))
        a11, a12, a22 = p.vent_faces == 0 ? film_lengthwise(h, h1, h2, chi, p) :
                                            film_vented(h, h1, h2, chi, p)
        d11 += a11; d12 += a12; d22 += a22
    end
    return d11, d12, d22
end
 
# --------------------------------- force functions ---------------------------------
# ALL forces below are generalized forces of the WHOLE array (N/2 beams, both walls) on
# the coordinates x1 (shuttle) and x2 (compliant tips, same frame as x1).
 
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
 
# Electrode bending + tip contact for one offset class (weight w, offset eps).
#   q  = bending of the compliant electrodes: x2 - x1 (orient 0) | x2 (orient 1)
#   tp = tip relative to the facing stiff face: x2 (orient 0) | x2 - x1 (orient 1), plus eps
# Returns Fc1 (beam on x1), Fc2 (beam on x2), Fw1 (contact reaction on x1), Fw (contact on x2),
# the state label and Pw >= 0, the exact contact loss including the zero-force unloading branch.
function collision(x1, x1dot, x2, x2dot, p; eps = 0.0, w = 1.0)
    q   = p.orient == 0 ? x2 - x1 : x2
    tp  = (p.orient == 0 ? x2 : x2 - x1) + eps
    tv  = p.orient == 0 ? x2dot : x2dot - x1dot
    Fc1 = p.orient == 0 ? w*p.nb*p.ke*q : 0.0
    Fw  = 0.0; Pw = 0.0
    for r in (-1.0, 1.0)
        zw = r*tp - p.gc; vr = r*tv
        s  = softpos(zw, p.epsw)
        A  = p.kw*s^p.pw*dsoftpos(zw, p.epsw)
        gate = max(1 + p.cw*vr, 0.0)
        Fw -= w*p.nb*r*A*gate
        Pw += w*p.nb*A*vr*(gate - 1)
    end
    Fw1 = p.orient == 0 ? 0.0 : -Fw               # the shuttle's stiff faces take the reaction
    collision_state = abs(tp) > p.gc ? "contact" : "translational"
    return Fc1, -w*p.nb*p.ke*q, Fw1, Fw, collision_state, Pw
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
 
# States: [x1, x1dot, x2, x2dot, Vout, (7 work integrals), x2_2, x2dot_2, ..., x2_K, x2dot_K].
#   Class 1 (the median offset class) keeps states 3-4, so n_cls = 1 is the v5.1 layout exactly.
#   M*[x1ddot, x2ddot_j] = [f1, f2_j]: arrowhead mass matrix (M11; w_j*M12; w_j*M22), base -beta*a(t)
#   dVout/dt = -Vout/(R*Ct) + ((Vbias - Vout)/Ct)*(dC1*x1dot + sum_j dC2_j*x2dot_j)
# Optional work integrals (states 6-12) are not physics: Wbase, Wbias, ER, Dfilm, Dstruct, Dwall,
#   throughput. Each class is evaluated with the one-class kernels on shifted coordinates:
#   (x1 + eps, x2 + eps) for orient 0 and (x1 - eps, x2) for orient 1.
nextra(p) = 2*(p.n_cls - 1)
function classtip(z, p, j)
    j == 1 && return z[3], z[4]
    b = length(z) - nextra(p) + 2*(j - 2)
    return z[b + 1], z[b + 2]
end
shifted(x1, x2, e, p) = p.orient == 0 ? (x1 + e, x2 + e) : (x1 - e, x2)

function CoupledSystem!(dz, z, p, t, current_acceleration)
    z1, z2, z5 = z[1], z[2], z[5]
    Fext = current_acceleration
    Vc   = p.Vbias - z5
    nc   = p.n_cls
    f1   = spring(z1, p) - p.c1*z2 - p.beta[1]*Fext
    Ctotal = p.cp; dC1 = 0.0; d11 = 0.0; s2 = 0.0
    Pw = 0.0; Ps = p.c1*z2*z2; Pb = -Fext*p.beta[1]*z2
    v2 = zeros(nc); f2 = zeros(nc); d12 = zeros(nc); d22 = zeros(nc)
    for j in 1:nc
        x2, v2[j] = classtip(z, p, j)
        w = p.wcls[j]; e = p.ecls[j]
        xa, xb = shifted(z1, x2, e, p)
        a11, a12, a22 = film(xa, xb, p)
        ct, _, _, c1, c2 = electrostatic(xa, xb, z5, p)
        d11 += w*a11; d12[j] = w*a12; d22[j] = w*a22
        Ctotal += w*(ct - p.cp); dC1 += w*c1; s2 += w*c2*v2[j]
        fc1, fc2, fw1, fw, _, pw = collision(z1, z2, x2, v2[j], p; eps = e, w = w)
        qd  = p.orient == 0 ? v2[j] - z2 : v2[j]
        fb2 = -w*p.nb*p.ce*qd                      # beam damping on the class tips
        f1 += fc1 + fw1 + (p.orient == 0 ? -fb2 : 0.0)
        f2[j] = fc2 + fb2 + fw + 0.5*Vc*Vc*w*c2 - w*p.beta[2]*Fext
        Pw += pw; Ps += w*p.nb*p.ce*qd^2; Pb += -Fext*w*p.beta[2]*v2[j]
    end
    f1 += -d11*z2 + 0.5*Vc*Vc*dC1
    Pf = d11*z2*z2
    for j in 1:nc
        f1    -= d12[j]*v2[j]
        f2[j] -= d12[j]*z2 + d22[j]*v2[j]
        Pf    += 2*d12[j]*z2*v2[j] + d22[j]*v2[j]^2
    end
    # arrowhead solve: the Schur complement M11 - M12^2/M22 does not depend on n_cls
    a1 = (f1 - (p.M[1,2]/p.M[2,2])*sum(f2))/(p.M[1,1] - p.M[1,2]^2/p.M[2,2])
    dz[1] = z2
    dz[2] = a1
    for j in 1:nc
        ix = j == 1 ? 3 : length(z) - nextra(p) + 2*(j - 2) + 1
        dz[ix]     = v2[j]
        dz[ix + 1] = (f2[j] - p.wcls[j]*p.M[1,2]*a1)/(p.wcls[j]*p.M[2,2])
    end
    dz[5] = -z5/(p.Rload*Ctotal) + (Vc/Ctotal)*(dC1*z2 + s2)
    if length(dz) - nextra(p) >= 12
        Pe = p.Vbias*z5/p.Rload
        PR = z5*z5/p.Rload
        dz[6] = Pb; dz[7] = Pe; dz[8] = PR; dz[9] = Pf; dz[10] = Ps; dz[11] = Pw
        dz[12] = abs(Pb) + abs(Pe) + PR + Pf + Ps + Pw
    end
    return nothing
end

# Class-summed forces, capacitance and powers at one state (same algebra as CoupledSystem!).
# Forces "on x2" are the totals on all tips; Fc is the beam force on x1 (zero for orient 1).
function forces(z, p, a)
    z1, z2, z5 = z[1], z[2], z[5]; Vc = p.Vbias - z5
    Fs = spring(z1, p) - p.c1*z2
    Fc = 0.0; Fct = 0.0; Fb = 0.0; Fw = 0.0; Fd1 = 0.0; Fd2 = 0.0; Fe1 = 0.0; Fe2 = 0.0
    Ct = p.cp; Pw = 0.0; Ps = p.c1*z2*z2; Pb = -a*p.beta[1]*z2; Pf = 0.0
    for j in 1:p.n_cls
        x2, v2 = classtip(z, p, j); w = p.wcls[j]; e = p.ecls[j]
        xa, xb = shifted(z1, x2, e, p)
        a11, a12, a22 = film(xa, xb, p)
        ct, _, _, c1, c2 = electrostatic(xa, xb, z5, p)
        fc1, fc2, _, fw, _, pw = collision(z1, z2, x2, v2, p; eps = e, w = w)
        qd = p.orient == 0 ? v2 - z2 : v2
        Fc += fc1; Fct += fc2; Fw += fw; Fb += -w*p.nb*p.ce*qd
        Fd1 += -w*(a11*z2 + a12*v2); Fd2 += -w*(a12*z2 + a22*v2)
        Fe1 += 0.5*Vc*Vc*w*c1; Fe2 += 0.5*Vc*Vc*w*c2
        Ct += w*(ct - p.cp); Pw += pw; Ps += w*p.nb*p.ce*qd^2; Pb += -a*w*p.beta[2]*v2
        Pf += w*(a11*z2*z2 + 2*a12*z2*v2 + a22*v2*v2)
    end
    return (; Fs, Fc, Fct, Fb, Fw, Fd1, Fd2, Fe1, Fe2, Ct, Pb, Pe = p.Vbias*z5/p.Rload,
              PR = z5*z5/p.Rload, Pf, Ps, Pw)
end
 
# Stored energy by reservoir, in base-relative coordinates:
#   Tk  kinetic (consistent mass matrix)      Usp suspension + soft/hard stoppers
#   Ube electrode bending                     Uw  tip-wall contact potential
#   Ue  electrical field energy, Ct*Vc^2/2
# and the exact residual
#   ledger = E - E0 - Wbase - Wbias + ER + Dfilm + Dstruct + Dwall  (= 0 analytically)
function energy_parts(z, p, Ctotal)
    x1, v1, Vout = z[1], z[2], z[5]
    Tk  = 0.5*p.M[1,1]*v1*v1
    Usp = p.k1*x1^2/2 + p.k3*x1^4/4
    Ube = 0.0
    Uw  = 0.0
    for r in (-1.0, 1.0)
        Usp += p.kss*softpos(r*x1 - p.gss, p.epss)^2/2
        isfinite(p.ghs) && (Usp += p.khs*softpos(r*x1 - p.ghs, p.epss)^(p.phs + 1)/(p.phs + 1))
    end
    for j in 1:p.n_cls
        x2, v2 = classtip(z, p, j); w = p.wcls[j]
        q   = p.orient == 0 ? x2 - x1 : x2
        tp  = (p.orient == 0 ? x2 : x2 - x1) + p.ecls[j]
        Tk  += w*(p.M[1,2]*v1*v2 + 0.5*p.M[2,2]*v2*v2)
        Ube += w*p.nb*p.ke*q^2/2
        for r in (-1.0, 1.0)
            Uw += w*p.nb*p.kw*softpos(r*tp - p.gc, p.epsw)^(p.pw + 1)/(p.pw + 1)
        end
    end
    Ue = 0.5*Ctotal*(p.Vbias - Vout)^2
    return (; Tk, Usp, Ube, Uw, Ue, E = Tk + Usp + Ube + Uw + Ue)
end
# Total capacitance, summed over the offset classes
function capacitance(z, p)
    Ct = p.cp
    for j in 1:p.n_cls
        x2, _ = classtip(z, p, j)
        xa, xb = shifted(z[1], x2, p.ecls[j], p)
        Ct += p.wcls[j]*(electrostatic(xa, xb, z[5], p)[1] - p.cp)
    end
    return Ct
end
energy_parts(z, p) = energy_parts(z, p, capacitance(z, p))
energy(z, p) = energy_parts(z, p).E
ledger(z, E0, p) = energy(z, p) - E0 - z[6] - z[7] + z[8] + z[9] + z[10] + z[11]
 
# Derived quantities of the current parameter set, printed for inspection (values only, no
# stored references). The film line includes one independent correctness number: the vented
# solver against the classical rectangular-plate solution (uniform gap, no slip, open edges).
function report(p)
    nbke = p.nb*p.ke
    k11 = p.k1 + nbke; k12 = -nbke
    K  = p.orient == 0 ? [k11 k12; k12 nbke] : [p.k1 0.0; 0.0 nbke]
    fr = sqrt.(eigvals(Symmetric(K), Symmetric(p.M)))./(2*pi)
    Mc = p.orient == 0 ? p.M[1,1] : p.M[1,1] + p.M[2,2]      # tips move with the shuttle in contact
    xc = p.orient == 0 ? (p.gc, p.gc) : (-p.gc, 0.0)          # rigid closure of every gap to contact
    bsum(d) = p.orient == 0 ? d[1] + 2*d[2] + d[3] : d[1]
    d0 = film(0.0, 0.0, p); dc = film(xc[1], xc[2], p)
    Cc, Fe1c, Fe2c, _, _ = electrostatic(xc[1], xc[2], 0.0, p)
    Fsc = p.k1*p.gc + p.k3*p.gc^3
    Fh  = p.orient == 0 ? Fe1c + Fe2c : abs(Fe1c)
    qr = create_params(Params{Float64}(panels = p.panels, gap_slope = 0.0, sigmap = 0.0, seal = false,
                                       vent_faces = 2, orient = 0, kinetic = false, exit_factor = false,
                                       chi_rest = 0.0); verbose = false)
    hr = qr.gc + qr.h_eff
    bt = 1 - 192/pi^5*(qr.Tf/qr.Leff)*sum(tanh(k*pi*qr.Leff/(2*qr.Tf))/k^5 for k in 1:2:199)
    dq = film(0.0, 0.0, qr)
    plate = (dq[1] + 2*dq[2] + dq[3])/(qr.N*qr.eta*qr.Leff*qr.Tf^3/hr^3*bt)
    println("\n--- Orientation: ", p.orient == 0 ? "compliant electrodes on the shuttle (design intent)" :
            "compliant electrodes on the substrate (as built)", " ---")
    println("\n--- Springs ---")
    @printf("ke = %.6f N/m per electrode   k1 = %.6f N/m   k3 = %.4e N/m^3   kss = %.4f N/m\n",
            p.ke, p.k1, p.k3, p.kss)
    println("\n--- Geometry ---")
    @printf("travel to contact gc = %.4f um   gap slope = %.5f   h_eff = %.0f nm   coating gap hd = %.1f nm\n",
            p.gc*1e6, p.a, p.h_eff*1e9, p.hd*1e9)
    println("\n--- Mass matrix (x1, x2) ---")
    @printf("M11 = %.6e   M12 = %.6e   M22 = %.6e   mtot = %.6e kg\n", p.M[1,1], p.M[1,2], p.M[2,2], p.mtot)
    println("\n--- Modes ---")
    @printf("shuttle f1 = %.1f Hz   free electrode mode = %.2f kHz   tips-pinned contact mode = %.0f Hz\n",
            fr[1], fr[2]/1e3, sqrt((p.k1 + nbke)/Mc)/(2*pi))
    @printf("zeta_c of ce in the contact mode = %.3g   shuttle Q of c1 = %s\n",
            p.nb*p.ce/(2*sqrt((p.k1 + nbke)*Mc)), p.c1 > 0 ? @sprintf("%.3g", sqrt(p.k1*p.mtot)/p.c1) : "Inf")
    println("\n--- Film, rigid closure of the whole array (vent_faces = ", p.vent_faces, ", kinetic = ", p.kinetic,
            ", exit factor = ", p.exit_factor, ", sealing ", p.seal_shift ? "after nominal contact" : "centred on contact", ") ---")
    @printf("tip pocket conductance = %.4e m^2 (chi_rest = %.2f)   free path = %.2f nm\n", p.kpocket, p.chi_rest, p.lambda*1e9)
    @printf("b(rest) = %.4e N s/m   b(contact) = %.4e N s/m\n", bsum(d0), bsum(dc))
    @printf("vented film / rectangular-plate solution (uniform gap, %d panels) = %.5f\n", p.panels, plate)
    println("\n--- Electrostatics at nominal contact, Vbias = ", p.Vbias, " V ---")
    @printf("Ct = %.4f pF   Fe on x1 = %.3f uN   Fe on tips = %.3f uN\n", Cc*1e12, Fe1c*1e6, Fe2c*1e6)
    @printf("suspension force at gc = %.2f uN   static hold voltage = %.2f V\n", Fsc*1e6, p.Vbias*sqrt(Fsc/Fh))
    println("\n--- Offset classes: n_cls = ", p.n_cls, ", offsets (nm) = ", round.(p.ecls .* 1e9; digits = 1))
    return nothing
end

# Energy accounting at five representative states for the parameters actually being run: free
# flight, approach, and contact on both walls. dE/dt (central differences of E along the vector
# field) must equal the power balance Pb + Pbias - PR - Pf - Ps - Pw, which holds only if every
# force is the exact gradient of its energy and every loss is booked. No stored references.
# The states are written for orient 0 and mapped to orient 1 with overlap and bending preserved.
function energy_check(p; tol = 1e-4, verbose = true)
    base = [(0.0, 0.0), (p.gc - 1e-6, p.gc - 1e-6), (p.gc + 0.1e-6, p.gc - 2e-9),
            (p.gc + 0.2e-6, p.gc + 5e-9), (-p.gc - 0.1e-6, -p.gc - 2e-9)]
    worst = 0.0
    for (x1d, x2d) in base
        x1, x2, v1, v2 = p.orient == 0 ? (x1d, x2d, 0.003, -0.005) : (-x1d, x2d - x1d, -0.003, -0.008)
        z = [x1, v1, x2, v2, 0.2, zeros(7)...]
        for j in 2:p.n_cls
            append!(z, [x2 + 1e-9*(j - 1), v2*(1 - 0.1*(j - 1))])
        end
        dz = zero(z)
        CoupledSystem!(dz, z, p, 0.0, 2.0)
        iphys = [1:5; 13:length(z)]; g = zeros(length(iphys))
        for (k, j) in enumerate(iphys)
            hj = (j == 5 ? 3.0 : (j in (2, 4) || (j > 12 && iseven(j - 12))) ? 0.02 : p.gc)*1e-6
            zp = copy(z); zm = copy(z); zp[j] += hj; zm[j] -= hj
            g[k] = (energy(zp, p) - energy(zm, p))/(2*hj)
        end
        exact = dz[6] + dz[7] - (dz[8] + dz[9] + dz[10] + dz[11])
        worst = max(worst, abs(dot(g, dz[iphys]) - exact)/max(dz[12], 1e-25))
    end
    verbose && @printf("\nenergy accounting at 5 states (flight, approach, contact on both walls): worst |dE/dt - power balance| / throughput = %.2e  %s\n",
                       worst, worst <= tol ? "(OK)" : "(TOO LARGE: a force and its energy disagree)")
    return worst
end

# Initialize a default Params instance and calculate dependent parameters
p = Params{Float64}()
p = create_params(p; verbose = false)
 
end # module AnalyticalModel
 
import .AnalyticalModel
 
# --------------------------------------- External Force ------------------------------------
 
# Sine Wave External Force
f = 20.0        # Frequency (Hz)
alpha = 1.0    # Applied acceleration constant (g). 4.95 -> panel (d); 2.7 -> panel (e).
                # Quasi-static contact threshold at 3 V is between 2.0 and 2.1
g = 9.80665     # Gravitational constant (m/s^2)
A = alpha*g
n_ramp = 4      # Ramp-up duration in drive cycles (C1 cosine ramp, zero end slopes)
ramp(t) = 0.5*(1 - cos(pi*min(t*f/n_ramp, 1.0)))
Fext_sine = t -> A*ramp(t)*sin(2*pi*f*t)

# ---------------------------------- Experimental Input (optional) ----------------------------------
# use_experiment = true drives the model with a measured base acceleration instead of the sine, and
# compares model and measured output when the output file is present. Files are the two-column
# LabVIEW text exports kept in ./data beside this script:
#   acceleration: "Samples" or "Time (s)", then "Acceleration (g)"   (sample index -> exp_fs_accel)
#   output      : "Time (s)", then "Voltage (V)"                     (optional: enables the comparison)
# Both records are assumed to start at the same instant (exp_align = :none); :xcorr shifts the measured
# output by the best-matching lag. The lag is printed either way. In experiment mode alpha (peak |a|
# in the window, g) and f (drive frequency) come from the record, so the ramp, the last-two-cycle
# windows, the animation and the run folder work unchanged. Runtime scales with the window length:
# the full 9.9 s record of accelT1 is ~1000 drive cycles; (2.0, 3.5) brackets the jump.
use_experiment   = false
exp_dir          = joinpath(@__DIR__, "data")
exp_accel_file   = joinpath(exp_dir, "accelT1.tmp.txt")
exp_output_file  = joinpath(exp_dir, "voltageT1.tmp.txt")   # "" or a missing file: model only
exp_fs_accel     = 100e3          # accel sampling rate when the file stores sample indices [S/s]
exp_window       = (1.0, 4.0)     # part of the record to simulate [s]
exp_lowpass      = 2000.0         # zero-phase low-pass on the measured acceleration [Hz]; 0 disables
exp_accel_scale  = 1.0            # sign / sensitivity correction on the acceleration
exp_output_scale = 1.0            # readout gain correction on the measured voltage
exp_align        = :none          # :none (common start) or :xcorr (best-matching lag)

# Two-column text export (header line, whitespace separated, CRLF tolerated)
function read_two_column(path)
    lines = readlines(path)
    c1 = Float64[]; c2 = Float64[]
    sizehint!(c1, length(lines)); sizehint!(c2, length(lines))
    for ln in @view lines[2:end]
        sp = split(strip(ln))
        length(sp) >= 2 || continue
        push!(c1, parse(Float64, sp[1])); push!(c2, parse(Float64, sp[2]))
    end
    return c1, c2, strip(lines[1])
end

# Linear interpolation on a uniform grid (t0, dt); holds the end values outside the record
function uniform_interp(y, t0, dt, t)
    x = (t - t0)/dt
    x <= 0 && return y[1]
    i = floor(Int, x) + 1
    i >= length(y) && return y[end]
    w = x - (i - 1)
    return (1 - w)*y[i] + w*y[i + 1]
end

# Zero-phase 2nd-order Butterworth low-pass (biquad run forward and backward); fc <= 0 disables
function lowpass_zero_phase(x, fs, fc)
    (fc <= 0 || fc >= fs/2) && return copy(x)
    K = tan(pi*fc/fs); nrm = 1/(1 + sqrt(2)*K + K^2)
    b0 = K^2*nrm; b1 = 2*b0; b2 = b0; a1 = 2*(K^2 - 1)*nrm; a2 = (1 - sqrt(2)*K + K^2)*nrm
    function onepass(u)
        y = similar(u); y[1] = u[1]; y[2] = u[2]
        for n in 3:length(u)
            y[n] = b0*u[n] + b1*u[n-1] + b2*u[n-2] - a1*y[n-1] - a2*y[n-2]
        end
        return y
    end
    return reverse(onepass(reverse(onepass(x))))
end

# Drive frequency from upward zero crossings with hysteresis (robust to noise and amplitude sweeps)
function drive_frequency(t, a)
    h = 0.2*maximum(abs, a); armed = false; tc = Float64[]
    for i in 2:length(a)
        a[i] < -h && (armed = true)
        if armed && a[i-1] < 0 && a[i] >= 0
            push!(tc, t[i-1] + (t[i] - t[i-1])*(-a[i-1])/(a[i] - a[i-1])); armed = false
        end
    end
    length(tc) >= 3 || error("drive_frequency: fewer than three drive cycles in the acceleration window")
    return 1/sort(diff(tc))[div(length(tc), 2)]
end

Fext_exp = t -> 0.0
exp_has_output = false
if use_experiment
    c1col, acol, ahdr = read_two_column(exp_accel_file)
    ta_exp    = occursin("Time", ahdr) ? c1col : c1col ./ exp_fs_accel
    dta       = ta_exp[2] - ta_exp[1]
    a_exp_raw = exp_accel_scale .* acol                                   # measured, g
    a_exp     = lowpass_zero_phase(a_exp_raw, 1/dta, exp_lowpass)         # model input before the ramp, g
    exp_has_output = !isempty(exp_output_file) && isfile(exp_output_file)
    if exp_has_output
        tv_exp, vcol, _ = read_two_column(exp_output_file)
        v_exp = exp_output_scale .* vcol                                   # measured output, V
        dtv   = tv_exp[2] - tv_exp[1]
    end
    t_rec_end = exp_has_output ? min(ta_exp[end], tv_exp[end]) : ta_exp[end]
    exp_t0 = max(exp_window[1], ta_exp[1]); exp_t1 = min(exp_window[2], t_rec_end)
    exp_t1 > exp_t0 || error("exp_window lies outside the record")
    exp_T  = exp_t1 - exp_t0
    iw     = findall(t -> exp_t0 <= t <= exp_t1, ta_exp)
    f_exp  = drive_frequency(ta_exp[iw], a_exp[iw])
    f      = round(f_exp; digits = 1)                     # used by the ramp, the windows and the animation
    alpha  = round(maximum(abs, a_exp[iw]); digits = 2)    # peak |a| in the window (g), for the run folder
    exp_tag = replace(replace(basename(exp_accel_file), r"(\.tmp)?\.(txt|csv)$" => ""), "accel" => "")
    isempty(exp_tag) && (exp_tag = "data")
    Fext_exp = let a_in = a_exp, t0a = ta_exp[1], dt_in = dta, toff = exp_t0
        t -> g*ramp(t)*uniform_interp(a_in, t0a, dt_in, t + toff)
    end
    @printf("experimental input: %s (%d samples, %.1f kS/s), window %.3f-%.3f s, drive %.4f Hz, peak %.3f g%s\n",
            basename(exp_accel_file), length(a_exp), 1e-3/dta, exp_t0, exp_t1, f_exp, alpha,
            exp_has_output ? string("; output: ", basename(exp_output_file)) : "; no output file (model only)")
end
 
# ------------------------------------- Set Input Force ------------------------------------
 
# Set to `true` to use sine forcing, `false` for a near-contact displaced IC
# (free evolution: one contact episode probe, no external force)
use_sine = true
Fext_input = use_experiment ? Fext_exp : use_sine ? Fext_sine : (t -> 0.0)
 
# ------------------------------------ Initialize Parameters --------------------------------
 
p_new = deepcopy(AnalyticalModel.p)

# Electrode orientation and gap disorder (v5.2). Defaults: fabricated orientation, no disorder.
#   p_new.orient = 0                           # design intent: compliant electrodes on the shuttle
#   p_new.n_cls = 5; p_new.sig_off = 90e-9     # staggered engagement; class tips appended to z
#   AnalyticalModel.create_params(p_new; verbose = false)
 
# To change parameters, set the fields and REBUILD the derived quantities, e.g.
#   p_new.Vbias = 5.0;  AnalyticalModel.create_params(p_new; verbose = false)
# Submitted corrected model (lengthwise film, no structural loss):
#   p_new.vent_faces = 0; p_new.c1 = 0.0; p_new.ce = 0.0; p_new.panels = 512
#   AnalyticalModel.create_params(p_new; verbose = false)
AnalyticalModel.create_params(p_new; verbose = true)
 
run_energy_check = true  # energy accounting at five representative states (about a second)
run_energy_check && AnalyticalModel.energy_check(p_new)
 
# Initial conditions
if use_sine || use_experiment
    x10, x10dot, x20, x20dot = 0.0, 0.0, 0.0, 0.0
else
    # Near-contact probe IC (episode-scale evaluation without forcing): the median tip starts
    # 30 nm short of its stiff face and closes at 4 mm/s in either orientation
    if p_new.orient == 0
        x20, x20dot = p_new.gc - 30e-9, 4e-3
        x10, x10dot = p_new.gc + 0.15e-6, 4e-3
    else
        x20, x20dot = 0.0, 0.0
        x10, x10dot = -(p_new.gc - 30e-9), -4e-3
    end
end
 
# Equilibrated start: q = Vbias*Ct  <=>  Vout = 0 exactly.
# use_ledger appends the seven work integrals (states 6-12): an exact energy acceptance
# test for every run, at roughly twice the Jacobian cost. Set false for speed.
use_ledger = true
z0 = [x10, x10dot, x20, x20dot, 0.0]
use_ledger && (z0 = vcat(z0, zeros(7)))
z0 = vcat(z0, repeat([x20, x20dot], p_new.n_cls - 1))   # class tips 2..n_cls, after the work states
n_cycles = 10
tspan  = use_experiment ? (0.0, exp_T) : use_sine ? (0.0, n_cycles/f) : (0.0, 600e-6)
 
# State scaling (corrected-model numerics). The finite-difference Jacobian perturbs each
# state by ~1.5e-8*max(|z|, 1): in SI metres that is 15 nm, far wider than the contact
# physics, so the solver integrates O(1) scaled states and abstol is one scalar.
xs = p_new.gc                                  # displacement scale
vs = xs*sqrt(p_new.k1/p_new.mtot)              # velocity scale
Es = p_new.k1*xs^2                             # energy scale
zscale = [xs, vs, xs, vs, max(abs(p_new.Vbias), 1.0)]
use_ledger && (zscale = vcat(zscale, fill(Es, 7)))
zscale = vcat(zscale, repeat([xs, vs], p_new.n_cls - 1))
abstol = 1e-10                                 # on the scaled states
reltol = 1e-7
dtmax  = (use_sine || use_experiment) ? 2e-5 : 1e-6                # brackets the ~100 us contact events
 
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

# Gas-film linearity: the film stays linear while the tip closes or opens slower than about 0.08 m/s
# (peak film pressure about 129 kPa per m/s at contact, one open face)
let Z = reduce(hcat, sol.u) .* zscale, q = p_new
    tp  = q.orient == 0 ? Z[3, :] : Z[3, :] .- Z[1, :]          # median tip relative to its stiff face
    tpd = q.orient == 0 ? Z[4, :] : Z[4, :] .- Z[2, :]
    d   = abs.(tp) .- q.gc
    vin  = maximum([abs(tpd[i-1]) for i in 2:length(d) if d[i-1] < 0 && d[i] >= 0]; init = 0.0)
    vout = maximum([abs(tpd[i])   for i in 2:length(d) if d[i-1] >= 0 && d[i] < 0]; init = 0.0)
    println("gas-film linearity: max tip speed at impact ", round(1e3*vin; digits = 2), " mm/s, at release ",
            round(1e3*vout; digits = 2), " mm/s (limit 80 mm/s; peak film pressure ~",
            round(100*129e3*max(vin, vout)/101325; digits = 2), "% of 1 atm)")
end
 
println(">>> collision_model version v7 (gas-film paper closure: kinetic = ", p_new.kinetic, ", exit factor = ",
        p_new.exit_factor, ", pocket chi_rest = ", p_new.chi_rest, "; orient = ", p_new.orient,
        ", n_cls = ", p_new.n_cls, ", sig_off = ", p_new.sig_off, ", vent_faces = ",
        p_new.vent_faces, ", c1 = ", p_new.c1, ", ce = ", p_new.ce, ") <<<")
println("Type of sol.u: ", typeof(sol.u))
println("Size of sol.u: ", size(sol.u))
println("Solver status: ", sol.retcode)
println("Solver stats:  ", sol.stats)
 
# Energy acceptance test: max |ledger| / throughput should be < 1e-5 (submitted model: ~1e-9..1e-7)
led_rel = NaN
if use_ledger
    E0     = AnalyticalModel.energy(z0, p_new)
    stride = max(1, div(length(sol.u), 20000))
    maxres = maximum(abs(AnalyticalModel.ledger(u .* zscale, E0, p_new)) for u in sol.u[1:stride:end])
    zend   = sol.u[end] .* zscale
    led_rel = maxres/max(zend[12], eps()*Es)
    @printf("energy ledger over the run: max |E - E0 - W_in + losses| / throughput = %.3e  %s\n", led_rel,
            led_rel < 1e-5 ? "(OK)" : "(LARGE: tighten reltol/abstol or dtmax before trusting the run)")
    @printf("load energy ER(t_end) = %.6e J\n", zend[8])
    if use_sine && n_cycles > 4
        ER4 = (sol(tspan[2] - 4/f) .* zscale)[8]
        @printf("last-four-cycle mean load power = %.6e W   (submitted model at 4.95 g: 9.8853e-13)\n",
                (zend[8] - ER4)/(4/f))
    end
end
 
use_ledger || println("Run-level energy check skipped: set use_ledger = true to carry the work integrals.")

# ------------------------------------ Run Folder & Saving ------------------------------------
# Every figure (PDF), the output data (XLSX) and the animation are written to RUN_Xg_YV_ZHz beside
# this file (X = alpha, Y = Vbias, Z = f; whole numbers print without decimals, e.g. RUN_2g_3V_200Hz).
numtag(x) = (r = round(Float64(x); digits = 4); isinteger(r) ? string(Int(r)) : string(r))
run_dir = joinpath(@__DIR__, string(use_experiment ? string("RUN_EXP_", exp_tag, "_") : "RUN_", numtag(alpha), "g_", numtag(p_new.Vbias), "V_", numtag(f), "Hz"))
mkpath(run_dir)
println("Saving figures, data and animation to: ", run_dir)

# Display a figure and save it as a PDF in the run folder
function showsave(pl, name)
    display(pl)
    savefig(pl, joinpath(run_dir, string(name, ".pdf")))
    return pl
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
    Ct = zeros(n); Fs = zeros(n); Fc = zeros(n); Fct = zeros(n); Fw = zeros(n); Fb = zeros(n)
    Fe1 = zeros(n); Fe2 = zeros(n); Fd1 = zeros(n); Fd2 = zeros(n)
    Tk = zeros(n); Usp = zeros(n); Ube = zeros(n); Uw = zeros(n); Ue = zeros(n); E = zeros(n)
    Pb = zeros(n); Pe = zeros(n); PR = zeros(n); Pf = zeros(n); Ps = zeros(n); Pw = zeros(n)
    for i in 1:n
        z  = U[:, i]
        ai = Fext(tg[i])
        F  = AM.forces(z, p, ai)                      # class-summed, same algebra as the ODE
        Fs[i] = F.Fs; Fc[i] = F.Fc; Fct[i] = F.Fct; Fw[i] = F.Fw; Fb[i] = F.Fb
        Fd1[i] = F.Fd1; Fd2[i] = F.Fd2; Fe1[i] = F.Fe1; Fe2[i] = F.Fe2; Ct[i] = F.Ct
        en = AM.energy_parts(z, p, F.Ct)
        Tk[i] = en.Tk; Usp[i] = en.Usp; Ube[i] = en.Ube; Uw[i] = en.Uw; Ue[i] = en.Ue; E[i] = en.E
        # instantaneous powers: dE/dt = Pb + Pe - PR - Pf - Ps - Pw
        Pb[i] = F.Pb                                  # base excitation (signed)
        Pe[i] = F.Pe                                  # bias source (signed)
        PR[i] = F.PR                                  # load resistor
        Pf[i] = F.Pf                                  # squeeze film
        Ps[i] = F.Ps                                  # structure: shuttle c1 + beam ce
        Pw[i] = F.Pw                                  # tip contact (Hunt-Crossley)
    end
    V    = U[5,:]                                  # Vout is the state
    Q    = [(p.Vbias - V[i])*Ct[i] for i in 1:n]   # charge as observable
    th   = p.orient                                # 1 when the shuttle carries the stiff faces
    tau  = U[3,:] .- th .* U[1,:] .+ p.ecls[1]      # median-class tip relative to its stiff face
    taud = U[4,:] .- th .* U[2,:]
    q    = U[3,:] .- (1 - th) .* U[1,:]            # bending of the compliant electrodes
    pen  = abs.(tau) .- p.gc                       # nominal tip overlap
    htip = [p.h_eff + AM.softpos(p.gc - abs(tau[i]), p.epsg) for i in 1:n]   # near-face tip air gap
    ae   = [Fext(t) for t in tg]
    led  = size(U, 1) >= 12                        # work-integral states present?
    col(k) = led ? U[k,:] : fill(NaN, n)
    return (; t = tg, x1 = U[1,:], x1dot = U[2,:], x2 = U[3,:], x2dot = U[4,:], tau, taudot = taud, q,
              Q, V, Ct, pen, htip, Fs, Fc, Fct, Fb, Fw, Fd1, Fd2, Fe1, Fe2, ae,
              Tk, Usp, Ube, Uw, Ue, E, Pb, Pe, PR, Pf, Ps, Pw,
              Wb = col(6), Wbias = col(7), ER = col(8), Df = col(9), Ds = col(10),
              Dw = col(11), thr = col(12), led)
end
 
# Uniform 10 us overview grid for the state and force plots
Wover = sample_window(sol, p_new, sol.t[1], sol.t[end]; dt = max(1e-5, (sol.t[end] - sol.t[1])/2e5))   # <= 2e5 samples
to = Wover.t
 
p3  = plot(to, Wover.x1,    xlabel = "Time (s)", ylabel = "x1 (m)",     title = "Shuttle Mass Displacement (x1)", label = "");    showsave(p3, "01_x1_displacement")
p4  = plot(to, Wover.x1dot, xlabel = "Time (s)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot)", label = "");    showsave(p4, "02_x1_velocity")
p5  = plot(to, Wover.tau,   xlabel = "Time (s)", ylabel = "tau (m)",    title = "Compliant Tip Relative to its Stiff Face (tau)", label = "")
hline!(p5, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); showsave(p5, "03_tip_relative_tau")
p6  = plot(to, Wover.taudot, xlabel = "Time (s)", ylabel = "taudot (m/s)", title = "Tip Closing Velocity (taudot)", label = ""); showsave(p6, "04_tip_closing_velocity")
p7  = plot(to, Wover.Q,     xlabel = "Time (s)", ylabel = "Q (C)",      title = "Charge (observable)", label = "");               showsave(p7, "05_charge")
p8  = plot(to, Wover.V,     xlabel = "Time (s)", ylabel = "Vout (V)",   title = "Output Voltage (state)", label = "");            showsave(p8, "06_output_voltage")
 
# Diagnostics: penetration and total capacitance
p8b = plot(to, Wover.pen .* 1e9, xlabel = "Time (s)", ylabel = "|tau|-gc (nm)",
           title = "Penetration (contact when > 0)", label = "")
hline!(p8b, [0.0]; ls = :dash, lc = :gray, label = ""); showsave(p8b, "07_penetration")
p8c = plot(to, Wover.Ct .* 1e12, xlabel = "Time (s)", ylabel = "Ctotal (pF)",
           title = "Total Capacitance", label = ""); showsave(p8c, "08_total_capacitance")
 
p9   = plot(to, Wover.Fs, xlabel = "Time (s)", ylabel = "Fs (N)", title = "Suspension + Stopper Force on x1 (incl. -c1*x1dot)", label = ""); showsave(p9, "09_suspension_force")
p10  = plot(to, Wover.Fc, xlabel = "Time (s)", ylabel = "Fc (N)", title = "Beam Force on x1 (zero when orient = 1)", label = ""); showsave(p10, "10_beam_force_x1")
p10a = plot(to, Wover.Fb, xlabel = "Time (s)", ylabel = "Fb (N)", title = "Beam Damping Force on the Tips", label = ""); showsave(p10a, "11_beam_damping_force")
p10b = plot(to, Wover.Fw, xlabel = "Time (s)", ylabel = "Fw (N)", title = "Tip Contact Force on the Tips (reaction on x1 when orient = 1)", label = ""); showsave(p10b, "12_tip_contact_force")
p11  = plot(to, [Wover.Fd1 Wover.Fd2], xlabel = "Time (s)", ylabel = "Fd (N)", title = "Squeeze-Film Force (projected Reynolds film)",
            label = ["on x1" "on x2"], legend = :topright); showsave(p11, "13_squeeze_film_force")
p12  = plot(to, [Wover.Fe1 Wover.Fe2], xlabel = "Time (s)", ylabel = "Fe (N)", title = "Electrostatic Force (attractive)",
            label = ["on x1" "on x2"], legend = :topright); showsave(p12, "14_electrostatic_force")
p13  = plot(to, Wover.ae, xlabel = "Time (s)", ylabel = "a_ext (m/s^2)", title = "Applied Base Acceleration", label = ""); showsave(p13, "15_base_acceleration")
 
# ==================== (2) LAST-TWO-CYCLE TWIN OF EVERY STATE / FORCE PLOT ====================
Tdrive = 1/f
Wzoom  = sample_window(sol, p_new, sol.t[end] - 2*Tdrive, sol.t[end]; dt = 2e-6)
# NOTE: 2 us sampling resolves the chatter structure and the ~100 us contacts; the alias-free
# view of a single impact is the 0.1 us close-up in section (4).
tzs = Wzoom.t .* 1e3
 
p3z  = plot(tzs, Wzoom.x1,    xlabel = "t (ms)", ylabel = "x1 (m)",      title = "Shuttle Mass Displacement (x1) - last 2 cycles", label = "")
hline!(p3z, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); showsave(p3z, "01_x1_displacement_last2cycles")
p4z  = plot(tzs, Wzoom.x1dot, xlabel = "t (ms)", ylabel = "x1dot (m/s)", title = "Shuttle Mass Velocity (x1dot) - last 2 cycles", label = ""); showsave(p4z, "02_x1_velocity_last2cycles")
p5z  = plot(tzs, Wzoom.tau,   xlabel = "t (ms)", ylabel = "tau (m)",     title = "Compliant Tip Relative to its Stiff Face - last 2 cycles", label = "")
hline!(p5z, [p_new.gc, -p_new.gc]; ls = :dash, lc = :gray, label = ""); showsave(p5z, "03_tip_relative_tau_last2cycles")
p6z  = plot(tzs, Wzoom.taudot, xlabel = "t (ms)", ylabel = "taudot (m/s)", title = "Tip Closing Velocity - last 2 cycles", label = ""); showsave(p6z, "04_tip_closing_velocity_last2cycles")
p7z  = plot(tzs, Wzoom.Q,     xlabel = "t (ms)", ylabel = "Q (C)",       title = "Charge (observable) - last 2 cycles", label = ""); showsave(p7z, "05_charge_last2cycles")
p8z  = plot(tzs, Wzoom.V,     xlabel = "t (ms)", ylabel = "Vout (V)",    title = "Output Voltage (state) - last 2 cycles", label = ""); showsave(p8z, "06_output_voltage_last2cycles")
p8bz = plot(tzs, Wzoom.pen .* 1e9, xlabel = "t (ms)", ylabel = "|tau|-gc (nm)", title = "Penetration (contact when > 0) - last 2 cycles", label = "")
hline!(p8bz, [0.0]; ls = :dash, lc = :gray, label = ""); showsave(p8bz, "07_penetration_last2cycles")
p8cz = plot(tzs, Wzoom.Ct .* 1e12, xlabel = "t (ms)", ylabel = "Ctotal (pF)", title = "Total Capacitance - last 2 cycles", label = ""); showsave(p8cz, "08_total_capacitance_last2cycles")
p9z   = plot(tzs, Wzoom.Fs, xlabel = "t (ms)", ylabel = "Fs (N)", title = "Suspension + Stopper Force on x1 - last 2 cycles", label = ""); showsave(p9z, "09_suspension_force_last2cycles")
p10z  = plot(tzs, Wzoom.Fc, xlabel = "t (ms)", ylabel = "Fc (N)", title = "Beam Force on x1 - last 2 cycles", label = ""); showsave(p10z, "10_beam_force_x1_last2cycles")
p10az = plot(tzs, Wzoom.Fb, xlabel = "t (ms)", ylabel = "Fb (N)", title = "Beam Damping Force on the Tips - last 2 cycles", label = ""); showsave(p10az, "11_beam_damping_force_last2cycles")
p10bz = plot(tzs, Wzoom.Fw, xlabel = "t (ms)", ylabel = "Fw (N)", title = "Tip Contact (Hunt-Crossley wall) Force - last 2 cycles", label = ""); showsave(p10bz, "12_tip_contact_force_last2cycles")
p11z  = plot(tzs, hcat(Wzoom.Fd1, Wzoom.Fd2), xlabel = "t (ms)", ylabel = "Fd (N)", title = "Squeeze-Film Force - last 2 cycles",
             label = ["on x1" "on x2"], legend = :topright); showsave(p11z, "13_squeeze_film_force_last2cycles")
p12z  = plot(tzs, hcat(Wzoom.Fe1, Wzoom.Fe2), xlabel = "t (ms)", ylabel = "Fe (N)", title = "Electrostatic Force - last 2 cycles",
             label = ["on x1" "on x2"], legend = :topright); showsave(p12z, "14_electrostatic_force_last2cycles")
p13z  = plot(tzs, Wzoom.ae, xlabel = "t (ms)", ylabel = "a_ext (m/s^2)", title = "Applied Base Acceleration - last 2 cycles", label = ""); showsave(p13z, "15_base_acceleration_last2cycles")
 
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
    showsave(e1, "20_stored_energy_full"); showsave(e2, "21_work_and_losses_full"); showsave(e3, "22_conservation_check_full"); showsave(e4, "23_power_flows_full")
    e1z, e2z, e3z, e4z = energy_plots(Wzoom, tzs, "t (ms)", "last 2 cycles")
    showsave(e1z, "20_stored_energy_last2cycles"); showsave(e2z, "21_work_and_losses_last2cycles"); showsave(e3z, "22_conservation_check_last2cycles"); showsave(e4z, "23_power_flows_last2cycles")
 
    # Energy budget over the last two cycles: inputs = sinks (+ change of stored energy)
    dlt(v) = v[end] - v[1]
    budget = [dlt(Wzoom.Wb), dlt(Wzoom.Wbias), -dlt(Wzoom.E), dlt(Wzoom.ER), dlt(Wzoom.Df), dlt(Wzoom.Ds), dlt(Wzoom.Dw)]
    e5 = bar(["W base", "W bias", "-dE stored", "E load", "D film", "D struct", "D wall"], budget .* 1e12;
             legend = false, ylabel = "Energy over the last 2 cycles (pJ)", xrotation = 30,
             title = "Energy Budget: first three bars = last four")
    showsave(e5, "24_energy_budget_last2cycles")
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
thS   = p_new.orient                                   # 1 when the shuttle carries the stiff faces
tauS  = Zs[3,:] .- thS .* Zs[1,:] .+ p_new.ecls[1]    # median-class tip relative to its stiff face
taudS = Zs[4,:] .- thS .* Zs[2,:]
qS    = Zs[3,:] .- (1 - thS) .* Zs[1,:]              # bending of the compliant electrodes
qdS   = Zs[4,:] .- (1 - thS) .* Zs[2,:]
 
ph1  = plot(Zs[1,:] .* 1e6, Zs[2,:] .* 1e3, xlabel = "x1 (um)", ylabel = "x1dot (mm/s)",
            title = "Shuttle Phase Plane - full trajectory", label = "", lw = 0.5)
vline!(ph1, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); showsave(ph1, "30_shuttle_phase_full")
ph1z = plot(Zs[1,izz] .* 1e6, Zs[2,izz] .* 1e3, xlabel = "x1 (um)", ylabel = "x1dot (mm/s)",
            title = "Shuttle Phase Plane - last 2 cycles", label = "", lw = 0.8)
vline!(ph1z, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); showsave(ph1z, "30_shuttle_phase_last2cycles")
 
ph2  = plot(tauS .* 1e6, taudS .* 1e3, xlabel = "tau (um)", ylabel = "taudot (mm/s)",
            title = "Tip Phase Plane - full trajectory", label = "", lw = 0.5)
vline!(ph2, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); showsave(ph2, "31_tip_phase_full")
ph2z = plot(tauS[izz] .* 1e6, taudS[izz] .* 1e3, xlabel = "tau (um)", ylabel = "taudot (mm/s)",
            title = "Tip Phase Plane - last 2 cycles", label = "", lw = 0.8)
vline!(ph2z, [gcu, -gcu]; ls = :dash, lc = :gray, label = ""); showsave(ph2z, "31_tip_phase_last2cycles")
 
# Electrode bending plane: loops appear only while the tips are on a wall (contact mode)
ph3  = plot(qS .* 1e9, qdS .* 1e3, xlabel = "bending q (nm)",
            ylabel = "qdot (mm/s)", title = "Electrode Bending Plane - full trajectory", label = "", lw = 0.5); showsave(ph3, "32_bending_plane_full")
ph3z = plot(qS[izz] .* 1e9, qdS[izz] .* 1e3, xlabel = "bending q (nm)",
            ylabel = "qdot (mm/s)", title = "Electrode Bending Plane - last 2 cycles", label = "", lw = 0.8); showsave(ph3z, "32_bending_plane_last2cycles")
 
# Near-wall portrait: nominal overlap vs wall-normal tip velocity, both walls folded together
nearw(idx) = [i for i in idx if abs(tauS[i]) - p_new.gc > -200e-9]
inear = nearw(1:size(Zs, 2)); inearz = nearw(izz)
ph4  = scatter((abs.(tauS[inear]) .- p_new.gc) .* 1e9, sign.(tauS[inear]) .* taudS[inear] .* 1e3,
               xlabel = "|tau| - gc (nm)", ylabel = "wall-normal tip velocity (mm/s)", ms = 1.2, msw = 0,
               title = "Phase Portrait at the Contact Boundary - full trajectory", label = "")
vline!(ph4, [0.0]; ls = :dash, lc = :gray, label = ""); showsave(ph4, "33_contact_boundary_portrait_full")
ph4z = scatter((abs.(tauS[inearz]) .- p_new.gc) .* 1e9, sign.(tauS[inearz]) .* taudS[inearz] .* 1e3,
               xlabel = "|tau| - gc (nm)", ylabel = "wall-normal tip velocity (mm/s)", ms = 1.5, msw = 0,
               title = "Phase Portrait at the Contact Boundary - last 2 cycles", label = "")
vline!(ph4z, [0.0]; ls = :dash, lc = :gray, label = ""); showsave(ph4z, "33_contact_boundary_portrait_last2cycles")
 
# Electrical plane: the area enclosed per cycle is the energy delivered to the load
ph5  = plot(Wover.Q .* 1e12, Wover.V .* 1e3, xlabel = "Q (pC)", ylabel = "Vout (mV)",
            title = "Electrical Plane (Q, Vout) - full trajectory", label = "", lw = 0.5); showsave(ph5, "34_electrical_plane_full")
ph5z = plot(Wzoom.Q .* 1e12, Wzoom.V .* 1e3, xlabel = "Q (pC)", ylabel = "Vout (mV)",
            title = "Electrical Plane (Q, Vout) - last 2 cycles", label = "", lw = 0.8); showsave(ph5z, "34_electrical_plane_last2cycles")
 
# ========================== (4) LABELLED COLLISION CLOSE-UPS =============================
# Same-side contact sequences from the accepted steps: entries closer together than `gap`
# belong to one sequence (one wall, one half cycle); a change of wall always starts a new one.
function contact_sequences(tt, x2, x1dot, gc; gap = 8e-3)
    dls  = abs.(x2) .- gc
    ient = [i for i in 2:length(dls) if dls[i-1] < 0 && dls[i] >= 0]
    iext = [i for i in 2:length(dls) if dls[i-1] >= 0 && dls[i] < 0]
    seqs = NamedTuple[]
    isempty(ient) && return seqs, ient, iext
    k0 = 1
    for k in 2:length(ient)+1
        if k > length(ient) || tt[ient[k]] - tt[ient[k-1]] > gap || sign(x2[ient[k]]) != sign(x2[ient[k-1]])
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
    r    = sign(Wd.tau[argmax(abs.(Wd.tau))])
    dl   = r .* Wd.tau .- p.gc                      # nominal overlap delta; contact when > 0
    sg   = 1 - 2*p.orient                           # shuttle approach sign (-1 when it carries the stiff faces)
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
    s1 = (r .* sg .* Wd.x1 .- p.gc) .* 1e9; s2 = dl .* 1e9
    pa = plot(tu, hcat(s1, s2); label = ["shuttle approach - gc" "tip  r*tau - gc"], ylabel = "nm", legend = :bottomright,
              title = "Shuttle and tip relative to the contact boundary")
    hlabel!(pa, 0.0, "x = gc: nominal contact (tip air gap = h_eff)", span(s1)); marks!(pa)
 
    # (b) the physical tip air gap on a log axis, with the three lengths that matter there
    hn = Wd.htip .* 1e9
    pb = plot(tu, hn; yscale = :log10, label = "", ylabel = "tip air gap (nm)",
              title = "Gap closure: h(tip) = h_eff + softpos(gc - r*tau)")
    hline!(pb, [p.h_eff, p.h_eff + p.ls, p.hd] .* 1e9; ls = :dash, lc = :gray, label = "")
    annotate!(pb, xR, p.h_eff*1e9*1.12, text(string("h_eff = ", nm(p.h_eff), " nm: residual gap (floor)"), 7, :right))
    annotate!(pb, xL, (p.h_eff + p.ls)*1e9*1.12, text(string("h_eff + ls = ", nm(p.h_eff + p.ls), " nm: tip vent starts to seal"), 7, :left))
    annotate!(pb, xR, p.hd*1e9*1.30, text(string("hd = 2Tp/ep = ", nm(p.hd), " nm: coating-equivalent gap"), 7, :right))
    marks!(pb)
 
    # (c) nominal overlap with the sealing window
    top = 1.4*max(maximum(dl)*1e9, 1.0)
    pc = plot(tu, dl .* 1e9; label = "", ylabel = "delta = r*tau - gc (nm)", ylims = (-3*p.ls*1e9, top),
              title = "Nominal overlap and the tip-vent sealing window")
    hlabel!(pc, 0.0, "delta = 0", top; side = :right)
    hlabel!(pc, p.ls*1e9, string("+ls = ", nm(p.ls), " nm: vent sealed (chi = 1)"), top)
    hlabel!(pc, -p.ls*1e9, string("-ls: vent open (chi = 0)"), top)
    marks!(pc)
 
    # (d) electrode bending
    bn = r .* Wd.q .* 1e9
    pd = plot(tu, bn; label = "", ylabel = "r*q (nm)", title = "Bending of the compliant electrodes")
    hline!(pd, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pd)
 
    # (e) wall-normal velocities
    pe = plot(tu, hcat(r .* sg .* Wd.x1dot, r .* Wd.taudot) .* 1e3; label = ["shuttle approach" "tip  r*taudot"], ylabel = "mm/s",
              legend = :topright, title = "Velocities normal to the wall")
    hline!(pe, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pe)
 
    # (f) forces on the tip coordinate x2, projected on the wall normal
    Ft = hcat(r .* Wd.Fct, r .* Wd.Fw, r .* Wd.Fe2, r .* Wd.Fd2, r .* Wd.Fb) .* 1e6
    pf = plot(tu, Ft; label = ["bending" "wall contact" "electrostatic" "squeeze film" "beam damping"], ylabel = "uN",
              legend = :bottomright, title = "Forces on the tip coordinate (+ = into the wall)")
    hline!(pf, [0.0]; ls = :dash, lc = :gray, label = ""); marks!(pf)
 
    # (g) output voltage
    pg = plot(tu, Wd.V .* 1e3; label = "", ylabel = "Vout (mV)", xlabel = xlab, title = "Output voltage")
    marks!(pg)
 
    # (h) phase portrait at the contact boundary
    ph = plot(dl .* 1e9, r .* Wd.taudot .* 1e3; label = "", xlabel = "delta = r*tau - gc (nm)", ylabel = "r*taudot (mm/s)",
              xlims = (-3*p.ls*1e9, top), title = "Phase portrait at the contact boundary")
    vline!(ph, [0.0]; ls = :dash, lc = :gray, label = "")
 
    # metrics of the FIRST contact in the window (shuttle-level restitution, contact time, peaks)
    met = (; side = r, n_contacts = length(tin), t_in = NaN, t_contact = NaN, v_in = NaN, v_out = NaN,
             e_shuttle = NaN, overlap_max = NaN, bending_max = NaN, Fw_max = NaN)
    jo = isempty(tin) ? nothing : findfirst(>(tin[1]), tout)
    if jo !== nothing
        i1 = findfirst(>=(tin[1]), Wd.t); i2 = findfirst(>=(tout[jo]), Wd.t)
        vin = r*sg*Wd.x1dot[i1]; vout = r*sg*Wd.x1dot[i2]
        met = (; side = r, n_contacts = length(tin), t_in = tin[1], t_contact = tout[jo] - tin[1], v_in = vin, v_out = vout,
                 e_shuttle = -vout/vin, overlap_max = maximum(dl[i1:i2]), bending_max = maximum(abs, bn[i1:i2])*1e-9,
                 Fw_max = maximum(abs, Wd.Fw[i1:i2]))
        annotate!(pe, xL, minimum(r .* sg .* Wd.x1dot)*1e3 + 0.08*span(r .* sg .* Wd.x1dot)*1e3,
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
 
seqs, ient, iext = contact_sequences(sol.t, tauS, taudS, p_new.gc; gap = min(8e-3, 0.25/f))
if isempty(seqs)
    println("\nNo contact this run (closest approach ",
            round(-maximum(abs.(tauS) .- p_new.gc)*1e9; digits = 1), " nm). Collision close-ups skipped.")
else
    println("\n================ CONTACT SEQUENCES (one line per wall visit) ================")
    println("   start (s)    wall   contacts   hardest tip entry (mm/s)       duration (ms)")
    for s in seqs
        @printf("   %9.5f     %+d     %5d        %8.3f                     %8.3f\n",
                s.t0, Int(s.side), s.n, s.vmax*1e3, (s.t1 - s.t0)*1e3)
    end
    @printf("wall +1: %d visits, %d contacts     wall -1: %d visits, %d contacts\n",
            count(s -> s.side > 0, seqs), sum([s.n for s in seqs if s.side > 0]; init = 0),
            count(s -> s.side < 0, seqs), sum([s.n for s in seqs if s.side < 0]; init = 0))

    # Chatter-interval ratio: successive impact-to-impact intervals within one wall visit, the
    # calibration target for ce (bench at 20 Hz: successive ratios 0.68-0.75, about 0.73)
    let tin = crossings(sol.t, abs.(tauS) .- p_new.gc)[1], sg = sign.(tauS)
        wall(tq) = sg[clamp(searchsortedfirst(sol.t, tq), 1, length(sol.t))]
        vis = Vector{Vector{Float64}}(); cur = Float64[]
        for tk in tin
            if !isempty(cur) && (wall(tk) != wall(cur[end]) || tk - cur[end] > min(8e-3, 0.25/f))
                push!(vis, cur); cur = Float64[]
            end
            push!(cur, tk)
        end
        isempty(cur) || push!(vis, cur)
        rr = Float64[]
        for v in vis
            length(v) >= 3 || continue
            dv = diff(v); append!(rr, dv[2:end] ./ dv[1:end-1])
        end
        if isempty(rr)
            println("chatter-interval ratio: no wall visit with three or more impacts")
        else
            rs = sort(rr); nr = length(rs)
            @printf("chatter-interval ratio: %d ratios from %d visits, median %.3f, interquartile %.3f-%.3f (bench ~0.73)\n",
                    nr, count(v -> length(v) >= 3, vis), rs[cld(nr, 2)], rs[cld(nr, 4)], rs[cld(3*nr, 4)])
        end
    end
 
    recent = [s for s in seqs if s.t0 >= sol.t[end] - 2*Tdrive]
    isempty(recent) && (recent = seqs)
 
    # (4a) the hardest single impact of the last two cycles, 0.1 us sampling
    sA  = recent[argmax([s.vmax for s in recent])]
    jA  = findfirst(>(sA.ihard), iext)
    tA1 = jA === nothing ? sol.t[end] : sol.t[iext[jA]]
    Wimp = sample_window(sol, p_new, sol.t[sA.ihard] - 40e-6, tA1 + 80e-6; dt = 1e-7)
    figA, mA = collision_figure(Wimp, p_new; unit = :us,
                                tag = string("hardest impact, t = ", round(sol.t[sA.ihard]; digits = 5), " s, wall ", Int(sA.side)))
    showsave(figA, "40_collision_closeup_hardest_impact")
 
    # (4b) the longest chatter sequence of the last two cycles: approach - chatter - dwell - release
    sB   = recent[argmax([s.n for s in recent])]
    Wep  = sample_window(sol, p_new, sB.t0 - 0.5e-3, sB.t1 + 0.5e-3; dt = max(1e-7, (sB.t1 - sB.t0 + 1e-3)/40000))
    figB, mB = collision_figure(Wep, p_new; unit = :ms,
                                tag = string("longest sequence, ", sB.n, " contacts from t = ", round(sB.t0; digits = 5), " s, wall ", Int(sB.side)))
    showsave(figB, "41_collision_closeup_longest_sequence")
 
    println("\n================ CLOSE-UP METRICS (first contact of each window) ================")
    for (nm_, m) in (("hardest impact  ", mA), ("longest sequence", mB))
        @printf("%s: wall %+d, %d contacts in window; first contact %.1f us; shuttle %.3f -> %.3f mm/s (e = %.3f); peak overlap %.2f nm; peak bending %.1f nm; peak wall force %.1f uN\n",
                nm_, Int(m.side), m.n_contacts, m.t_contact*1e6, m.v_in*1e3, m.v_out*1e3, m.e_shuttle,
                m.overlap_max*1e9, m.bending_max*1e9, m.Fw_max*1e6)
    end
end




# ============================ (5) MODEL VS EXPERIMENT (measured input) ============================
# Runs with use_experiment = true and a measured output file. Long records are drawn as min/max
# envelopes so impact spikes stay visible; the last two drive cycles are shown at full resolution,
# and the output peak-to-peak of every drive cycle is plotted against the drive amplitude in it.
function envelope(t, y, nb)
    n = length(t); nb = min(nb, n)
    edges = round.(Int, range(1, n + 1; length = nb + 1))
    tc = zeros(nb); lo = zeros(nb); hi = zeros(nb)
    for k in 1:nb
        r = edges[k]:(edges[k+1] - 1)
        tc[k] = t[r[1]]; lo[k] = minimum(@view y[r]); hi[k] = maximum(@view y[r])
    end
    return tc, lo, hi
end

function compare_experiment(W, gacc, t0, fdrive, nramp, a_raw, ta0, dta, v_meas, tv0, dtv, align, title)
    tm   = W.t; td = tm .+ t0; dtg = tm[2] - tm[1]
    a_in = W.ae ./ gacc                                  # model input (g): filtered, ramped
    a_ms = [uniform_interp(a_raw, ta0, dta, t) for t in td]
    vmod = W.V .* 1e3                                    # mV
    vat(lag) = [uniform_interp(v_meas, tv0, dtv, t + lag)*1e3 for t in td]
    post = tm .>= nramp/fdrive                           # after the input ramp
    vm0  = vmod[post] .- sum(vmod[post])/count(post)
    nl   = ceil(Int, 1/(fdrive*dtg)); cbest = -Inf; lbest = 0.0
    for k in -nl:nl                                      # lags within one drive period
        ve = vat(k*dtg)[post]; ve .-= sum(ve)/length(ve)
        c  = sum(vm0 .* ve)/sqrt(sum(abs2, vm0)*sum(abs2, ve) + 1e-300)
        if c > cbest
            cbest = c; lbest = k*dtg
        end
    end
    vexp = vat(align == :xcorr ? lbest : 0.0)
    il = findall(>=(tm[end] - min(0.1, 0.25*(tm[end] - tm[1]))), tm)
    ws(x) = (maximum(x), minimum(x), maximum(x) - minimum(x), sqrt(sum(abs2, x .- sum(x)/length(x))/length(x)))
    sm = ws(vmod[il]); se = ws(vexp[il])
    @printf("\nmodel vs experiment over the last %.3f s (model | measured): peak %.2f | %.2f mV, dip %.2f | %.2f mV, p-p %.2f | %.2f mV, rms %.2f | %.2f mV\n",
            tm[end] - tm[il[1]], sm[1], se[1], sm[2], se[2], sm[3], se[3], sm[4], se[4])
    @printf("best-matching lag of the measured output: %+.3f ms (correlation %.3f), %s\n", lbest*1e3, cbest,
            align == :xcorr ? "applied" : "not applied (exp_align = :none, common start assumed)")
    # output peak-to-peak of every drive cycle against the drive amplitude in that cycle
    ci = floor.(Int, tm .* fdrive) .+ 1; nc = maximum(ci)
    lom = fill(Inf, nc); him = fill(-Inf, nc); loe = fill(Inf, nc); hie = fill(-Inf, nc); am = zeros(nc)
    for i in eachindex(tm)
        c = ci[i]
        lom[c] = min(lom[c], vmod[i]); him[c] = max(him[c], vmod[i])
        loe[c] = min(loe[c], vexp[i]); hie[c] = max(hie[c], vexp[i])
        am[c]  = max(am[c], abs(a_in[i]))
    end
    cy = 1:max(nc - 1, 1)                                # drop the partial last cycle
    ok = "#000000"; om = "#D55E00"; oa = "#0072B2"; og = "#999999"
    te, alo, ahi = envelope(td, a_ms, 4000); _, ilo, ihi = envelope(td, a_in, 4000)
    _, elo, ehi = envelope(td, vexp, 4000);  _, mlo, mhi = envelope(td, vmod, 4000)
    pa = plot(te, alo; fillrange = ahi, c = og, fillcolor = og, fillalpha = 0.6, lw = 0, label = "measured (raw)",
              ylabel = "a (g)", title = "Base acceleration")
    plot!(pa, te, ilo; fillrange = ihi, c = oa, fillcolor = oa, fillalpha = 0.4, lw = 0, label = "model input (filtered)")
    pb = plot(te, elo; fillrange = ehi, c = ok, fillcolor = ok, fillalpha = 0.45, lw = 0, label = "measured",
              ylabel = "Vout (mV)", title = "Output voltage")
    plot!(pb, te, mlo; fillrange = mhi, c = om, fillcolor = om, fillalpha = 0.45, lw = 0, label = "model")
    iz = findall(>=(tm[end] - 2/fdrive), tm)
    pc = plot(td[iz] .* 1e3, vexp[iz]; c = ok, label = "measured", xlabel = "record time (ms)", ylabel = "Vout (mV)",
              title = "Last two drive cycles")
    plot!(pc, td[iz] .* 1e3, vmod[iz]; c = om, label = "model")
    pd = scatter(am[cy], (hie .- loe)[cy]; c = ok, ms = 2.5, msw = 0, label = "measured",
                 xlabel = "drive amplitude in the cycle (g)", ylabel = "output p-p (mV)", title = "Output per drive cycle")
    scatter!(pd, am[cy], (him .- lom)[cy]; c = om, ms = 2.5, msw = 0, label = "model")
    fx = plot(pa, pb, pc, pd; layout = (4, 1), size = (1000, 1500), plot_title = title,
              left_margin = 5*Plots.mm, titlefontsize = 10)
    cols  = Any[tm, td, a_ms, a_in, vmod, vexp]
    names = ["t model (s)", "t record (s)", "a measured (g)", "a model input (g)", "Vout model (mV)", "Vout measured (mV)"]
    return fx, cols, names
end

exp_cols = Any[]; exp_names = String[]
if use_experiment && exp_has_output
    fx, exp_cols, exp_names = compare_experiment(Wover, g, exp_t0, f_exp, n_ramp, a_exp_raw, ta_exp[1], dta,
                                                 v_exp, tv_exp[1], dtv, exp_align,
                                                 string("Model vs experiment: ", basename(exp_accel_file),
                                                        ", Vbias = ", p_new.Vbias, " V"))
    showsave(fx, "50_model_vs_experiment")
elseif use_experiment
    exp_cols  = Any[Wover.t, Wover.t .+ exp_t0, Wover.ae ./ g]
    exp_names = ["t model (s)", "t record (s)", "a model input (g)"]
end

# ------------------------------------- Output Data (XLSX) -------------------------------------
# Model states and forces on the uniform 10 us overview grid (Wover), one sheet each, plus the
# stored energies / powers and the run settings. Units are in the column headers.
Uall   = Array(sol(Wover.t)) .* zscale             # every state: class tips and work integrals included
snames = ["x1 (m)", "x1dot (m/s)", "x2 (m)", "x2dot (m/s)", "Vout (V)"]
use_ledger && append!(snames, ["W_base (J)", "W_bias (J)", "E_load (J)", "D_film (J)", "D_struct (J)",
                               "D_wall (J)", "throughput (J)"])
for j in 2:p_new.n_cls
    append!(snames, [string("x2 class ", j, " (m)"), string("x2dot class ", j, " (m/s)")])
end
state_cols = Any[Wover.t]
for k in 1:size(Uall, 1)
    push!(state_cols, Uall[k, :])
end
push!(state_cols, Wover.tau, Wover.taudot, Wover.q, Wover.Q, Wover.Ct)
state_names = vcat(["t (s)"], snames, ["tau, median tip rel. stiff face (m)", "taudot (m/s)", "q, bending (m)",
                                       "Q (C)", "Ctotal (F)"])
force_cols  = Any[Wover.t, Wover.ae, Wover.Fs, Wover.Fc, Wover.Fct, Wover.Fb, Wover.Fw,
                  Wover.Fd1, Wover.Fd2, Wover.Fe1, Wover.Fe2]
force_names = ["t (s)", "a_ext (m/s^2)", "Fs suspension + stoppers on x1 (N)", "Fc beam on x1 (N)",
               "Fct beam on tips (N)", "Fb beam damping on tips (N)", "Fw contact on tips (N)",
               "Fd1 film on x1 (N)", "Fd2 film on tips (N)", "Fe1 electrostatic on x1 (N)",
               "Fe2 electrostatic on tips (N)"]
energy_cols  = Any[Wover.t, Wover.Tk, Wover.Usp, Wover.Ube, Wover.Uw, Wover.Ue, Wover.E,
                   Wover.Pb, Wover.Pe, Wover.PR, Wover.Pf, Wover.Ps, Wover.Pw]
energy_names = ["t (s)", "kinetic (J)", "suspension + stoppers (J)", "bending (J)", "contact (J)",
                "electrical (J)", "total E (J)", "P base (W)", "P bias (W)", "P load (W)", "P film (W)",
                "P struct (W)", "P contact (W)"]
info = [("alpha (g)", alpha), ("f (Hz)", f), ("Vbias (V)", p_new.Vbias), ("Rload (Ohm)", p_new.Rload),
        ("orient", p_new.orient), ("n_cls", p_new.n_cls), ("mu_off (m)", p_new.mu_off),
        ("sig_off (m)", p_new.sig_off), ("k1 (N/m)", p_new.k1), ("ke (N/m)", p_new.ke), ("c1 (N s/m)", p_new.c1),
        ("ce (N s/m)", p_new.ce), ("vent_faces", p_new.vent_faces), ("n_cycles", n_cycles),
        ("use_sine", use_sine), ("input", use_experiment ? basename(exp_accel_file) : use_sine ? "sine" : "free probe"),
        ("solver retcode", sol.retcode), ("energy residual / throughput", led_rel)]
info_cols = Any[[first(x) for x in info], [string(last(x)) for x in info]]
xlsx_path = joinpath(run_dir, "output_data.xlsx")
XLSX.openxlsx(xlsx_path, mode = "w") do xf
    sh = xf[1]
    XLSX.rename!(sh, "States")
    XLSX.writetable!(sh, state_cols, state_names)
    XLSX.writetable!(XLSX.addsheet!(xf, "Forces"), force_cols, force_names)
    XLSX.writetable!(XLSX.addsheet!(xf, "Energy"), energy_cols, energy_names)
    XLSX.writetable!(XLSX.addsheet!(xf, "Run info"), info_cols, ["setting", "value"])
    isempty(exp_cols) || XLSX.writetable!(XLSX.addsheet!(xf, "Experiment"), exp_cols, exp_names)
end
println("Saved output data: ", xlsx_path)


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------

using Plots, Printf

"""
    animate_electrode_contact(sol, p, zscale; frequency, outfile, kwargs...)

Animate the last two forcing cycles of the existing SCALED solution. No new
simulation is run. Root = state 1, tip = state 3; both are restored to SI units.
Beam shape: w(s,t) = (1-orient)*(1-phi(s))*x1(t) + phi(s)*x2(t), using AnalyticalModel.bshape.
With orient = 1 the compliant beam stands on the fixed base and the stiff electrodes ride
on the shuttle (top bar); contact uses tau, the tip relative to its stiff face.

`window=(ta,tb)` selects an explicit interval in seconds, useful for slow motion.
`seconds=10, fps=30` makes 300 frames. The movie is an overview, not a count of
contact crossings; very fast rebounds need a narrow window / more frames.
`frame_times` optionally supplies increasing physical times for variable-speed
playback. Use `animate_electrode_cycle` below for both wall visits of one cycle.
The diagram compresses the vertical dimension; horizontal geometry is in µm.
Orange means nominal contact, r*tau >= gc. The tiny residual film is reported
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
    th = p.orient; bend = (U[3,:] .- (1-th).*U[1,:]).*1e9
    TAU = (U[3,:] .- th.*U[1,:] .+ p.ecls[1]).*1e6
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
    xlim = 1.12max(maximum(abs, X1), maximum(abs, TAU), p.gc*1e6)
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
            tau=z[3]-th*z[1]+p.ecls[1]; taud=z[4]-th*z[2]
            w=(1-th).*(1 .- phi).*x1 .+ phi.*x2
            sx=th*x1; xb0=(1-th)*x1
            onleft=-tau>=p.gc; onright=tau>=p.gc
            side=tau<0 ? -1 : tau>0 ? 1 : (taud<0 ? -1 : 1)
            wall=side<0 ? "LEFT" : "RIGHT"
            j=isempty(visits) ? nothing : argmin([abs(t-(v.entry+v.exit)/2) for v in visits])
            visit=j===nothing ? nothing : visits[j]
            state=onleft ? "LEFT CONTACT" : onright ? "RIGHT CONTACT" :
                visit!==nothing && visit.entry<=t<=visit.exit ? "BETWEEN REBOUNDS - "*(visit.side<0 ? "LEFT" : "RIGHT") :
                side*taud>0 ? "APPROACH "*wall : "RELEASE / AWAY FROM "*wall
            hm=p.h_eff+AnalyticalModel.softpos(p.gc+tau,p.epsg)
            hp=p.h_eff+AnalyticalModel.softpos(p.gc-tau,p.epsg)
            a=plot(; common..., xlims=(-extent,extent), ylims=(-48,135),
                framestyle=:none, ticks=false, legend=false,
                title=@sprintf("2D electrode schematic  |  t = %.6f s",t))
            plot!(a,box(-extent+3+sx,extent-3+sx,118,130); c=th==0 ? gray : teal,lc=:black,label="")
            annotate!(a,sx,124,text(th==0 ? "Fixed support" : "Shuttle  x1(t)  (stiff electrodes)",10,th==0 ? :black : :white))
            for (r,active) in ((-1,onleft),(1,onright))
                ix=r.*[inner_bottom,inner_top,inner_top] .+ sx
                ox=r.*[2fixed_center-inner_bottom,2fixed_center-inner_top,2fixed_center-inner_top] .+ sx
                sy=[ybottom,100.,118.]
                plot!(a,Shape(vcat(ix,reverse(ox)),vcat(sy,reverse(sy)));
                    c=active ? orange : gray,lc=:black,label="")
            end
            plot!(a,box(xb0-shuttle_half,xb0+shuttle_half,-9,0);c=th==0 ? teal : gray,lc=:black,label="")
            plot!(a,Shape(vcat(w.-halfwidth,reverse(w.+halfwidth)),vcat(yy,reverse(yy)));
                c=teal,lc=:black,label="")
            plot!(a,w,yy;c=:white,ls=:dash,lw=1,label="")
            scatter!(a,[xb0,x2],[0,100];c=red,ms=4,msc=red,label="")
            annotate!(a,xb0,-4.5,text(th==0 ? "Shuttle  x1(t)" : "Substrate (fixed)",10,th==0 ? :white : :black))
            annotate!(a,x2,105,text("x2(t)",10,:black))
            # Tip-gap brackets; the numerical readout resolves nanometre clearances.
            for (xa,xb) in ((-inner_top+sx,x2-p.wb*1e6/2),(x2+p.wb*1e6/2,inner_top+sx))
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
            annotate!(a,0,-38,text(@sprintf("Nominal overlap |tau| − gc = %.2f nm",(abs(tau)-p.gc)*1e9),9,:black))
            annotate!(a,0,-45,text("Horizontal scale retained; vertical dimension compressed",8,:gray35))
            b=plot(tx,hcat(X1,TAU); common...,c=[teal red],label=["Shuttle x1" "Tip rel. face tau"],
                xlims=(0,(tb-ta)*1e3),ylims=(-xlim,xlim),ylabel="Displacement (µm)",
                xlabel="",title=variable_speed ? "One cycle; playback slower around wall visits" : "Motion and nominal contact limits",
                legend=:outertop,legend_column=2)
            for v in visits
                vspan!(b,([v.lo,v.hi].-ta).*1e3;c=orange,alpha=0.09,label="")
            end
            hline!(b,[-p.gc,p.gc].*1e6;c=:gray,ls=:dash,lw=1,label="")
            vline!(b,[tm];c=:black,lw=1,label="")
            scatter!(b,[tm,tm],[x1,tau*1e6];c=[teal,red],ms=4,msc=:white,label="")
            c=plot(tx,bend;common...,c=teal,label="",xlims=(0,(tb-ta)*1e3),ylims=blim,
                xlabel="Time from window start (ms)",ylabel="bending q (nm)",title="Bending of the compliant electrodes",legend=false)
            vline!(c,[tm];c=:black,lw=1,label="")
            scatter!(c,[tm],[(z[3]-(1-th)*z[1])*1e9];c=red,ms=4,msc=:white,label="")
            if visit!==nothing
                # Wall-relative coordinates resolve nanometre contact on the same
                # trajectory. This close-up switches between the two wall visits.
                v=visit; sel=findall(q->v.lo<=q<=v.hi,tt)
                td=(tt[sel].-v.entry).*1e3
                d1=(v.side.*(1-2th).*U[1,sel].-p.gc).*1e9
                d2=(v.side.*(U[3,sel].-th.*U[1,sel].+p.ecls[1]).-p.gc).*1e9
                d=plot(td,hcat(d1,d2);common...,c=[teal red],label=["Shuttle" "Tip"],
                    legend=false,xlims=((v.lo-v.entry)*1e3,(v.hi-v.entry)*1e3),
                    xlabel="Time from first contact in this visit (ms)",ylabel="approach - gc (nm)",
                    title="Visit $j/$(length(visits)): "*(v.side<0 ? "left" : "right")*" approach and release")
                hline!(d,[0.];c=:gray,ls=:dash,label="")
                vline!(d,([v.entry,v.exit].-v.entry).*1e3;c=orange,ls=:dot,label="")
                if v.lo<=t<=v.hi
                    vline!(d,[(t-v.entry)*1e3];c=:black,label="")
                    scatter!(d,fill((t-v.entry)*1e3,2),(v.side.*[(1-2th)*z[1],tau].-p.gc).*1e9;
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
    tauof(t)=(s=sol(t).*zscale; s[3]-p.orient*s[1]+p.ecls[1])   # median tip relative to its stiff face
    x2=[tauof(t) for t in tt]; delta=abs.(x2).-p.gc
    hit=delta.>=0; visits=NamedTuple[]; peaks=Float64[]
    # Refine observed brackets in the dense solution; this is for movie timing,
    # not an independently converged physical event-detection calculation.
    function crossing(a,b,r)
        ga=r*tauof(a)-p.gc
        for _ in 1:35
            m=(a+b)/2; gm=r*tauof(m)-p.gc
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
    outfile=joinpath(run_dir, "electrode_one_cycle_two_hits.gif"))