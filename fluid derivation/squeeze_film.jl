#=
squeeze_film_fd.jl — frequency-domain squeeze-film impedance of a laterally drained film
(thickness-mode theory with first-order slip, gas inertia (Womersley), compressibility (squeeze number)
and a lateral exit factor f(h, ω) calibrated from Stokes cells).  Fluid solution only.

Companion to "Squeeze-film impedance of laterally drained films closing to line contact".
Time dependence exp(iωt).  The generalized force on the structure is Q = -Z(ω) v, where v holds the
velocity amplitudes of the generalized coordinates; Re Z is damping, Im Z > 0 is mass-like, Im Z < 0 spring-like.

Usage (Julia ≥ 1.9):
    ] add Plots
    include("squeeze_film_fd.jl")
    selftest()          # compares with reference values from the verified Python reference implementation
    run_case_study()    # prints the case-study numbers quoted in the paper
    make_plots("figs")  # writes the figures (PDF and PNG)

Model summary (paper, Eq. numbers in brackets refer to the paper):
  field equation   ∂y(G ∂y P) + ∂ζ((G/f) ∂ζ P) − iω (12η h/(n pa)) P = 12η Σ_j h_j v_j
  mobility         G(h, ω) = (3h³/λ²)[1 − tanh λ / (λ(1 + (2b/h) λ tanh λ))],  λ = (h/2)√(iω/ν)
  thickness modes  k_n = (2n+1)π/(2ℓd),  w_n = 8W/((2n+1)²π²)   (ℓd = W one open face, W/2 two)
  impedance        Z_ij = 12η Σ_n w_n ∫ h_i X_n^[j] dy + strip tail,   L_n X_n^[j] = h_j
SI units throughout.
=#

using LinearAlgebra, Printf
using Plots

# ----------------------------------------------------------------------------- fluid and mobility
Base.@kwdef struct Fluid
    eta::Float64 = 1.849e-5      # dynamic viscosity, Pa s
    pa::Float64 = 101325.0       # ambient pressure, Pa
    rho::Float64 = 1.204         # density, kg/m^3
    lam::Float64 = 70e-9         # mean free path, m
    sig_p::Float64 = 1.016       # slip coefficient
    npoly::Float64 = 1.0         # polytropic index (1 = isothermal)
end
nu(fl::Fluid) = fl.eta / fl.rho
bslip(fl::Fluid) = fl.sig_p * fl.lam            # Navier slip length b
kp(fl::Fluid) = 6 * bslip(fl)                    # k_p = 6 σ_p λ, G0 = h²(h + k_p)

"Womersley mobility with first-order slip (complex). q = −G ∇p/(12η)."
function G_mob(h::Real, w::Real, fl::Fluid; slip::Bool=true, wom::Bool=true)
    b = slip ? bslip(fl) : 0.0
    G0 = h^2 * (h + 6b)
    (!wom || w == 0) && return complex(G0)
    lam = 0.5 * h * sqrt(im * w / nu(fl))
    if abs(lam) < 1e-3
        return G0 - im * w / nu(fl) * h^5 * (0.1 + b / h + 3b^2 / h^2)
    end
    t = tanh(lam)
    beta = 2b / h
    return 3h^3 / lam^2 * (1 - t / (lam * (1 + beta * lam * t)))
end

# ----------------------------------------------------------------------------- configuration
Base.@kwdef struct Config
    W::Float64 = 25e-6           # lateral width of the film (device-layer thickness)
    faces::Symbol = :one         # :one (one open face, other blocked) or :two (both open)
    slip::Bool = true
    wom::Bool = true             # gas inertia across the gap
    comp::Bool = true            # compressibility
    exit::Bool = true            # lateral exit factor f(h, ω)
    exit_omega::Bool = true      # frequency dependence of f
    N::Int = 24                  # exact thickness modes
    M::Int = 2000                # explicit tail terms (local balance)
end
ld(cfg::Config) = cfg.faces == :one ? cfg.W : cfg.W / 2

# ----------------------------------------------------------------------------- lateral exit factor
# f0(h) = (1 + 2 c(X) X)^3,  c(X) = c0 + c1 X/(X + x0),  X = h/(2 ℓd)   (fits to 2-D Stokes cells, < 0.21 %)
const F0_COEF = Dict(:one => (0.4207, 0.2449, 1.0015), :two => (0.4066, 0.1279, 0.2733))
# frequency correction f/f0 − 1 relative to the Womersley strip, one open face, comb geometry of the case study
const H_TAB = [0.5, 1.0, 2.0, 4.0, 8.0, 13.76, 18.426666666666666, 25.42666666666667, 32.42666666666667]   # μm
const W_TAB = 2π .* [0.0, 300.0, 1000.0, 3000.0, 10000.0]                                                # rad/s
const DF_RE = [
    0.000000e+00  1.029372e-09  3.529997e-09  2.074261e-08  1.997910e-07;
    0.000000e+00  1.925826e-08  4.894635e-08  2.043413e-07  1.692664e-06;
    0.000000e+00  3.199716e-07  7.003564e-07  2.219530e-06  1.511431e-05;
    0.000000e+00  4.828726e-06  9.860482e-06  2.611164e-05  1.460743e-04;
    0.000000e+00  6.333553e-05  1.260565e-04  3.076258e-04  1.512227e-03;
    0.000000e+00  4.165356e-04  8.279783e-04  1.998578e-03  9.354829e-03;
    0.000000e+00  1.081416e-03  2.129399e-03  4.945629e-03  2.136432e-02;
    0.000000e+00  2.932301e-03  5.715869e-03  1.260737e-02  4.640307e-02;
    0.000000e+00  6.023225e-03  1.159824e-02  2.384944e-02  6.579589e-02]
const DF_IM = [
    0.000000e+00  3.201702e-07  1.065880e-06  3.194973e-06  1.061631e-05;
    0.000000e+00  2.132756e-06  7.082087e-06  2.119825e-05  7.032852e-05;
    0.000000e+00  1.360231e-05  4.487508e-05  1.338531e-04  4.425611e-04;
    0.000000e+00  8.150571e-05  2.645503e-04  7.821472e-04  2.564419e-03;
    0.000000e+00  4.468457e-04  1.395322e-03  4.034791e-03  1.292767e-02;
    0.000000e+00  1.565250e-03  4.596279e-03  1.277487e-02  3.882903e-02;
    0.000000e+00  2.738017e-03  7.506939e-03  1.986747e-02  5.657455e-02;
    0.000000e+00  4.734163e-03  1.135626e-02  2.670880e-02  6.212154e-02;
    0.000000e+00  6.722794e-03  1.323891e-02  2.422856e-02  2.829901e-02]

function f0_exit(h::Real, cfg::Config)
    c0, c1, x0 = F0_COEF[cfg.faces]
    X = h / (2 * ld(cfg))
    return (1 + 2 * (c0 + c1 * X / (X + x0)) * X)^3
end

"Lateral exit factor f(h, ω) = f0(h)[1 + Δf(h, ω)] (Δf tabulated for one open face; zero otherwise)."
function f_exit(h::Real, w::Real, cfg::Config)
    f0 = f0_exit(h, cfg)
    (w == 0 || cfg.faces != :one || !cfg.exit_omega) && return complex(f0)
    lH = log.(H_TAB)
    lh = log(clamp(h * 1e6, H_TAB[1], H_TAB[end]))
    i = clamp(searchsortedlast(lH, lh), 1, length(lH) - 1)
    th = (lh - lH[i]) / (lH[i+1] - lH[i])
    wf = min(w, W_TAB[end])
    j = clamp(searchsortedlast(W_TAB, wf), 1, length(W_TAB) - 1)
    tw = (wf - W_TAB[j]) / (W_TAB[j+1] - W_TAB[j])
    DF(a, b) = complex(DF_RE[a, b], DF_IM[a, b])
    df = (1 - th) * (1 - tw) * DF(i, j) + th * (1 - tw) * DF(i + 1, j) +
         (1 - th) * tw * DF(i, j + 1) + th * tw * DF(i + 1, j + 1)
    return f0 * (1 + df)
end

# ----------------------------------------------------------------------------- grid
"Graded grid towards y = 0 with Simpson midpoints; returns nodes y, Simpson weights ws, control volumes V."
function grid(Np::Int, lg::Float64, L::Float64)
    z = collect(0:Np) ./ Np .* log1p(L / lg)
    ye = lg .* expm1.(z)
    ye[end] = L
    y = zeros(2Np + 1)
    y[1:2:end] .= ye
    y[2:2:end] .= 0.5 .* (ye[1:end-1] .+ ye[2:end])
    d = diff(ye)
    ws = zeros(2Np + 1)
    ws[1:2:end-1] .+= d ./ 6
    ws[2:2:end] .+= 4 .* d ./ 6
    ws[3:2:end] .+= d ./ 6
    V = zeros(2Np + 1)
    dy = diff(y)
    V[1:end-1] .+= dy ./ 2
    V[2:end] .+= dy ./ 2
    return y, ws, V
end

# ----------------------------------------------------------------------------- modal solver for one film
"""
    film_impedance(y, V, h, H, w, fl, cfg; kap0=Inf, kapL=Inf, return_P=false)

Impedance matrix of one film. `H` is m×K (sensitivities h_i at the nodes). `kap0`, `kapL` are the end
(Robin) conductances at y = 0 and y = L: Inf = ambient pressure, 0 = sealed.
"""
function film_impedance(y, V, h, H::AbstractMatrix, w::Real, fl::Fluid, cfg::Config;
                        kap0::Real=Inf, kapL::Real=Inf, return_P::Bool=false)
    m = size(H, 1)
    K = length(y)
    G = [G_mob(hh, w, fl; slip=cfg.slip, wom=cfg.wom) for hh in h]
    f = cfg.exit ? [f_exit(hh, w, cfg) for hh in h] : ones(ComplexF64, K)
    c = (cfg.comp && w > 0) ? [12 * fl.eta * im * w * hh / (fl.npoly * fl.pa) for hh in h] : zeros(ComplexF64, K)
    l = ld(cfg)
    q = 2 .* collect(0:cfg.N-1) .+ 1
    k = q .* π ./ (2l)
    wn = 8 * cfg.W ./ (q .* π) .^ 2
    g = 2 .* G[1:end-1] .* G[2:end] ./ ((G[1:end-1] .+ G[2:end]) .* diff(y))   # harmonic-mean conductances
    B = H .* transpose(V)
    lo = isinf(kap0) ? 2 : 1
    hi = isinf(kapL) ? K - 1 : K
    Bs = B[:, lo:hi]
    RHS = ComplexF64.(Matrix(transpose(Bs)))
    Z = zeros(ComplexF64, m, m)
    Ps = Matrix{ComplexF64}[]
    for n in 1:cfg.N
        d = zeros(ComplexF64, K)
        d[1:end-1] .+= g
        d[2:end] .+= g
        d .+= (G .* k[n]^2 ./ f .+ c) .* V
        isinf(kap0) || (d[1] += kap0)
        isinf(kapL) || (d[end] += kapL)
        off = -g[lo:hi-1]
        T = Tridiagonal(copy(off), d[lo:hi], copy(off))
        X = T \ RHS
        Z .+= 12 * fl.eta * wn[n] .* (Bs * X)
        if return_P
            P = zeros(ComplexF64, K, m)
            P[lo:hi, :] .= -12 * fl.eta .* X
            push!(Ps, P)
        end
    end
    # tail: modes n ≥ N by local balance, explicit to N+M plus a quartic remainder
    qt = 2 .* collect(cfg.N:cfg.N+cfg.M-1) .+ 1
    kt = qt .* π ./ (2l)
    wt = 8 * cfg.W ./ (qt .* π) .^ 2
    S = [sum(wt ./ (G[i] .* kt .^ 2 ./ f[i] .+ c[i])) for i in 1:K]
    Q = 2 * (cfg.N + cfg.M) + 1
    S .+= f ./ G .* (8 * cfg.W * (2l)^2 / π^4) ./ (6.0 * (Q - 1)^3)
    Z .+= 12 * fl.eta .* ((H .* transpose(V .* S)) * transpose(H))
    return return_P ? (Z, Ps, k) : Z
end

# ----------------------------------------------------------------------------- MEMS case study (comb with sidewall taper)
Base.@kwdef struct Device
    L::Float64 = 400e-6          # overlap
    Lf::Float64 = 450e-6         # electrode length
    wt::Float64 = 9e-6           # compliant electrode width at root
    wb::Float64 = 30e-6          # ... at tip
    g0::Float64 = 14e-6          # bare tip gap
    Tp::Float64 = 120e-9         # coating thickness
    heff::Float64 = 50e-9        # residual air gap at contact
    eps::Float64 = 2e-9          # regularization of the positive part
    ls::Float64 = 25e-9          # sealing width
    nb::Int = 80                 # electrodes in parallel
    E::Float64 = 170e9
    Tf::Float64 = 25e-6
    theta::Int = 1               # 1: compliant set on the shuttle (design), 0: anchored (as built)
end
alpha(d::Device) = (d.wb - d.wt) / d.Lf
gcl(d::Device) = d.g0 - 2d.Tp - d.heff

"Normalized deflection shape of the tapered cantilever (Castigliano), tabulated."
function phi_table(d::Device; n::Int=4001)
    sg = collect(range(0, d.Lf, length=n))
    bw = d.wt .+ (d.wb - d.wt) .* sg ./ d.Lf
    EI = d.E * d.Tf .* bw .^ 3 ./ 12
    cum(fv) = vcat(0.0, cumsum(0.5 .* (fv[2:end] .+ fv[1:end-1]) .* diff(sg)))
    trap(fv) = sum(0.5 .* (fv[2:end] .+ fv[1:end-1]) .* diff(sg))
    I1 = cum((d.Lf .- sg) ./ EI)
    I2 = cum(sg .* (d.Lf .- sg) ./ EI)
    ke = 1 / trap((d.Lf .- sg) .^ 2 ./ EI)
    return sg, ke .* (sg .* I1 .- I2)
end

function interp1(x::Real, xp::AbstractVector, fp::AbstractVector)
    x <= xp[1] && return fp[1]
    x >= xp[end] && return fp[end]
    i = searchsortedlast(xp, x)
    t = (x - xp[i]) / (xp[i+1] - xp[i])
    return (1 - t) * fp[i] + t * fp[i+1]
end

sig(z, e) = z >= 0 ? 0.5 * (z + sqrt(z^2 + e^2)) : e^2 / (2 * (sqrt(z^2 + e^2) - z))
function dsig(z, e)
    R = sqrt(z^2 + e^2)
    return z >= 0 ? 0.5 * (1 + z / R) : e^2 / (2R * (R - z))
end
smoother(t) = (s = clamp(t, 0.0, 1.0); s^3 * (10 - 15s + 6s^2))

"Gap and sensitivities on side r of a compliant electrode (two generalized coordinates x1, x2)."
function gapfield(d::Device, phit, y, x1, x2, r)
    sg, ph = phit
    B2 = [interp1(d.Lf - yi, sg, ph) for yi in y]
    B1 = (2d.theta - 1) .- d.theta .* B2
    dr = gcl(d) .+ alpha(d) .* y .- r .* (B1 .* x1 .+ B2 .* x2)
    h = d.heff .+ sig.(dr, d.eps)
    sp = dsig.(dr, d.eps)
    H = vcat(transpose(-r .* B1 .* sp), transpose(-r .* B2 .* sp))
    return h, H
end

"Tip-end conductance: sealing law (contact) in series with the tip pocket."
function tip_kappa(d::Device, x1, x2, r, I0, kpocket)
    tau = x2 - (1 - d.theta) * x1
    delta = r * tau - gcl(d)
    chi_s = smoother((delta + d.ls) / (2d.ls))
    kseal = chi_s <= 0 ? Inf : (chi_s >= 1 ? 0.0 : (1 - chi_s) / (chi_s * I0))
    (kseal == 0 || kpocket == 0) && return 0.0
    isinf(kseal) && return kpocket
    isinf(kpocket) && return kseal
    return 1 / (1 / kseal + 1 / kpocket)
end

"Pocket conductance calibrated so the tip end at rest has sealing fraction chi_rest (3-D check: 0.75)."
function pocket_from_rest(d::Device, fl::Fluid; Np::Int=512, chi_rest::Float64=0.75)
    phit = phi_table(d)
    y, ws, V = grid(Np, d.heff / alpha(d), d.L)
    h, _ = gapfield(d, phit, y, 0.0, 0.0, 1)
    I0 = sum(ws ./ (h .^ 2 .* (h .+ kp(fl))))
    return chi_rest < 1 ? (1 - chi_rest) / (chi_rest * I0) : 0.0
end

"2×2 impedance of the whole comb (both gaps of every electrode) at state (x1, x2) and frequency w."
function device_impedance(d::Device, fl::Fluid, cfg::Config, x1, x2, w, kpocket;
                          Np::Int=512, return_fields::Bool=false, phit=phi_table(d))
    y, ws, V = grid(Np, d.heff / alpha(d), d.L)
    Z = zeros(ComplexF64, 2, 2)
    fields = []
    for r in (1, -1)
        h, H = gapfield(d, phit, y, x1, x2, r)
        I0 = sum(ws ./ (h .^ 2 .* (h .+ kp(fl))))
        k0 = tip_kappa(d, x1, x2, r, I0, kpocket)
        if return_fields
            Zr, Ps, kk = film_impedance(y, V, h, H, w, fl, cfg; kap0=k0, return_P=true)
            push!(fields, (y=y, h=h, Ps=Ps, k=kk))
        else
            Zr = film_impedance(y, V, h, H, w, fl, cfg; kap0=k0)
        end
        Z .+= Zr
    end
    Z .*= d.nb
    return return_fields ? (Z, fields) : Z
end

"Velocity pattern of rigid relative translation and the state at a given tip gap (theta = 1 convention)."
rigid(d::Device) = d.theta == 1 ? [1.0, 1.0] : [1.0, 0.0]
state_at_gap(d::Device, gap::Real) = gap <= d.heff + d.eps / 2 + 1e-12 ? (gcl(d), gcl(d)) : (gcl(d) + d.heff - gap, gcl(d) + d.heff - gap)
zrt(Z::AbstractMatrix, v) = transpose(v) * Z * v

"Pressure field p(y, ζ) (per unit velocity pattern v) on one wall from the modal solution (exact modes only)."
function pressure_field(field, cfg::Config, v::AbstractVector; nz::Int=61)
    y, Ps, k = field.y, field.Ps, field.k
    W = cfg.W
    zeta = collect(range(0, W, length=nz))
    p = zeros(ComplexF64, length(y), nz)
    for n in eachindex(Ps)
        if cfg.faces == :one            # ζ = 0 blocked, ζ = W open
            cn = 2 * (-1)^(n - 1) / (k[n] * W)
            psi = cos.(k[n] .* zeta)
        else                            # both open
            cn = 4 / (k[n] * W)
            psi = sin.(k[n] .* zeta)
        end
        Pn = Ps[n] * v
        p .+= cn .* (Pn .* transpose(psi))
    end
    return y, zeta, p
end

"Compact form of the saturated contact function Φsat(κ) (exact limits, < 0.5 % error)."
function phi_sat(kap::Real)
    a = 7 * 1.2020569031595942 / 4
    dd = 12 * (0.390725221631 - log(a))
    return 12 * log(1 + a / kap) + dd / (1 + 1.70135 * kap + kap^2 / 2.91591^2)
end

# ----------------------------------------------------------------------------- self-test
function selftest()
    fl = Fluid()
    npass = 0
    nfail = 0
    function check(name, val, ref; rtol=1e-7)
        err = maximum(abs.(val .- ref)) / maximum(abs.(ref))
        ok = err < rtol
        ok ? (npass += 1) : (nfail += 1)
        @printf("  %-44s rel. error %.2e  %s\n", name, err, ok ? "pass" : "FAIL")
    end
    println("selftest: values from the verified Python reference implementation")
    check("mobility, h = 1 μm, 1 kHz, no slip", G_mob(1e-6, 2π * 1e3, fl; slip=false),
          complex(9.999999983086157e-19, -4.091375926323227e-23))
    check("mobility, h = 5 μm, 300 kHz, slip", G_mob(5e-6, 2π * 3e5, fl),
          complex(1.2261397856574372e-16, -3.975331877170366e-17))
    y, ws, V = grid(800, 20e-6, 400e-6)
    h = fill(0.3e-6, length(y))
    Zu = film_impedance(y, V, h, ones(1, length(y)), 2π * 5e4, fl, Config(faces=:one, exit=false))
    check("uniform gap, one face, 50 kHz", Zu[1, 1], complex(0.004147517180803255, -0.0030798891626986594))
    d = Device()
    cfg = Config()
    kq = pocket_from_rest(d, fl)
    check("tip-pocket conductance", kq, 7.35124854117732e-12)
    M(a, b, c, e, f, g, i, j) = [complex(a, b) complex(c, e); complex(f, g) complex(i, j)]
    Z = device_impedance(d, fl, cfg, 0.0, 0.0, 0.0, kq)
    check("device, rest, quasi-static", Z, M(4.079354614915459e-06, 0, 3.4796419238902166e-06, 0,
          3.4796419238902154e-06, 0, 9.999189644220241e-06, 0))
    Z = device_impedance(d, fl, cfg, 0.0, 0.0, 2π * 1e3, kq)
    check("device, rest, 1 kHz", Z, M(4.101108280425428e-06, 1.5510232031831499e-07, 3.492371267034713e-06,
          9.963966119847769e-08, 3.4923712670347164e-06, 9.963966119847775e-08, 1.0018394754110667e-05, 1.9547900074755584e-07))
    Z = device_impedance(d, fl, cfg, gcl(d), gcl(d), 0.0, kq)
    check("device, contact, quasi-static", Z, M(2.2769025269078227e-05, 0, 0.0001297644243724945, 0,
          0.0001297644243724945, 0, 0.0034105606903015568, 0))
    Z = device_impedance(d, fl, cfg, gcl(d), gcl(d), 2π * 1e5, kq)
    check("device, contact, 100 kHz", Z, M(2.3025302069834133e-05, 1.2875063920648716e-05, 0.0001277517129643285,
          -5.80299639217988e-06, 0.0001277517129643287, -5.802996392179888e-06, 0.0033272595072590264, -0.000493911662795534))
    Z = device_impedance(d, fl, Config(exit=false), gcl(d), gcl(d), 0.0, Inf)
    check("device, contact, original Reynolds model", Z, M(1.615251412363449e-05, 0, 0.00011851523462518454, 0,
          0.00011851523462518462, 0, 0.003336466955055775, 0))
    @printf("selftest: %d passed, %d failed\n", npass, nfail)
    return nfail == 0
end

# ----------------------------------------------------------------------------- case-study summary
function run_case_study()
    fl = Fluid(); d = Device(); cfg = Config(); kq = pocket_from_rest(d, fl); v = rigid(d)
    println("rigid-translation impedance of the comb (one open face), Z = c + iX")
    @printf("  %-10s %12s %12s %12s %14s\n", "state", "c(0) N s/m", "c(1k)/c(0)", "X/c at 1k", "X/c at 100k")
    for (name, gap) in (("rest", 13.76e-6), ("gap 5 um", 5e-6), ("gap 1 um", 1e-6), ("contact", 0.0))
        x1, x2 = state_at_gap(d, gap)
        z0 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 0.0, kq), v))
        z1 = zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e3, kq), v)
        z5 = zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e5, kq), v)
        @printf("  %-10s %12.4e %12.5f %12.5f %14.5f\n", name, z0, real(z1) / z0, imag(z1) / real(z1), imag(z5) / real(z5))
    end
    Dstar = d.nb * fl.eta * cfg.W / alpha(d)^3
    kap = π * kp(fl) / (2 * alpha(d) * ld(cfg))
    @printf("contact groups: Λ = %.4f, κ = %.4f; saturated contact damping D⋆Φsat(κ) = %.3f mN s/m\n",
            π * (d.heff + d.eps / 2) / (2 * alpha(d) * ld(cfg)), kap, 1e3 * Dstar * phi_sat(kap))
end

# ----------------------------------------------------------------------------- plots (Okabe–Ito, Computer Modern)
const OI = (black=RGB(0, 0, 0), orange=RGB(230 / 255, 159 / 255, 0), sky=RGB(86 / 255, 180 / 255, 233 / 255),
            green=RGB(0, 158 / 255, 115 / 255), blue=RGB(0, 114 / 255, 178 / 255),
            vermillion=RGB(213 / 255, 94 / 255, 0), purple=RGB(204 / 255, 121 / 255, 167 / 255))

function setup_style()
    default(fontfamily="Computer Modern", framestyle=:box, grid=false, linewidth=1.6,
            guidefontsize=10, tickfontsize=8, legendfontsize=8, titlefontsize=10, legend_background_color=:transparent)
end

function plot_mobility(fl::Fluid=Fluid())
    lam = 10 .^ range(-2, 2, length=161)
    h = 10e-6
    p = plot(xscale=:log10, xlabel="Womersley parameter |λ| = (h/2)√(ω/ν)", ylabel="G(ω)/G(0)",
             title="Mobility with gas inertia and slip", legend=:bottomleft)
    for (bh, col) in ((0.0, OI.black), (0.05, OI.blue), (0.2, OI.vermillion))
        flb = bh > 0 ? Fluid(lam=bh * h / 1.016) : fl
        w = (2 .* lam ./ h) .^ 2 .* nu(flb)
        r = [G_mob(h, wi, flb; slip=bh > 0) / G_mob(h, 0.0, flb; slip=bh > 0) for wi in w]
        plot!(p, lam, real.(r), color=col, label="Re, b/h = $(bh)")
        plot!(p, lam, -imag.(r), color=col, ls=:dash, label="−Im, b/h = $(bh)")
    end
    return p
end

function plot_exit()
    hs = 10 .^ range(log10(0.1), log10(40), length=200)
    p1 = plot(xscale=:log10, xlabel="gap h (μm)", ylabel="f₀(h)", title="Lateral exit factor (W = 25 μm)", legend=:topleft)
    for (faces, col) in ((:one, OI.blue), (:two, OI.vermillion))
        cfg = Config(faces=faces)
        plot!(p1, hs, [f0_exit(h * 1e-6, cfg) for h in hs], color=col, label=faces == :one ? "one open face (fit)" : "two open faces (fit)")
    end
    one_h = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 13.76, 16.093, 18.427, 20.76, 23.093, 25.427, 27.76, 30.093, 32.427]
    one_f = [1.0125, 1.0253, 1.0513, 1.1058, 1.2247, 1.5037, 2.0051, 2.2457, 2.5092, 2.7969, 3.1105, 3.4518, 3.8231, 4.2278, 4.6701]
    two_h = [0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 13.76, 16.09, 18.43, 20.76, 23.09, 25.43, 27.76, 30.09, 32.43]
    two_f = [1.0253, 1.0514, 1.1059, 1.2243, 1.5007, 1.831, 2.219, 2.664, 3.170, 3.666, 4.407, 5.242, 6.179, 7.222, 8.381, 9.664, 11.083, 12.657]
    scatter!(p1, one_h, one_f, color=OI.blue, ms=3, label="Stokes, one face")
    scatter!(p1, two_h, two_f, color=OI.vermillion, ms=3, marker=:diamond, label="Stokes, two faces")
    p2 = plot(xscale=:log10, yscale=:log10, xlabel="gap h (μm)", ylabel="|f/f₀ − 1|", title="Frequency dependence (one face)", legend=:topleft)
    for (j, col, lab) in ((3, OI.green, "1 kHz"), (5, OI.orange, "10 kHz"))
        plot!(p2, H_TAB, abs.(complex.(DF_RE[:, j], DF_IM[:, j])), color=col, marker=:circle, ms=3, label=lab)
    end
    return plot(p1, p2, layout=(1, 2), size=(900, 360))
end

function plot_frequency(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    kq = pocket_from_rest(d, fl); v = rigid(d)
    F = 10 .^ range(0, 5, length=41)
    p1 = plot(xscale=:log10, xlabel="frequency (Hz)", ylabel="c(f)/c(0)", title="Damping", legend=:bottomleft)
    p2 = plot(xscale=:log10, xlabel="frequency (Hz)", ylabel="Im Z / Re Z", title="Reactance (+ mass-like, − spring-like)", legend=:topleft)
    for (name, gap, col) in (("rest", 13.76e-6, OI.black), ("gap 5 μm", 5e-6, OI.blue), ("gap 1 μm", 1e-6, OI.green), ("contact", 0.0, OI.vermillion))
        x1, x2 = state_at_gap(d, gap)
        z0 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 0.0, kq), v))
        zs = [zrt(device_impedance(d, fl, cfg, x1, x2, 2π * fh, kq), v) for fh in F]
        plot!(p1, F, real.(zs) ./ z0, color=col, label=name)
        plot!(p2, F, imag.(zs) ./ real.(zs), color=col, label=name)
    end
    vline!(p1, [1000], color=:gray, ls=:dot, label="1 kHz")
    vline!(p2, [1000], color=:gray, ls=:dot, label="")
    return plot(p1, p2, layout=(1, 2), size=(900, 360))
end

function plot_gap(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    kq = pocket_from_rest(d, fl); v = rigid(d)
    gaps = 10 .^ range(log10(0.051), log10(13.76), length=28)
    rows = map(gaps) do g
        x1, x2 = state_at_gap(d, g * 1e-6)
        o = real(zrt(device_impedance(d, fl, Config(cfg; exit=false), x1, x2, 0.0, Inf), v))
        c0 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 0.0, kq), v))
        c10 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e4, kq), v))
        c100 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e5, kq), v))
        (o, c0, c10, c100)
    end
    p = plot(xscale=:log10, yscale=:log10, xlabel="tip gap (μm)", ylabel="damping Re Z (mN s/m)",
             title="Damping across the stroke", legend=:topright)
    plot!(p, gaps, 1e3 .* getindex.(rows, 1), color=OI.black, ls=:dash, label="Reynolds, quasi-static (original)")
    plot!(p, gaps, 1e3 .* getindex.(rows, 2), color=OI.blue, label="full model, quasi-static")
    plot!(p, gaps, 1e3 .* getindex.(rows, 3), color=OI.green, ls=:dot, label="full model, 10 kHz")
    plot!(p, gaps, 1e3 .* getindex.(rows, 4), color=OI.vermillion, ls=:dashdot, label="full model, 100 kHz")
    return p
end

Config(c::Config; kw...) = Config(; W=c.W, faces=c.faces, slip=c.slip, wom=c.wom, comp=c.comp, exit=c.exit,
                                   exit_omega=c.exit_omega, N=c.N, M=c.M, kw...)

function plot_contact(fl::Fluid=Fluid(), cfg::Config=Config())
    he = 10 .^ range(log10(0.5e-9), log10(1e-6), length=20)
    res = map(he) do h_
        d = Device(heff=h_); kq = pocket_from_rest(d, fl); v = rigid(d)
        [zrt(device_impedance(d, fl, cfg, gcl(d), gcl(d), 2π * fh, kq), v) for fh in (0.0, 1e4, 1e5)]
    end
    hmin = (he .+ 1e-9) .* 1e9
    d0 = Device()
    Dstar = d0.nb * fl.eta * cfg.W / alpha(d0)^3
    kap = π * kp(fl) / (2 * alpha(d0) * ld(cfg))
    p = plot(xscale=:log10, xlabel="residual gap h_min (nm)", ylabel="contact damping Re Z (mN s/m)",
             title="Contact limit", legend=:topright)
    for (j, col, lab) in ((1, OI.blue, "quasi-static"), (2, OI.green, "10 kHz"), (3, OI.vermillion, "100 kHz"))
        plot!(p, hmin, 1e3 .* real.(getindex.(res, j)), color=col, label=lab)
    end
    hline!(p, [1e3 * Dstar * phi_sat(kap)], color=OI.black, ls=:dash, label="saturated limit D⋆Φsat(κ)")
    return p
end

function plot_pressure(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config(); w::Real=0.0)
    kq = pocket_from_rest(d, fl); v = rigid(d)
    plots_ = []
    for (name, gap, ymax) in (("contact", 0.0, 40e-6), ("rest", 13.76e-6, 400e-6))
        x1, x2 = state_at_gap(d, gap)
        _, fields = device_impedance(d, fl, cfg, x1, x2, w, kq; return_fields=true)
        y, zeta, p = pressure_field(fields[1], cfg, v)
        keep = y .<= ymax
        hm = heatmap(1e6 .* y[keep], 1e6 .* zeta, transpose(real.(p[keep, :])), color=:viridis,
                     xlabel="y (μm)", ylabel="ζ (μm)  (ζ = W open)", title="Re p per unit closing speed, $(name)",
                     colorbar_title="Pa s/m")
        push!(plots_, hm)
    end
    return plot(plots_..., layout=(1, 2), size=(950, 360))
end

function make_plots(outdir::AbstractString="figs")
    setup_style()
    mkpath(outdir)
    for (name, f) in (("mobility", plot_mobility), ("exit_factor", plot_exit), ("frequency", plot_frequency),
                      ("damping_vs_gap", plot_gap), ("contact_limit", plot_contact), ("pressure_field", plot_pressure))
        p = f()
        savefig(p, joinpath(outdir, name * ".pdf"))
        savefig(p, joinpath(outdir, name * ".png"))
        println("wrote ", joinpath(outdir, name), ".{pdf,png}")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    selftest()
    run_case_study()
    make_plots()
end