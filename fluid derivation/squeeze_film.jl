
using LinearAlgebra, Printf
using Plots

# ----------------------------------------------------------------------------- fluid and mobility
Base.@kwdef struct Fluid
    eta::Float64 = 1.849e-5      # dynamic viscosity, Pa s
    pa::Float64 = 101325.0       # ambient pressure, Pa
    rho::Float64 = 1.204         # density, kg/m^3
    lam::Float64 = eta * sqrt(2 / (rho * pa))   # equivalent free path η v0/p, v0 = sqrt(2kT/m) = sqrt(2p/ρ), m
    sig_p::Float64 = 1.016191    # viscous slip coefficient, BGK, diffuse walls (referred to lam)
    npoly::Float64 = 1.0         # polytropic index (1 = isothermal)
end
nu(fl::Fluid) = fl.eta / fl.rho
bslip(fl::Fluid) = fl.sig_p * fl.lam            # Navier slip length b
kp(fl::Fluid) = 6 * bslip(fl)                    # k_p = 6 σ_p λ, G0 = h²(h + k_p)

"Womersley mobility with first-order slip (complex). q = −G ∇p/(12η)."
function G_mob(h::Real, w::Real, fl::Fluid; slip::Bool=true, wom::Bool=true, kinetic::Bool=false)
    # kinetic: times the BGK flow factor (exact in the steady rarefied and inertial continuum limits; within 1e-4 of
    # oscillatory BGK for the case study at 100 kHz)
    kinetic && slip && return G_mob(h, w, fl; slip=true, wom=wom, kinetic=false) * kin_factor(h, fl)
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
    kinetic::Bool = true         # linearized-BGK flow factor in the mobility (false: first-order slip)
    N::Int = 24                  # exact thickness modes
    M::Int = 2000                # explicit tail terms (local balance)
end
ld(cfg::Config) = cfg.faces == :one ? cfg.W : cfg.W / 2

# ----------------------------------------------------------------------------- lateral exit factor
# f0(h) = (1 + 2 c(X) X)^3,  c(X) = c0 + c1 X/(X + x0),  X = h/(2 ℓd)   (fits to 2-D Stokes cells, < 0.21 %)
const F0_COEF = Dict(:one => (0.4207, 0.2449, 1.0015), :two => (0.4066, 0.1279, 0.2733))
# frequency correction f/f0 − 1 relative to the Womersley strip, one open face, comb geometry of the case study:
# oscillatory 2-D Stokes cells at half-decade frequencies, 100 Hz – 100 kHz (rows: gaps H_TAB; columns: F_TAB)
const H_TAB = [0.5, 1.0, 2.0, 4.0, 8.0, 13.76, 18.426666666666666, 25.42666666666667, 32.42666666666667]   # μm
const F_TAB = [100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0]   # Hz
const W_TAB = 2π .* F_TAB                                                                    # rad/s
const DF_RE = [
    4.930871266850545e-10 1.029372365124459e-09 3.5299974054936456e-09 2.0742613626723028e-08 1.9979098264677475e-07 1.4269531876109909e-06 7.3759677263751655e-06
    1.0413612194781763e-08 1.9258259875698513e-08 4.894635430297001e-08 2.043413191987753e-07 1.6926643149339782e-06 1.1651598524675677e-05 5.977400871715055e-05
    1.790121824107871e-07 3.1997157923235875e-07 7.00356361571508e-07 2.2195297963989447e-06 1.5114307648111946e-05 9.804848466243854e-05 0.000504068697839033
    2.7369657242815038e-06 4.828726461569843e-06 9.860481912093988e-06 2.6111644112036814e-05 0.00014607430318291925 0.0008756076765399357 0.004641897835421016
    3.6056393674765985e-05 6.33355290051707e-05 0.000126056536631447 0.00030762576564113964 0.0015122274879417752 0.00833796987444102 0.044819882720136706
    0.00023705671681972795 0.00041653564806010124 0.000827978266419338 0.001998578119734029 0.009354829394971453 0.046446508629397254 0.17996862181414475
    0.0006160889087654109 0.001081416403281521 0.002129399465274018 0.0049456293118139705 0.02136431863589694 0.09390560927073333 0.21402616208106107
    0.00167139478758771 0.0029323008487662783 0.005715868832835147 0.012607365081189092 0.046403072578815596 0.13656911961710705 0.11772021269121158
    0.0034347901702922456 0.006023224541664485 0.011598235905109444 0.023849443228960432 0.06579588861893337 0.08947758568071351 -0.02649471189698449]
const DF_IM = [
    1.06774053210707e-07 3.201701608259986e-07 1.0658800767759546e-06 3.1949725155076796e-06 1.0616312090545994e-05 3.129325628652446e-05 9.879225953666778e-05
    7.151856559015418e-07 2.132755570093563e-06 7.0820872542659904e-06 2.1198247416837e-05 7.0328523749068e-05 0.00020635418578850384 0.0006438570670016288
    4.606923714058222e-06 1.360230535472074e-05 4.487508449518803e-05 0.0001338530861915424 0.0004425611483669063 0.001287742918615706 0.0039416406649351865
    2.8285442185999945e-05 8.150570747944809e-05 0.000264550336139876 0.0007821471739812263 0.002564418714168723 0.007329010820110048 0.021557869455570193
    0.00016368205696471038 0.00044684566058424434 0.0013953224729276025 0.004034790724196143 0.012927669386713745 0.03513281921049477 0.08902571934034353
    0.0006186902519569605 0.0015652497585513202 0.004596278545178677 0.012774870801075214 0.03882902699020296 0.0915547713259698 0.12528119467806492
    0.0011649050054218322 0.00273801681061391 0.007506939295488353 0.019867471195868517 0.0565745539333774 0.10687201811248953 0.04383347165639794
    0.0022637071471180415 0.004734163238234336 0.011356256560514031 0.026708801870476126 0.06212154223533083 0.04843101599949223 -0.06289446562152184
    0.0036535746160106367 0.006722794430620192 0.013238911501491412 0.02422855696261024 0.028299013802735626 -0.054816164160152316 -0.10947153965117401]

# log|Δf| and arg Δf: monotone cubic (Fritsch–Carlson) in log ω, linear in log h. Exact for a power law in either
# variable, so it carries both the iω scaling of thin gaps and the √(iω) Stokes-layer scaling of wide gaps.
"Fritsch–Carlson slopes of each row of Y on nodes x (zero at interior extrema, secant slopes at the ends)."
function pchip_slopes(x::Vector{Float64}, Y::Matrix{Float64})
    hx = diff(x)
    S = diff(Y, dims=2) ./ hx'
    D = zeros(size(Y))
    D[:, 1] = S[:, 1]
    D[:, end] = S[:, end]
    for k in 2:length(x)-1
        w1 = 2hx[k] + hx[k-1]
        w2 = hx[k] + 2hx[k-1]
        for r in axes(Y, 1)
            a, b = S[r, k-1], S[r, k]
            D[r, k] = a * b > 0 ? (w1 + w2) / (w1 / a + w2 / b) : 0.0
        end
    end
    return D
end
const LOG_H = log.(H_TAB)
const LOG_W = log.(W_TAB)
const LOG_ABS_DF = log.(abs.(complex.(DF_RE, DF_IM)))
const ARG_DF = angle.(complex.(DF_RE, DF_IM))          # continuous over the table, no unwrapping needed
const D_LA = pchip_slopes(LOG_W, LOG_ABS_DF)
const D_PH = pchip_slopes(LOG_W, ARG_DF)

# the same for TWO open faces (cells with the device layer symmetric about its mid-plane; same gaps and frequencies)
const DF2_RE = [
    5.194200847213892e-09 9.38320177112928e-09 2.125398168573156e-08 7.28872135979941e-08 5.320932772168163e-07 3.5800958475551425e-06 1.8116484451491388e-05
    8.948479668369202e-08 1.578109933841887e-07 3.2127889793009956e-07 8.441542942438929e-07 4.7151545945478546e-06 2.908574366133898e-05 0.00014392382249428515
    1.3682706179185544e-06 2.3915590212642e-06 4.632238136847988e-06 1.0329275158760254e-05 4.4179492744600424e-05 0.00024029874802788953 0.0011603733542888683
    1.822864367873045e-05 3.17446829332102e-05 6.0139498614164566e-05 0.00012300684010035923 0.0004291101395894614 0.0020251386310452535 0.009699264040461175
    0.0001962008228715284 0.0003413572938490983 0.000642202594707042 0.0012721707001277505 0.004000433864122632 0.016802435646541936 0.07837845629451712
    0.0010370326793915918 0.001806071293123157 0.0034074125445431314 0.006777746921763805 0.020904626770354318 0.07918718193035446 0.2575571344287444
    0.0023607829503646816 0.0041136777601114005 0.007757042273612846 0.015265316150101738 0.044587509021164706 0.1470501928552157 0.26899768782966027
    0.005539169759354756 0.009666669035275799 0.018255978674144524 0.03549863799037256 0.09324488125996866 0.2024060047434355 0.11107801672015172
    0.010227955983772574 0.017883503798902245 0.03381826819089517 0.06440146887108145 0.13933580085345998 0.14139513329181352 -0.0774430530872946]
const DF2_IM = [
    2.2658304703743533e-07 6.738513883654653e-07 2.2325764409021173e-06 6.675122411096662e-06 2.2134332402046252e-05 6.497531171465839e-05 0.00020282670294215175
    1.5356782158240231e-06 4.497288984158307e-06 1.4757714552001688e-05 4.38994246748213e-05 0.0001449960027166423 0.0004229510152095991 0.001301767011702222
    1.013326668806818e-05 2.8721217253778222e-05 9.216883965839367e-05 0.00027089608552531504 0.000886735770350801 0.0025556074999033424 0.007683178587911727
    6.489571290118013e-05 0.00017230363558315344 0.0005267726676038446 0.0015060818652694466 0.004826750313604147 0.013549665579451952 0.038783894873441055
    0.0003951257586097278 0.0009443264836266707 0.0026350775445153317 0.007104419878622172 0.021640616463783034 0.0566121765703413 0.1356273820742361
    0.0015249156554319003 0.0032999244760524805 0.008282108448878651 0.020554480166799626 0.056966202461020445 0.12376388532105237 0.13918176042068284
    0.0029564355031239348 0.00596486492226665 0.013676825039259807 0.03116126701348497 0.07711542958469847 0.12426905948488967 0.001408234289127941
    0.005976146266773338 0.011103134828167473 0.02235341975125629 0.04334235564850303 0.0780709567705293 0.019542876668365575 -0.14462411328648458
    0.009977808540296323 0.01730431983288824 0.03037820517711916 0.045995641625171224 0.02874789999724247 -0.12886729033250122 -0.19755826430575496]
const LOG_ABS_DF2 = log.(abs.(complex.(DF2_RE, DF2_IM)))
const ARG_DF2 = angle.(complex.(DF2_RE, DF2_IM))
const D_LA2 = pchip_slopes(LOG_W, LOG_ABS_DF2)
const D_PH2 = pchip_slopes(LOG_W, ARG_DF2)

"Cubic Hermite in log ω along table row i. Below the band: log|Δf| extended linearly (power law of the first interval,
Δf → 0) with the phase held. Above it: Δf held, since cell and strip both turn inertial and their ratio saturates;
this keeps Re(G/f) > 0 (passivity) up to 2.8 MHz for the case study."
function along_w(T::Matrix{Float64}, DT::Matrix{Float64}, i::Int, lw::Float64, hold::Bool)
    n = length(LOG_W)
    lw <= LOG_W[1] && return hold ? T[i, 1] : T[i, 1] + DT[i, 1] * (lw - LOG_W[1])
    lw >= LOG_W[n] && return T[i, n]
    j = clamp(searchsortedlast(LOG_W, lw), 1, n - 1)
    dx = LOG_W[j+1] - LOG_W[j]
    t = (lw - LOG_W[j]) / dx
    return (2t^3 - 3t^2 + 1) * T[i, j] + (t^3 - 2t^2 + t) * dx * DT[i, j] +
           (-2t^3 + 3t^2) * T[i, j+1] + (t^3 - t^2) * dx * DT[i, j+1]
end

"Steady exit factor f0(h) = (1 + 2 c(X) X)^3 with c(X) = c0 + c1 X/(X + x0), X = h/(2ℓd)."
function f0_exit(h::Real, cfg::Config)
    c0, c1, x0 = F0_COEF[cfg.faces]
    X = h / (2 * ld(cfg))
    return (1 + 2 * (c0 + c1 * X / (X + x0)) * X)^3
end

"Lateral exit factor f(h, ω) = f0(h)[1 + Δf(h, ω)], one or two open faces; h clamped to the table range."
function f_exit(h::Real, w::Real, cfg::Config)
    f0 = f0_exit(h, cfg)
    (w == 0 || !cfg.exit_omega) && return complex(f0)
    LA, DLA, PH, DPH = cfg.faces == :one ? (LOG_ABS_DF, D_LA, ARG_DF, D_PH) : (LOG_ABS_DF2, D_LA2, ARG_DF2, D_PH2)
    lh = log(clamp(h * 1e6, H_TAB[1], H_TAB[end]))
    i = clamp(searchsortedlast(LOG_H, lh), 1, length(LOG_H) - 1)
    th = (lh - LOG_H[i]) / (LOG_H[i+1] - LOG_H[i])
    lw = log(Float64(w))
    la = (1 - th) * along_w(LA, DLA, i, lw, false) + th * along_w(LA, DLA, i + 1, lw, false)
    ph = (1 - th) * along_w(PH, DPH, i, lw, true) + th * along_w(PH, DPH, i + 1, lw, true)
    return f0 * (1 + exp(la + im * ph))
end

# ----------------------------------------------------------------------------- kinetic flow factor
# Linearized BGK plane Poiseuille flow, diffuse walls (validated against Barichello et al. 2001 to seven digits):
# reduced flow rate G_P(δ), δ = h/lam, tabulated for 1e-5 ≤ δ ≤ 60; ln(δ G_P) (monotone) by PCHIP in ln δ.
const KIN_D = [9.999999999999999e-06, 1.1388973796922415e-05, 1.2970872414698539e-05, 1.4772492605422543e-05, 1.682435311983875e-05, 1.916121168320134e-05, 2.1822653777726372e-05, 2.4853763205383563e-05, 2.8305885790102787e-05, 3.2237499156215916e-05, 3.671520331684516e-05, 4.181484885242285e-05, 4.7622821790251516e-05, 5.423750695046804e-05, 6.177095454692779e-05, 7.03507782745846e-05, 8.012231703623429e-05, 9.125109692743828e-05, 0.0001039256351847022, 0.00011836063359470917, 0.00013480061545972778, 0.00015352406772798543, 0.00017484815845509683, 0.00019913410950852365, 0.00022679331552660548, 0.0002582943127849667, 0.00029417071602020684, 0.0003350302576576041, 0.0003815650825638619, 0.0004345634727140361, 0.0004949232003839766, 0.0005636667360662092, 0.0006419585687254839, 0.0007311249317924354, 0.0008326762690460736, 0.0009483328209484852, 0.0010800537648543815, 0.0012300704027193954, 0.0014009239584940997, 0.0015955086254770129, 0.00181712059283214, 0.002069513881761337, 0.002356963937174706, 0.002684340052077382, 0.0030571878515138653, 0.003481823233316095, 0.003965439356975268, 0.004516228492987621, 0.0051435207967550425, 0.005857942357816869, 0.006671595201705823, 0.007598262293590093, 0.008653641016384118, 0.00985560907835718, 0.011224527354612058, 0.012783584792451562, 0.014559191223196672, 0.01658142473453697, 0.01888454118172828, 0.021507554468560564, 0.024494897427831785, 0.02789717449638785, 0.03177201893475335, 0.036185069112322873, 0.04121108039600719, 0.046935191477298896, 0.05345436658884934, 0.06087903804114902, 0.06933497690324891, 0.07896542351613227, 0.08993351392881116, 0.10242504336003874, 0.11665161349761234, 0.13285421694930283, 0.15130731956462556, 0.1723235097804087, 0.19625879374827782, 0.2235186259414737, 0.25456477739715466, 0.28992315793955825, 0.3301927248894628, 0.37605562917005037, 0.42828877068028764, 0.4877769586793909, 0.5555279001142092, 0.632689269786006, 0.72056815151868, 0.8206531796543067, 0.9346397559443963, 1.0644587690012692, 1.2123093028059746, 1.3806958883422527, 1.5724709293848431, 1.7908830211186217, 2.0396319800873237, 2.3229315176579513, 2.6455806186651625, 3.013044834362333, 3.431548866750505, 3.908182012628031, 4.451018253542416, 5.069253025921794, 5.773358988219298, 6.57526342370561, 7.488550284044557, 8.528690296191936, 9.713303030539645, 11.062455369638311, 12.599001433443439, 14.348969719287528, 16.342004014579885, 18.611865551125124, 21.1970049073607, 24.141213346316686, 27.494364622711448, 31.31325982510913, 35.66258956443913, 40.6160298079796, 46.25748992180995, 52.68253406308962, 60.0]
const KIN_G = [6.853741530642733, 6.780417036395427, 6.70709864625514, 6.633787027380521, 6.560482917199392, 6.48718713044266, 6.41390056682413, 6.340624219433568, 6.267359183879465, 6.194106668252924, 6.120868003978412, 6.047644657601969, 5.974438243602629, 5.901250538281386, 5.8280834948168305, 5.7549392595535895, 5.6818201896005185, 5.608728871830426, 5.535668143342408, 5.462641113487963, 5.389651187520524, 5.316702091962343, 5.243797901751522, 5.170943069249882, 5.0981424551656565, 5.0254013614606245, 4.9527255662769125, 4.8801213609322325, 4.8075955889972075, 4.735155687468187, 4.662809730019578, 4.590566472306283, 4.518435399255351, 4.446426774265302, 4.3745516901962125, 4.302822122003805, 4.231250980833252, 4.159852169350721, 4.088640638050669, 4.017632442232379, 3.946844799300257, 3.8762961459934604, 3.806006195114114, 3.7359959912760687, 3.6662879651664415, 3.5969059857725214, 3.527875410009863, 3.459223129163532, 3.3909776115586183, 3.323168940877811, 3.2558288495786902, 3.1889907469017262, 3.122689741036616, 3.056962655101532, 2.9918480367196945, 2.927386161125368, 2.863619027925489, 2.8005903518623176, 2.7383455481908214, 2.6769317135839907, 2.6163976038294536, 2.556793609967119, 2.498171734955004, 2.440585573428212, 2.3840902976452587, 2.328742653289812, 2.2746009694234677, 2.221725187563767, 2.1701769156012727, 2.1200195130720405, 2.0713182151818263, 2.0241403039436197, 1.9785553358634986, 1.9346354368070882, 1.8924556760352085, 1.8520945329388177, 1.8136344717810189, 1.7771626418067477, 1.742771722478073, 1.7105609363879222, 1.6806372556853115, 1.6531168316761178, 1.628126681746306, 1.605806672958995, 1.5863118477019629, 1.5698151436621617, 1.556510568243731, 1.5466168963235072, 1.5403819699286554, 1.538087688900252, 1.5400557926824134, 1.5466545447246993, 1.5583064421830306, 1.5754970841082905, 1.5987853405392412, 1.628814972327681, 1.6663278568362379, 1.7121789780994436, 1.7673533428030044, 1.8329849880382287, 1.9103782575637671, 2.0010315465136164, 2.1066637580670333, 2.2292437879428237, 2.371023623157809, 2.534574541208286, 2.7228299396244884, 2.9391323933134523, 3.1872874929151367, 3.471626484771623, 3.7970783108401807, 4.169252301257675, 4.594533428667247, 5.080191711431437, 5.634507196738212, 6.266911469811236, 6.98814994844532, 7.810465596892156, 8.747807880198096, 9.816070350514547, 11.033360902484304]
const LKD = log.(KIN_D)
const LKG = log.(KIN_D .* KIN_G)
const DKG = pchip_slopes(LKD, reshape(LKG, 1, :))[1, :]

"Reduced Poiseuille flow rate G_P(δ) (linearized BGK, diffuse walls); → δ/6 + σ_p in the slip limit."
function gp_bgk(dl::Real)
    dl < KIN_D[1] && return 0.35887 + 0.56410 * log(1 / dl)
    dl > KIN_D[end] && return dl / 6 + 1.016191 + 1.0650 / dl - 2.1246 / dl^2
    lw = log(dl)
    j = clamp(searchsortedlast(LKD, lw), 1, length(LKD) - 1)
    dx = LKD[j+1] - LKD[j]
    t = (lw - LKD[j]) / dx
    v = (2t^3 - 3t^2 + 1) * LKG[j] + (t^3 - 2t^2 + t) * dx * DKG[j] +
        (-2t^3 + 3t^2) * LKG[j+1] + (t^3 - t^2) * dx * DKG[j+1]
    return exp(v) / dl
end

"BGK flow rate over the first-order-slip flow rate, Q(δ)/(1 + 6σ_p/δ) ≥ 1; → 1 for thick films."
kin_factor(h::Real, fl::Fluid) = (dl = h / fl.lam; 6 * gp_bgk(dl) / dl / (1 + 6 * fl.sig_p / dl))

"Steady mobility with slip and (optionally) the BGK factor, used for the end conductances."
G_steady(h::Real, fl::Fluid; kinetic::Bool=true) = real(G_mob(h, 0.0, fl; slip=true, wom=true, kinetic=kinetic))

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
    G = [G_mob(hh, w, fl; slip=cfg.slip, wom=cfg.wom, kinetic=cfg.kinetic) for hh in h]
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
    seal_centred::Bool = false   # false: sealing starts at contact (tip vented at nominal contact); true: original law centred on contact
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

"Tip-end conductance: contact sealing law in series with the tip pocket. Sealing starts at nominal contact and completes
at an overlap of 2ls (seal_centred = true restores the law centred on contact, χ = 1/2 at nominal contact)."
function tip_kappa(d::Device, x1, x2, r, I0, kpocket)
    tau = x2 - (1 - d.theta) * x1
    delta = r * tau - gcl(d)
    chi_s = smoother(d.seal_centred ? (delta + d.ls) / (2 * d.ls) : delta / (2 * d.ls))
    kseal = chi_s <= 0 ? Inf : (chi_s >= 1 ? 0.0 : (1 - chi_s) / (chi_s * I0))
    (kseal == 0 || kpocket == 0) && return 0.0
    isinf(kseal) && return kpocket
    isinf(kpocket) && return kseal
    return 1 / (1 / kseal + 1 / kpocket)
end

"Pocket conductance calibrated so the tip end at rest has sealing fraction chi_rest (3-D check: 0.75)."
function pocket_from_rest(d::Device, fl::Fluid; Np::Int=512, chi_rest::Float64=0.75, kinetic::Bool=true)
    phit = phi_table(d)
    y, ws, V = grid(Np, d.heff / alpha(d), d.L)
    h, _ = gapfield(d, phit, y, 0.0, 0.0, 1)
    I0 = sum(ws ./ G_steady.(h, Ref(fl); kinetic=kinetic))
    return chi_rest < 1 ? (1 - chi_rest) / (chi_rest * I0) : 0.0
end

pocket_for(d::Device, fl::Fluid, cfg::Config) = pocket_from_rest(d, fl; kinetic=cfg.kinetic) * (cfg.faces == :two ? 2.0 : 1.0)

function device_impedance(d::Device, fl::Fluid, cfg::Config, x1, x2, w, kpocket;
                          Np::Int=512, return_fields::Bool=false, phit=phi_table(d))
    y, ws, V = grid(Np, d.heff / alpha(d), d.L)
    Z = zeros(ComplexF64, 2, 2)
    fields = []
    for r in (1, -1)
        h, H = gapfield(d, phit, y, x1, x2, r)
        I0 = sum(ws ./ G_steady.(h, Ref(fl); kinetic=cfg.kinetic))
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
          complex(1.2319407624094012e-16, -4.0096701647538716e-17))
    check("BGK flow rate, δ = 0.68 (table)", gp_bgk(0.68), 1.5620295501648154)
    check("BGK flow rate, δ = 3e-4 (table)", gp_bgk(3e-4), 4.941766682114505)
    check("BGK flow rate, δ = 150 (slip tail)", gp_bgk(150.0), 26.023196573333333)
    check("kinetic factor, h = 51 nm", kin_factor(51e-9, fl), 1.3824988249729968)
    y, ws, V = grid(800, 20e-6, 400e-6)
    h = fill(0.3e-6, length(y))
    Zu = film_impedance(y, V, h, ones(1, length(y)), 2π * 5e4, fl, Config(faces=:one, exit=false))
    check("uniform gap, one face, 50 kHz", Zu[1, 1], complex(0.003963414055252005, -0.0025863083418439714))
    d = Device()
    cfg = Config()
    check("exit factor, h = 10 μm, 2 kHz (off-node)", f_exit(10e-6, 2π * 2e3, cfg), complex(1.6630860440759023, 0.007319744271161413))
    check("exit factor, h = 30 μm, 50 Hz (below band)", f_exit(30e-6, 2π * 50.0, cfg), complex(4.219635810452772, 0.008963189930125516))
    check("exit factor, h = 5 μm, 200 kHz (above band)", f_exit(5e-6, 2π * 2e5, cfg), complex(1.3033931869186992, 0.04424849377169374))
    kq = pocket_from_rest(d, fl)
    check("tip-pocket conductance", kq, 7.363676717037202e-12)
    M(a, b, c, e, f, g, i, j) = [complex(a, b) complex(c, e); complex(f, g) complex(i, j)]
    Z = device_impedance(d, fl, cfg, 0.0, 0.0, 0.0, kq)
    check("device, rest, quasi-static", Z, M(4.074353435986358e-06, 0.0, 3.4744223824113605e-06, 0.0, 3.474422382411361e-06, 0.0, 9.98098995317206e-06, 0.0))
    Z = device_impedance(d, fl, cfg, 0.0, 0.0, 2π * 1e3, kq)
    check("device, rest, 1 kHz", Z, M(4.09532157111549e-06, 1.5495856279814438e-07, 3.486551183329197e-06, 9.928415690749089e-08, 3.4865511833291954e-06, 9.928415690749088e-08, 9.999083817469267e-06, 1.9424948942095587e-07))
    Z = device_impedance(d, fl, cfg, 0.0, 0.0, 2π * 2e4, kq)
    check("device, rest, 20 kHz (off-node)", Z, M(4.352965265027708e-06, 2.6860155799830123e-06, 3.6515463671379053e-06, 1.7422837515963682e-06, 3.651546367137907e-06, 1.7422837515963686e-06, 1.028716103934848e-05, 3.5042283161038e-06))
    Z = device_impedance(d, fl, cfg, gcl(d), gcl(d), 0.0, kq)
    check("device, contact, quasi-static", Z, M(2.2480337164960124e-05, 0.0, 0.00012042028764719671, 0.0, 0.00012042028764719669, 0.0, 0.00271665899588216, 0.0))
    Z = device_impedance(d, fl, cfg, gcl(d), gcl(d), 2π * 1e5, kq)
    check("device, contact, 100 kHz", Z, M(2.409754579480602e-05, 1.3976131499144832e-05, 0.00012006571742329255, -1.9028536884916797e-06, 0.00012006571742329269, -1.9028536884916914e-06, 0.0026707190005156086, -0.00032752467863648066))
    Z = device_impedance(Device(seal_centred=true), fl, cfg, gcl(d), gcl(d), 0.0, kq)
    check("device, contact, centred sealing law", Z, M(2.2529155279676396e-05, 0.0, 0.00012485948436139068, 0.0, 0.0001248594843613907, 0.0, 0.0031286352030933714, 0.0))
    Z = device_impedance(d, fl, Config(exit=false), gcl(d), gcl(d), 0.0, Inf)
    check("device, contact, original Reynolds model", Z, M(1.589424904839425e-05, 0.0, 0.00010941304010825427, 0.0, 0.00010941304010825434, 0.0, 0.0026473377327829355, 0.0))
    c2 = Config(faces=:two)
    check("exit factor, two faces, h = 10 μm, 2 kHz", f_exit(10e-6, 2π * 2e3, c2), complex(2.6681111804714766, 0.020554019532760776))
    kq2 = pocket_for(d, fl, c2)
    check("tip pocket, two faces (doubled)", kq2, 1.4727353434074404e-11)
    Z = device_impedance(d, fl, c2, 0.0, 0.0, 2π * 1e3, kq2)
    check("device, two faces, rest, 1 kHz", Z, M(2.611543189311072e-06, 1.2813926985160897e-07, 1.954384327930891e-06, 7.338781538153949e-08, 1.954384327930892e-06, 7.338781538153949e-08, 5.179018672612699e-06, 1.3105968152557517e-07))
    Z = device_impedance(d, fl, c2, gcl(d), gcl(d), 0.0, kq2)
    check("device, two faces, contact, quasi-static", Z, M(9.032623939079388e-06, 0.0, 4.889543534621846e-05, 0.0, 4.8895435346218465e-05, 0.0, 0.0015428252356469215, 0.0))
    @printf("selftest: %d passed, %d failed\n", npass, nfail)
    return nfail == 0
end

# ----------------------------------------------------------------------------- case-study summary
function run_case_study(; faces::Symbol=:one)
    fl = Fluid(); d = Device(); cfg = Config(faces=faces); kq = pocket_for(d, fl, cfg); v = rigid(d)
    println("rigid-translation impedance of the comb (", faces == :one ? "one open face" : "two open faces", "), Z = c + iX")
    @printf("  %-10s %12s %12s %14s %12s %14s\n", "state", "c(0) N s/m", "c(1k)/c(0)", "c(100k)/c(0)", "X/c at 1k", "X/c at 100k")
    for (name, gap) in (("rest", 13.76e-6), ("gap 5 um", 5e-6), ("gap 1 um", 1e-6), ("contact", 0.0))
        x1, x2 = state_at_gap(d, gap)
        z0 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 0.0, kq), v))
        z1 = zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e3, kq), v)
        z5 = zrt(device_impedance(d, fl, cfg, x1, x2, 2π * 1e5, kq), v)
        @printf("  %-10s %12.4e %12.5f %14.5f %12.5f %14.5f\n", name, z0, real(z1) / z0, real(z5) / z0,
                imag(z1) / real(z1), imag(z5) / real(z5))
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
            guidefontsize=10, tickfontsize=8, legendfontsize=8, titlefontsize=10, legend_background_color=:transparent,
            margin=4Plots.mm)
end

function plot_mobility(fl::Fluid=Fluid())
    lam = 10 .^ range(-2, 2, length=161)
    h = 10e-6
    p = plot(xscale=:log10, xlabel="Womersley parameter |λ| = (h/2)√(ω/ν)", ylabel="G(ω)/G(0)",
             title="Mobility with gas inertia and slip", legend=:bottomleft)
    for (bh, col) in ((0.0, OI.black), (0.05, OI.blue), (0.2, OI.vermillion))
        flb = bh > 0 ? Fluid(lam=bh * h / 1.016191) : fl
        w = (2 .* lam ./ h) .^ 2 .* nu(flb)
        r = [G_mob(h, wi, flb; slip=bh > 0) / G_mob(h, 0.0, flb; slip=bh > 0) for wi in w]
        plot!(p, lam, real.(r), color=col, label="Re, b/h = $(bh)")
        plot!(p, lam, -imag.(r), color=col, ls=:dash, label="−Im, b/h = $(bh)")
    end
    return p
end

function plot_exit(cfg::Config=Config())
    hs = 10 .^ range(log10(0.1), log10(40), length=200)
    p1 = plot(xscale=:log10, xlabel="gap h (μm)", ylabel="f0(h)", title="Lateral exit factor (W = 25 μm)", legend=:topleft)
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
    RE, IM = cfg.faces == :one ? (DF_RE, DF_IM) : (DF2_RE, DF2_IM)
    p2 = plot(xscale=:log10, yscale=:log10, xlabel="gap h (μm)", ylabel="|f/f0 − 1|",
              title=cfg.faces == :one ? "Frequency dependence (one face)" : "Frequency dependence (two faces)", legend=:topleft)
    for (j, col, lab) in ((3, OI.green, "1 kHz"), (5, OI.orange, "10 kHz"), (7, OI.vermillion, "100 kHz"))
        plot!(p2, H_TAB, abs.(complex.(RE[:, j], IM[:, j])), color=col, marker=:circle, ms=3, label=lab)
    end
    return plot(p1, p2, layout=(1, 2), size=(900, 360), left_margin=6Plots.mm, bottom_margin=6Plots.mm)
end

function plot_frequency(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    kq = pocket_for(d, fl, cfg); v = rigid(d)
    F = 10 .^ range(0, 5, length=41)
    p1 = plot(xscale=:log10, yscale=:log10, xlabel="frequency (Hz)", ylabel="|c(f)/c(0) − 1|",
              title="Damping change (dashed: decrease, extrapolated < 100 Hz)", legend=:topleft, ylims=(1e-7, 1))
    p2 = plot(xscale=:log10, xlabel="frequency (Hz)", ylabel="Im Z / Re Z", title="Reactance (+ mass-like, − spring-like)", legend=:topleft)
    for (name, gap, col) in (("rest", 13.76e-6, OI.black), ("gap 5 μm", 5e-6, OI.blue), ("gap 1 μm", 1e-6, OI.green), ("contact", 0.0, OI.vermillion))
        x1, x2 = state_at_gap(d, gap)
        z0 = real(zrt(device_impedance(d, fl, cfg, x1, x2, 0.0, kq), v))
        zs = [zrt(device_impedance(d, fl, cfg, x1, x2, 2π * fh, kq), v) for fh in F]
        dc = real.(zs) ./ z0 .- 1
        plot!(p1, F, [x > 0 ? x : NaN for x in dc], color=col, label=name)
        any(dc .< 0) && plot!(p1, F, [x < 0 ? -x : NaN for x in dc], color=col, ls=:dash, label="")
        plot!(p2, F, imag.(zs) ./ real.(zs), color=col, label=name)
    end
    vline!(p1, [1000], color=:gray, ls=:dot, label="1 kHz")
    vline!(p1, [100], color=:gray, ls=:dashdot, label="100 Hz")
    vline!(p2, [1000], color=:gray, ls=:dot, label="")
    Fa, Fb = F[F .<= 1e4], F[F .>= 3e3]           # reference slopes: Stokes layer at the exit, analytic expansion
    plot!(p1, Fa, 5e-5 .* (Fa ./ 10) .^ 0.5, color=:gray, lw=0.8, label="slope 1/2")
    plot!(p1, Fb, 3e-5 .* (Fb ./ 1e4) .^ 2, color=:gray, lw=0.8, ls=:dashdot, label="slope 2")
    return plot(p1, p2, layout=(1, 2), size=(900, 360), left_margin=6Plots.mm, bottom_margin=6Plots.mm)
end

function plot_gap(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    kq = pocket_for(d, fl, cfg); v = rigid(d)
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
             title="Damping across the stroke" * facetag(cfg), legend=:topright)
    plot!(p, gaps, 1e3 .* getindex.(rows, 1), color=OI.black, ls=:dash, label="Reynolds, quasi-static (original)")
    plot!(p, gaps, 1e3 .* getindex.(rows, 2), color=OI.blue, label="full model, quasi-static")
    plot!(p, gaps, 1e3 .* getindex.(rows, 3), color=OI.green, ls=:dot, label="full model, 10 kHz")
    plot!(p, gaps, 1e3 .* getindex.(rows, 4), color=OI.vermillion, ls=:dashdot, label="full model, 100 kHz")
    return p
end

Config(c::Config; kw...) = Config(; W=c.W, faces=c.faces, slip=c.slip, wom=c.wom, comp=c.comp, exit=c.exit,
                                   exit_omega=c.exit_omega, N=c.N, M=c.M, kinetic=c.kinetic, kw...)

function plot_contact(fl::Fluid=Fluid(), cfg::Config=Config())
    he = 10 .^ range(log10(0.5e-9), log10(1e-6), length=20)
    res = map(he) do h_
        d = Device(heff=h_); kq = pocket_for(d, fl, cfg); v = rigid(d)
        [zrt(device_impedance(d, fl, cfg, gcl(d), gcl(d), 2π * fh, kq), v) for fh in (0.0, 1e4, 1e5)]
    end
    hmin = (he .+ 1e-9) .* 1e9
    d0 = Device()
    Dstar = d0.nb * fl.eta * cfg.W / alpha(d0)^3
    kap = π * kp(fl) / (2 * alpha(d0) * ld(cfg))
    p = plot(xscale=:log10, xlabel="residual gap h_min (nm)", ylabel="contact damping Re Z (mN s/m)",
             title="Contact limit" * facetag(cfg), legend=:topright)
    for (j, col, lab, sty) in ((1, OI.blue, "quasi-static", :solid), (2, OI.green, "10 kHz", :dash), (3, OI.vermillion, "100 kHz", :solid))
        plot!(p, hmin, 1e3 .* real.(getindex.(res, j)), color=col, ls=sty, label=lab)
    end
    hline!(p, [1e3 * Dstar * phi_sat(kap)], color=OI.black, ls=:dash, label="first-order-slip limit D⋆Φsat(κ)")
    return p
end

function plot_pressure(d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config(); w::Real=0.0)
    kq = pocket_for(d, fl, cfg); v = rigid(d)
    plots_ = []
    for (name, gap, ymax) in (("contact", 0.0, 40e-6), ("rest", 13.76e-6, 400e-6))
        x1, x2 = state_at_gap(d, gap)
        _, fields = device_impedance(d, fl, cfg, x1, x2, w, kq; return_fields=true)
        y, zeta, p = pressure_field(fields[1], cfg, v; nz=201)
        keep = y .<= ymax
        hm = heatmap(1e6 .* y[keep], 1e6 .* zeta, transpose(real.(p[keep, :])), color=:viridis,
                     xlabel="y (μm)", ylabel=zlabel(cfg), title="Re p per unit closing speed, $(name)",
                     colorbar_title="Pa s/m")
        push!(plots_, hm)
    end
    return plot(plots_..., layout=(1, 2), size=(950, 360), left_margin=6Plots.mm, bottom_margin=6Plots.mm)
end

# ----------------------------------------------------------------------------- flow-field visualizations
facetag(cfg::Config) = cfg.faces == :one ? "" : " (two open faces)"
zlabel(cfg::Config) = cfg.faces == :one ? "ζ (μm), open at top" : "ζ (μm), both faces open"

"Film pressure, depth-averaged gas flux and dissipation density of the closing film, per unit closing speed."
function film_flow(gap::Real; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config(), nz::Int=81)
    kq = pocket_for(d, fl, cfg)
    v = rigid(d)
    x1, x2 = state_at_gap(d, gap)
    Z, fields = device_impedance(d, fl, cfg, x1, x2, 0.0, kq; return_fields=true)
    F = fields[1]
    y, zeta, pc = pressure_field(F, cfg, v; nz=nz)
    p = real.(pc)
    p .*= sign(p[argmax(abs.(p))])
    G = [real(G_mob(hh, 0.0, fl; kinetic=cfg.kinetic)) for hh in F.h]
    f = [real(f_exit(hh, 0.0, cfg)) for hh in F.h]
    ny = length(y)
    py = zeros(ny, nz)
    pz = zeros(ny, nz)
    for i in 1:ny
        i1, i2 = max(i - 1, 1), min(i + 1, ny)
        py[i, :] = (p[i2, :] .- p[i1, :]) ./ (y[i2] - y[i1])
    end
    for j in 1:nz
        j1, j2 = max(j - 1, 1), min(j + 1, nz)
        pz[:, j] = (p[:, j2] .- p[:, j1]) ./ (zeta[j2] - zeta[j1])
    end
    a = G ./ (12 * fl.eta)
    b = G ./ f ./ (12 * fl.eta)
    return (y=y, zeta=zeta, p=p, qy=-a .* py, qz=-b .* pz, e=a .* py .^ 2 .+ b .* pz .^ 2, c=real(zrt(Z, v)))
end

"Distances from the tip within which the given fractions of the film dissipation occur."
function dissipation_quantiles(F, fr)
    dz = F.zeta[2] - F.zeta[1]
    Ey = vec(sum(F.e, dims=2)) .* dz
    E = zeros(length(F.y))
    for i in 2:length(F.y)
        E[i] = E[i-1] + 0.5 * (Ey[i] + Ey[i-1]) * (F.y[i] - F.y[i-1])
    end
    E ./= E[end]
    return [F.y[min(searchsortedfirst(E, q), length(E))] for q in fr]
end

"Pressure with flux arrows, and dissipation density, at nominal contact and at rest."
function plot_film_flow(; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    panels = Any[]
    for (name, gap, ymax) in (("nominal contact", 0.0, 60e-6), ("rest", 13.76e-6, 400e-6))
        F = film_flow(gap; d=d, fl=fl, cfg=cfg)
        keep = F.y .<= ymax
        yk = 1e6 .* F.y[keep]
        zk = 1e6 .* F.zeta
        P = F.p[keep, :]
        pm = maximum(P)
        h1 = heatmap(yk, zk, transpose(P ./ pm), color=:viridis, clims=(0, 1), legend=false,
                     xlabel="y from the tip (μm)", ylabel=zlabel(cfg), title="pressure / peak and gas flux, $(name)",
                     colorbar_title=@sprintf("p/pmax, pmax = %.3g Pa per m/s", pm), titlefontsize=9)
        Ly, Lz = 1e6 * ymax, 1e6 * cfg.W
        iy = unique([min(searchsortedfirst(yk, t), length(yk)) for t in range(0.04 * Ly, 0.96 * Ly, length=14)])
        iz = round.(Int, range(3, length(zk) - 2, length=8))
        QY = F.qy[keep, :]
        QZ = F.qz[keep, :]
        xs, zs, us, ws = Float64[], Float64[], Float64[], Float64[]
        for i in iy, j in iz
            u, w = QY[i, j] / Ly, QZ[i, j] / Lz
            m = hypot(u, w)
            m > 0 || continue
            push!(xs, yk[i]); push!(zs, zk[j]); push!(us, 0.045 * Ly * u / m); push!(ws, 0.045 * Lz * w / m)
        end
        quiver!(h1, xs, zs, quiver=(us, ws), color=:white, lw=0.8, label="")
        E = F.e[keep, :]
        L = log10.(max.(E ./ maximum(E), 1e-6))
        y50, y90 = 1e6 .* dissipation_quantiles(F, (0.5, 0.9))
        h2 = heatmap(yk, zk, transpose(L), color=:magma, clims=(-4, 0), legend=false, xlabel="y from the tip (μm)",
                     title=@sprintf("dissipation, %s: half within %.3g μm, 90%% within %.3g μm", name, y50, y90),
                     colorbar_title="log10 of dissipation / max", titlefontsize=9)
        for (yy, st) in ((y50, :dash), (y90, :dot))
            yy < Ly && vline!(h2, [yy], color=:white, ls=st, label="")
        end
        push!(panels, h1, h2)
    end
    return plot(panels..., layout=(2, 2), size=(1050, 640), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

"Rigid-translation damping across the stroke: full model and with each ingredient removed (columns:
gap, full, ideal tip vent, no exit factor and ideal vent, first-order slip, full at 100 kHz)."
function stroke_table(; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config(), n::Int=25)
    kq = pocket_for(d, fl, cfg)
    kqs = pocket_for(d, fl, Config(cfg; kinetic=false))
    v = rigid(d)
    gaps = 10 .^ range(log10(0.0515e-6), log10(13.76e-6), length=n)
    T = zeros(n, 6)
    for (i, g) in enumerate(gaps)
        x1, x2 = state_at_gap(d, i == 1 ? 0.0 : g)
        cz(cf, k, w) = real(zrt(device_impedance(d, fl, cf, x1, x2, w, k), v))
        T[i, :] = [g, cz(cfg, kq, 0.0), cz(cfg, Inf, 0.0), cz(Config(cfg; exit=false), Inf, 0.0),
                   cz(Config(cfg; kinetic=false), kqs, 0.0), cz(cfg, kq, 2π * 1e5)]
    end
    return T
end

"What each ingredient does to the damping across the stroke."
function plot_correction_budget(; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    T = stroke_table(d=d, fl=fl, cfg=cfg)
    g = 1e6 .* T[:, 1]
    p1 = plot(g, T[:, 3] ./ T[:, 4], xscale=:log10, color=OI.blue, lw=2, label="exit factor (Stokes cells)",
              xlabel="tip gap (μm), contact at 0.051", ylabel="damping ratio (with / without)",
              title="Large corrections" * facetag(cfg), legend=:topleft)
    plot!(p1, g, T[:, 2] ./ T[:, 3], color=OI.green, lw=2, label="tip pocket (3-D calibrated)")
    plot!(p1, g, T[:, 6] ./ T[:, 2], color=OI.vermillion, lw=2, ls=:dash, label="100 kHz vs quasi-static")
    hline!(p1, [1.0], color=:gray, lw=0.6, label="")
    p2 = plot(g, T[:, 2] ./ T[:, 5], xscale=:log10, color=OI.purple, lw=2, label="kinetic (BGK) vs first-order slip",
              xlabel="tip gap (μm), contact at 0.051", title="Rarefaction", ylims=(0.92, 1.01), legend=:topleft)
    hline!(p2, [1.0], color=:gray, lw=0.6, label="")
    return plot(p1, p2, layout=(1, 2), size=(1000, 380), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

"Contact damping against residual gap (nm) for the four combinations of rarefaction model and sealing law."
function contact_table(; fl::Fluid=Fluid(), cfg::Config=Config(), n::Int=26)
    he = 10 .^ range(log10(0.5e-9), log10(1e-6), length=n)
    T = zeros(n, 5)
    for (i, h_) in enumerate(he)
        dv = Device(heff=h_)
        dc = Device(heff=h_, seal_centred=true)
        k1 = pocket_for(dv, fl, cfg)
        k0 = pocket_for(dv, fl, Config(cfg; kinetic=false))
        v = rigid(dv)
        zc(dd, cf, k) = 1e3 * real(zrt(device_impedance(dd, fl, cf, gcl(dd), gcl(dd), 0.0, k), v))
        T[i, :] = [1e9 * (h_ + dv.eps / 2), zc(dv, cfg, k1), zc(dv, Config(cfg; kinetic=false), k0),
                   zc(dc, cfg, k1), zc(dc, Config(cfg; kinetic=false), k0)]
    end
    return T
end

"Contact damping: kinetic factor against first-order slip, vented tip against the centred sealing law."
function plot_contact_comparison(; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config())
    T = contact_table(fl=fl, cfg=cfg)
    Dstar = d.nb * fl.eta * cfg.W / alpha(d)^3
    kap = π * kp(fl) / (2 * alpha(d) * ld(cfg))
    sat = 1e3 * Dstar * phi_sat(kap)
    p = plot(T[:, 1], T[:, 2], xscale=:log10, color=OI.blue, lw=2.2, label="kinetic (BGK), tip vented [model]",
             xlabel="residual gap h_min (nm)", ylabel="contact damping (mN s/m)",
             title="Contact damping: rarefaction and tip end condition" * facetag(cfg), legend=:topright, size=(760, 440))
    plot!(p, T[:, 1], T[:, 3], color=OI.purple, lw=1.6, label="first-order slip, tip vented")
    plot!(p, T[:, 1], T[:, 4], color=OI.blue, lw=1.2, ls=:dash, label="kinetic, sealing law centred on contact")
    plot!(p, T[:, 1], T[:, 5], color=OI.purple, lw=1.2, ls=:dash, label="first-order slip, centred sealing")
    hline!(p, [sat], color=:black, ls=:dot, label=@sprintf("first-order-slip limit %.2f mN s/m", sat))
    vline!(p, [51.0], color=:gray, lw=0.7, label="device, 51 nm")
    return p
end

"Linearized-BGK Poiseuille flow rate and the kinetic flow factor."
function plot_kinetic(; fl::Fluid=Fluid())
    dd = 10 .^ range(-4, 3, length=400)
    GP = gp_bgk.(dd)
    p1 = plot(dd, GP, xscale=:log10, yscale=:log10, color=OI.blue, lw=2, label="linearized BGK",
              xlabel="δ (gap / free path)", ylabel="reduced flow rate G_P", ylims=(1, 300), legend=:topleft,
              title="Plane Poiseuille flow of a rarefied gas")
    plot!(p1, dd, dd ./ 6 .+ fl.sig_p, color=:black, ls=:dash, label="first-order slip, δ/6 + σp")
    sm = dd .< 0.3
    plot!(p1, dd[sm], 0.35887 .+ log.(1 ./ dd[sm]) ./ sqrt(π), color=OI.orange, ls=:dashdot, label="free molecular")
    i = argmin(GP)
    scatter!(p1, [dd[i]], [GP[i]], color=OI.vermillion, ms=4, label="Knudsen minimum")
    p2 = plot(dd, 6 .* GP ./ (dd .+ 6 * fl.sig_p), xscale=:log10, color=OI.purple, lw=2, ylims=(0.9, 2.7), label="",
              xlabel="δ (gap / free path)", ylabel="R = BGK / first-order slip flow", title="Kinetic flow factor")
    hline!(p2, [1.0], color=:gray, lw=0.6, label="")
    vline!(p2, [51e-9, 1e-6, 13.76e-6] ./ fl.lam, color=:gray, ls=:dot, label="contact, 1 μm, rest")
    return plot(p1, p2, layout=(1, 2), size=(1000, 380), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

"Animation: oscillatory channel-flow (Womersley) profiles across the gap over one cycle."
function animate_womersley(outfile::AbstractString; lams=(0.5, 2.0, 6.0), nframes::Int=36, fps::Int=12)
    zs = collect(range(-0.5, 0.5, length=201))
    prof = [2 .* (1 .- cosh.(2 * L * cis(π / 4) .* zs) ./ cosh(L * cis(π / 4))) ./ (L * cis(π / 4))^2 for L in lams]
    anim = @animate for i in 0:nframes-1
        ph = 2π * i / nframes
        ps = Any[]
        for (k, (L, u)) in enumerate(zip(lams, prof))
            env = abs.(u)
            ttl = @sprintf("|λ| = %g%s", L, L < 1 ? ", quasi-steady" : ", inertial")
            k == 1 && (ttl *= @sprintf("   (phase %.0f deg)", rad2deg(ph)))
            p = plot(Shape(vcat(-env, reverse(env)), vcat(zs, reverse(zs))), fillcolor=:gray90, linecolor=:gray90,
                     label="", xlims=(-1.05, 1.05), ylims=(-0.5, 0.5), xlabel="u / u_Poiseuille",
                     ylabel=(k == 1 ? "z/h across the gap" : ""), title=ttl, titlefontsize=9)
            plot!(p, 1 .- 4 .* zs .^ 2, zs, color=:gray, ls=:dot, label="")
            plot!(p, real.(u .* cis(ph)), zs, color=OI.blue, lw=2, label="")
            dl = 1 / (sqrt(2) * L)
            dl < 0.5 && hline!(p, [0.5 - dl, dl - 0.5], color=OI.orange, ls=:dash, label="")
            push!(ps, p)
        end
        plot(ps..., layout=(1, 3), size=(1000, 360), left_margin=4Plots.mm, bottom_margin=5Plots.mm)
    end
    gif(anim, outfile, fps=fps)
end

"Animation: film pressure as the electrode closes from rest to contact, with the damping curve."
function animate_closing(outfile::AbstractString; d::Device=Device(), fl::Fluid=Fluid(), cfg::Config=Config(),
                         nframes::Int=34, fps::Int=6)
    gaps = vcat(collect(10 .^ range(log10(13.76e-6), log10(0.0515e-6), length=nframes)), zeros(6))
    T = stroke_table(d=d, fl=fl, cfg=cfg)
    ticks = (collect(-1.0:1.0:2.0), ["0.1", "1", "10", "100"])
    anim = @animate for g in gaps
        F = film_flow(g; d=d, fl=fl, cfg=cfg, nz=61)
        ly = log10.(1e6 .* F.y[2:end])
        P = F.p[2:end, :]
        pm = maximum(P)
        lab = g == 0 ? "contact (51 nm)" : @sprintf("%.3g μm", 1e6 * g)
        h1 = heatmap(ly, 1e6 .* F.zeta, transpose(P ./ pm), color=:viridis, clims=(0, 1), legend=false,
                     xlims=(log10(0.05), log10(400.0)), xticks=ticks, xlabel="y from the tip (μm)",
                     ylabel=zlabel(cfg), titlefontsize=9,
                     title=@sprintf("film pressure / peak, tip gap %s, peak %.3g Pa per m/s", lab, pm))
        h2 = plot(1e6 .* T[:, 1], 1e3 .* T[:, 2], xscale=:log10, yscale=:log10, color=OI.blue, lw=1.6, label="",
                  xlabel="tip gap (μm)", ylabel="damping (mN s/m)", title="rigid translation", titlefontsize=9)
        scatter!(h2, [1e6 * max(g, 0.0515e-6)], [1e3 * F.c], color=OI.vermillion, ms=6, label="")
        plot(h1, h2, layout=@layout([a{0.68w} b]), size=(1050, 400), left_margin=5Plots.mm, bottom_margin=6Plots.mm)
    end
    gif(anim, outfile, fps=fps)
end

"Write all figures and animations: the as-built configuration (bottom face on the PCB, top face open) to outdir, and,
with two_faces = true, the configuration with both faces vented to joinpath(outdir, \"two_faces\")."
function make_plots(outdir::AbstractString="figs"; two_faces::Bool=true)
    setup_style()
    d, fl = Device(), Fluid()
    sets = Any[(outdir, Config(), true)]
    two_faces && push!(sets, (joinpath(outdir, "two_faces"), Config(faces=:two), false))
    for (dir, cfg, full) in sets
        mkpath(dir)
        figs = Any[("exit_factor", () -> plot_exit(cfg)), ("frequency", () -> plot_frequency(d, fl, cfg)),
                   ("damping_vs_gap", () -> plot_gap(d, fl, cfg)), ("contact_limit", () -> plot_contact(fl, cfg)),
                   ("pressure_field", () -> plot_pressure(d, fl, cfg)),
                   ("film_flow_dissipation", () -> plot_film_flow(d=d, fl=fl, cfg=cfg)),
                   ("correction_budget", () -> plot_correction_budget(d=d, fl=fl, cfg=cfg)),
                   ("contact_comparison", () -> plot_contact_comparison(d=d, fl=fl, cfg=cfg))]
        full && append!(figs, Any[("mobility", () -> plot_mobility(fl)), ("kinetic_flow_factor", () -> plot_kinetic(fl=fl))])
        for (name, fn) in figs
            p = fn()
            savefig(p, joinpath(dir, name * ".pdf"))
            savefig(p, joinpath(dir, name * ".png"))
            println("wrote ", joinpath(dir, name), ".{pdf,png}")
        end
        full && animate_womersley(joinpath(dir, "womersley_cycle.gif"))
        animate_closing(joinpath(dir, "closing_stroke.gif"); d=d, fl=fl, cfg=cfg)
        println("wrote the animations in ", dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    selftest()
    run_case_study()
    run_case_study(faces=:two)
    make_plots()
end