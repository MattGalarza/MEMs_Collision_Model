/-
  SqueezeFilm.Core — algebraic skeleton of
  "Squeeze films between doubly tapered plates closing to contact".

  Core Lean only (no Mathlib): every statement holds in any field of characteristic zero,
  so it applies to ℝ (see `SqueezeFilm.Analysis.RealInstances`).
  Proved with `grind`'s commutative-ring/field normalizer. Calculus steps (derivatives,
  integrals, series) are formalized separately in `SqueezeFilm.Analysis` (needs Mathlib).
-/
namespace SqueezeFilm.Core

variable {α : Type} [Lean.Grind.Field α] [Lean.Grind.IsCharP α 0]

/-! ## Field helpers (grind does not infer these on its own) -/

theorem mul_ne_zero' {x y : α} (hx : x ≠ 0) (hy : y ≠ 0) : x * y ≠ 0 := by
  intro h
  have : x * y * y⁻¹ = x := by grind
  grind

theorem inv_ne_zero' {x : α} (hx : x ≠ 0) : x⁻¹ ≠ 0 := by
  intro h
  have : x * x⁻¹ = 1 := by grind
  grind

theorem div_ne_zero' {x y : α} (hx : x ≠ 0) (hy : y ≠ 0) : x / y ≠ 0 := by
  have := mul_ne_zero' hx (inv_ne_zero' hy)
  grind

theorem div_eq_of_eq_mul' {x y d : α} (hd : d ≠ 0) (e : x = y * d) : x / d = y := by grind

theorem two_ne_zero' : (2 : α) ≠ 0 := by grind

theorem pow2_ne_zero' {x : α} (hx : x ≠ 0) : x ^ 2 ≠ 0 := by
  have := mul_ne_zero' hx hx; grind

theorem pow3_ne_zero' {x : α} (hx : x ≠ 0) : x ^ 3 ≠ 0 := by
  have := mul_ne_zero' (mul_ne_zero' hx hx) hx; grind

/-! ## §2.4 Scaling: every coefficient of the reduced momentum and mass balances -/

/-- z-momentum in-plane viscous term over the pressure-gradient scale: `ε²`. -/
theorem zmom_viscous (μ V ld hs : α) (hl2 : ld ^ 2 ≠ 0) (hs2 : hs ^ 2 ≠ 0) (hs3 : hs ^ 3 ≠ 0)
    (hz : μ * V * ld ^ 2 / hs ^ 3 / hs ≠ 0) :
    (μ * V / hs ^ 2) / (μ * V * ld ^ 2 / hs ^ 3 / hs) = hs ^ 2 / ld ^ 2 := by
  have e : μ * V * ld ^ 2 / hs ^ 3 / hs = μ * V * ld ^ 2 / hs ^ 4 := by grind
  rw [e]
  have h4 : hs ^ 4 ≠ 0 := by grind
  exact (show (μ * V / hs ^ 2) / (μ * V * ld ^ 2 / hs ^ 4) = hs ^ 2 / ld ^ 2 by grind)

/-- With `U = V ℓ_d / h*` and `P* = μ V ℓ_d² / h*³`, dividing each x-momentum term by `μU/h*²`
gives the coefficients claimed in (2.7); dividing the z-momentum terms by `P*/h*` gives (2.11);
the storage/drainage ratio of continuity is `σ/12` (2.12). -/
theorem scaling_coefficients (ρ μ ω V ld hs pa : α) (hρ : ρ ≠ 0)
    (hμ : μ ≠ 0) (hV : V ≠ 0) (hld : ld ≠ 0) (hhs : hs ≠ 0) (hpa : pa ≠ 0) :
    let U := V * ld / hs
    let Ps := μ * V * ld ^ 2 / hs ^ 3
    -- x-momentum, divided by μU/h*²
    (ρ * ω * U) / (μ * U / hs ^ 2) = ρ * ω * hs ^ 2 / μ ∧            -- unsteady: Re_ω
    (ρ * U ^ 2 / ld) / (μ * U / hs ^ 2) = ρ * V * hs / μ ∧          -- u ∂ₓu: Re*
    (ρ * V * U / hs) / (μ * U / hs ^ 2) = ρ * V * hs / μ ∧          -- w ∂_z u: the same Re*
    (Ps / ld) / (μ * U / hs ^ 2) = 1 ∧                              -- pressure gradient
    (μ * U / ld ^ 2) / (μ * U / hs ^ 2) = (hs ^ 2 / ld ^ 2) ∧           -- in-plane viscous: ε²
    -- z-momentum, divided by P*/h*
    (μ * V / hs ^ 2) / (Ps / hs) = (hs ^ 2 / ld ^ 2) ∧                  -- ε²
    (ρ * ω * V) / (Ps / hs) = (hs ^ 2 / ld ^ 2) * (ρ * ω * hs ^ 2 / μ) ∧ -- ε² Re_ω
    -- continuity: storage over drainage
    (ω * ρ * Ps / pa) / (ρ * V / hs) = μ * ω * ld ^ 2 / (pa * hs ^ 2) := by
  intro U Ps
  have hs2 := pow2_ne_zero' hhs
  have hs3 := pow3_ne_zero' hhs
  have hl2 := pow2_ne_zero' hld
  have hU : U ≠ 0 := div_ne_zero' (mul_ne_zero' hV hld) hhs
  have hPs : Ps ≠ 0 := div_ne_zero' (mul_ne_zero' (mul_ne_zero' hμ hV) hl2) hs3
  have hx : μ * U / hs ^ 2 ≠ 0 := div_ne_zero' (mul_ne_zero' hμ hU) hs2
  have hz : Ps / hs ≠ 0 := div_ne_zero' hPs hhs
  have hc : ρ * V / hs ≠ 0 := div_ne_zero' (mul_ne_zero' hρ hV) hhs
  have hpa2 : pa * hs ^ 2 ≠ 0 := mul_ne_zero' hpa hs2
  refine ⟨by grind, by grind, by grind, by grind, by grind, ?_, by grind, by grind⟩
  -- z-momentum viscous term, proved separately (grind needs the intermediate form)
  exact zmom_viscous μ V ld hs hl2 hs2 hs3 hz

/-- `Re*/Re_ω` is the amplitude-to-gap ratio (2.10). -/
theorem reynolds_ratio (ρ μ ω V hs : α) (hμ : μ ≠ 0) (hω : ω ≠ 0) (hhs : hs ≠ 0) :
    (ρ * V * hs / μ) = (ρ * ω * hs ^ 2 / μ) * (V / (ω * hs)) := by grind

/-- Stokes' hypothesis `λ_b = -2μ/3` leaves the compressive coefficient `μ + λ_b = μ/3` (2.4). -/
theorem stokes_hypothesis (μ : α) : μ + (-(2 : α) / 3) * μ = μ / 3 := by grind

/-! ## §2.5 Flow across the gap -/

/-- The profile (2.14) satisfies the first-order Maxwell slip conditions at both walls.
`du` is its z-derivative (power rule; proved formally in `Analysis.Profile`). -/
theorem profile_slip_conditions (P μ b hb ha Ua Ub : α) (hd : (ha - hb) + 2 * b ≠ 0) :
    let h := ha - hb
    let u := fun z : α => P / (2 * μ) * ((z - hb) * (z - ha) - b * h) + Ub + (Ua - Ub) * (z - hb + b) / (h + 2 * b)
    let du := fun z : α => P / (2 * μ) * (2 * z - ha - hb) + (Ua - Ub) / (h + 2 * b)
    u hb - Ub = b * du hb ∧ u ha - Ua = -b * du ha := by
  intro h u du
  constructor <;> grind

/-- The antiderivative `Ianti` of the profile (in `s = z - h_b`) has derivative equal to the profile. -/
theorem flux_antiderivative_matches (P μ b hb h Ua Ub s : α) :
    let u := fun z : α => P / (2 * μ) * ((z - hb) * (z - (hb + h)) - b * h) + Ub + (Ua - Ub) * (z - hb + b) / (h + 2 * b)
    -- derivative of  P/(2μ)(s³/3 - h s²/2 - b h s) + Ub s + (Ua-Ub)(s²/2 + b s)/(h+2b)
    P / (2 * μ) * (s ^ 2 - h * s - b * h) + Ub + (Ua - Ub) * (s + b) / (h + 2 * b) = u (hb + s) := by
  intro u
  have e1 : (hb + s - hb) * (hb + s - (hb + h)) = s ^ 2 - h * s := by grind
  have e2 : hb + s - hb + b = s + b := by grind
  simp only [u, e1, e2]

/-- Flux law (2.15): the antiderivative evaluated across the gap gives
`q = -G/(12μ) ∂ₓp + h (U_a + U_b)/2` with `G = h²(h + 6b)`; slip leaves the Couette flux unchanged. -/
theorem poiseuille_part (P μ b h : α) (hμ : μ ≠ 0) :
    P / (2 * μ) * (h ^ 3 / 3 - h * h ^ 2 / 2 - b * h * h) = -(h ^ 2 * (h + 6 * b)) / (12 * μ) * P := by grind

theorem couette_part (Ua Ub h b : α) (hd : h + 2 * b ≠ 0) :
    Ub * h + (Ua - Ub) * (h ^ 2 / 2 + b * h) / (h + 2 * b) = h * (Ua + Ub) / 2 := by grind

theorem flux_law (P μ b h Ua Ub : α) (hμ : μ ≠ 0) (hd : h + 2 * b ≠ 0) :
    let Ianti := fun s : α => P / (2 * μ) * (s ^ 3 / 3 - h * s ^ 2 / 2 - b * h * s) + Ub * s + (Ua - Ub) * (s ^ 2 / 2 + b * s) / (h + 2 * b)
    Ianti h - Ianti 0 = -(h ^ 2 * (h + 6 * b)) / (12 * μ) * P + h * (Ua + Ub) / 2 := by
  intro Ianti
  have a1 := poiseuille_part P μ b h hμ
  have a2 := couette_part Ua Ub h b hd
  have z0 : Ianti 0 = 0 := by simp only [Ianti] <;> grind
  have zh : Ianti h = P / (2 * μ) * (h ^ 3 / 3 - h * h ^ 2 / 2 - b * h * h) + Ub * h
      + (Ua - Ub) * (h ^ 2 / 2 + b * h) / (h + 2 * b) := by simp only [Ianti] <;> grind
  rw [z0, zh]; grind

/-- Local dissipation (2.16): viscous part `μ∫(∂_z u)²` (antiderivative `(2s-h)³/6` of `(2s-h)²`)
plus wall-slip part `(μ/b)(u_{s,a}² + u_{s,b}²)` equals `G (∂ₓp)²/(12μ)`. -/
theorem dissipation_identity (P μ b h : α) (hμ : μ ≠ 0) (hb : b ≠ 0) :
    let viscous := μ * (P / (2 * μ)) ^ 2 * ((2 * h - h) ^ 3 / 6 - (2 * 0 - h) ^ 3 / 6)
    let us := P / (2 * μ) * (-(b * h))          -- slip velocity at either wall
    viscous + μ / b * (us ^ 2 + us ^ 2) = h ^ 2 * (h + 6 * b) * P ^ 2 / (12 * μ) := by
  intro viscous us; grind

/-! ## §2.6 Mass balance across the gap -/

/-- Wall kinematics: a material wall surface `z = h_a(x,y,t)` moving with the solid, plus impermeability
`(u - u_wall)·n = 0`, give the kinematic condition used in (2.17) for the gas, slip included. -/
theorem kinematic_condition (w W u U v V ht hx hy : α)
    (hwall : W = ht + U * hx + V * hy) (himperm : (w - W) - (u - U) * hx - (v - V) * hy = 0) :
    w - ht - u * hx - v * hy = 0 := by grind

/-- Collecting the Leibniz terms (2.17): with both wall brackets zero the integrated continuity
equation reduces to `∂ₜ(ρh) + ∇·(ρq) = 0` (2.18). -/
theorem leibniz_collection (dtRhoH dxRhoQx dyRhoQy ρ wa wb hta htb ua va ub vb hxa hya hxb hyb : α)
    (ha : wa - hta - ua * hxa - va * hya = 0) (hb : wb - htb - ub * hxb - vb * hyb = 0) :
    (dtRhoH - ρ * hta + ρ * htb) + (dxRhoQx - ρ * ua * hxa + ρ * ub * hxb)
      + (dyRhoQy - ρ * va * hya + ρ * vb * hyb) + (ρ * wa - ρ * wb)
      = dtRhoH + dxRhoQx + dyRhoQy := by grind

/-! ## §4 Strip laws -/

/-- Moments behind the small-taper expansion (4.2): `∫₀¹ η² = 1/3`, `∫₀¹ η²(η-½) = 1/12`,
`∫₀¹ η²(η-½)² = 1/30` (antiderivatives evaluated), and the resulting one-vent coefficients. -/
theorem small_taper_moments (δ : α) :
    ((1 : α) ^ 3 / 3 - 0) = 1 / 3 ∧
    ((1 : α) ^ 4 / 4 - 1 ^ 3 / 6) = 1 / 12 ∧
    ((1 : α) ^ 5 / 5 - 1 ^ 4 / 4 + 1 ^ 3 / 12) = 1 / 30 ∧
    3 * ((1 : α) / 3 - 3 * δ * (1 / 12) + 6 * δ ^ 2 * (1 / 30)) = 1 - 3 / 4 * δ + 3 / 5 * δ ^ 2 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> grind

/-- One-vent strip law for a uniform gap: `12μ ∫₀^W y²/G dy = 4μW³/G`. -/
theorem strip_uniform (μ W G : α) (hG : G ≠ 0) : 12 * μ * (W ^ 3 / 3) / G = 4 * μ * W ^ 3 / G := by grind

/-- Outer integral of the class-B strip limit: `∫₀¹ η² · 1/(2η²) dη = 1/2`, so `Φ_B → 6(α/β)²` (4.4). -/
theorem strip_limit_B (η : α) (hη : η ≠ 0) : η ^ 2 * (1 / (2 * η ^ 2)) = 1 / 2 := by grind

/-! ## §5 Corner solutions -/

/-- Algebraic core of the corner ODE `L[f] = (cos³ψ f')' - 2cos³ψ f` (Appendix B), with
`c = cos ψ`, `s = sin ψ`, `s² + c² = 1`; the derivatives used are formalized in `Analysis.Corner`.
Particular solution `1/c`: `cos³ f' = s c`, `(s c)' = c² - s²`, so `L = -1`.
Homogeneous `s/c²`: `cos³ f' = 1 + s²`, derivative `2 s c`, so `L = 0`.
Homogeneous `1/c²`: `cos³ f' = 2 s`, derivative `2 c`, so `L = 0`. -/
theorem corner_ode_algebra (s c : α) (hc : c ≠ 0) (hpyth : s ^ 2 + c ^ 2 = 1) :
    (c ^ 2 - s ^ 2) - 2 * c ^ 3 * (1 / c) = -1 ∧
    (2 * s * c) - 2 * c ^ 3 * (s / c ^ 2) = 0 ∧
    (2 * c) - 2 * c ^ 3 * (1 / c ^ 2) = 0 := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- The derivative formulas behind `corner_ode_algebra`: `c³·(s/c²) = s c`, `c³·((c² + 2s²)/c³) = 1 + s²`,
`c³·(2s/c³) = 2s`. -/
theorem corner_flux_forms (s c : α) (hc : c ≠ 0) (hpyth : s ^ 2 + c ^ 2 = 1) :
    c ^ 3 * (s / c ^ 2) = s * c ∧ c ^ 3 * ((c ^ 2 + 2 * s ^ 2) / c ^ 3) = 1 + s ^ 2 ∧
    c ^ 3 * (2 * s / c ^ 3) = 2 * s := by
  refine ⟨?_, ?_, ?_⟩ <;> grind

/-- Proposition 5.1, sealed corner: the coefficients (5.2) satisfy both edge conditions
(no flux through the sealed tip and through the blocked narrow face). -/
theorem corner_sealed (A al be R : α) (hR : R ^ 2 = al ^ 2 + be ^ 2)
    (h1 : al + be ≠ 0) (h2 : R ^ 2 + al * be ≠ 0) (h3 : al ≠ 0) (h4 : R ≠ 0) :
    let c1 := A * al * be * (al - be) / ((al + be) * (R ^ 2 + al * be))
    let c2 := -(A * al * be + c1 * (R ^ 2 + al ^ 2)) / (2 * al * R)
    A * al * be + c1 * (R ^ 2 + al ^ 2) + 2 * c2 * al * R = 0 ∧           -- sealed tip
    -A * al * be + c1 * (R ^ 2 + be ^ 2) - 2 * c2 * be * R = 0 := by        -- blocked face
  intro c1 c2
  have hd1 : (al + be) * (R ^ 2 + al * be) ≠ 0 := mul_ne_zero' h1 h2
  have hd2 : 2 * al * R ≠ 0 := mul_ne_zero' (mul_ne_zero' two_ne_zero' h3) h4
  have e1 : c1 * ((al + be) * (R ^ 2 + al * be)) = A * al * be * (al - be) := by simp only [c1]; grind
  have e2 : c2 * (2 * al * R) = -(A * al * be + c1 * (R ^ 2 + al ^ 2)) := by simp only [c2]; grind
  clear_value c1 c2
  constructor <;> grind

/-- Proposition 5.1, open corner: the coefficients (5.3) satisfy the sealed-tip condition and `p = 0`
on the open narrow face. -/
theorem corner_open (A al be R : α) (hR : R ^ 2 = al ^ 2 + be ^ 2)
    (h1 : R ^ 2 + al ^ 2 + 2 * al * be ≠ 0) (h4 : R ≠ 0) :
    let c1 := A * al * (2 * al - be) / (R ^ 2 + al ^ 2 + 2 * al * be)
    let c2 := (c1 * be - A * al) / R
    A * al * be + c1 * (R ^ 2 + al ^ 2) + 2 * c2 * al * R = 0 ∧             -- sealed tip
    A * al - c1 * be + c2 * R = 0 := by                                      -- open face
  intro c1 c2
  have e1 : c1 * (R ^ 2 + al ^ 2 + 2 * al * be) = A * al * (2 * al - be) := by simp only [c1]; grind
  have e2 : c2 * R = c1 * be - A * al := by simp only [c2]; grind
  clear_value c1 c2
  constructor <;> grind

/-- Proposition 5.2: solvability of the slip corner, `2 A_s k_p ∫₀^{π/2} s² dθ = 12μ|ḣ| π/2` with
`∫ s² = πR²/4 + αβ` (formalized in `Analysis.Corner`), gives `A_s = 12πμ|ḣ|/(k_p(πR² + 4αβ))`. -/
theorem slip_corner_constant (As kp μ hd π R al be : α) (hk : kp ≠ 0) (hden : π * R ^ 2 + 4 * al * be ≠ 0)
    (hsolv : 2 * As * kp * (π * R ^ 2 / 4 + al * be) = 12 * μ * hd * (π / 2)) :
    As = 12 * π * μ * hd / (kp * (π * R ^ 2 + 4 * al * be)) := by grind

/-! ## §5.4 Finiteness -/

/-- The bracket of the finiteness proof (Appendix C.3) minus `1/2` is `-(t-2)²/2`
(with `t = g/(αX+g)`), hence at most zero. -/
theorem finiteness_bracket (t : α) :
    -(3 : α) / 2 + 2 * t - t ^ 2 / 2 - 1 / 2 = -((t - 2) ^ 2) / 2 := by grind

/-! ## §6 Contact damping -/

/-- Proposition 6.5: at contact `G = (αx+βy)²(αx+βy+k_p)` is homogeneous of degree three in `(α, β, k_p)`. -/
theorem homogeneity (c al be kp x y : α) :
    (c * al * x + c * be * y) ^ 2 * ((c * al * x + c * be * y) + c * kp)
      = c ^ 3 * ((al * x + be * y) ^ 2 * ((al * x + be * y) + kp)) := by grind

/-- Partial fractions behind the strip function `F(a,k)` (6.7), with `w = u + a`:
`(w-a)²/(w²(w+k)) = A/w + B/w² + C/(w+k)` with `A = -a(2k+a)/k²`, `B = a²/k`, `C = (a+k)²/k²`,
stated with the common denominator `w²(w+k)` cleared (equivalent for `w ≠ 0`, `w + k ≠ 0`);
`A + C = 1` makes the coefficient of `ln U` equal to one. -/
theorem slice_partial_fractions (a k w : α) (hk : k ≠ 0) :
    let A := -(a * (2 * k + a)) / k ^ 2
    let B := a ^ 2 / k
    let C := (a + k) ^ 2 / k ^ 2
    A * w * (w + k) + B * w ^ 0 * (w + k) + C * w ^ 2 = (w - a) ^ 2 ∧ A + C = 1 := by
  intro A B C
  have hk2 : k ^ 2 ≠ 0 := by have := mul_ne_zero' hk hk; grind
  have poly : -(a * (2 * k + a)) * w * (w + k) + a ^ 2 * k * (w + k) + (a + k) ^ 2 * w ^ 2 = k ^ 2 * (w - a) ^ 2 := by grind
  have eA : A * w * (w + k) = -(a * (2 * k + a)) * w * (w + k) / k ^ 2 := by simp only [A]; grind
  have eB : B * w ^ 0 * (w + k) = a ^ 2 * k * (w + k) / k ^ 2 := by simp only [B]; grind
  have eC : C * w ^ 2 = (a + k) ^ 2 * w ^ 2 / k ^ 2 := by simp only [C]; grind
  have sumC : A + C = 1 := by simp only [A, C]; grind
  refine ⟨?_, sumC⟩
  rw [eA, eB, eC]
  have comb : -(a * (2 * k + a)) * w * (w + k) / k ^ 2 + a ^ 2 * k * (w + k) / k ^ 2 + (a + k) ^ 2 * w ^ 2 / k ^ 2
      = (-(a * (2 * k + a)) * w * (w + k) + a ^ 2 * k * (w + k) + (a + k) ^ 2 * w ^ 2) / k ^ 2 := by grind
  rw [comb, poly]
  exact div_eq_of_eq_mul' hk2 (by grind)

/-- Small-taper constant of (6.9): with the class-A constant `c₀`, the geometric-mean residual gap
`βW/e` gives `C_X = ln(2eℓ_d/(πW)) + c₀`; in logarithmic variables this is linear bookkeeping:
`ln(2αℓ_d/(π h_eff)) = ln(α/β) + ln(2eℓ_d/(πW))` when `ln h_eff = ln β + ln W - 1`. -/
theorem small_taper_constant (Lα Lβ Lld LW L2 Lπ : α) :
    let LhEff := Lβ + LW - 1                       -- ln(βW/e)
    (L2 + Lα + Lld - Lπ - LhEff) = (Lα - Lβ) + (L2 + 1 + Lld - Lπ - LW) := by
  intro LhEff; grind

/-! ## Appendix D — rigid tilted plates -/

/-- Moy et al.'s parametrization in mean-gap coordinates: `h̄ = h₀(1 + f_T/2)`, and the half
fractional tilt `τ = (h_max - h_min)/(h_max + h_min)` equals `f_T/(2 + f_T)`. -/
theorem moy_mapping (h0 fT : α) (h0ne : h0 ≠ 0) (hden : 2 + fT ≠ 0) :
    let hL := h0 * (1 + fT)
    ((hL - h0) / (hL + h0) = fT / (2 + fT)) ∧ ((h0 + hL) / 2 = h0 * (1 + fT / 2)) := by
  intro hL
  constructor <;> grind

end SqueezeFilm.Core
