import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.Linarith
import SqueezeFilm.Analysis.Tactics

/-!
# The strip law (§4, Proposition 4.1)

For a film whose gap does not vary along `x`, the width-wise problem is `(G p')' = -12μ` with no
flux at the blocked face `y = 0` and `p(W) = 0`. Its solution is `p(y) = 12μ ∫_y^W s/G(s) ds`, and
the damping per unit length is `∫₀^W p = 12μ ∫₀^W y²/G` — the order exchange is integration by parts.
-/

namespace SqueezeFilm.Strip

/-- Order exchange: `∫₀^W (∫_y^W f) dy = ∫₀^W y f(y) dy` for continuous `f`. -/
theorem integral_tail_eq (f : ℝ → ℝ) (hf : Continuous f) (W : ℝ) :
    ∫ y in (0:ℝ)..W, (∫ s in y..W, f s) = ∫ y in (0:ℝ)..W, y * f y := by
  have hF : ∀ y ∈ Set.uIcc (0:ℝ) W, HasDerivAt (fun y => ∫ s in y..W, f s) (-f y) y := fun y _ =>
    intervalIntegral.integral_hasDerivAt_left (hf.intervalIntegrable _ _)
      (hf.stronglyMeasurableAtFilter _ _) hf.continuousAt
  have hid : ∀ y ∈ Set.uIcc (0:ℝ) W, HasDerivAt (fun y : ℝ => y) 1 y := fun y _ => hasDerivAt_id'
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul hid hF
    intervalIntegrable_const (hf.neg.intervalIntegrable _ _)
  simp only [intervalIntegral.integral_same, mul_zero, zero_mul, sub_zero, one_mul, zero_sub] at hparts
  have hneg : ∫ y in (0:ℝ)..W, y * -f y = -∫ y in (0:ℝ)..W, y * f y := by
    simp only [mul_neg, intervalIntegral.integral_neg]
  linarith

/-- **Strip law** (one vent): the pressure `p(y) = 12μ ∫_y^W s/G(s) ds` has flux `G p' = -12μ y`
(so `(G p')' = -12μ`, and no flux at `y = 0`), vanishes at the vent, and integrates to
`12μ ∫₀^W y²/G`. -/
theorem strip_law (μ W : ℝ) (G : ℝ → ℝ) (hG : Continuous G) (hpos : ∀ y, 0 < G y) :
    let p := fun y : ℝ => 12 * μ * ∫ s in y..W, s / G s
    (∀ y, HasDerivAt p (12 * μ * -(y / G y)) y) ∧
    (∀ y, G y * (12 * μ * -(y / G y)) = -12 * μ * y) ∧
    p W = 0 ∧
    ∫ y in (0:ℝ)..W, p y = 12 * μ * ∫ y in (0:ℝ)..W, y ^ 2 / G y := by
  intro p
  have hf : Continuous (fun s => s / G s) := continuous_id.div hG (fun y => (hpos y).ne')
  refine ⟨fun y => ?_, fun y => ?_, ?_, ?_⟩
  · exact (intervalIntegral.integral_hasDerivAt_left (hf.intervalIntegrable _ _)
      (hf.stronglyMeasurableAtFilter _ _) hf.continuousAt).const_mul (12 * μ)
  · have := (hpos y).ne'
    field_ring
  · simp [p]
  · simp only [p]
    rw [intervalIntegral.integral_const_mul, integral_tail_eq _ hf W]
    congr 1
    apply intervalIntegral.integral_congr
    intro y _
    have := (hpos y).ne'
    simp only
    field_ring

end SqueezeFilm.Strip
