import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.DefEqTransformations
import SqueezeFilm.Analysis.Polynomial

/-!
# Linear taper and the finiteness bound (§4.2 eq. (4.3); Appendix C.3)
-/

namespace SqueezeFilm.Taper
open SqueezeFilm.Calculus

/-- Closed form (4.3): `∫₀^W y²/(h₀ + βy)³ dy = β⁻³[ln(h_W/h₀) - 3/2 + 2h₀/h_W - h₀²/(2h_W²)]`. -/
theorem taper_integral (h0 β W : ℝ) (h0pos : 0 < h0) (βpos : 0 < β) (Wnn : 0 ≤ W) :
    ∫ y in (0:ℝ)..W, y ^ 2 / (h0 + β * y) ^ 3
      = 1 / β ^ 3 * (Real.log ((h0 + β * W) / h0) - 3 / 2 + 2 * h0 / (h0 + β * W)
          - h0 ^ 2 / (2 * (h0 + β * W) ^ 2)) := by
  have hpos : ∀ y ∈ Set.uIcc (0:ℝ) W, 0 < h0 + β * y := by
    intro y hy
    rw [Set.uIcc_of_le Wnn] at hy
    have : 0 ≤ β * y := mul_nonneg βpos.le hy.1
    linarith
  have hβ : β ≠ 0 := βpos.ne'
  have hF : ∀ y ∈ Set.uIcc (0:ℝ) W, HasDerivAt
      (fun y => 1 / β ^ 3 * (Real.log (h0 + β * y) + 2 * h0 / (h0 + β * y)
        - h0 ^ 2 / (2 * ((h0 + β * y) * (h0 + β * y)))))
      (y ^ 2 / (h0 + β * y) ^ 3) y := by
    intro y hy
    have hy0 := hpos y hy
    have hH : HasDerivAt (fun y => h0 + β * y) β y := by
      convert hasDerivAt_cubic 0 0 β h0 y using 1
      · funext x; ring
      · ring
    have hlog := hH.log hy0.ne'
    have hA := (hasDerivAt_const (x := y) (c := 2 * h0)).div hH hy0.ne'
    have hHH := (hH.mul hH).const_mul 2
    have hB := (hasDerivAt_const (x := y) (c := h0 ^ 2)).div hHH
      (mul_pos two_pos (mul_pos hy0 hy0)).ne'
    have h := ((hlog.add hA).sub hB).const_mul (1 / β ^ 3)
    have hy0' := hy0.ne'
    convert h using 1
    field_ring
  have hW := hpos W Set.right_mem_uIcc
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF (ContinuousOn.intervalIntegrable ?_)]
  · beta_reduce
    simp only [mul_zero, add_zero]
    rw [Real.log_div hW.ne' h0pos.ne']
    have h0' := h0pos.ne'
    have hW' := hW.ne'
    field_ring
  · apply ContinuousOn.div (by fun_prop) (by fun_prop)
    intro y hy
    exact (pow_pos (hpos y hy) 3).ne'

/-- **Finiteness bound** (Appendix C.3): the bracket of (4.3) is at most `ln(h_W/h₀) + 1/2`, since
`-3/2 + 2t - t²/2 ≤ 1/2` for every real `t` (here `t = h₀/h_W`). Hence a contact on a set where
`∫ ln g > -∞` gives finite damping. -/
theorem finiteness_bound (h0 β W : ℝ) (h0pos : 0 < h0) (βpos : 0 < β) (Wnn : 0 ≤ W) :
    ∫ y in (0:ℝ)..W, y ^ 2 / (h0 + β * y) ^ 3
      ≤ 1 / β ^ 3 * (Real.log ((h0 + β * W) / h0) + 1 / 2) := by
  rw [taper_integral h0 β W h0pos βpos Wnn]
  have hHW : 0 < h0 + β * W := by positivity
  have key : -3 / 2 + 2 * h0 / (h0 + β * W) - h0 ^ 2 / (2 * (h0 + β * W) ^ 2) ≤ 1 / 2 := by
    have e1 : 2 * h0 / (h0 + β * W) = 2 * (h0 / (h0 + β * W)) := by ring
    have hH := hHW.ne'
    have e2 : h0 ^ 2 / (2 * (h0 + β * W) ^ 2) = (h0 / (h0 + β * W)) ^ 2 / 2 := by field_ring
    rw [e1, e2]
    nlinarith [sq_nonneg (h0 / (h0 + β * W) - 2)]
  have hb : 0 ≤ 1 / β ^ 3 := by positivity
  apply mul_le_mul_of_nonneg_left _ hb
  linarith

end SqueezeFilm.Taper
