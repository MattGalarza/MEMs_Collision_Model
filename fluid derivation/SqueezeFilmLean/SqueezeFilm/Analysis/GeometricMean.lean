import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Tactic.Positivity
import SqueezeFilm.Analysis.Tactics

/-!
# Geometric-mean residual gap of a linear taper (eq. (6.8)): `exp⟨ln(βy)⟩ = βW/e`
-/

namespace SqueezeFilm.GeometricMean

theorem geometric_mean_linear (β W : ℝ) (hβ : 0 < β) (hW : 0 < W) :
    1 / W * ∫ y in (0:ℝ)..W, Real.log (β * y) = Real.log (β * W / Real.exp 1) := by
  have hcongr : ∫ y in (0:ℝ)..W, Real.log (β * y) = ∫ y in (0:ℝ)..W, (Real.log β + Real.log y) := by
    apply intervalIntegral.integral_congr_ae
    refine Filter.Eventually.of_forall (fun y hy => ?_)
    rw [Set.uIoc_of_le hW.le] at hy
    rw [Real.log_mul hβ.ne' hy.1.ne']
  rw [hcongr, intervalIntegral.integral_add intervalIntegrable_const intervalIntegral.intervalIntegrable_log',
    intervalIntegral.integral_const, integral_log,
    Real.log_div (by positivity) (Real.exp_pos 1).ne', Real.log_mul hβ.ne' hW.ne', Real.log_exp]
  simp only [smul_eq_mul, Real.log_zero, mul_zero, sub_zero, zero_mul, add_zero]
  have hW' := hW.ne'
  field_ring

end SqueezeFilm.GeometricMean
