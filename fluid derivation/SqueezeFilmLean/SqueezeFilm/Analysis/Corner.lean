import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Tactic.FunProp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.DefEqTransformations
import SqueezeFilm.Analysis.Tactics

/-!
# Corner solutions (§5, Appendix B): the angular ODE and the slip-corner integral

With `ψ` the angle measured from the line of contact, a pressure `p = f(ψ)/r` solves the corner
equation iff `L[f] := (cos³ψ f')' - 2 cos³ψ f = -A`. We verify, by actual differentiation, that
`1/cos ψ` is a particular solution (`L = -1`) and that `sin ψ/cos²ψ` and `1/cos²ψ` are homogeneous.
-/

namespace SqueezeFilm.Corner
open Real

theorem hasDerivAt_sec (ψ : ℝ) (hc : cos ψ ≠ 0) :
    HasDerivAt (fun x => 1 / cos x) (sin ψ / cos ψ ^ 2) ψ := by
  have h := (hasDerivAt_const (x := ψ) (c := (1:ℝ))).div (Real.hasDerivAt_cos ψ) hc
  convert h using 1
  field_ring

theorem hasDerivAt_f2 (ψ : ℝ) (hc : cos ψ ≠ 0) :
    HasDerivAt (fun x => sin x / (cos x * cos x)) ((cos ψ ^ 2 + 2 * sin ψ ^ 2) / cos ψ ^ 3) ψ := by
  have hcc : cos ψ * cos ψ ≠ 0 := mul_ne_zero hc hc
  have h := (Real.hasDerivAt_sin ψ).div ((Real.hasDerivAt_cos ψ).mul (Real.hasDerivAt_cos ψ)) hcc
  convert h using 1
  field_ring

theorem hasDerivAt_f3 (ψ : ℝ) (hc : cos ψ ≠ 0) :
    HasDerivAt (fun x => 1 / (cos x * cos x)) (2 * sin ψ / cos ψ ^ 3) ψ := by
  have hcc : cos ψ * cos ψ ≠ 0 := mul_ne_zero hc hc
  have h := (hasDerivAt_const (x := ψ) (c := (1:ℝ))).div
    ((Real.hasDerivAt_cos ψ).mul (Real.hasDerivAt_cos ψ)) hcc
  convert h using 1
  field_ring

/-- Particular solution `f = 1/cos ψ`: the flux `cos³ψ f' = sin ψ cos ψ`, its derivative is
`cos²ψ - sin²ψ`, and `L[f] = -1`. -/
theorem corner_particular (ψ : ℝ) (hc : cos ψ ≠ 0) :
    cos ψ ^ 3 * (sin ψ / cos ψ ^ 2) = sin ψ * cos ψ ∧
    HasDerivAt (fun x => sin x * cos x) (cos ψ ^ 2 - sin ψ ^ 2) ψ ∧
    (cos ψ ^ 2 - sin ψ ^ 2) - 2 * cos ψ ^ 3 * (1 / cos ψ) = -1 := by
  have hs := Real.sin_sq_add_cos_sq (x := ψ)
  refine ⟨by field_ring, ?_, ?_⟩
  · have h := (Real.hasDerivAt_sin ψ).mul (Real.hasDerivAt_cos ψ)
    convert h using 1
    ring
  · have e : 2 * cos ψ ^ 3 * (1 / cos ψ) = 2 * cos ψ ^ 2 := by field_ring
    rw [e]
    linarith

/-- Homogeneous solution `f = sin ψ/cos²ψ`: flux `cos²ψ + 2 sin²ψ`, derivative `2 sin ψ cos ψ`, `L[f] = 0`. -/
theorem corner_homogeneous_odd (ψ : ℝ) (hc : cos ψ ≠ 0) :
    cos ψ ^ 3 * ((cos ψ ^ 2 + 2 * sin ψ ^ 2) / cos ψ ^ 3) = cos ψ ^ 2 + 2 * sin ψ ^ 2 ∧
    HasDerivAt (fun x => cos x * cos x + 2 * (sin x * sin x)) (2 * sin ψ * cos ψ) ψ ∧
    2 * sin ψ * cos ψ - 2 * cos ψ ^ 3 * (sin ψ / (cos ψ * cos ψ)) = 0 := by
  refine ⟨by field_ring, ?_, by field_ring⟩
  have h := ((Real.hasDerivAt_cos ψ).mul (Real.hasDerivAt_cos ψ)).add
    (((Real.hasDerivAt_sin ψ).mul (Real.hasDerivAt_sin ψ)).const_mul 2)
  convert h using 1
  ring

/-- Homogeneous solution `f = 1/cos²ψ`: flux `2 sin ψ`, derivative `2 cos ψ`, `L[f] = 0`. -/
theorem corner_homogeneous_even (ψ : ℝ) (hc : cos ψ ≠ 0) :
    cos ψ ^ 3 * (2 * sin ψ / cos ψ ^ 3) = 2 * sin ψ ∧
    HasDerivAt (fun x => 2 * sin x) (2 * cos ψ) ψ ∧
    2 * cos ψ - 2 * cos ψ ^ 3 * (1 / (cos ψ * cos ψ)) = 0 :=
  ⟨by field_ring, (Real.hasDerivAt_sin ψ).const_mul 2, by field_ring⟩

/-- **Slip-corner solvability integral** (Proposition 5.2):
`∫₀^{π/2} (α cos θ + β sin θ)² dθ = π(α² + β²)/4 + αβ`. -/
theorem slip_corner_integral (α β : ℝ) :
    ∫ θ in (0:ℝ)..(π / 2), (α * cos θ + β * sin θ) ^ 2 = π * (α ^ 2 + β ^ 2) / 4 + α * β := by
  have hF : ∀ x ∈ Set.uIcc (0:ℝ) (π / 2), HasDerivAt
      (fun θ => α ^ 2 * (θ / 2 + sin θ * cos θ / 2) + β ^ 2 * (θ / 2 - sin θ * cos θ / 2)
        + α * β * (sin θ * sin θ))
      ((α * cos x + β * sin x) ^ 2) x := by
    intro x _
    have hsc := (Real.hasDerivAt_sin x).mul (Real.hasDerivAt_cos x)
    have hss := (Real.hasDerivAt_sin x).mul (Real.hasDerivAt_sin x)
    have hid : HasDerivAt (fun θ : ℝ => θ) 1 x := hasDerivAt_id' (x := x)
    have h := ((((hid.div_const 2).add (hsc.div_const 2)).const_mul (α ^ 2)).add
      (((hid.div_const 2).sub (hsc.div_const 2)).const_mul (β ^ 2))).add (hss.const_mul (α * β))
    convert h using 1
    have hs := Real.sin_sq_add_cos_sq (x := x)
    linear_combination ((α ^ 2 + β ^ 2) / 2) * hs
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF
    (by apply Continuous.intervalIntegrable; fun_prop)]
  beta_reduce
  simp only [Real.sin_pi_div_two, Real.cos_pi_div_two, Real.sin_zero, Real.cos_zero]
  ring

end SqueezeFilm.Corner
