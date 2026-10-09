import SqueezeFilm.Analysis.Polynomial
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# Flow across the gap (§2.5): profile, flux law and dissipation, by calculus

Coordinates across the gap are `s = z - h_b ∈ [0, h]`. The core layer checks the algebra; here the
derivatives and integrals themselves are proved.
-/

namespace SqueezeFilm.Profile
open SqueezeFilm.Calculus

/-- The profile (2.14) has `∂_z u = du` and `∂²_z u = P/μ`: it solves `μ ∂²_z u = ∂ₓ p`. -/
theorem profile_derivatives (P μ b h Ua Ub s : ℝ) (hμ : μ ≠ 0) (hd : h + 2 * b ≠ 0) :
    HasDerivAt (fun s : ℝ => P / (2 * μ) * (s * (s - h) - b * h) + Ub + (Ua - Ub) * (s + b) / (h + 2 * b))
      (P / (2 * μ) * (2 * s - h) + (Ua - Ub) / (h + 2 * b)) s ∧
    HasDerivAt (fun s : ℝ => P / (2 * μ) * (2 * s - h) + (Ua - Ub) / (h + 2 * b)) (P / μ) s := by
  constructor
  · have h' := hasDerivAt_cubic 0 (P / (2 * μ)) (-(P / (2 * μ)) * h + (Ua - Ub) / (h + 2 * b))
        (-(P / (2 * μ)) * b * h + Ub + (Ua - Ub) * b / (h + 2 * b)) s
    convert h' using 1
    · funext x; field_ring
    · field_ring
  · have h' := hasDerivAt_cubic 0 0 (P / μ) (-(P / (2 * μ)) * h + (Ua - Ub) / (h + 2 * b)) s
    convert h' using 1
    · funext x; field_ring
    · field_ring

/-- **Flux law (2.15)** by the fundamental theorem of calculus:
`∫₀ʰ u ds = -h²(h + 6b)/(12μ) ∂ₓp + h (U_a + U_b)/2`. -/
theorem flux_integral (P μ b h Ua Ub : ℝ) (hμ : μ ≠ 0) (hd : h + 2 * b ≠ 0) :
    ∫ s in (0:ℝ)..h, (P / (2 * μ) * (s * (s - h) - b * h) + Ub + (Ua - Ub) * (s + b) / (h + 2 * b))
      = -(h ^ 2 * (h + 6 * b)) / (12 * μ) * P + h * (Ua + Ub) / 2 := by
  have hF : ∀ x ∈ Set.uIcc (0:ℝ) h, HasDerivAt
      (fun s : ℝ => P / (6 * μ) * s ^ 3 + (-(P / (4 * μ)) * h + (Ua - Ub) / (2 * (h + 2 * b))) * s ^ 2
          + (-(P / (2 * μ)) * b * h + Ub + (Ua - Ub) * b / (h + 2 * b)) * s + 0)
      (P / (2 * μ) * (x * (x - h) - b * h) + Ub + (Ua - Ub) * (x + b) / (h + 2 * b)) x := by
    intro x _
    have h' := hasDerivAt_cubic (P / (6 * μ)) (-(P / (4 * μ)) * h + (Ua - Ub) / (2 * (h + 2 * b)))
      (-(P / (2 * μ)) * b * h + Ub + (Ua - Ub) * b / (h + 2 * b)) 0 x
    convert h' using 1
    field_ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF
    (by apply Continuous.intervalIntegrable; fun_prop)]
  beta_reduce
  field_ring

/-- **Dissipation identity (2.16)** by calculus: viscous dissipation of the Poiseuille part
`μ ∫₀ʰ (∂_s u)² ds = h³P²/(12μ)`, plus slip dissipation `(μ/b)·2 u_s²`, equals `G P²/(12μ)`. -/
theorem dissipation_integral (P μ b h : ℝ) (hμ : μ ≠ 0) (hb : b ≠ 0) :
    (∫ s in (0:ℝ)..h, μ * (P / (2 * μ) * (2 * s - h)) ^ 2)
      + μ / b * (2 * (P / (2 * μ) * (-(b * h))) ^ 2) = h ^ 2 * (h + 6 * b) * P ^ 2 / (12 * μ) := by
  have hF : ∀ x ∈ Set.uIcc (0:ℝ) h, HasDerivAt
      (fun s : ℝ => (μ * (P / (2 * μ)) ^ 2 * (4 / 3)) * s ^ 3 + (-(μ * (P / (2 * μ)) ^ 2 * 2 * h)) * s ^ 2
          + (μ * (P / (2 * μ)) ^ 2 * h ^ 2) * s + 0)
      (μ * (P / (2 * μ) * (2 * x - h)) ^ 2) x := by
    intro x _
    have h' := hasDerivAt_cubic (μ * (P / (2 * μ)) ^ 2 * (4 / 3)) (-(μ * (P / (2 * μ)) ^ 2 * 2 * h))
      (μ * (P / (2 * μ)) ^ 2 * h ^ 2) 0 x
    convert h' using 1
    field_ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hF
    (by apply Continuous.intervalIntegrable; fun_prop)]
  beta_reduce
  field_ring

end SqueezeFilm.Profile
