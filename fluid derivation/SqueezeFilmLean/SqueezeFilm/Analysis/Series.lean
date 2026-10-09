import Mathlib.NumberTheory.ZetaValues
import Mathlib.Topology.Algebra.InfiniteSum.NatInt
import Mathlib.Tactic.Linarith

/-!
# Modal sum rules (3.5): sums over odd integers

`Σ_{q odd} q⁻² = π²/8` and `Σ_{q odd} q⁻⁴ = π⁴/96`, from `ζ(2) = π²/6`, `ζ(4) = π⁴/90` by
splitting into even and odd terms.
-/

namespace SqueezeFilm.Series
open Real

theorem hasSum_odd_inv_sq :
    HasSum (fun k : ℕ => 1 / ((2 * k + 1 : ℕ) : ℝ) ^ 2) (π ^ 2 / 8) := by
  set f : ℕ → ℝ := fun n => 1 / (n : ℝ) ^ 2 with hf
  have hz : HasSum f (π ^ 2 / 6) := hasSum_zeta_two
  have he : HasSum (fun k => f (2 * k)) (π ^ 2 / 24) := by
    have := hz.mul_left (1 / 4)
    convert this using 1
    · funext k; simp only [hf]; push_cast; ring
    · ring
  have hos : Summable (fun k => f (2 * k + 1)) :=
    hz.summable.comp_injective (fun a b hab => by simp only at hab; omega)
  have huniq := hz.unique (he.even_add_odd hos.hasSum)
  have hval : ∑' k, f (2 * k + 1) = π ^ 2 / 8 := by linarith
  have h := hos.hasSum
  rw [hval] at h
  exact h

theorem hasSum_odd_inv_fourth :
    HasSum (fun k : ℕ => 1 / ((2 * k + 1 : ℕ) : ℝ) ^ 4) (π ^ 4 / 96) := by
  set f : ℕ → ℝ := fun n => 1 / (n : ℝ) ^ 4 with hf
  have hz : HasSum f (π ^ 4 / 90) := hasSum_zeta_four
  have he : HasSum (fun k => f (2 * k)) (π ^ 4 / 1440) := by
    have := hz.mul_left (1 / 16)
    convert this using 1
    · funext k; simp only [hf]; push_cast; ring
    · ring
  have hos : Summable (fun k => f (2 * k + 1)) :=
    hz.summable.comp_injective (fun a b hab => by simp only at hab; omega)
  have huniq := hz.unique (he.even_add_odd hos.hasSum)
  have hval : ∑' k, f (2 * k + 1) = π ^ 4 / 96 := by linarith
  have h := hos.hasSum
  rw [hval] at h
  exact h

end SqueezeFilm.Series
