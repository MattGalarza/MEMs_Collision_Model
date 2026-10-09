import SqueezeFilm.Core
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.CharP.Basic

/-!
# The core (field-generic) theorems apply to ℝ

Mathlib supplies `Lean.Grind.Field ℝ` (from `Field ℝ`) and `Lean.Grind.IsCharP ℝ 0` (from `CharP ℝ 0`),
so every theorem of `SqueezeFilm.Core` specializes to real numbers. These checks fail to compile if
either instance is missing.
-/

namespace SqueezeFilm.RealInstances

example : Lean.Grind.Field ℝ := inferInstance
example : Lean.Grind.IsCharP ℝ 0 := inferInstance

noncomputable example := @SqueezeFilm.Core.corner_sealed ℝ _ _
noncomputable example := @SqueezeFilm.Core.flux_law ℝ _ _
noncomputable example := @SqueezeFilm.Core.scaling_coefficients ℝ _ _

end SqueezeFilm.RealInstances
