import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Add
import SqueezeFilm.Analysis.Tactics

namespace SqueezeFilm.Calculus

/-- Derivative of a cubic polynomial; every polynomial antiderivative below reduces to it. -/
theorem hasDerivAt_cubic (c3 c2 c1 c0 s : ℝ) :
    HasDerivAt (fun x : ℝ => c3 * x ^ 3 + c2 * x ^ 2 + c1 * x + c0)
      (3 * c3 * s ^ 2 + 2 * c2 * s + c1) s := by
  have h3 := (hasDerivAt_pow 3 s).const_mul c3
  have h2 := (hasDerivAt_pow 2 s).const_mul c2
  have h1 := (hasDerivAt_id' (x := s)).const_mul c1
  have h0 := hasDerivAt_const (x := s) (c := c0)
  have h := ((h3.add h2).add h1).add h0
  convert h using 1 <;> (try norm_num) <;> (try ring)

end SqueezeFilm.Calculus
