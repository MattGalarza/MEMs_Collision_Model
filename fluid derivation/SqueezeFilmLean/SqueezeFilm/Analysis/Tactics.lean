import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-- Close an identity in a field: plain `ring`, or `field_simp` (using nonzero hypotheses in
context) followed by `ring`. Robust to `field_simp` closing the goal by itself. -/
macro "field_ring" : tactic =>
  `(tactic| first | ring | (field_simp; ring) | field_simp | (norm_num; ring) | (field_simp; ring_nf))
