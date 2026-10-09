import Mathlib.LinearAlgebra.BilinearMap
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Variational principles, monotonicity and the damping matrix (§2.7, §3.3, Props. 3.1–3.2)

The weak form of the film equation is: find `p*` with `k p* v = ℓ v` for all admissible `v`, where
`k p q = ∫ G ∇p·∇q / (12μ)` (symmetric, nonnegative) and `ℓ v = ∫ ḣ v`. The damping is `ℓ p*`.
Everything below is stated for an abstract real vector space `V`, so it covers the continuous
problem (V a Sobolev space) and every conforming discretization (V finite-dimensional) at once.
-/

namespace SqueezeFilm.Variational

variable {V : Type*} [AddCommGroup V] [Module ℝ V]

/-- Expansion of the quadratic form of a symmetric bilinear form. -/
theorem quadratic_expansion (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (hsymm : ∀ u v, k u v = k v u) (p q : V) :
    k (p - q) (p - q) = k p p - 2 * k q p + k q q := by
  rw [LinearMap.map_sub₂, map_sub, map_sub, hsymm p q]
  ring

/-- **Pressure principle** (§3.3): every trial pressure gives a lower bound on the damping,
`2 ℓ p - k p p ≤ ℓ p*`, with equality at `p = p*`. -/
theorem pressure_principle (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (ℓ : V →ₗ[ℝ] ℝ)
    (hsymm : ∀ u v, k u v = k v u) (hpos : ∀ v, 0 ≤ k v v)
    (pstar : V) (hstar : ∀ v, k pstar v = ℓ v) (p : V) :
    2 * ℓ p - k p p ≤ ℓ pstar := by
  have h := hpos (p - pstar)
  rw [quadratic_expansion k hsymm p pstar, hstar p, hstar pstar] at h
  linarith

/-- The bound is attained by the solution itself: `ℓ p* = 2 ℓ p* - k p* p*`. -/
theorem pressure_principle_attained (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (ℓ : V →ₗ[ℝ] ℝ)
    (pstar : V) (hstar : ∀ v, k pstar v = ℓ v) :
    2 * ℓ pstar - k pstar pstar = ℓ pstar := by
  rw [hstar pstar]; ring

/-- **Monotonicity in the mobility** (Proposition 3.2): if the film operator grows
(`k₁ ≤ k₂` as quadratic forms, i.e. more gap or more slip everywhere: `G₁ ≤ G₂`),
the damping can only fall. -/
theorem damping_monotone (k₁ k₂ : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (ℓ : V →ₗ[ℝ] ℝ)
    (hsymm₁ : ∀ u v, k₁ u v = k₁ v u) (hpos₁ : ∀ v, 0 ≤ k₁ v v)
    (hle : ∀ v, k₁ v v ≤ k₂ v v)
    (p₁ p₂ : V) (h₁ : ∀ v, k₁ p₁ v = ℓ v) (h₂ : ∀ v, k₂ p₂ v = ℓ v) :
    ℓ p₂ ≤ ℓ p₁ := by
  have hp := pressure_principle k₁ ℓ hsymm₁ hpos₁ p₁ h₁ p₂
  have e : k₂ p₂ p₂ = ℓ p₂ := h₂ p₂
  have hl := hle p₂
  linarith

/-- **Galerkin lower bound, and monotonicity in venting** (Proposition 3.2): the solution of the
same problem on a subspace `S` (a finite-element space; or the pressures that also vanish on an
added vent) has damping no larger than the full solution. -/
theorem subspace_lower_bound (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (ℓ : V →ₗ[ℝ] ℝ)
    (hsymm : ∀ u v, k u v = k v u) (hpos : ∀ v, 0 ≤ k v v)
    (pstar : V) (hstar : ∀ v, k pstar v = ℓ v)
    (S : Submodule ℝ V) (pS : V) (hmem : pS ∈ S) (hS : ∀ w ∈ S, k pS w = ℓ w) :
    ℓ pS ≤ ℓ pstar := by
  have hp := pressure_principle k ℓ hsymm hpos pstar hstar pS
  have e : k pS pS = ℓ pS := hS pS hmem
  linarith

/-- Galerkin solutions are monotone under refinement: a larger subspace never gives less damping. -/
theorem subspace_monotone (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (ℓ : V →ₗ[ℝ] ℝ)
    (hsymm : ∀ u v, k u v = k v u) (hpos : ∀ v, 0 ≤ k v v)
    (S T : Submodule ℝ V) (hST : S ≤ T)
    (pS pT : V) (hmS : pS ∈ S) (hmT : pT ∈ T)
    (hS : ∀ w ∈ S, k pS w = ℓ w) (hT : ∀ w ∈ T, k pT w = ℓ w) :
    ℓ pS ≤ ℓ pT := by
  -- p_S is a trial function for the problem on T; apply the pressure principle inside T.
  have h := hpos (pS - pT)
  rw [quadratic_expansion k hsymm pS pT, hT pS (hST hmS)] at h
  have e1 : k pS pS = ℓ pS := hS pS hmS
  have e2 : k pT pT = ℓ pT := hT pT hmT
  linarith

/-- **Damping matrix is a Gram matrix** (2.23): `D_ij = k(pᵢ, pⱼ)` is positive semidefinite. -/
theorem gram_psd {m : ℕ} (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (hpos : ∀ v, 0 ≤ k v v)
    (φ : Fin m → V) (c : Fin m → ℝ) :
    0 ≤ ∑ i, ∑ j, c i * (c j * k (φ i) (φ j)) := by
  have h := hpos (∑ i, c i • φ i)
  have e : k (∑ i, c i • φ i) (∑ j, c j • φ j) = ∑ i, ∑ j, c i * (c j * k (φ i) (φ j)) := by
    simp [map_sum, LinearMap.sum_apply, map_smul, smul_eq_mul, Finset.mul_sum]
  rw [e] at h
  exact h

/-- Symmetry of the damping matrix follows from symmetry of the form. -/
theorem damping_symmetric (k : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (hsymm : ∀ u v, k u v = k v u) (p q : V) :
    k p q = k q p := hsymm p q

end SqueezeFilm.Variational
