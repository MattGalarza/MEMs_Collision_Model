# Lean 4 formalization of "Squeeze films between doubly tapered plates closing to contact"

Toolchain `leanprover/lean4:v4.35.0-rc4`; Mathlib pinned to commit `49fd9dee36e7bbf345fe7981d2de6aeeb591f75b` (the commit that uses this toolchain).

    curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh   # once
    cd SqueezeFilmLean
    lake exe cache get        # prebuilt Mathlib for the pinned commit
    lake build                # builds everything; success = no errors and no `sorry`
    lean SqueezeFilm/Core.lean   # the core layer alone needs no Mathlib at all

## Status

| layer | files | status |
|---|---|---|
| Core (algebraic skeleton) | `SqueezeFilm/Core.lean` | **Compiled and checked**: 32 theorems (7 are field helpers), each depending only on `propext`, `Classical.choice`, `Quot.sound`; no `sorry`. See `core_build.log`. |
| Analysis (calculus, integrals, series, variational logic) | `SqueezeFilm/Analysis/*.lean` | **Written, not compiled.** Mathlib could not be built where this was written (1 core, 3 GB). Every lemma name and namespace used was checked against the pinned Mathlib source; proof-level fixes may still be needed. |

Core theorems are stated for any field of characteristic zero (`Lean.Grind.Field`, `Lean.Grind.IsCharP α 0`); Mathlib supplies both instances for ℝ, and `Analysis/RealInstances.lean` fails to compile if it does not.

## Map to the paper

| paper | Lean |
|---|---|
| §2.3 Stokes hypothesis, (2.4) | `Core.stokes_hypothesis` |
| §2.4 scaling (2.7)–(2.12), Re*/Re_ω | `Core.scaling_coefficients`, `Core.zmom_viscous`, `Core.reynolds_ratio` |
| (2.14) profile: slip conditions; ODE | `Core.profile_slip_conditions`; `Profile.profile_derivatives` |
| (2.15) flux law | `Core.flux_law` (+ `poiseuille_part`, `couette_part`, `flux_antiderivative_matches`); `Profile.flux_integral` (by FTC) |
| (2.16) dissipation | `Core.dissipation_identity`; `Profile.dissipation_integral` (by FTC) |
| Prop 2.1 compatibility: kinematics, Leibniz collection | `Core.kinematic_condition`, `Core.leibniz_collection` |
| (2.23) damping matrix is a Gram matrix; symmetry | `Variational.gram_psd`, `Variational.damping_symmetric` |
| §3.3 pressure principle (lower bounds) | `Variational.pressure_principle`, `pressure_principle_attained` |
| Prop 3.2 monotonicity in G (gap, slip) | `Variational.damping_monotone` |
| Prop 3.2 monotonicity in venting; Galerkin lower bound | `Variational.subspace_lower_bound`, `Variational.subspace_monotone` |
| (3.5) sum rules | `Series.hasSum_odd_inv_sq`, `Series.hasSum_odd_inv_fourth` |
| Prop 4.1 strip law | `Strip.strip_law` (with `Strip.integral_tail_eq`); `Core.strip_uniform` |
| (4.2) small-taper moments | `Core.small_taper_moments` |
| (4.3) linear-taper closed form | `Taper.taper_integral` |
| (4.4) class-B strip limit (outer integral) | `Core.strip_limit_B` |
| Prop 5.1 corner ODE solutions | `Corner.corner_particular`, `corner_homogeneous_odd`, `corner_homogeneous_even` (by differentiation); `Core.corner_ode_algebra`, `Core.corner_flux_forms` |
| Prop 5.1 coefficients (5.2), (5.3) | `Core.corner_sealed`, `Core.corner_open` |
| Prop 5.2 slip corner | `Corner.slip_corner_integral`; `Core.slip_corner_constant` |
| Prop 5.3 finiteness bound | `Taper.finiteness_bound`; `Core.finiteness_bracket` |
| Prop 6.5 exact scaling | `Core.homogeneity` |
| (6.7) slice function partial fractions | `Core.slice_partial_fractions` |
| (6.8) geometric mean βW/e; (6.9) constant | `GeometricMean.geometric_mean_linear`; `Core.small_taper_constant` |
| App. D Moy mapping | `Core.moy_mapping` |

The variational theorems are proved for an arbitrary symmetric nonnegative bilinear form on a real vector space, so one proof covers the continuous weak problem and every conforming discretization. Monotonicity in venting and the finite-element lower bound turn out to be the same subspace lemma. The monotonicity proof runs through the pressure principle rather than the flux principle used in the paper; the two routes are equivalent.

## Not formalized, and why

The asymptotic thin-film limit of the Navier–Stokes equations as a theorem about solutions; existence and regularity of weak solutions of the film equation where it degenerates at contact (the variational theorems take a solution as a hypothesis); the formal matched asymptotics of Proposition 6.6; the numerical contact maps; $J_C = 1 - \ln 2$ and the small-taper series beyond their moments (both are checked symbolically in `lean_verification/` and numerically in the Julia script). Mathlib has none of the PDE infrastructure the first three would need.

## If an Analysis proof fails

The usual causes are normal-form drift in `field_simp` (replace `field_ring` by `field_simp` then `ring_nf`, or vice versa), a `convert ... using 1` that leaves its goals in a different order (swap the bullets), and simp-set differences in `gram_psd` (finish with `Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => by ring`). `exact?` finds renamed lemmas. A statement should never need to change; if one does, that is a finding about the paper and worth flagging.
