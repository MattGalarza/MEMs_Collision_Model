# Understanding and reproducing the MEMS plots

The corrected Julia model has now been run natively. The visualization package copies the corrected model without changing its equations or physical defaults. The new files add detailed sampling, diagnostics, and figure rendering; they do not introduce a new force law merely to improve the appearance of a plot.

## 1. What was run

The driven case uses a 20 Hz sinusoidal base acceleration with a 4.95 g peak amplitude and a 3 V bias. The amplitude rises through a four-cycle cosine ramp, then remains fixed for four more cycles: 0.4 seconds total. The solver is Rodas5P, with relative tolerance 1e-7 and absolute tolerance 1e-10 on scaled states. Both its completion status and energy checks passed.

A separate 600-microsecond unforced contact probe was also run. It starts with prescribed mechanical displacement and velocity close to the wall. Its voltage is a transient response to that initialization and is not a steady harvesting-voltage result.

Measured from the numerical run (these are model predictions, not device measurements):

| Quantity | Result | Interpretation |
|---|---:|---|
| Mean resistor power, last four cycles | 0.9885 pW | Computed from accumulated resistor energy, not coarse voltage samples |
| Largest sampled absolute voltage, final cycle | 8.612 mV | Sampled on accepted solver states; an optimized continuous-time peak was not calculated |
| Largest sampled tip overlap, final cycle | 16.17 nm | Allowed by the compliant contact law |
| Largest sampled shuttle–tip displacement difference | 83.09 nm | Electrode bending coordinate magnitude |
| Largest sampled absolute shuttle displacement | 13.809 μm | Compare with the nominal tip-contact travel of 13.71 μm |
| Integrated fluid loss, entire run | 2.281 nJ | Dominant modeled dissipative channel for these defaults |
| Integrated contact loss, entire run | 0.7820 pJ | Includes the unloading-cutoff contribution |
| Integrated resistor energy, entire run | 0.2668 pJ | Includes the initial forcing ramp |
| Maximum energy residual / integrated power throughput | 6.66e-11 | Numerical energy consistency; not a physical-accuracy estimate |

The final-cycle state nearly recurs after one drive period, but that observation alone does not establish stability or exclude other attractors. A driven tolerance/grid study, Floquet analysis, and experimental comparison remain separate tasks.

## 2. How to read the four figures

### Figure 1: driven overview

Panel (a) shows the shuttle and tip following nearly the same micrometer-scale motion. Once the forcing amplitude is large enough, the tip reaches the contact region and the response develops plateaus near the walls. The curves can overlap visually even though the electrode is bending by tens of nanometers.

Panel (b) shows short voltage pulses associated with fast changes in capacitance. This overview is sampled every 50 microseconds, so its visible pulse heights must not be treated as the resolved maxima. The dense final-cycle plot and contact close-up resolve those features.

Panel (c) shows the total capacitance, including the parasitic component. Both sides contribute, so large capacitance can occur near either wall. Panel (d) shows step-like increases in resistor energy when voltage pulses deposit energy in the load. The vertical line marks the end of the forcing ramp, not a bifurcation.

### Figure 2: contact close-up

The first positive-wall entry in the final cycle occurs at approximately 0.37972376 s. The plot uses 0.1-microsecond sampling, with time zero defined by x2 = gc.

The tip slows sharply as contact develops while the shuttle continues moving. Their separation changes, storing bending energy. The wall reaction and bending force then become large and oppose one another. The electrode subsequently unloads and rebounds; the circuit responds to the associated change in capacitance with a voltage pulse of the opposite sign.

The force panel shows forces conjugate to the tip coordinate, summed over the mobile beams. It is not a complete Newtonian force balance for an isolated tip mass: the mechanical equations also contain base excitation and off-diagonal inertia.

### Figure 3: final cycle

The phase portrait displays velocity against shuttle position. The contact regions distort the orbit away from a sinusoidal ellipse. The voltage and instantaneous power panels expose short pulses; the sampled voltage maximum is approximately 8.61 mV even though the full-run overview can show a smaller apparent maximum.

Instantaneous power is Vout^2/R and is always nonnegative. Its peak and its cycle average are different quantities. The reported 0.9885 pW average uses the integrated energy state over four cycles and therefore does not depend on the displayed sample spacing.

The capacitance-versus-shuttle-position curve can have a finite width because capacitance depends on both x1 and x2. A projection onto x1 alone does not uniquely specify the beam deformation. This should not automatically be called material hysteresis or evidence of chaos.

### Figure 4: energy and separate probe

Panels (a)–(c) show source work, stored energy, dissipated energy, and the energy residual for the driven run. Fluid loss is much larger than resistor dissipation with this particular gas closure and parameter set. That is a useful design/model-validation observation, not an experimentally established device efficiency.

The bias-source work is tracked explicitly. It can exchange energy during a transient; excluding it can lead to an incorrect interpretation of shaker-to-load energy conversion.

Panel (d) is the separate unforced contact probe. Its larger voltage follows its chosen near-contact initial state and must not be compared with the driven operating voltage as though the excitation were the same.

## 3. Make the plots with Julia

Extract `MEMS_Visualization.zip`, then open a terminal in the resulting `MEMS_Visualization` directory.

```sh
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. run_visualization.jl
julia --project=. make_figures.jl
```

The first command installs the dependencies in the included environment. Julia 1.11.7 was used for the supplied results. The second command runs the model and writes the raw CSV data and numerical summaries. The third command prepares display data and creates the four figures.

If you only want to redraw the supplied results, run the third command alone. The supplied CSV files already contain the numerical results. Expect the full driven simulation to take several minutes or longer depending on hardware because rapid contact motion is resolved within much slower forcing cycles.

The simulation and data preparation are native Julia. The Julia plotting script calls pdfLaTeX with PGFPlots to produce vector PDFs. Install a LaTeX distribution containing PGFPlots to regenerate these PDFs. If Poppler's `pdftoppm` is also on your PATH, the script additionally produces 300 dpi PNGs. The supplied PDFs and PNGs can be viewed without installing any of those programs. No Python code is required.

### File responsibilities

| File | Purpose |
|---|---|
| `collision_model_corrected.jl` | Mechanical/electrical equations, capacitance, fluid damping, contact, energy accounting, and solver |
| `run_visualization.jl` | Executes the reference case, locates a contact entry, exports fine contact/final-cycle data, and runs the unforced probe |
| `make_figures.jl` | Prepares manageable display tables, preserving local extrema and neighboring samples, then renders the figures |
| `figure_style.tex` | Shared figure size, fonts, colors, and axis styling |
| `01_*.tex` through `04_*.tex` | The four figure layouts |
| `results/*.csv` | Numerical data, including the full dense final-cycle table |
| `results/*summary.txt`, `visualization_metrics.txt` | Solver checks and numerical measurements |
| `figures/*.pdf`, `figures/*.png` | Final figures |

The figure data reduction does not rerun or change the dynamics. The full data remain in `last_cycle.csv`; `last_cycle_plot.csv` is the smaller table used for display. `energy_plot.csv` similarly reduces rendering cost while retaining per-bin extrema and adjacent samples. It avoids the memory limit encountered when repeatedly loading the full multicolumn data table into pdfLaTeX.

### Modify the physical model

Create a new immutable parameter record rather than editing derived quantities in place:

```julia
include("collision_model_corrected.jl")
using .CorrectedMEMS

p = Params(Vbias=2.0, seal_at_contact=false)
r = simulate(p=p, kind=:drive, cycles=8, freq=20.0,
             tag="open_tip_2V", outdir="my_results")

# simulate returns scaled solver states:
u_physical = r.sol(0.35) .* r.scale
```

For a customized set of all four plots, edit the parameter record passed to `simulate` in `run_visualization.jl`, and update the figure annotations as appropriate. The provided visualization script is a fixed reference-case example: its final-cycle interval, ramp annotation, selected contact wall, and plotted gc are currently written for 20 Hz, four ramp cycles, and gc = 13.71 μm. Do not change frequency or geometry and leave those annotations/selection settings untouched. A case without positive-wall contact will intentionally stop at the contact-selection check; use the general `simulate` outputs for that case.

## 4. What changed from your original model

The central improvement is consistency: one displacement field now determines inertia, gaps, capacitance derivatives, and pressure-force projection. This removes mismatches between constitutive laws that could otherwise change the predicted dynamics for numerical rather than physical reasons.

| Aspect | Earlier formulation | Corrected formulation and benefit |
|---|---|---|
| Beam geometry through contact | Separate translation and rotation branches with a blending region | One assumed beam shape and two independent coordinates throughout approach/contact. A pinned tip and moving shuttle naturally produce bending. |
| Inertia and base excitation | Empirical tip mass combined with other assumed kinematics | Full consistent mass matrix, including its off-diagonal term, and matching base-force vector. Mass and base-input work follow the same shape. |
| Dielectric coating | Lumped dielectric correction after integrating air-gap capacitance | Coatings enter each local strip before parallel integration. This gives a consistent capacitance and coordinate gradient. |
| Electrostatic force | Branch-dependent force construction | Both generalized forces come from the same capacitance gradient: Qi = 0.5 Vc^2 ∂Ct/∂xi. Mechanical work and capacitor energy are compatible. |
| Squeeze-film force | Inconsistent pressure quadrature and an unweighted rotational pressure resultant | Graded quadrature and virtual-work projection produce a symmetric positive-semidefinite damping matrix. The film cannot generate mechanical energy in this closure. |
| Contact loss | The unloading cutoff was omitted from part of the energy ledger | The exact loss associated with the clipped force and contact potential is included, even when unloading force is zero. |
| Structural coefficients | A stored cubic coefficient with an extra force factor; different stopper idealization | The actual cubic force coefficient is explicit, and a two-cantilever tip-loaded stopper is used under stated geometric assumptions. |
| Parameter updates | Mutable records could retain stale derived quantities | New immutable parameters trigger reconstruction of stiffness, inertia, shape, and quadrature. |
| Numerical verification | Successful integration did not by itself establish a consistent energy budget | Scaled states, solver-completion checks, an exact ledger, independent resistor-energy quadrature, and explicit convergence tests. |

Several changes are modeling choices that need validation, rather than proven physical improvements:

- The contact-travel datum is now gc = g0 − 2Tp − h_eff. This is appropriate if g0 denotes a bare mean gap; a different metrology datum requires a corresponding change.
- The default gap slope is based on the full beam length, not the overlap length, under the stated mirrored taper geometry. A measured gap profile should take precedence.
- Beam count and shuttle-only mass must be checked against the device layout to avoid double-counting mass or forces.
- The cubic coefficient assumes ideal guided spans sharing displacement equally. A serpentine suspension may require a measured or finite-element force law.
- Structural damping is set to zero as an uncalibrated reference, with fluid/contact loss retained. This does not establish that the real structural damping is zero or that the new power prediction is more accurate.
- The hard stop is disabled until its position is measured. Its removal from the reference calculation is not evidence that the fabricated structure has no hard stop.
- The sealing closure remains a candidate model. First-order slip at the narrowest reference gap, gas compressibility, transverse drainage, higher beam modes, and distributed contact remain unresolved physical-validity questions.

The changes therefore make the model more internally consistent, testable, and reproducible. Whether it predicts the fabricated device better must be established with measurements and comparisons against more detailed mechanics/fluid models.
