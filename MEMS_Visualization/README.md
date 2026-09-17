# MEMS model visualization: native Julia run

These figures come from a fresh run of the delivered corrected Julia model, with its physical equations and default parameters unchanged.

## Reference experiment

- 20 Hz sinusoidal base acceleration, 4.95 g peak, 3 V bias.
- Eight forcing cycles (0.4 s), including a four-cycle cosine amplitude ramp.
- 512 graded film panels; relative tolerance 1e-7 and scaled absolute tolerance 1e-10.
- A separate unforced 600 microsecond near-contact initial-condition probe.
- No experimentally calibrated parameter adjustment is made for these plots.

## Reproduce

Use Julia 1.11.7 and the included environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run_visualization.jl
julia --project=. make_figures.jl
```

The first command installs the simulation dependencies. The driven calculation resolves fast contact motion inside slow forcing cycles, so allow several minutes or longer depending on your computer. The second script prepares plot data in Julia and invokes pdfLaTeX/PGFPlots to produce vector PDFs. A LaTeX distribution including PGFPlots is required. Poppler (`pdftoppm`) is optional and adds 300 dpi PNG previews. Existing figures can be opened without these tools. No Python is needed.

## Figures

1. `01_driven_overview`: shuttle/tip motion, output voltage, capacitance, and integrated resistor energy over the full ramp and driven response.
2. `02_contact_closeup`: 0.1-microsecond sampling around positive-tip contact entry in the last forcing cycle; displacements, velocities, voltage, and signed generalized tip forces.
3. `03_last_cycle`: phase portrait, voltage, instantaneous load power, and capacitance over the last cycle. Data combine uniform times with every accepted solver time. The plotting script keeps local extrema and their neighboring samples to reduce file size; the full data are retained separately.
4. `04_energy_and_probe`: stored energy/source work, cumulative losses, numerical energy residual, and the separate initialized contact-probe voltage.

Each is supplied as a vector PDF and a 300 dpi PNG. The CSV data and Julia scripts are included.

## Read the outputs correctly

- The nominal contact travel is 13.71 micrometers. Positive tip overlap is allowed by compliant contact; it is not rigid-wall penetration in a constraint solver.
- Shuttle and tip displacements can look identical at the full-device scale while their nanometer-scale difference drives beam bending.
- The broad overview samples at 50 microseconds and can miss brief peaks. Use the dense final-cycle data and contact close-up for those features. Reported peak metrics use accepted solver states and are sampled maxima, not optimized continuous-time extrema.
- Late mean load power is calculated from integrated resistor energy over the last four cycles. Do not estimate it from coarse plotted voltage samples.
- `Wbias` is electrical source work, not shaker work. The energy ledger includes both sources.
- The unforced probe starts with a specified displacement and velocity close to contact. Its larger transient voltage is not the operating voltage of the driven response.
- Numerical agreement and apparent recurrence do not establish experimental accuracy, Floquet stability, or uniqueness of the attractor. Contact/gas calibration and a driven convergence study remain necessary for publication-level predictions.

See `results/visualization_metrics.txt` and `results/visualization_summary.txt` for achieved quantities and solver checks. The original corrected source file is copied unchanged into this package.
