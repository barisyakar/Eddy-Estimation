# A Self-Exciting Model of Eddy Formation at Submesoscale

MATLAB code and figure assets accompanying the paper **"A Self-Exciting Model
of Eddy Formation at Submesoscale."** The model replaces
independent Poisson eddy births with a marked spatio-temporal Hawkes process.
Each event is described by its formation time, center, amplitude, radius, and
lifetime.

## Model and workflow

The constant-background model is fitted first. Its triggering kernel combines
the parent strain mark `|a_i|/b_i`, exponential temporal decay, Gaussian spatial
decay, and truncation at the parent lifetime. The code then:

1. estimates `(lambda0, K0, omega, rho)` by EM;
2. simulates catalogs with the immigration-birth representation;
3. compares simulated rates with the Volterra expected rate;
4. replaces the scalar background by a piecewise-constant `lambda0(t)` using
   4-hour and 8-hour bins; and
5. reconstructs velocity fields by superposing active eddies with the compactly
   supported basic-eddy profile defined in the paper.

The simulated eddy amplitude decays linearly during its lifetime:

```text
a_i(t) = a_i * (1 - (t - t_i)/l_i),  t_i <= t < t_i + l_i.
```

The two simulated velocity panels in the manuscript use consecutive 15-minute
blocks 252 and 253 of the first-window synthetic catalog. Their repository
filenames retain those snapshot identifiers:

```text
figures/simulated_191-252.pdf
figures/simulated_191-253.pdf
```

## Public repository contents

```text
.
|-- exp_kernel/                     Constant-background EM and simulation
|   `-- results/                    Two included seed-42 simulated catalogs
|-- time_varying_intensity/         Time-varying EM, simulation, and summaries
|-- figure_code/matlab/             Figure and velocity-field routines
|-- figures/                        PDFs referenced by the manuscript
|-- make_overleaf_figures.m         Figure-generation entry point
|-- CITATION.cff                    GitHub citation metadata
|-- .gitignore                      Public-release exclusions
`-- README.md
```

The included constant-background catalogs are:

```text
exp_kernel/results/Lifetime-first14days/validation/sim_catalog.mat
exp_kernel/results/Lifetime-last14days/validation/sim_catalog.mat
```

They contain 6,210 and 8,386 events, respectively, and store event times,
locations, marks, lifetimes, generations, parent indices, and simulation
parameters. Load one with:

```matlab
s = load(fullfile("exp_kernel", "results", ...
    "Lifetime-first14days", "validation", "sim_catalog.mat"));
catalog = s.catalog;
```

The numerical summaries used by the paper are provided in:

```text
time_varying_intensity/time_varying_fit_summary.csv
time_varying_intensity/simulation_summary.csv
```

The constant-background representative estimates are:

| Period | lambda0 | K0 | omega | rho |
|---|---:|---:|---:|---:|
| First 14 days | 4.9833e-4 | 2.1154 | 0.77657 | 0.020512 |
| Last 14 days | 7.1145e-4 | 0.79015 | 0.83823 | 0.019279 |

## Source data

The VHF-radar fields, processed observational event catalogs, and Mathematica
files are not included because they are subject to separate data-sharing
conditions. Re-estimation and regeneration of observation-derived figures
require authorized copies of:

```text
data/Lifetime-first14days.txt
data/Lifetime-last14days.txt
```

The repository includes the final figure PDFs so the paper can be compiled
without those source observations. The two simulated velocity panels can be
rebuilt directly from the included synthetic catalog.

## Running the code

Use MATLAB R2025a or a compatible recent release with the Statistics and
Machine Learning Toolbox. The Signal Processing Toolbox is also required for
the autocorrelation figures.

From the repository root, rebuild supported paper figures with:

```matlab
make_overleaf_figures
```

After placing authorized event catalogs in `data/`, run the constant-background
multi-start estimation from `exp_kernel/`:

```matlab
run_sweep
```

Run the 4-hour and 8-hour time-varying fits, simulations, Volterra solutions,
and figures from `time_varying_intensity/`:

```matlab
run_timevarying_bin_analysis
```

Generated fit workspaces, intermediate validation runs, and source observations
are excluded by `.gitignore`. The final PDFs, compact result summaries, MATLAB
source, and the two synthetic catalogs form the public release.
