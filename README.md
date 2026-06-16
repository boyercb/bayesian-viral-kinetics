# Bayesian Viral Kinetics

Bayesian joint model for within-host SARS-CoV-2 viral kinetics and infectiousness inference, integrating viral RNA (qPCR), infectious virus culture (PFU/TCID50/FFA), antigen (LFD), and symptom onset across five cohorts.

This repository contains the analysis pipeline, Stan models, and figure/table generation code used for the manuscript:

**Inferring infectiousness: a joint model of the within-host viral kinetics of SARS-CoV-2**

## Overview

Direct infectiousness is hard to observe. This project models a latent infectious-virus trajectory jointly with multiple observed proxies to:

- infer infectious virus dynamics in individuals with only PCR data,
- estimate population-level infectiousness curves and durations,
- evaluate policy-relevant quantities (for example, residual risk by isolation day),
- produce individualized Bayesian filtering updates as new test history accumulates.

The implementation uses:

- `targets` for reproducible pipelines,
- `cmdstanr` + Stan (NUTS) for Bayesian inference,
- modular R scripts in `R/` for data processing, diagnostics, prediction, policy analysis, and publication figures.

## Repository structure

- `_targets.R`: end-to-end pipeline definition and main execution entrypoint.
- `R/`: analysis functions and plotting modules.
- `stan/`: core Stan model files (`kinetics_model.stan`, `kinetics_model_gq.stan`, `sim.stan`).
- `data/`: harmonized cohort input data used by the analysis.
- `output/`: generated artifacts (figures, tables, posterior draw objects, CmdStan CSVs).
- `sample_trajectories.qmd`: worked example for sampling realistic trajectories from fitted posterior draws.
- `doc/`: local manuscript and submission materials (ignored by git in the public-facing repository by default).
- `_archive/`: legacy project material retained for local provenance.

## Data included

The pipeline reads harmonized cohort-level files from `data/`:

- `nba_dat.csv`
- `ataccc_dat.csv`, `ataccc_sx_dat.csv`
- `uiuc_dat.csv`, `uiuc_sx_dat.csv`
- `hct_dat.csv`, `hct_sx_dat.csv`
- `legacy_dat.csv`

Raw source files and metadata snapshots are present under `data/raw/` for provenance.

Important:

- Respect original study data-use agreements and licenses.
- Before public redistribution of raw files, confirm all permissions and de-identification requirements.

For release and reuse guidance, see `DATA_POLICY.md`.

## Requirements

- R (recommended >= 4.2)
- CmdStan toolchain (via `cmdstanr`)
- C++ toolchain compatible with CmdStan

Core R packages used by the pipeline include:

- `targets`, `tarchetypes`, `stantargets`
- `tidyverse`, `purrr`
- `cmdstanr`, `posterior`, `loo`, `bayesplot`
- `deSolve`, `mvtnorm`, `truncnorm`
- `knitr`, `kableExtra`
- `patchwork`, `colorspace`, `viridis`
- `qs2` (used for direct object reads in some scripts)

## Setup

From an R session at repository root:

```r
install.packages(c(
  "targets", "tarchetypes", "stantargets", "tidyverse", "purrr",
  "cmdstanr", "posterior", "loo", "bayesplot", "deSolve",
  "mvtnorm", "truncnorm", "knitr", "kableExtra", "patchwork",
  "colorspace", "viridis", "qs2"
))

# Install CmdStan if needed (one-time)
cmdstanr::install_cmdstan()
```

## Run the full pipeline

```r
# In R, from repository root
library(targets)
tar_make()
```

Useful targets commands:

```r
tar_visnetwork()   # visualize DAG
tar_outdated()     # list targets needing rebuild
tar_read(stan_data)
```

## Pipeline outputs

The full pipeline produces:

- fitted posterior objects (`kinetics_mcmc`, `gq_fit`),
- diagnostics (`convergence`, `loo_result`, `waic_result`, PPC objects),
- summary tables in `output/tables/`,
- figures in `output/figures/`,
- posterior draw snapshots (for example `output/pop_draws_200.rds`),
- persistent CmdStan CSV files in `output/stan_csv/` and `output/stan_csv_gq/`.

## Reproduce manuscript and policy figures

Most figures are generated through `_targets.R` targets, including:

- trajectory fit figures,
- diagnostics panels,
- publication figures,
- policy analysis figures,
- Bayesian filtering figures,
- site-specific comparison figure,
- prior predictive and ODE illustration figures.

You can also run selected scripts directly when targets outputs already exist. Example:

```bash
Rscript R/pub_figures.R all
```

This script writes journal-style figure files to `output/figures/{style}/`.

## Sampling trajectories for simulation

A demonstration is provided in `sample_trajectories.qmd`.

Render it with Quarto:

```bash
quarto render sample_trajectories.qmd
```

The document shows how to:

- load fitted posterior objects,
- sample individual-level RNA/PFU/LFD/symptom trajectories,
- compare trajectory distributions across covariate profiles.

## Runtime and compute notes

- Full model fitting is computationally intensive.
- The manuscript reports a representative runtime of about 16.5 hours wall time on a 16-core workstation for the main model configuration.
- A generated-quantities-only Stan model is used post hoc to compute log-likelihood and posterior predictive quantities without inflating sampling output size.

## Reproducibility notes

- The canonical workflow is the `targets` pipeline in `_targets.R`.
- Generated artifacts in `output/` and `_targets/` are expected to be rebuildable.
- CmdStan CSV outputs are deliberately persisted to make fitted objects robust across serialization and path changes.
- Parameter recovery simulation defaults to full-size data (`recovery_frac = 1.0`) and supports optional per-source stratified subsampling for faster debug runs.
- Recovery coverage should be interpreted from replicate-aggregated summaries (`recovery_coverage_summary`), not a single simulated dataset.

## Computational efficiency

The pipeline uses a two-stage Stan workflow to balance inference quality with memory efficiency:

1. **Sampling phase** (`kinetics_model.stan`): Runs MCMC to draw from the posterior. The model computes the joint log-likelihood during sampling but does not store per-observation contributions to keep MCMC output lean (~4,000 iterations per chain × 4 chains = 16,000 draws across all model parameters).

2. **Post-hoc diagnostics phase** (`kinetics_model_gq.stan`): After sampling, a lightweight generated-quantities-only model recomputes diagnostics (per-observation log-likelihoods, posterior predictive replicates, fitted values) from the saved posterior draws. This avoids storing ~250,000 extra columns per iteration during sampling, reducing output size ~5–10×.

This approach is standard practice in production Bayesian workflows and enables efficient leave-one-out cross-validation (LOO-CV via `loo::loo`) and posterior predictive checks without compromising the sampling phase.

<!-- ## Data policy

This repository contains harmonized analysis inputs and raw-source files used for data provenance.

- Do not assume all raw files can be redistributed outside this repository context.
- Before rehosting or republishing datasets, verify study-level permissions and legal requirements.
- See `DATA_POLICY.md` for details and recommended practice.

## Contributing

Contributions are welcome. See `CONTRIBUTING.md` for:

- setup expectations,
- minimal reproducibility checks,
- pull request scope and review notes.

## Citation

If you use this codebase, cite the manuscript and repository metadata in `CITATION.cff`.

Current manuscript citation target:

Boyer CB, Kissler SM, Hakki S, Jonnerby J, Lalvani A, Lipsitch M.
*Inferring infectiousness: a joint model of the within-host viral kinetics of SARS-CoV-2*. -->

## License

This repository is distributed under the license in `LICENSE`.
