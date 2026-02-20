# Comprehensive Status Report: Bayesian Joint Model of Viral Kinetics

*Generated: February 10, 2026*

## 1. Project Overview

**Research objective:** Build a Bayesian generative model for the joint distribution of within-host viral shedding markers during acute SARS-CoV-2 infection:

$$f(V_t, R_t, Y_t, \mathbf{Z}_t \mid \mathbf{X}, t, S; \theta)$$

factorized into four components: (1) infectious virus trajectory, (2) viral RNA trajectory conditional on infectious virus, (3) symptom onset conditional on both, and (4) observation models linking latent trajectories to measured biomarkers (qPCR Ct values, viral culture, LFD rapid tests, symptom diaries).

**Data sources:** Five prospective, longitudinally-sampled cohorts (~2,000 infections):

| Cohort | N infections | Ct/RNA | Viral culture | LFD | Symptoms |
|--------|-------------|--------|---------------|-----|----------|
| NBA | ~1,989 | Yes | No | No | No |
| ATACCC | ~57 | Yes | Yes | Yes | Yes |
| UIUC | ~60 | Yes | Yes (TCID50) | Yes | Yes |
| HCT | 18 | Yes | Yes (PFU+FFA) | Yes | Yes |
| Legacy | ~157 | Yes | No | No | Yes |

---

## 2. What Has Been Done

### Theory & Manuscript

The manuscript (`4_manuscripts/main.tex`) has a **well-developed theory section** covering:

- Introduction and motivation (complete)
- Background on proxies for infectiousness (complete for culture, RNA, antigens; **Symptoms subsection is empty**)
- The full model specification including:
  - Piecewise exponential model for infectious virus (V_t)
  - Linked piecewise exponential for viral RNA (R_t)
  - Symptom onset model (discrete-time hazard)
  - Observation models for all assay types (culture, TCID50, PFU/FFA, qPCR, LFD)
  - Covariate effects and individual random effects
  - Missing data framework

### Two Stan Models Implemented

1. **`kinetics_model.stan`** (357 lines) — an earlier/simpler version:
   - Single piecewise exponential for RNA with covariate effects
   - RNA-to-PFU transformation via multiplicative scaling (ρ parameters)
   - LFD modeled as logistic function of RNA_hat and PFU_hat
   - Test error mixture model for RNA
   - No symptom modeling
   - No dedicated individual effects for PFU trajectory
   - A saved `.rds` fit object exists for this model

2. **`kinetics_model.stan`** (~670 lines) — the **current/evolved version**:
   - Separate individual random effects for both RNA and PFU trajectories
   - RNA-to-PFU transformation via **log-affine form**: log δ' = a₀ + a₁ · log δ (equivalently δ' = e^{a₀} · δ^{a₁}), **consistent across code, manuscript, and presentation**
   - Symptom onset via **discrete-time cloglog hazard** with individual random effect (NCP): P(Y=1) = 1 - exp(-exp(η₀ + η₁ log V + η₂ log R + u_i)), **consistent across code, manuscript, and presentation**
   - Source-level random effects (toggleable) for PFU, RNA, LFD, and symptoms
   - Multiple viral culture assay types (PFU, TCID50, simple culture)
   - Test error mixture for both RNA and PFU
   - Full set of model specification flags, all implemented:
     - `test_error` — mixture model for false positives/negatives ✅
     - `ind_effects` — individual REs on peak, proliferation, clearance ✅
     - `adj_pfu` — covariate effects on PFU kinetics ✅ (toggled off)
     - `adj_rna` — covariate effects on RNA kinetics ✅ (toggled on)
     - `adj_lfd` — covariate effects on LFD positivity ✅ (toggled off)
     - `adj_sym` — covariate effects on symptom hazard ✅ (toggled off)
     - `source_pfu` — source-level REs on PFU ✅ (toggled off)
     - `source_rna` — source-level REs on RNA ✅ (toggled off)
     - `source_lfd` — source-level REs on LFD ✅ (toggled off)
     - `source_sym` — source-level REs on symptoms ✅ (toggled off)
   - Planned extension: correlated individual random effects (`ind_corr`) — see `DESIGN_NOTE_IND_CORR.md`

### Code Pipeline

The master script (`__master_run.R`) orchestrates:

1. **`functions.R`** — helper functions including `ct_to_rna()` calibrations, `pefun()`, `predict_kinetics()`, `prior_predictive()`
2. **`clean_data.R`** — data harmonization for all 5 cohorts into a stacked dataset + Stan-ready list
3. **`prior_predictive.R`** — prior predictive checks (plots generated)
4. **`fit_model.R`** — compiles `kinetics_model.stan`, runs MCMC (4 chains, 1000 warmup, 2000 sampling), saves to `.rds`
5. **`check_model.R`** — **essentially empty** (has placeholder comments only)
6. **`model_summaries.R`** — extracts posterior summaries and creates LaTeX tables (references `kinetics_model.stan` variables like `rho_dp`, not the `kinetics_model` variables like `tau_dp`)
7. **`prediction.R`** — generates posterior predictive trajectories + individual-level fit plots (also references **old** model variables)

### Generated Outputs

- **Figures (`3_figures/`):** Prior predictive plots, trace plots, individual fit plots for ATACCC/NBA/UIUC/HCT, TCL ODE examples, testing policy probability plots
- **Supplement tables:** Posterior estimates for peak, proliferation, and clearance (appear to come from the **older** `kinetics_model` fit)
- **Presentation (`demo.tex`):** ~90 slides with full model description, data tables, prior checks, select results, and application analysis (testing/isolation policy)

### Side Projects

- **Simulation study (`1_code/simstudy/`):** Compares target-cell limited ODE model vs. piecewise exponential vs. smoothed piecewise — appears complete as a standalone investigation
- **`run_ct_model.R`:** Fits a Ct-only model to NBA data — earlier standalone analysis
- **`chain_binomial.R`:** Toy chain binomial transmission model — appears to be a sketch
- **`paper2_meta/`:** A separate manuscript framed as a meta-analysis/IPD analysis — partially started (data/methods sections with content recycled from paper 1; empty results/discussion)
- **`paper3_application/`:** Empty folder

---

## 3. Critical Gaps & Mismatches

### A. Model Code vs. Manuscript Misalignment

| Component | Manuscript | Generative Stan Model | Status |
|---|---|---|---|
| RNA-to-PFU transformation | log δ' = a₀ + a₁ · log δ (log-affine) | log δ' = a₀ + a₁ · log δ (log-affine) | **Resolved** — now consistent |
| Symptom onset | Discrete-time cloglog hazard with individual RE: P(Y=1) = 1-exp(-exp(η₀ + η₁ log V + η₂ log R + u_i)) | Discrete-time cloglog hazard with NCP individual RE | **Resolved** — now consistent |
| Covariate effects on PFU | Manuscript implies β' effects on RNA-derived PFU params | Code has `adj_pfu = 0` (disabled but **fully implemented**) | Ready to toggle on |
| Covariate effects on symptoms | Can be added to cloglog linear predictor | Code has `adj_sym = 0` (disabled but **fully implemented**) | Ready to toggle on |
| Source random effects (RNA/PFU) | Described in manuscript | Toggled off but **fully implemented** (`source_pfu = 0`, `source_rna = 0`) | Ready to toggle on |
| Source random effects (LFD) | Implied | Toggled off but **fully implemented** (`source_lfd = 0`) | Ready to toggle on |
| Individual RE correlation | Mentioned conceptually | Commented-out; full design plan in `DESIGN_NOTE_IND_CORR.md` | Planned — not yet implemented |

### B. Post-processing Code Uses Wrong Model

`1_code/model_summaries.R` and `1_code/prediction.R` reference parameter names from `kinetics_model.stan` (e.g., `rho_dp`, `rho_wp`, `theta_tp`, `log_dpi`, `beta_dp`, `gamma`, `gamma0`) rather than from `kinetics_model.stan` (which uses `tau_dp`, `tau_wp`, `tau_tp`, `dp_i_rna`, `beta_dp_rna`, `tau_lfd`, `tau0_lfd`). This means:

- **If you re-fit with `kinetics_model.stan`, all summary/prediction code will break.**
- The existing results in the supplement/presentation likely come from the **older** `kinetics_model` fit.

### C. Incomplete Manuscript Sections

- **Section 2.4 (Symptoms):** Empty background text (the model specification itself is complete)
- **Section 5 (Computation):** Empty
- **Section 6 (Model checking and inference):** Empty
- **Section 7 (Results):** Only contains two figure references (probability plots) and subsection headers
- **Section 7.3 (Population parameters):** Empty
- **Section 7.4 (Transmission models):** Empty
- **Section 8 (Discussion):** Empty
- **Bibliography:** No `.bib` file referenced (commented out)
- **Figure/table captions:** All empty in supplement
- Various placeholder text ("Figure X", "between X and Y%")

### D. Incomplete Code

- ~~`1_code/check_model.R` is **essentially empty**~~ **DONE.** Replaced by `R/diagnostics.R` with `check_convergence()`, `compute_loo()`, `posterior_predictive_check()`, `plot_traces()`
- ~~No `generated quantities` block in either Stan model for posterior predictive checks or LOO-CV~~ **DONE.** `generated quantities` block added to `kinetics_model.stan` with observation-level `log_lik` for LOO-CV and posterior predictive draws (`rna_rep`, `pfu_rep`, `lfd_rep`, `sym_rep`)
- ~~The `predict_kinetics()` function in `1_code/functions.R` only works with `kinetics_model.stan` parameters~~ **DONE.** `predict_kinetics()` in `R/model.R` updated to use `kinetics_model.stan` parameter names
- Legacy dataset calibration (`ct_to_rna` type "nba" with a -2 offset) is a workaround that should be documented/justified
- HCT calibration functions (`hct-cn`, `hct-ct`) are **empty** — the raw qPCR values are used directly (already in copies/ml)

### E. Data Issues

- `hct-cn` and `hct-ct` calibration functions in `ct_to_rna()` have no implementation — just comments
- The ATACCC data joins symptom data from a separate file but `sym_exist` is set to 0 initially and only updated after the join — there may be edge cases
- Some hardcoded exclusions (ATACCC IDs 12, 18, 23, 25, 41, 56) lack documentation
- ~~The NBA plot code references source 2 (ATACCC) for the NBA section — likely a copy-paste error~~ **DONE.** `plot_source_trajectories()` takes `source_id` as parameter

---

## 4. Recommendations

### Modeling & Theory

1. ~~**Reconcile model specifications.**~~ **DONE.** RNA-to-PFU transformation unified to log-affine form across Stan model, manuscript, functions.R, and presentation.
2. ~~**Upgrade symptom model.**~~ **DONE.** Replaced random-effect offset with discrete-time cloglog hazard with individual RE, consistent across all files.
3. ~~**Add `generated quantities` block**~~ **DONE.** Block added to `kinetics_model.stan` with: (a) observation-level `log_lik` array that accumulates across RNA/PFU/LFD/symptom components for LOO-CV, (b) posterior predictive draws (`rna_rep`, `pfu_rep`, `lfd_rep`, `sym_rep`) faithful to the observation model including LOD censoring, assay-specific PFU types, and cloglog symptom hazard.
4. ~~**Consider correlated individual random effects.**~~ **Designed.** Full implementation plan in `DESIGN_NOTE_IND_CORR.md`. Recommends starting with 4×4 RNA-only correlation via Cholesky NCP.
5. ~~**Source-level effects should be evaluated.**~~ **All implemented.** `source_lfd` was the last gap — now fully implemented. All source flags (`source_pfu`, `source_rna`, `source_lfd`, `source_sym`) are ready to toggle on.
6. ~~**Covariate effects on symptoms.**~~ **Implemented.** `adj_sym` flag added with `beta_sym` coefficients on the cloglog linear predictor.

### Code Quality & Organization

1. ~~**Consolidate Stan models.**~~ **DONE.** Project restructured to `targets`-idiomatic layout. Active Stan model (`kinetics_model.stan`) lives in top-level `stan/`. All old models (`kinetics_model.stan`, test files, archived variants) moved to `stan/_archive/` and `_archive/code/`.
2. ~~**Update `model_summaries.R` and `prediction.R`**~~ **DONE.** `predict_kinetics()` in `R/model.R` now uses correct parameter names from `kinetics_model.stan` (e.g., `tau0_dp`/`tau_dp` instead of `rho_dp`/`theta_tp`). Parameter summaries in `R/summaries.R` updated to match.
3. ~~**Implement `check_model.R`**~~ **DONE.** `R/diagnostics.R` implements `check_convergence()` (R-hat, ESS, divergences), `compute_loo()` (PSIS-LOO via `loo` package on `log_lik`), `posterior_predictive_check()` (extracts rna_rep/pfu_rep/lfd_rep), and `plot_traces()`.
4. ~~**Add documentation.**~~ **DONE.** All functions in `R/` have roxygen-style headers with `@param`/`@return`. Each file has a purpose header.
5. ~~**Remove dead code.**~~ **DONE.** Old imperative scripts archived to `_archive/code/`. New `R/` files contain only clean, active function definitions.
6. ~~**Fix the copy-paste error**~~ **DONE.** `plot_source_trajectories()` in `R/predictions.R` now takes `source_id` as a parameter — no hardcoded source filtering.
7. **Document data exclusions and calibration decisions** — especially the ATACCC exclusions and the Legacy offset. *(Still TODO — requires domain knowledge writeup)*

### Project Structure (migrated to `targets`)

The project has been restructured from numbered folders + `__master_run.R` to a `targets`-based pipeline:

```
_targets.R              ← Pipeline definition (replaces __master_run.R)
R/                      ← Function definitions (auto-loaded by targets)
  utils.R               ← Utility helpers (pefun, ct_to_rna, calc_corners, etc.)
  data.R                ← Data cleaning & Stan data construction
  model.R               ← Model fitting, predictions, prior predictive
  diagnostics.R         ← Convergence, LOO-CV, posterior predictive checks
  summaries.R           ← Parameter & data summary tables
  predictions.R         ← Posterior prediction computation & plotting
  plots.R               ← Prior predictive & ODE comparison plots
  simulation.R          ← Simulation study (piecewise vs target-cell)
stan/                   ← Active Stan models only
  kinetics_model.stan
  sim.stan
  _archive/             ← Old/inactive models
data/                   ← Datasets (read-only)
  raw/                  ← Original source files
output/                 ← Generated outputs (gitignored, reproduced by pipeline)
  figures/
  tables/
doc/                    ← Manuscripts and presentations
  manuscripts/
  presentations/
_archive/               ← Old code (gitignored, preserved for reference)
```

Usage: `targets::tar_make()` runs the full pipeline; `targets::tar_visnetwork()` shows the DAG.

### Manuscript

1. **Write the empty sections** (Symptoms background, Computation, Model checking, Results, Discussion).
2. **Set up bibliography** — currently no `.bib` file is linked.
3. **Add figure/table captions** in the supplement.
4. **Resolve the multi-paper structure.** There are three paper folders (`paper1_model`, `paper2_meta`, `paper3_application`) but only `paper2_meta` has content, and it largely duplicates the main manuscript. Decide whether this is one paper or a series.
5. **Update the presentation** to reflect the final model specification once reconciled.

### Priority Roadmap (suggested)

1. ~~Reconcile Stan model with manuscript theory~~ **DONE**
2. ~~Update all post-processing code to work with the active model~~ **DONE** (refactored to `R/` functions with correct param names)
3. ~~Implement model checking and diagnostics~~ **DONE** (`R/diagnostics.R` + `generated quantities` block)
4. Complete results analysis (run `targets::tar_make()`)
5. Write remaining manuscript sections
6. ~~Clean up code documentation and organization~~ **DONE** (`targets` pipeline + roxygen headers)
7. Evaluate toggling on extensions (source effects, covariate effects, correlated REs)
