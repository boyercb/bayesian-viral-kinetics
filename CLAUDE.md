# Bayesian Joint Model of Viral Kinetics — Project Status

*Last updated: March 1, 2026*

---

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

## 2. Candidate Model Fit (Run 3)

### Fit Metadata

| Field | Value |
|-------|-------|
| **Git commit** | `ec42f02` (main) |
| **Fit timestamp** | 2026-02-28 10:18:59 UTC |
| **Stan model** | `stan/kinetics_model.stan` |
| **GQ model** | `stan/kinetics_model_gq.stan` |
| **MCMC config** | 4 chains × 1000 warmup × 4000 sampling, adapt_delta=0.95, max_treedepth=12, 4 threads/chain |
| **Recovery config** | 4 chains × 1000 warmup × 2000 sampling, same tuning, 50% subsample |
| **Targets store** | `_targets/` (all 34 targets completed) |
| **Stan CSV dir** | `output/stan_csv/` (on remote: `christopherboyer/Dropbox/...`) |

### Key Targets and Timestamps (all 2026-02-28)

| Target | Timestamp | Runtime |
|--------|-----------|---------|
| `kinetics_mcmc` | 10:18:59 | 16.5 hr |
| `gq_fit` | 10:31:37 | 3.6 min |
| `convergence` | 10:27:17 | 6.9 min |
| `ppc` | 10:41:54 | 8.4 min |
| `loo_result` | 11:00:15 | 10.0 min |
| `waic_result` | 10:50:17 | 8.4 min |
| `predictions` | 10:19:59 | 1.0 min |
| `param_summary` | 10:33:07 | 0.7 min |
| `recovery_mcmc` | 17:15:34 | 6.3 hr |
| **Total pipeline** | — | **23.5 hr** |

### Convergence Summary

| Metric | Value |
|--------|-------|
| **Divergences** | 0 |
| **Max treedepth hits** | 0 |
| **E-BFMI** | 0.684, 0.703, 0.696, 0.688 (all healthy) |
| **Max Rhat** | 1.029 (`L_Omega_rna[4,3]`) |
| **Params Rhat > 1.01** | 17 (all RNA correlation matrix; none > 1.05) |
| **PFU params Rhat > 1.01** | **0** |
| **Min ESS bulk** | 183 (`L_Omega_rna[4,1]`) |
| **Min ESS tail** | 328 (`L_Omega_rna[4,3]`) |

**Verdict:** Production-quality fit. All PFU convergence issues resolved by NCP. Mild remaining ESS shortfall in RNA correlation matrix only.

### Model Comparison (LOO/WAIC)

| Metric | Value | SE |
|--------|-------|----|
| **LOO elpd** | -45,743 | 208 |
| **p_loo** | 4,736 | 56 |
| **WAIC elpd** | -45,320 | 209 |
| **p_waic** | 4,313 | 62 |
| **Pareto k > 1** | 658 (2.7%) | — |
| **Pareto k > 0.7** | 1,707 (6.9%) | — |

### Population Parameter Estimates

| Parameter | Estimate | 90% CI |
|-----------|----------|--------|
| Peak log-RNA ($d_p$) | 16.18 | (15.65, 16.73) |
| Proliferation rate ($w_p$) | 5.34 | (4.64, 6.18) |
| Clearance rate ($w_r$) | 13.74 | (12.42, 15.19) |
| RNA observation SD | 2.24 | (2.21, 2.27) |
| PFU observation SD | 2.22 | (2.01, 2.47) |
| False positive rate | 0.019 | (0.016, 0.023) |
| False negative rate | ~0 | (~0, ~0) |

**RNA-to-PFU transformation (log-affine):**

| Parameter | Estimate | 90% CI |
|-----------|----------|--------|
| PFU peak intercept ($a_0^{dp}$) | -0.93 | (-2.19, 0.30) |
| PFU peak elasticity ($a_1^{dp}$) | 1.07 | (0.65, 1.50) |
| PFU prolif intercept ($a_0^{wp}$) | -0.64 | (-1.33, 0.00) |
| PFU prolif elasticity ($a_1^{wp}$) | 0.77 | (0.45, 1.12) |
| PFU clear intercept ($a_0^{wr}$) | -0.02 | (-0.63, 0.56) |
| PFU clear elasticity ($a_1^{wr}$) | 0.66 | (0.44, 0.88) |

**Symptom model (cloglog hazard):**

| Parameter | Estimate | 90% CI |
|-----------|----------|--------|
| Intercept ($\eta_0$) | -1.47 | (-1.90, -1.02) |
| log-PFU effect ($\eta_1$) | 0.71 | (0.23, 1.31) |
| log-RNA effect ($\eta_2$) | 1.22 | (0.59, 1.85) |
| Symptom RE SD | 0.83 | (0.58, 1.13) |

**Individual RE standard deviations (RNA):**

| Parameter | Estimate | 90% CI |
|-----------|----------|--------|
| RE SD: peak time ($\tau_p$) | 0.239 | (0.193, 0.287) |
| RE SD: peak height ($d_p$) | 0.161 | (0.154, 0.168) |
| RE SD: proliferation ($w_p$) | 0.631 | (0.594, 0.669) |
| RE SD: clearance ($w_r$) | 0.508 | (0.486, 0.531) |

**RNA individual RE correlations:**

| Pair | Estimate | 90% CI |
|------|----------|--------|
| Corr($\tau_p$, $d_p$) | -0.22 | (-0.40, -0.04) |
| Corr($\tau_p$, $w_p$) | -0.54 | (-0.72, -0.35) |
| Corr($\tau_p$, $w_r$) | 0.28 | (0.12, 0.44) |
| Corr($d_p$, $w_p$) | -0.25 | (-0.33, -0.17) |
| Corr($d_p$, $w_r$) | 0.21 | (0.15, 0.28) |
| Corr($w_p$, $w_r$) | -0.38 | (-0.44, -0.30) |

### Parameter Recovery (Simulation Study)

- **Coverage:** 26/34 parameters (76.5%) covered by 95% CI
- **Non-covered (8):** `wp_mean_rna` (narrow miss), `sigma_rna` (tight CI), `sigma_pfu` (upward bias +0.48), `eta_sym_intercept` (bias toward 0), `fp` (overestimated 2×), `fn` (prior-dominated), `alpha_tcid50_log_b` (shrinkage), `alpha_cult_1` (not identified — too few observations)

### Covariate Effects (selected highlights)

- **Delta variant** increases peak RNA by 17% (1.12, 1.22)
- **Vaccination (boosted)** reduces peak RNA by 13% (0.84, 0.90) and prolongs proliferation by 45% (1.20, 1.75)
- **Recurrence** reduces peak by 4% (0.94, 0.99), shortens proliferation by 17% (0.73, 0.95), and shortens clearance by 24% (0.70, 0.83)
- **Age ≥50** increases clearance time by 19% (1.11, 1.29)
- **Unvaccinated unboosted** shows 32% slower clearance (0.60, 0.78)

### Generated Outputs (from this fit)

**Figures** (`output/figures/`):
- `ataccc_fit.pdf`, `hct_fit.pdf`, `legacy_fit.pdf`, `nba_fit.pdf`, `uiuc_fit.pdf` — individual trajectory fits by cohort
- `corr_densities.pdf`, `corr_matrix.pdf` — RNA RE correlation posteriors
- `forest_covariates.pdf` — covariate effect forest plots
- `param_recovery.pdf` — recovery simulation coverage
- `ppc_lfd_calibration.pdf`, `ppc_pfu.pdf`, `ppc_rna.pdf` — posterior predictive checks
- `prior_pe.pdf`, `prior_sym.pdf`, `prior_trajectories.pdf`, `prior_trans_pe.pdf` — prior predictive checks
- `tcl_example_1.pdf`, `tcl_example_2.pdf` — ODE comparison
- `trace_plots.pdf` — MCMC trace plots

**Tables** (`output/tables/`):
- `convergence.tex` — convergence diagnostics
- `params.tex` — parameter estimates
- `table1.tex` — data summary (Table 1)

---

## 3. Development History

### Convergence Journey

#### Run 1 (commit `4b9863`, ~Feb 2026)
- **Configuration:** Centered parameterization, PFU REs for all 2,261 individuals (9,044 params), wide priors
- **Result:** Catastrophic. E-BFMI < 0.07, max Rhat 2.46, 19% of params with Rhat > 1.01
- **Diagnosis:** Too many weakly-informed PFU individual effects; Neal's funnel geometry

#### Fixes Applied After Run 1
1. **Restrict PFU REs to PFU-informed individuals** (`b9548f6`): 9,044 → ~1,100 PFU RE params (N_pfu_ind=275)
2. **Tight prior on PFU RE SD** (`9ed8b70`): half-normal(0, 0.3) on sigma_ind_pfu
3. **Restore fmin constraint** (`72c94e5`): PFU ≤ RNA at all times
4. **Fix prediction indexing** (`2b15b93`): subscript out-of-bounds with PFU RE restriction

#### Run 2 (~Feb 2026)
- **Configuration:** PFU REs restricted, tight prior, centered parameterization, 2000 iterations, 1 thread/chain
- **Result:** Major improvement. E-BFMI 0.13–0.58, max Rhat 1.44, 570 params with Rhat > 1.01
- **Diagnosis:** Residual funnel in sigma_ind_pfu[1] (Rhat 1.44, ESS 8) — all 275 `tp_i_pfu` partially stuck
- **Runtime:** 11.0 hr

#### Fixes Applied After Run 2
1. **Non-centered parameterization (NCP) for all 4 PFU REs** (`2c46dca`): `tp_i_pfu ~ N(0, σ)` → `z_tp_pfu ~ N(0,1)`, reconstruct in transformed parameters. Breaks funnel coupling.
2. **Doubled sampling iterations** (`2c46dca`): 2000 → 4000 iter_sampling
3. **Threading** (`ec42f02`): threads_per_chain = 4 (exploits reduce_sum parallelism in Stan model)

#### Run 3 — Candidate Fit (2026-02-28, commit `ec42f02`)
- **Result:** Production-quality. E-BFMI 0.68–0.70, max Rhat 1.029, **17** params with Rhat > 1.01 (all RNA correlation), **zero PFU convergence issues**, min ESS bulk 183, min ESS tail 328
- **Runtime:** 16.5 hr main fit, 23.5 hr total pipeline

### Git History (chronological)

```
02e3070 Initial commit
f9046a7 initial setup
8b1fe7f updates
71c7f62 Restructure project to targets pipeline layout
0853101 Fix Stan to_int error: mark pfu_obs as data vector in partial_sum_ll
e46b16a fix: disable flat-top (use_wf=0), LKJ(2)->LKJ(4) for regularization
22fe215 feat: mechanistic interval-censored normal TCID50 model, fix sim bugs
1942404 docs: update TCID50 specification to mechanistic interval-censored normal
8387013 feat: sample_trajectories() for ABM agent trajectory generation
1e022a5 Add sigma_ind_pfu as learned hierarchical SD for PFU individual effects
31feaf7 Fix PFU trajectory scale: variance decomposition + RNA ceiling
72c94e5 Mode 2: residual PFU REs with propagated RNA variation + fmin constraint
b9548f6 Restrict PFU individual effects to PFU-informed individuals
9ed8b70 Add tighter half-normal(0, 0.3) prior on PFU RE hyperparameters
2b15b93 Fix subscript out of bounds in predictions with PFU RE restriction
2c46dca NCP for PFU individual effects; increase sampling to 4000 iterations
ec42f02 Set threads_per_chain = 4 for main and recovery runs
```

---

## 4. Current Model Specification

### Stan Model Flags (current settings in `_targets.R`)

| Flag | Value | Description |
|------|-------|-------------|
| `ind_effects` | 1 | Individual REs on peak time, peak height, proliferation, clearance |
| `test_error` | 1 | Mixture model for false positives/negatives |
| `adj_rna` | 1 | Covariate effects on RNA kinetics |
| `ind_corr` | 1 | Correlated RNA individual REs (4×4 Cholesky LKJ) |
| `use_smooth` | 1 | Smoothed piecewise exponential (vs. sharp corners) |
| `use_wf` | 0 | Flat-top/plateau period (disabled) |
| `adj_pfu` | — | **Not toggled on** (fully implemented in Stan) |
| `adj_lfd` | — | **Not toggled on** (fully implemented in Stan) |
| `adj_sym` | — | **Not toggled on** (fully implemented in Stan) |
| `source_pfu` | — | **Not toggled on** (fully implemented in Stan) |
| `source_rna` | — | **Not toggled on** (fully implemented in Stan) |
| `source_lfd` | — | **Not toggled on** (fully implemented in Stan) |
| `source_sym` | — | **Not toggled on** (fully implemented in Stan) |

### Parameterization Details

- **RNA trajectory:** Smoothed piecewise exponential with individual REs (correlated 4×4 via Cholesky NCP) + covariate effects on peak height, proliferation rate, clearance rate
- **PFU trajectory:** Derived from RNA via log-affine transformation + independent PFU individual REs (4 params × N_pfu_ind=275 individuals) using **non-centered parameterization** (`z_*_pfu ~ std_normal()`, reconstructed in transformed parameters)
- **PFU constraint:** PFU ≤ RNA enforced via `fmin(pfu_hat, rna_hat)`
- **Symptom onset:** Discrete-time cloglog hazard: $P(Y=1) = 1 - \exp(-\exp(\eta_0 + \eta_1 \log V + \eta_2 \log R + u_i))$
- **Observation models:** Normal(rna_hat, σ_rna) for qPCR; Normal(pfu_hat, σ_pfu) for PFU/TCID50; Bernoulli-logistic for LFD; test error mixture for RNA and PFU
- **Prior on PFU RE SD:** half-normal(0, 0.3) — `prior_pfu_i_sd = 0.3` set in `R/data.R`

---

## 5. Project Structure

```
_targets.R              ← Pipeline definition (34 targets)
R/                      ← Function definitions (auto-loaded by targets)
  utils.R               ← Utility helpers (pefun, ct_to_rna, calc_corners, etc.)
  data.R                ← Data cleaning & Stan data construction
  model.R               ← Model fitting, initialization (NCP), predictions
  diagnostics.R         ← Convergence, LOO-CV, WAIC, PPC
  summaries.R           ← Parameter & data summary tables
  predictions.R         ← Posterior prediction computation & plotting
  trajectories.R        ← ABM trajectory generation (sample_trajectories)
  plots.R               ← Prior predictive & ODE comparison plots
  simulation.R          ← Parameter recovery simulation study
stan/                   ← Active Stan models
  kinetics_model.stan   ← Main MCMC model (~670 lines, NCP for PFU REs)
  kinetics_model_gq.stan ← Generated quantities (log_lik, posterior predictive)
  sim.stan              ← Simulation study model
  _archive/             ← Old/inactive model variants
data/                   ← Datasets (read-only, not gitignored)
  raw/                  ← Original source files
output/                 ← Generated outputs (gitignored, reproduced by pipeline)
  figures/              ← 19 PDF figures
  tables/               ← 3 LaTeX tables
  stan_csv/             ← MCMC CSV output (on remote machine)
  stan_csv_gq/          ← GQ CSV output (on remote machine)
doc/                    ← Manuscripts and presentations
  manuscripts/
    main.tex            ← Primary manuscript
    setup.tex           ← LaTeX preamble
    supplement.tex      ← Supplementary material
    paper1_model/       ← Empty
    paper2_meta/        ← Partially started meta-analysis manuscript
    paper3_application/ ← Empty
  presentations/
    demo.tex            ← ~90-slide presentation
    demo.bib            ← Presentation bibliography
_archive/               ← Old imperative code (preserved for reference)
CLAUDE.md               ← This file
DESIGN_NOTE_IND_CORR.md ← Design plan for correlated PFU individual REs
DIAGNOSIS_PATHFINDER.md ← Notes on Pathfinder initialization issues
```

**Usage:** `targets::tar_make()` runs the full pipeline; `targets::tar_visnetwork()` shows the DAG; `targets::tar_read(target_name)` reads cached results.

---

## 6. What Has Been Completed

### Model Development
- [x] Piecewise exponential model for RNA with covariate effects
- [x] Log-affine RNA-to-PFU transformation (consistent across code, manuscript, presentation)
- [x] Discrete-time cloglog symptom hazard with individual RE
- [x] Observation models for all assay types (qPCR, PFU, TCID50, culture, LFD)
- [x] Test error mixture model (false positives/negatives)
- [x] Correlated RNA individual REs (4×4 Cholesky LKJ)
- [x] PFU individual REs restricted to PFU-informed individuals (N=275)
- [x] NCP for PFU individual effects (breaks funnel geometry)
- [x] All source/covariate toggle flags implemented in Stan model
- [x] Generated quantities model for log-lik and posterior predictive draws
- [x] reduce_sum parallelism (threading)

### Infrastructure
- [x] `targets` pipeline (34 targets, full DAG)
- [x] Pathfinder MAP initialization with fallback
- [x] Prior predictive checks
- [x] Posterior predictive checks (RNA, PFU, LFD)
- [x] LOO-CV and WAIC
- [x] Parameter recovery simulation study
- [x] Trace plots, correlation plots, forest plots, trajectory fits
- [x] LaTeX table generation (Table 1, parameters, convergence)
- [x] ABM trajectory generation (`sample_trajectories()`)

### Documentation
- [x] All R functions have roxygen-style headers
- [x] Design note for correlated PFU REs (`DESIGN_NOTE_IND_CORR.md`)
- [x] Pathfinder diagnosis notes (`DIAGNOSIS_PATHFINDER.md`)
- [x] This status document

---

## 7. Remaining Work — Phased Implementation Plan

### Phase 1: Immediate Post-Fit Tasks

**Goal:** Finalize the candidate fit, fix minor gaps, prepare for sensitivity analysis.

#### 1a. Investigate High Pareto-k Observations
- 658 observations (2.7%) have Pareto k > 1, indicating influential/poorly-fit points
- Extract and characterize these: which individuals, which data types (RNA/PFU/LFD/symptom), which time points, which cohorts
- Determine if they represent outliers, data quality issues, or model misspecification
- Consider moment-matching or integrated IS corrections (`loo::loo_moment_match()`)
- **Output:** Table/plot of problematic observations; decision on whether to address in model or document as limitation

#### 1b. Add `sigma_ind_pfu` to Parameter Summary
- `R/summaries.R` `summarize_parameters()` currently only extracts RNA RE SDs, not PFU RE SDs
- Add PFU RE SD estimates (`sigma_ind_pfu[1:4]`) to the `corr_params` output table
- Also add PFU RE summary to `output/tables/params.tex`

#### 1c. Document Data Exclusions and Calibration Decisions
- ATACCC hardcoded exclusions (IDs 12, 18, 23, 25, 41, 56) — document rationale in `R/data.R`
- Legacy `ct_to_rna` offset (-2) — document justification in `R/utils.R`
- HCT `hct-cn`/`hct-ct` calibration stubs — confirm and document that raw values are already copies/ml
- **Output:** Comments in code + a Data Notes section added to this document

#### 1d. Clean Up Stale Files
- Remove `_review_output.txt` from workspace root
- Remove `CLAUDE_old.md` if still present
- Verify `.gitignore` covers `output/`, `_targets/`, `_archive/`

### Phase 2: Manuscript Completion

**Goal:** Write all remaining manuscript sections using candidate fit results.

#### 2a. Empty Sections to Write

| Section | Content Needed |
|---------|---------------|
| **2.4 Symptoms** | Background on symptom diaries as infectiousness proxies, prior literature |
| **5 Computation** | Stan implementation details, NCP, threading, Pathfinder init, convergence criteria, software versions |
| **6 Model checking** | Prior/posterior predictive checks, LOO/WAIC, parameter recovery results |
| **7.1–7.2** | Individual-level trajectory fits, covariate effects interpretation |
| **7.3 Population parameters** | Population-level estimates table and interpretation |
| **7.4 Transmission models** | Application to testing/isolation policy (probability curves) |
| **8 Discussion** | Summary of findings, limitations, comparison to prior modeling approaches, future directions |

#### 2b. Supporting Manuscript Tasks
- Set up bibliography (`.bib` file — currently commented out in `main.tex`)
- Add figure/table captions in supplement
- Incorporate parameter recovery results as a figure/table
- Cross-reference all generated figures from `output/figures/`
- Update presentation (`doc/presentations/demo.tex`) with final results

#### 2c. Resolve Multi-Paper Structure
- `paper1_model/` — empty
- `paper2_meta/` — partially started, content mostly duplicates main manuscript
- `paper3_application/` — empty
- **Decision needed:** Is this one paper or a series? If one paper, consolidate and remove empty folders. If series, define scope boundaries.

### Phase 3: Sensitivity Analysis (Post-Submission)

**Goal:** Formally compare model variants via WAIC/LOO as a follow-up analysis. Can be added after initial submission; some variants may be redundant and the grid should be refined based on reviewer feedback or scientific priorities.

#### 3a. Build Sensitivity Analysis Infrastructure
- Create `R/sensitivity.R` with functions:
  - `run_sensitivity_model(label, flag_overrides, stan_data, ...)` — fits a model variant with specified flag overrides, returns fit + diagnostics
  - `compare_models(results_list)` — computes WAIC/LOO comparison table using `loo::loo_compare()`, returns formatted output
  - `save_sensitivity_result(label, fit, loo, waic, convergence)` — stores results with metadata to `output/sensitivity/`
- Add sensitivity targets to `_targets.R` (or a separate `_targets_sensitivity.R` pipeline)

#### 3b. Candidate Sensitivity Analysis Grid

Each variant toggles one or more flags from the base model. Not all may be necessary — prioritize based on scientific questions and reviewer feedback:

| Label | Change from base | Flags | Priority |
|-------|-----------------|-------|----------|
| `base` | Candidate fit (current) | — | — |
| `source_all` | All source REs | all source flags = 1 | High — tests cross-cohort heterogeneity |
| `adj_pfu` | Covariate effects on PFU | `adj_pfu = 1` | Medium — do covariates affect infectiousness directly? |
| `adj_sym` | Covariate effects on symptoms | `adj_sym = 1` | Medium — do covariates affect symptom onset? |
| `no_corr` | Independent RNA REs | `ind_corr = 0` | High — does correlation structure matter? |
| `no_smooth` | Sharp piecewise exponential | `use_smooth = 0` | Low — smoothing is mechanistically motivated |
| `source_rna` | Source REs on RNA only | `source_rna = 1` | Low — subsumed by `source_all` |
| `source_pfu` | Source REs on PFU only | `source_pfu = 1` | Low — subsumed by `source_all` |
| `source_lfd` | Source REs on LFD only | `source_lfd = 1` | Low — subsumed by `source_all` |
| `source_sym` | Source REs on symptoms only | `source_sym = 1` | Low — subsumed by `source_all` |
| `adj_lfd` | Covariate effects on LFD | `adj_lfd = 1` | Low — LFD data limited |
| `adj_all` | All covariate effects | all adj flags = 1 | Medium |
| `full` | All toggles on | all source + adj flags = 1 | Medium — kitchen sink model |

#### 3c. Run and Compare
- Run each variant (~16 hr each; parallelizable on cluster)
- For each: check convergence (reject if max Rhat > 1.1 or divergences > 0)
- GQ pass + LOO/WAIC for each converged fit
- Compare via `loo::loo_compare()` on WAIC and PSIS-LOO elpd
- **Output:** Model comparison table, elpd difference forest plot

#### 3d. Interpret and Select Final Model
- If any variant improves WAIC by >4 SE, consider adopting
- If source effects are non-negligible, document heterogeneity across cohorts
- Incorporate results into manuscript revision if needed

### Phase 4: Extended Model Features (Optional/Future)

#### 4a. Correlated PFU Individual REs
- Full design plan in `DESIGN_NOTE_IND_CORR.md`
- Would add 4×4 Cholesky correlation for PFU REs (mirroring RNA)
- Only meaningful if PFU data is rich enough (N_pfu_ind=275, but each has few observations)

#### 4b. K-fold Cross-Validation
- Already scaffolded in `_targets.R` (commented out)
- Grouped by individual (leave-one-individual-out)
- Budget: ~K × 16 hr — likely needs cluster
- Provides more robust model comparison than PSIS-LOO for the 2.7% of observations with high Pareto k

#### 4c. Additional Applications
- Testing/isolation policy probability curves (partially done in presentation)
- ABM integration via `sample_trajectories()` (infrastructure exists in `R/trajectories.R`)
- Transmission model coupling (chain binomial sketch exists in `doc/manuscripts/chain_binomial.R`)
- Meta-analytic summary for paper 2 (`doc/manuscripts/paper2_meta/`)

---

## 8. Known Issues and Limitations

### Data Issues
- `hct-cn` and `hct-ct` calibration functions in `ct_to_rna()` are stubs — HCT raw values used directly (assumed already copies/ml)
- ATACCC data joins symptom data from separate file; `sym_exist` initially 0, updated post-join — potential edge cases
- ATACCC hardcoded exclusions (IDs 12, 18, 23, 25, 41, 56) lack documentation
- Legacy `ct_to_rna` uses NBA calibration with -2 offset — workaround needs justification

### Model Limitations
- 2.7% of observations have Pareto k > 1 (influential points)
- Recovery: `sigma_pfu` biased upward (+0.48); `fp` overestimated 2×; `alpha_cult_1` not identified (too few culture observations)
- PFU RE SDs not currently in `param_summary` output
- Source-level and extended covariate effects not yet evaluated (deferred to Phase 3, post-submission)
- Individual RE correlation for PFU not implemented (designed but deferred to Phase 4)

### Infrastructure
- Stan CSV files live on remote machine (`christopherboyer/Dropbox/...`); `fit$draws()` and `fit$summary()` calls fail locally — all diagnostics must use cached `targets` objects
- Recovery fit convergence not independently verifiable locally (same CSV path issue)

---

## 9. Technical Reference

### Running the Pipeline
```r
library(targets)
tar_make()              # Run full pipeline
tar_visnetwork()        # Visualize DAG
tar_read(kinetics_mcmc) # Read cached fit (draws only, not CSV-dependent methods)
tar_read(convergence)   # Convergence diagnostics
tar_read(param_summary) # Parameter estimates
tar_read(loo_result)    # LOO-CV
tar_read(waic_result)   # WAIC
tar_read(ppc)           # Posterior predictive check draws
tar_read(predictions)   # Posterior predictions (obs + grid)
tar_read(recovery_check) # Parameter recovery
```

### Key Parameter Name Mapping

| Concept | Stan parameter | R extraction |
|---------|---------------|--------------|
| RNA peak (population) | `dp_mean_rna` | `param_summary$pop_params` |
| RNA proliferation (population) | `wp_mean_rna` | `param_summary$pop_params` |
| RNA clearance (population) | `wr_mean_rna` | `param_summary$pop_params` |
| PFU transformation | `tau_dp[1:2]`, `tau_wp[1:2]`, `tau_wr[1:2]` | `param_summary$transformation_params` |
| PFU individual RE (NCP z-score) | `z_tp_pfu[j]` | draws (parameter) |
| PFU individual RE (reconstructed) | `tp_i_pfu[j]` | draws (transformed parameter) |
| RNA individual RE (NCP z-score) | `z_ind_rna[k,i]` | draws (parameter) |
| RNA individual RE (reconstructed) | `eff_rna_mat[k,i]` | draws (transformed parameter) |
| RNA RE SDs | `sigma_ind_rna[1:4]` | `param_summary$corr_params` |
| PFU RE SDs | `sigma_ind_pfu[1:4]` | *not yet in param_summary* |
| RNA RE correlations | `Omega_rna[i,j]` | `param_summary$corr_params` |
| Symptom model | `eta_sym_intercept`, `eta_sym_pfu`, `eta_sym_rna`, `sigma_sym` | `param_summary$symptom_params` |
| Covariate effects | `beta_dp_rna`, `beta_wp_rna`, `beta_wr_rna` | `param_summary$covariate_effects` |
| Observation error | `sigma_rna`, `sigma_pfu`, `fp`, `fn` | `param_summary$error_params` |
