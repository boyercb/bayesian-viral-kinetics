# Diagnosis: Pathfinder Initialization Failure

**Date:** February 19, 2026

## Problem

After adding several model features (smooth trajectories, flat-top `wf` parameter,
correlated individual random effects via Cholesky NCP, tighter test-error priors,
and MCMC speed optimizations including pathfinder initialization), the model
produced catastrophic MCMC results:

- **7,998 / 8,000 divergent transitions** (99.975%)
- **R-hat up to 4.54**, ESS ≈ 4 (one effective sample per chain)
- `dp_mean_rna ≈ 0.2` instead of the expected ~17
- `sigma_rna ≈ 6.7` (residual noise absorbing all signal)
- `total_ll ≈ -63,500` (far worse than the correct ~-42,000)

## Root Cause: Pathfinder

**Pathfinder completely fails on this model.** Diagnostic run confirmed:

- Only **2 of 8 pathfinder paths succeeded** (6 failed to initialize)
- Hundreds of `"Non-finite gradient"` and `"Non-finite function evaluation"` errors
- `lkj_corr_cholesky_lpdf: Random variable[k] is 0, but must be positive!` exceptions
- **0% of pathfinder draws have `dp_mean_rna > 5`** (should be ~17)
- Best pathfinder `lp__ = -95,193` (vs. correct mode `lp__ ≈ -42,000`)
- Pathfinder converges to degenerate mode: `dp_raw ≈ -5` (dp ≈ 0.5), `sigma_rna ≈ 11.5`

Because all 4 MCMC chains were initialized from pathfinder's best draws, they all
started in this degenerate mode. Warmup couldn't escape — it adapted the mass matrix
and step size around the wrong mode, producing frozen chains with near-zero mixing.

## Bisection Test

Ran 4 model configurations (1 chain × 50 warmup + 10 sampling each) using
`build_init()` (prior-mode initialization, bypassing pathfinder):

| Config | ind_corr | use_smooth | use_wf | dp_mean_rna | sigma_rna | total_ll | Divergences | Max TD hits | Time (s) |
|--------|----------|------------|--------|-------------|-----------|----------|-------------|------------|----------|
| baseline | 0 | 0 | 0 | **17.47** | 2.21 | -42,083 | 0/10 | 10/10 | 1,076 |
| +smooth | 0 | 1 | 0 | **15.63** | 2.24 | -42,254 | 0/10 | 10/10 | 1,482 |
| +smooth+wf | 0 | 1 | 1 | **18.55** | 2.25 | -42,308 | 0/10 | 10/10 | 1,518 |
| +smooth+wf+corr | 1 | 1 | 1 | **16.16** | 2.23 | -42,390 | 0/10 | 10/10 | 1,527 |

**All four configs work correctly** with `build_init()`:
- `dp_mean_rna` in the correct range (15–19)
- `sigma_rna ≈ 2.2` (reasonable)
- `total_ll ≈ -42,000` (correct mode)
- **Zero divergences**
- 100% max treedepth hits → need `max_treedepth ≥ 12`

### Conclusion

- **Smooth trajectories (`use_smooth=1`):** No issues.
- **Flat-top parameter (`use_wf=1`):** No issues.
- **Correlated individual effects (`ind_corr=1`):** No issues when properly initialized.
- **Pathfinder (`use_pathfinder=TRUE`):** **Root cause of failure.** Cannot navigate the
  posterior landscape of this model (non-finite gradients in large parameter regions,
  LKJ Cholesky boundary violations during optimization).

## Fix Applied

1. **Disabled pathfinder** (`use_pathfinder = FALSE` default in `fit_model()`)
2. **Increased `max_treedepth`** from 10 → 12 (all configs saturated at 10)
3. **Using `build_init()`** which starts chains at prior modes (dp_raw = 0 → dp = 17)

Files changed:
- `R/model.R`: `fit_model()` defaults to `use_pathfinder = FALSE`, `max_treedepth = 12`
- `_targets.R`: `kinetics_mcmc` and `recovery_mcmc` targets pass
  `max_treedepth = 12, use_pathfinder = FALSE`

## Archived Diagnostic Scripts

Located in `_archive/diagnostics/`:
- `_test_bisect.R` — Bisection test (4 configs × short MCMC runs)
- `_diag.R` — Per-chain parameter inspection
- `_diag_pathfinder.R` — Pathfinder mode analysis
- `_test_nopf.R` — MCMC test without pathfinder (full model)
- `_invalidate.R` — Target cache invalidation helper
