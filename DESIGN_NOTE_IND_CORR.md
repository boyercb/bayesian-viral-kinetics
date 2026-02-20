# Design Note: Correlated Individual Random Effects (`ind_corr`)

*Created: February 10, 2026*

## 1. Motivation

The model currently implements **independent** individual random effects for each kinetic parameter (peak δ, proliferation ωₚ, clearance ωᵣ) for both the RNA and PFU trajectories. That is, for individual *i*:

```
log δ_i = log Δ₀ + δ₀ᵢ,    δ₀ᵢ ~ N(0, σ_δ)
log ωₚᵢ = log Ωₚ₀ + ωₚ₀ᵢ,  ωₚ₀ᵢ ~ N(0, σ_ω_p)
log ωᵣᵢ = log Ωᵣ₀ + ωᵣ₀ᵢ,  ωᵣ₀ᵢ ~ N(0, σ_ω_r)
```

and similarly for the PFU parameters. These are modeled as uncorrelated. In reality, individuals who shed more virus at peak (high δ) likely also shed for longer (high ωᵣ), and proliferation and clearance dynamics may be correlated through shared immune function. Ignoring these correlations:

- Inflates posterior widths on derived quantities (e.g., total viral shedding = area under curve)
- Prevents borrowing of information across parameters within an individual
- Cannot answer questions like "conditional on a high peak, what is the expected clearance time?"

## 2. Random Effects Structure

There are **up to 9** individual random effects currently (when `ind_effects = 1`):

| Index | Parameter | Trajectory | Current variable |
|-------|-----------|------------|-----------------|
| 1 | Peak time | RNA | `tp_i_rna` |
| 2 | Peak height | RNA | `dp_i_rna` |
| 3 | Proliferation | RNA | `wp_i_rna` |
| 4 | Clearance | RNA | `wr_i_rna` |
| 5 | Peak time | PFU | `tp_i_pfu` |
| 6 | Peak height | PFU | `dp_i_pfu` |
| 7 | Proliferation | PFU | `wp_i_pfu` |
| 8 | Clearance | PFU | `wr_i_pfu` |
| 9 | Symptom susceptibility | Symptom | `z_sym` (currently NCP, σ_sym * z) |

A full 9×9 correlation matrix is likely overparameterized given the data. Practical options:

### Option A: Block-diagonal correlation

Separate correlation matrices for RNA (4×4), PFU (4×4), and symptoms (1×1, trivially diagonal). This assumes RNA and PFU individual effects are independent conditional on the RNA→PFU transformation, which is reasonable since the transformation already links them.

- **RNA block**: `(tp_i_rna, dp_i_rna, wp_i_rna, wr_i_rna) ~ MVN(0, Σ_RNA)`
- **PFU block**: `(tp_i_pfu, dp_i_pfu, wp_i_pfu, wr_i_pfu) ~ MVN(0, Σ_PFU)`
- **Total Cholesky parameters**: 2 × (4×3/2) = 12 correlation parameters + 8 scale parameters = 20 parameters

### Option B: RNA-only correlation (recommended starting point)

Only correlate the RNA trajectory random effects (4×4). PFU individual effects remain independent. Rationale: the RNA data is densest (all 5 cohorts contribute), and the PFU data is sparse (only 3 cohorts with viral culture).

- **RNA block**: `(tp_i_rna, dp_i_rna, wp_i_rna, wr_i_rna) ~ MVN(0, Σ_RNA)`  
- **Everything else**: independent as now
- **Total new parameters**: 4×3/2 = 6 correlation parameters + 4 scale parameters (replace existing σ's) = 10 parameters net

### Option C: Reduced-rank correlation (factor model)

Instead of estimating the full correlation matrix, use a factor model:
```
η_i = Λ f_i + ε_i,   f_i ~ N(0, I_q),   ε_i ~ N(0, diag(ψ))
```
where Λ is K×q with q << K (e.g., q=2 for a 2-factor model of 9 random effects). This is more scalable but harder to interpret.

**Recommendation: Start with Option B.**

## 3. Stan Implementation Plan

### 3.1 Data block

```stan
int<lower=0, upper=1> ind_corr; // should individual effects be correlated?
```

### 3.2 Parameters block

When `ind_corr = 1`, replace the independent `dp_i_rna`, `wp_i_rna`, `wr_i_rna`, `tp_i_rna` arrays with a matrix of standardized effects and a Cholesky factor:

```stan
// Option B: correlated RNA individual effects (non-centered)
// Dimension D_rna = 4 (tp, dp, wp, wr for RNA)
cholesky_factor_corr[ind_corr ? 4 : 0] L_Omega_rna;
vector<lower=0>[ind_corr ? 4 : 0] sigma_ind_rna;  // RE standard deviations
matrix[ind_corr ? 4 : 0, ind_corr ? sum(M) : 0] z_ind_rna;  // standardized effects
```

### 3.3 Transformed parameters block

```stan
if (ind_corr) {
  // Cholesky factor of the covariance matrix: L = diag(sigma) * L_Omega
  matrix[4, 4] L_Sigma_rna = diag_pre_multiply(sigma_ind_rna, L_Omega_rna);
  
  // Correlated individual effects (non-centered): eta = L * z
  matrix[4, sum(M)] eta_rna = L_Sigma_rna * z_ind_rna;
  
  // Unpack into existing variable names
  for (i in 1:sum(M)) {
    tp_i_rna[i] = eta_rna[1, i];
    dp_i_rna[i] = eta_rna[2, i];
    wp_i_rna[i] = eta_rna[3, i];
    wr_i_rna[i] = eta_rna[4, i];
  }
} 
// else: use the existing independent draws (no change needed)
```

### 3.4 Model block (priors)

```stan
if (ind_corr) {
  L_Omega_rna ~ lkj_corr_cholesky(2);  // LKJ(2): mild shrinkage toward identity
  sigma_ind_rna ~ normal(0, prior_i_sd) T[0, ];
  to_vector(z_ind_rna) ~ std_normal();  // NCP standardized effects
}
```

### 3.5 Generated quantities (optional)

Recover the correlation matrix for reporting:
```stan
generated quantities {
  corr_matrix[ind_corr ? 4 : 0] Omega_rna;
  if (ind_corr) {
    Omega_rna = multiply_lower_tri_self_transpose(L_Omega_rna);
  }
}
```

## 4. Integration with Existing Code

### 4.1 Parameter sizing challenge

Currently, `tp_i_rna`, `dp_i_rna` etc. are declared separately with different sizing rules:
- `tp_i_rna` is always size `sum(M)` (unconditional)
- `dp_i_rna` is size `sum(M)` only when `ind_effects = 1`

When `ind_corr = 1`, **all 4** RNA individual effects must exist (you need something to correlate). 
This implies `ind_corr` requires `ind_effects = 1`. Add a validation check:

```stan
if (ind_corr && !ind_effects) reject("ind_corr requires ind_effects = 1");
```

### 4.2 Interaction with `ind_effects`

When `ind_corr = 0, ind_effects = 1`: current behavior (independent REs)
When `ind_corr = 1, ind_effects = 1`: correlated RNA REs, independent PFU REs
When `ind_corr = 0, ind_effects = 0`: no individual effects (current behavior)
When `ind_corr = 1, ind_effects = 0`: invalid, reject

### 4.3 Changes to `functions.R`

The `prior_predictive()` function would need a correlated draw branch:
```r
if (data$ind_corr) {
  L_Omega_rna <- rethinking::rlkjcorr(draws, 4, eta = 2)  # or use direct Cholesky
  sigma_ind_rna <- truncnorm::rtruncnorm(draws, 0, Inf, 0, data$prior_i_sd)
  for (d in 1:draws) {
    L <- diag(sigma_ind_rna[d, ]) %*% t(chol(L_Omega_rna[d,,]))
    z <- matrix(rnorm(4 * sum(data$M)), nrow = 4)
    eta <- L %*% z
    tp_i_rna[, d] <- eta[1, ]
    dp_i_rna[, d] <- eta[2, ]
    wp_i_rna[, d] <- eta[3, ]
    wr_i_rna[, d] <- eta[4, ]
  }
}
```

### 4.4 Changes to `clean_data.R`

Add to Stan data list:
```r
ind_corr = 0,  # toggled off by default
```

### 4.5 Changes to `model_summaries.R` / `prediction.R`

Extract `Omega_rna` from posterior and report as correlation matrix table/heatmap.

## 5. Computational Considerations

- **LKJ prior**: `lkj_corr_cholesky(2)` with η=2 provides mild regularization toward the identity. For a 4×4 matrix this is computationally manageable.
- **Non-centered parameterization**: Essential. The centered version (directly sampling from MVN) creates strong posterior correlations between Σ and the random effects, causing divergent transitions. The NCP (`L * z` with `z ~ N(0,1)`) avoids this.
- **Warm-up**: May need longer warm-up (1500–2000) for the adaptation phase to learn the mass matrix for the 6 new correlation parameters.
- **Identifiability**: With ~2000 individuals contributing RNA data, a 4×4 correlation matrix should be well-identified. Monitor `R̂` and `n_eff` for the off-diagonal correlation terms.
- **Scaling**: The `z_ind_rna` matrix is 4 × sum(M), which is ~4 × 2000 = 8000 parameters. This is large but each column is iid, so HMC handles it well with NCP.

## 6. Staged Rollout Plan

1. **Stage 1**: Implement Option B (RNA-only, 4×4) with `ind_corr = 0` (toggled off). Verify the model compiles and produces identical results when off.
2. **Stage 2**: Toggle on (`ind_corr = 1`) and run on a subset of data (e.g., NBA only) to verify convergence.
3. **Stage 3**: Run on full dataset. Compare ELPD/LOO-CV with vs. without correlation.
4. **Stage 4** (optional): Extend to Option A (block-diagonal RNA + PFU) if the RNA correlation structure is informative and the model fits stably.
