#include "functions.stan"

data {
  // fixed parameters --------------------------------------------------------
  int<lower=0> K; // number of sources
  array[K] int<lower=0> M; // number of individuals per source
  array[K] int<lower=0> N; // number of observations per source
  array[K] real lod_rna; // limit of detection of qPCR by source
  array[K] real lod_pfu; // limit of detection for viral culture by source
  array[K] real fp_mean; // mean RNA of false positive
  
  // observed data -----------------------------------------------------------
  array[sum(N)] int id; // participant id
  array[sum(N)] real time; // time of observation
  vector[sum(N)] rna; // RNA copies per ml from PCR
  vector[sum(N)] pfu; // PFU of culturable virus 
  array[sum(N)] int<lower=0, upper=1> lfd; // rapid antigen test result
  array[sum(N)] int<lower=0, upper=1> sym; // daily symptom indicator (0/1)
  array[sum(N)] int<lower=0, upper=1> sym_ever; // ever symptomatic
  array[sum(N)] int<lower=0, upper=1> sym_at_risk; // at-risk for onset (no prior onset)

  array[sum(N)] int<lower=0> source; // source indicator
  array[sum(N)] int<lower=0> pfu_type; // source indicator
  
  int<lower=0> P; // number of covariates
  array[sum(M)] row_vector[P] x; // covariate vectors

  // grainsize for reduce_sum (within-chain threading).
  // Use 1 to let Stan auto-schedule, or tune to ~50-200 for throughput.
  // Set threads_per_chain > 1 in $sample() to activate parallelism.
  int<lower=1> grainsize;

  // PFU individual effects: restrict to individuals with PFU-informing data
  // (culture, LFD, or symptom observations).  Non-informed individuals
  // (e.g. NBA with RNA-only) get population-level RNA→PFU transform only.
  int<lower=0> N_pfu_ind;                        // number of PFU-informed individuals
  array[sum(M)] int<lower=0> pfu_ind_idx;        // mapping: 0=no RE, 1..N_pfu_ind=index

  // missingness indicators
  array[sum(N)] int<lower=0, upper=1> rna_exist; 
  array[sum(N)] int<lower=0, upper=1> pfu_exist; 
  array[sum(N)] int<lower=0, upper=1> lfd_exist; 
  array[sum(N)] int<lower=0, upper=1> sym_exist; 

  // model specification -----------------------------------------------------
  int<lower=0, upper=1> test_error; // should test error be included?
  int<lower=0, upper=1> ind_effects; // should individual effects be modeled?
  int<lower=0, upper=1> ind_corr; // should RNA individual effects be correlated (4×4 Cholesky)?
  
  int<lower=0, upper=1> adj_pfu; // should PFU be adjusted?
  int<lower=0, upper=1> adj_rna; // should RNA be adjusted?
  int<lower=0, upper=1> adj_lfd; // should antigen test positivity be adjusted?
  int<lower=0, upper=1> adj_sym; // should symptom onset be adjusted

  int<lower=0, upper=1> source_pfu; // should source effects be modeled?
  int<lower=0, upper=1> source_rna; // should source effects be modeled?
  int<lower=0, upper=1> source_lfd; // should source effects be modeled?
  int<lower=0, upper=1> source_sym; // should source effects be modeled?

  int<lower=0, upper=1> use_wf;     // include flat-top duration parameter?
  int<lower=0, upper=1> use_smooth;  // use smooth trajectory (vs piecewise)?

  // priors ------------------------------------------------------------------
  real<lower=0> prior_dp_mean;
  real<lower=0> prior_dp_cv;
  real<lower=0> prior_wp_mean;
  real<lower=0> prior_wp_cv;
  real<lower=0> prior_wr_mean;
  real<lower=0> prior_wr_cv;
  
  real<lower=0> prior_i_sd;
  real<lower=0> prior_pfu_i_sd;  // tighter prior SD for PFU RE hyperparameters
  real<lower=0> prior_k_sd;
  real<lower=0> prior_beta_sd;

  real<lower=0> prior_sigma_sd;
  real<lower=0, upper=1> prior_lfd_mean;  
  
  real prior_fn; // prior mean on false negative probability
  real prior_fp; // prior mean on false positive probability
  real<lower=0> prior_fp_kappa; // precision (concentration) for fp/fn Beta priors

  real<lower=0> prior_wf_mean;  // prior mean for flat-top duration (days)
  real<lower=0> prior_wf_cv;    // prior CV for flat-top duration (NCP scale)

  // Steepness for smooth post-peak sigmoid: inv_logit(kappa * (t - tp)).
  // Replaces the hard indicator I(t >= tp) with a differentiable
  // approximation so HMC can propagate gradients through tp.
  // kappa ~ 5 ≈ transition over ~1 day; kappa ~ 10 ≈ ~0.5 day.
  real<lower=0> kappa_postpeak;

}

transformed data {
  array[sum(N)] int<lower=0, upper=1> first_obs_per_id;
  first_obs_per_id = first_instance(id);

  // Scale factor for normalising viral-load covariates in the cloglog
  // symptom hazard.  Dividing by prior_dp_mean (~17 log copies/ml)
  // keeps the linear predictor O(1) and prevents exp(exp(eta)) overflow.
  real scale_vl = prior_dp_mean;
}

parameters {
  // observation noise for [RNA] and PFUs of culturable virus
  real<lower=0> sigma_rna;
  real<lower=0> sigma_pfu;
  
  // symptom onset discrete-time hazard (cloglog link)
  real zeta_sym_intercept; // baseline log-hazard
  real<lower=0> zeta_sym_pfu; // log V_t coefficient
  real<lower=0> zeta_sym_rna; // log R_t coefficient
  real zeta_sym_postpeak;     // post-peak indicator
  real zeta_sym_postpeak_rna; // post-peak × log R_t interaction
  real<lower=0> sigma_sym; // individual random effect SD
  array[sum(M)] real z_sym; // NCP random effects (std normal)

  // test characteristics 
  real<lower=0,upper=1> fp; // false positive probability
  real<lower=0,upper=1> fn; // false negative probability
  
  // population effects: RNA (non-centered)
  real dp_raw; // peak
  real wp_raw; // proliferation time
  real wr_raw; // clearance time
  array[use_wf ? 1 : 0] real wf_raw; // flat-top duration (NCP)
  
  // individual effects: RNA
  // When ind_corr = 1, tp/dp/wp/wr are replaced by z_ind_rna (NCP Cholesky)
  array[(1 - ind_corr) * sum(M)] real tp_i_rna; // onset (independent)
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real dp_i_rna; // peak (independent)
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real wp_i_rna; // proliferation (independent)
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real wr_i_rna; // clearance (independent)
  array[use_wf && ind_effects ? sum(M) : 0] real wf_i; // flat-top individual RE (always independent)
  
  // Correlated RNA individual effects (NCP Cholesky)
  // Active when ind_corr = 1; D_rna = 4: (tp, dp, wp, wr)
  cholesky_factor_corr[ind_corr ? 4 : 1] L_Omega_rna;
  vector<lower=0>[ind_corr ? 4 : 0] sigma_ind_rna;  // RE standard deviations
  matrix[ind_corr ? 4 : 0, ind_corr ? sum(M) : 0] z_ind_rna;  // standardized effects
  
  // source effects: RNA
  array[K && source_rna ? K : 0] real tp_k_rna; // onset
  array[K && source_rna ? K : 0] real dp_k_rna; // peak
  array[K && source_rna ? K : 0] real wp_k_rna; // proliferation time
  array[K && source_rna ? K : 0] real wr_k_rna; // clearance time
  array[use_wf && source_rna ? K : 0] real wf_k; // flat-top source RE
  
  // coefficient vectors for predictors of RNA peak, proliferation, and clearance 
  vector[P && adj_rna ? P : 0] beta_dp_rna;
  vector[P && adj_rna ? P : 0] beta_wp_rna;
  vector[P && adj_rna ? P : 0] beta_wr_rna;
  
  // transformation parameters
  vector[2] tau_tp; // onset 
  vector[2] tau_dp; // peak 
  vector[2] tau_wp; // proliferation time 
  vector[2] tau_wr; // clearance time 

  // individual effects: PFU (restricted to PFU-informed individuals)
  // NCP: sample z-scores, reconstruct actual effects in transformed parameters
  vector<lower=0>[ind_effects ? 4 : 1] sigma_ind_pfu;  // PFU RE SDs (learned from informed individuals)
  array[N_pfu_ind] real z_tp_pfu; // NCP z-score for onset
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_dp_pfu; // NCP z-score for peak
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_wp_pfu; // NCP z-score for proliferation
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_wr_pfu; // NCP z-score for clearance
  
  // source effects: PFU
  array[K && source_pfu ? K : 0] real tp_k_pfu; // onset
  array[K && source_pfu ? K : 0] real dp_k_pfu; // peak
  array[K && source_pfu ? K : 0] real wp_k_pfu; // proliferation time
  array[K && source_pfu ? K : 0] real wr_k_pfu; // clearance time
  
  // coefficient vectors for predictors of PFU peak, proliferation, and clearance 
  vector[P && adj_pfu ? P : 0] beta_dp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wr_pfu;
  
  // TCID50 interval-censored regression: [1]=a (intercept),
  // [2]=log(b) (inverse culture growth rate), [3]=log(sigma) (scale)
  vector[3] alpha_tcid50;
  
  // coefficient for binary culture result to PFUs 
  vector[2] alpha_cult;
  
  // non-centered parameterization for logistic intercept for LFD
  real tau0_lfd_raw;
  
  // coefficient vector for LFD positivity.
  // [1]=RNA, [2]=PFU, [3]=I(post-peak), [4]=I(post-peak)*RNA
  vector[4] tau_lfd;
  
  // coefficient vectors for predictors of LFD positivity
  vector[P && adj_lfd ? P : 0] beta_lfd;
  
  // source effects: LFD
  array[K && source_lfd ? K : 0] real lfd_k; // source-level LFD intercept shifts
  
  // coefficient vectors for predictors of symptom onset hazard
  vector[P && adj_sym ? P : 0] beta_sym;
  
  // source effects: symptom onset
  array[K && source_sym ? K : 0] real to_k_sym;
  
}

transformed parameters {
  // -----------------------------------------------------------------------
  // Population means — computed ONCE per iteration (not per observation).
  // Saved to output for posterior summaries.
  // -----------------------------------------------------------------------
  real<lower=0> dp_mean_rna = prior_dp_mean * exp(prior_dp_cv * dp_raw);
  real<lower=0> wp_mean_rna = prior_wp_mean * exp(prior_wp_cv * wp_raw);
  real<lower=0> wr_mean_rna = prior_wr_mean * exp(prior_wr_cv * wr_raw);

  // PFU population means (log-affine of clamped RNA means)
  real<lower=0> dp_mean_pfu = log_affine(tau_dp[1], tau_dp[2], safe_pos(dp_mean_rna));
  real<lower=0> wp_mean_pfu = log_affine(tau_wp[1], tau_wp[2], safe_pos(wp_mean_rna));
  real<lower=0> wr_mean_pfu = log_affine(tau_wr[1], tau_wr[2], safe_pos(wr_mean_rna));

  // NCP intercept for LFD
  real tau0_lfd = logit(prior_lfd_mean) + 1 * tau0_lfd_raw;

  // ── NCP reconstruction for PFU individual effects ─────────────────────────
  // tp_i_pfu[i] = sigma_ind_pfu[1] * z_tp_pfu[i], etc.
  array[N_pfu_ind] real tp_i_pfu;
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real dp_i_pfu;
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real wp_i_pfu;
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real wr_i_pfu;
  for (j in 1:N_pfu_ind) tp_i_pfu[j] = sigma_ind_pfu[1] * z_tp_pfu[j];
  if (ind_effects) {
    for (j in 1:N_pfu_ind) {
      dp_i_pfu[j] = sigma_ind_pfu[2] * z_dp_pfu[j];
      wp_i_pfu[j] = sigma_ind_pfu[3] * z_wp_pfu[j];
      wr_i_pfu[j] = sigma_ind_pfu[4] * z_wr_pfu[j];
    }
  }

  // -----------------------------------------------------------------------
  // Total log-likelihood accumulator.
  // Replaces per-observation log_lik arrays — avoids writing ~100 K
  // columns per MCMC iteration.  Per-observation log_lik and posterior
  // predictive draws are computed in a post-hoc generated-quantities pass.
  // -----------------------------------------------------------------------
  // ── Pre-compute effective individual RNA effects (4 × N_ind matrix) ────────
  // Rows = [tp, dp, wp, wr]; handles correlated (Cholesky NCP) and independent.
  matrix[4, sum(M)] eff_rna_mat = rep_matrix(0.0, 4, sum(M));
  if (ind_corr) {
    eff_rna_mat = diag_pre_multiply(sigma_ind_rna, L_Omega_rna) * z_ind_rna;
  } else {
    for (i in 1:sum(M)) eff_rna_mat[1, i] = tp_i_rna[i];
    if (ind_effects) {
      for (i in 1:sum(M)) {
        eff_rna_mat[2, i] = dp_i_rna[i];
        eff_rna_mat[3, i] = wp_i_rna[i];
        eff_rna_mat[4, i] = wr_i_rna[i];
      }
    }
  }

  // ── Pre-compute population flat-top duration ──────────────────────────────
  real wf_pop = use_wf ? prior_wf_mean * exp(prior_wf_cv * wf_raw[1]) : 0.0;

  // ── Log-likelihood via reduce_sum (within-chain thread parallelism) ────────
  // Each thread evaluates partial_sum_ll on a contiguous slice of observations.
  // Set threads_per_chain > 1 in $sample() and compile with stan_threads=TRUE.
  real total_ll = reduce_sum(
    partial_sum_ll, id, grainsize,
    // per-observation
    time, rna, pfu, lfd, sym, sym_at_risk, sym_exist,
    rna_exist, pfu_exist, lfd_exist, pfu_type, source,
    // individual-level
    eff_rna_mat,
    tp_i_pfu, dp_i_pfu, wp_i_pfu, wr_i_pfu,
    pfu_ind_idx,
    z_sym, wf_i, x,
    // source-level
    lod_rna, lod_pfu, fp_mean,
    tp_k_rna, dp_k_rna, wp_k_rna, wr_k_rna, wf_k,
    tp_k_pfu, dp_k_pfu, wp_k_pfu, wr_k_pfu,
    lfd_k, to_k_sym,
    // covariate betas
    beta_dp_rna, beta_wp_rna, beta_wr_rna,
    beta_dp_pfu, beta_wp_pfu, beta_wr_pfu,
    beta_lfd, beta_sym,
    // population / transformation parameters
    dp_mean_rna, wp_mean_rna, wr_mean_rna,
    tau_tp, tau_dp, tau_wp, tau_wr,
    tau0_lfd, tau_lfd,
    sigma_rna, sigma_pfu,
    fp, fn,
    zeta_sym_intercept, zeta_sym_pfu, zeta_sym_rna,
    zeta_sym_postpeak, zeta_sym_postpeak_rna,
    sigma_sym, alpha_tcid50, alpha_cult,
    wf_pop,
    // flags
    ind_effects,
    source_rna, source_pfu, source_lfd, source_sym,
    adj_rna, adj_pfu, adj_lfd, adj_sym,
    test_error, use_wf, use_smooth,
    scale_vl,
    kappa_postpeak
  );
}

model {
  // PRIORS
  fp ~ beta_proportion(prior_fp, prior_fp_kappa);
  fn ~ beta_proportion(prior_fn, prior_fp_kappa);

  sigma_rna ~ normal(0, prior_sigma_sd) T[0, ];
  sigma_pfu ~ normal(0, prior_sigma_sd) T[0, ];

  dp_raw ~ std_normal();
  wp_raw ~ std_normal();
  wr_raw ~ std_normal();

  if (use_wf) {
    wf_raw ~ std_normal();
    if (ind_effects) { wf_i ~ normal(0, prior_i_sd); }
    if (source_rna)  { wf_k ~ normal(0, prior_k_sd); }
  }

  // Validation: ind_corr requires ind_effects
  if (ind_corr && !ind_effects) reject("ind_corr requires ind_effects = 1");

  sigma_ind_pfu ~ normal(0, prior_pfu_i_sd) T[0, ];  // half-normal prior on PFU RE SDs (tight)
  z_tp_pfu ~ std_normal();  // NCP z-scores for PFU individual effects
  z_sym ~ std_normal();  // NCP for symptom random effects

  // RNA individual effects: correlated (Cholesky NCP) or independent
  if (ind_corr) {
    L_Omega_rna ~ lkj_corr_cholesky(4);  // LKJ(4): moderate shrinkage toward identity
    sigma_ind_rna ~ normal(0, prior_i_sd) T[0, ];
    to_vector(z_ind_rna) ~ std_normal();
  } else {
    tp_i_rna ~ normal(0, prior_i_sd);
    if (ind_effects) {
      dp_i_rna ~ normal(0, prior_i_sd);
      wp_i_rna ~ normal(0, prior_i_sd);
      wr_i_rna ~ normal(0, prior_i_sd);
    }
  }

  if (ind_effects) {
    z_dp_pfu ~ std_normal();
    z_wp_pfu ~ std_normal();
    z_wr_pfu ~ std_normal();
  }

  if (source_pfu) {
    tp_k_pfu ~ normal(0, prior_k_sd);
    dp_k_pfu ~ normal(0, prior_k_sd);
    wp_k_pfu ~ normal(0, prior_k_sd);
    wr_k_pfu ~ normal(0, prior_k_sd);
  }

  if (source_rna) {
    tp_k_rna ~ normal(0, prior_k_sd);
    dp_k_rna ~ normal(0, prior_k_sd);
    wp_k_rna ~ normal(0, prior_k_sd);
    wr_k_rna ~ normal(0, prior_k_sd);
  }

  if (source_sym) {
    to_k_sym ~ normal(0, prior_k_sd);
  }

  if (adj_pfu) {
    beta_dp_pfu ~ normal(0, prior_beta_sd);
    beta_wp_pfu ~ normal(0, prior_beta_sd);
    beta_wr_pfu ~ normal(0, prior_beta_sd);
  }

  if (adj_rna) {
    beta_dp_rna ~ normal(0, prior_beta_sd);
    beta_wp_rna ~ normal(0, prior_beta_sd);
    beta_wr_rna ~ normal(0, prior_beta_sd);
  }

  if (adj_lfd) {
    beta_lfd ~ normal(0, prior_beta_sd);
  }

  if (source_lfd) {
    lfd_k ~ normal(0, prior_k_sd);
  }

  if (adj_sym) {
    beta_sym ~ normal(0, prior_beta_sd);
  }

  tau_tp ~ std_normal();

  // symptom onset discrete-time hazard priors (cloglog link)
  zeta_sym_intercept ~ normal(-3, 1);  // low baseline daily hazard ~exp(-3)≈0.05
  zeta_sym_pfu ~ normal(0, 0.5) T[0, ];  // more virus -> more symptoms
  zeta_sym_rna ~ normal(0, 0.5) T[0, ];  // more RNA -> more symptoms
  zeta_sym_postpeak ~ std_normal();       // post-peak level shift
  zeta_sym_postpeak_rna ~ std_normal();   // post-peak × RNA interaction
  sigma_sym ~ normal(0, 1) T[0, ];  // individual heterogeneity in susceptibility

  tau_dp[1] ~ normal(-1, 1);  // log-affine intercept
  tau_wp[1] ~ normal(-1, 1);  // log-affine intercept
  tau_wr[1] ~ normal(-1, 1);  // log-affine intercept

  tau_dp[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)
  tau_wp[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)
  tau_wr[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)

  // TCID50 interval-censored regression priors
  alpha_tcid50[1] ~ normal(8, 3);       // intercept a: d*≈2 at high VL
  alpha_tcid50[2] ~ normal(-0.5, 1);    // log(b): b≈0.6, inverse culture growth rate
  alpha_tcid50[3] ~ normal(0, 1);       // log(sigma): sigma≈1 day
  alpha_cult ~ normal(0, prior_beta_sd);

  tau0_lfd_raw ~ std_normal();
  tau_lfd ~ std_normal();

  // LIKELIHOOD — single scalar accumulated in transformed parameters
  target += total_ll;
}

generated quantities {
  // -----------------------------------------------------------------------
  // Recovered 4x4 RNA individual-effect correlation matrix
  // (1x1 identity when ind_corr = 0)
  // -----------------------------------------------------------------------
  corr_matrix[ind_corr ? 4 : 1] Omega_rna;
  if (ind_corr) {
    Omega_rna = multiply_lower_tri_self_transpose(L_Omega_rna);
  } else {
    Omega_rna[1, 1] = 1.0;
  }
  // NOTE: log_lik, rna_rep, pfu_rep, lfd_rep, sym_rep are computed
  // in a separate post-hoc generated-quantities pass using
  // kinetics_model_gq.stan.  This avoids writing ~250 K columns
  // per MCMC iteration during sampling.
}
