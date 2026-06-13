// ==========================================================================
// kinetics_model_gq.stan — Post-hoc generated quantities
//
// This model is NOT for sampling.  It is used with CmdStanR's
// $generate_quantities() method to compute:
//   - Per-observation log-likelihoods (for LOO-CV via loo::loo)
//   - Posterior predictive draws (rna_rep, pfu_rep, lfd_rep, sym_rep)
//   - Fitted values (rna_hat, pfu_hat, lfd_hat)
//
// It reuses the same functions, data, parameters, and transformed data
// blocks as kinetics_model.stan.  The generated quantities block
// recomputes trajectories from the saved parameter draws and generates
// the diagnostic outputs.
// ==========================================================================
#include "functions.stan"

data {
  int<lower=0> K;
  array[K] int<lower=0> M;
  array[K] int<lower=0> N;
  array[K] real lod_rna;
  array[K] real lod_pfu;
  array[K] real fp_mean;

  array[sum(N)] int id;
  array[sum(N)] real time;
  vector[sum(N)] rna;
  vector[sum(N)] pfu;
  array[sum(N)] int<lower=0, upper=1> lfd;
  array[sum(N)] int<lower=0, upper=1> sym;
  array[sum(N)] int<lower=0, upper=1> sym_ever;
  array[sum(N)] int<lower=0, upper=1> sym_at_risk;

  array[sum(N)] int<lower=0> source;
  array[sum(N)] int<lower=0> pfu_type;

  int<lower=0> P;
  array[sum(M)] row_vector[P] x;

  array[sum(N)] int<lower=0, upper=1> rna_exist;
  array[sum(N)] int<lower=0, upper=1> pfu_exist;
  array[sum(N)] int<lower=0, upper=1> lfd_exist;
  array[sum(N)] int<lower=0, upper=1> sym_exist;

  int<lower=0, upper=1> test_error;
  int<lower=0, upper=1> ind_effects;
  int<lower=0, upper=1> ind_corr;

  int<lower=0, upper=1> adj_pfu;
  int<lower=0, upper=1> adj_rna;
  int<lower=0, upper=1> adj_lfd;
  int<lower=0, upper=1> adj_sym;

  int<lower=0, upper=1> source_pfu;
  int<lower=0, upper=1> source_rna;
  int<lower=0, upper=1> source_lfd;
  int<lower=0, upper=1> source_sym;

  int<lower=0, upper=1> use_wf;
  int<lower=0, upper=1> use_smooth;

  // PFU individual effects: restrict to PFU-informed individuals
  int<lower=0> N_pfu_ind;
  array[sum(M)] int<lower=0> pfu_ind_idx;

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

  real prior_fn;
  real prior_fp;
  real<lower=0> prior_fp_kappa;

  real<lower=0> prior_wf_mean;
  real<lower=0> prior_wf_cv;

  real<lower=0> kappa_postpeak;

}

transformed data {
  array[sum(N)] int<lower=0, upper=1> first_obs_per_id;
  first_obs_per_id = first_instance(id);
  real scale_vl = prior_dp_mean;
}

parameters {
  real<lower=0> sigma_rna;
  real<lower=0> sigma_pfu;

  real zeta_sym_intercept;
  real<lower=0> zeta_sym_pfu;
  real<lower=0> zeta_sym_rna;
  real zeta_sym_postpeak;
  real zeta_sym_postpeak_rna;
  real<lower=0> sigma_sym;
  array[sum(M)] real z_sym;

  real<lower=0,upper=1> fp;
  real<lower=0,upper=1> fn;

  real dp_raw;
  real wp_raw;
  real wr_raw;
  array[use_wf ? 1 : 0] real wf_raw;

  array[(1 - ind_corr) * sum(M)] real tp_i_rna;
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real dp_i_rna;
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real wp_i_rna;
  array[ind_effects * (1 - ind_corr) ? sum(M) : 0] real wr_i_rna;
  array[use_wf && ind_effects ? sum(M) : 0] real wf_i;

  cholesky_factor_corr[ind_corr ? 4 : 1] L_Omega_rna;
  vector<lower=0>[ind_corr ? 4 : 0] sigma_ind_rna;
  matrix[ind_corr ? 4 : 0, ind_corr ? sum(M) : 0] z_ind_rna;

  array[K && source_rna ? K : 0] real tp_k_rna;
  array[K && source_rna ? K : 0] real dp_k_rna;
  array[K && source_rna ? K : 0] real wp_k_rna;
  array[K && source_rna ? K : 0] real wr_k_rna;
  array[use_wf && source_rna ? K : 0] real wf_k;

  vector[P && adj_rna ? P : 0] beta_dp_rna;
  vector[P && adj_rna ? P : 0] beta_wp_rna;
  vector[P && adj_rna ? P : 0] beta_wr_rna;

  vector[2] tau_tp;
  vector[2] tau_dp;
  vector[2] tau_wp;
  vector[2] tau_wr;

  vector<lower=0>[ind_effects ? 4 : 1] sigma_ind_pfu;  // PFU RE SDs (learned from informed individuals)
  array[N_pfu_ind] real z_tp_pfu;  // NCP z-scores
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_dp_pfu;
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_wp_pfu;
  array[N_pfu_ind && ind_effects ? N_pfu_ind : 0] real z_wr_pfu;

  array[K && source_pfu ? K : 0] real tp_k_pfu;
  array[K && source_pfu ? K : 0] real dp_k_pfu;
  array[K && source_pfu ? K : 0] real wp_k_pfu;
  array[K && source_pfu ? K : 0] real wr_k_pfu;

  vector[P && adj_pfu ? P : 0] beta_dp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wr_pfu;

  vector[3] alpha_tcid50;  // [1]=a, [2]=log(b), [3]=log(sigma)
  vector[2] alpha_cult;

  real tau0_lfd_raw;
  vector[4] tau_lfd;
  vector[P && adj_lfd ? P : 0] beta_lfd;
  array[K && source_lfd ? K : 0] real lfd_k;

  vector[P && adj_sym ? P : 0] beta_sym;
  array[K && source_sym ? K : 0] real to_k_sym;
}

// No transformed parameters needed — the sampling model's total_ll
// is irrelevant here.  We only need the generated quantities.

generated quantities {
  // Correlation matrix
  corr_matrix[ind_corr ? 4 : 1] Omega_rna;
  if (ind_corr) {
    Omega_rna = multiply_lower_tri_self_transpose(L_Omega_rna);
  } else {
    Omega_rna[1, 1] = 1.0;
  }

  // Per-observation log-likelihoods for LOO-CV
  array[sum(N)] real log_lik;

  // Fitted values (for posterior predictive plots)
  array[sum(N)] real rna_hat;
  array[sum(N)] real pfu_hat;
  array[sum(N)] real lfd_hat;

  // Posterior predictive draws
  array[sum(N)] real rna_rep;
  array[sum(N)] real pfu_rep;
  array[sum(N)] int  lfd_rep;
  array[sum(N)] int  sym_rep;

  {
    // Population means
    real dp_mean_rna = prior_dp_mean * exp(prior_dp_cv * dp_raw);
    real wp_mean_rna = prior_wp_mean * exp(prior_wp_cv * wp_raw);
    real wr_mean_rna = prior_wr_mean * exp(prior_wr_cv * wr_raw);
    real tau0_lfd_val = logit(prior_lfd_mean) + 1 * tau0_lfd_raw;

    // Local variables
    real tp_rna;
    real dp_rna;
    real wp_rna;
    real wr_rna;
    real tp_pfu;
    real dp_pfu;
    real wp_pfu;
    real wr_pfu;
    real wf;
    real theta;

    // Effective RNA individual effects
    array[sum(M)] real eff_tp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_dp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_wp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_wr_i_rna = rep_array(0.0, sum(M));

    if (ind_corr) {
      matrix[4, sum(M)] eta_rna = diag_pre_multiply(sigma_ind_rna, L_Omega_rna) * z_ind_rna;
      for (i in 1:sum(M)) {
        eff_tp_i_rna[i] = eta_rna[1, i];
        eff_dp_i_rna[i] = eta_rna[2, i];
        eff_wp_i_rna[i] = eta_rna[3, i];
        eff_wr_i_rna[i] = eta_rna[4, i];
      }
    } else {
      for (i in 1:sum(M)) eff_tp_i_rna[i] = tp_i_rna[i];
      if (ind_effects) {
        for (i in 1:sum(M)) {
          eff_dp_i_rna[i] = dp_i_rna[i];
          eff_wp_i_rna[i] = wp_i_rna[i];
          eff_wr_i_rna[i] = wr_i_rna[i];
        }
      }
    }

    // NCP reconstruction for PFU individual effects
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

    for (n in 1 : sum(N)) {
      // Initialize log_lik to 0
      log_lik[n] = 0;

      // ---- RNA trajectory ----
      tp_rna = eff_tp_i_rna[id[n]];
      dp_rna = dp_mean_rna;
      wp_rna = wp_mean_rna;
      wr_rna = wr_mean_rna;

      if (ind_effects) {
        dp_rna = dp_rna * exp(eff_dp_i_rna[id[n]]);
        wp_rna = wp_rna * exp(eff_wp_i_rna[id[n]]);
        wr_rna = wr_rna * exp(eff_wr_i_rna[id[n]]);
      }
      if (source_rna) {
        tp_rna = tp_rna + tp_k_rna[source[n]];
        dp_rna = dp_rna * exp(dp_k_rna[source[n]]);
        wp_rna = wp_rna * exp(wp_k_rna[source[n]]);
        wr_rna = wr_rna * exp(wr_k_rna[source[n]]);
      }
      if (adj_rna) {
        dp_rna = dp_rna * exp(x[id[n], ] * beta_dp_rna);
        wp_rna = wp_rna * exp(x[id[n], ] * beta_wp_rna);
        wr_rna = wr_rna * exp(x[id[n], ] * beta_wr_rna);
      }
      dp_rna = safe_pos(dp_rna);
      wp_rna = safe_pos(wp_rna);
      wr_rna = safe_pos(wr_rna);

      // Flat-top
      if (use_wf) {
        wf = prior_wf_mean * exp(prior_wf_cv * wf_raw[1]);
        if (ind_effects) wf = wf * exp(wf_i[id[n]]);
        if (source_rna)  wf = wf * exp(wf_k[source[n]]);
      } else {
        wf = 0.0;
      }

      // RNA fitted value
      if (use_smooth) {
        rna_hat[n] = safe_vl(smooth(time[n], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      } else {
        rna_hat[n] = safe_vl(piecewise(time[n], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      }
      real post_peak_n = inv_logit(kappa_postpeak * (time[n] - tp_rna));

      // ---- PFU trajectory ----
      {
        int pidx = pfu_ind_idx[id[n]];
        tp_pfu = tau_tp[1] + tau_tp[2] * tp_rna;
        if (pidx > 0) tp_pfu += tp_i_pfu[pidx];
      }
      dp_pfu = log_affine(tau_dp[1], tau_dp[2], dp_rna);
      wp_pfu = log_affine(tau_wp[1], tau_wp[2], wp_rna);
      wr_pfu = log_affine(tau_wr[1], tau_wr[2], wr_rna);

      if (ind_effects && pfu_ind_idx[id[n]] > 0) {
        int pidx = pfu_ind_idx[id[n]];
        dp_pfu = dp_pfu * exp(dp_i_pfu[pidx]);
        wp_pfu = wp_pfu * exp(wp_i_pfu[pidx]);
        wr_pfu = wr_pfu * exp(wr_i_pfu[pidx]);
      }
      if (source_pfu) {
        tp_pfu = tp_pfu + tp_k_pfu[source[n]];
        dp_pfu = dp_pfu * exp(dp_k_pfu[source[n]]);
        wp_pfu = wp_pfu * exp(wp_k_pfu[source[n]]);
        wr_pfu = wr_pfu * exp(wr_k_pfu[source[n]]);
      }
      if (adj_pfu) {
        dp_pfu = dp_pfu * exp(x[id[n], ] * beta_dp_pfu);
        wp_pfu = wp_pfu * exp(x[id[n], ] * beta_wp_pfu);
        wr_pfu = wr_pfu * exp(x[id[n], ] * beta_wr_pfu);
      }
      dp_pfu = safe_pos(dp_pfu);
      wp_pfu = safe_pos(wp_pfu);
      wr_pfu = safe_pos(wr_pfu);

      // PFU fitted value
      if (use_smooth) {
        pfu_hat[n] = safe_vl(smooth(time[n], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      } else {
        pfu_hat[n] = safe_vl(piecewise(time[n], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      }
      // Biological constraint: infectious virus ≤ total viral RNA
      pfu_hat[n] = fmin(pfu_hat[n], rna_hat[n]);

      // ---- RNA log-likelihood + posterior predictive ----
      if (rna_exist[n] == 1) {
        real ll_rna;
        if (rna[n] <= lod_rna[source[n]]) {
          ll_rna = normal_lcdf(lod_rna[source[n]] | rna_hat[n], sigma_rna);
          if (test_error) ll_rna += log1m(fn);
        } else {
          ll_rna = normal_lpdf(rna[n] | rna_hat[n], sigma_rna);
          if (test_error) {
            ll_rna = log_sum_exp(
              log(fp) + exponential_lpdf(rna[n] - lod_rna[source[n]] | fp_mean[source[n]]),
              log1m(fp) + log1m(fn) + ll_rna
            );
          }
        }
        log_lik[n] += ll_rna;
        rna_rep[n] = normal_rng(rna_hat[n], sigma_rna);
        if (rna_rep[n] < lod_rna[source[n]]) rna_rep[n] = lod_rna[source[n]];
      } else {
        rna_rep[n] = lod_rna[source[n]];
      }

      // ---- PFU log-likelihood + posterior predictive ----
      if (pfu_exist[n] == 1) {
        real ll_pfu;
        if (pfu_type[n] == 1) {
          if (pfu[n] <= lod_pfu[source[n]]) {
            ll_pfu = normal_lcdf(lod_pfu[source[n]] | pfu_hat[n], sigma_pfu);
            if (test_error) ll_pfu += log1m(fn);
          } else {
            ll_pfu = normal_lpdf(pfu[n] | pfu_hat[n], sigma_pfu);
            if (test_error) {
              ll_pfu = log_sum_exp(
                log(fp) + exponential_lpdf(pfu[n] - lod_pfu[source[n]] | fp_mean[source[n]]),
                log1m(fp) + log1m(fn) + ll_pfu
              );
            }
          }
          pfu_rep[n] = normal_rng(pfu_hat[n], sigma_pfu);
          if (pfu_rep[n] < lod_pfu[source[n]]) pfu_rep[n] = lod_pfu[source[n]];
        } else if (pfu_type[n] == 2) {
          // Interval-censored normal for TCID50 days-to-positivity
          real mu_tcid = alpha_tcid50[1] - exp(alpha_tcid50[2]) * pfu_hat[n];
          real sig_tcid = exp(alpha_tcid50[3]);
          real d_obs = pfu[n];
          if (d_obs >= 5.5) {
            ll_pfu = normal_lccdf(5.0 | mu_tcid, sig_tcid);
          } else if (d_obs <= 2.5) {
            ll_pfu = normal_lcdf(2.0 | mu_tcid, sig_tcid);
          } else {
            real d_val = round(d_obs);
            ll_pfu = log_diff_exp(
              normal_lcdf(d_val | mu_tcid, sig_tcid),
              normal_lcdf(d_val - 1.0 | mu_tcid, sig_tcid)
            );
          }
          // Posterior predictive: draw d*, discretize to {2,3,4,5,6}
          {
            real d_star = normal_rng(mu_tcid, sig_tcid);
            if (d_star <= 2.0) pfu_rep[n] = 2.0;
            else if (d_star <= 3.0) pfu_rep[n] = 3.0;
            else if (d_star <= 4.0) pfu_rep[n] = 4.0;
            else if (d_star <= 5.0) pfu_rep[n] = 5.0;
            else pfu_rep[n] = 6.0;
          }
        } else if (pfu_type[n] == 3) {
          theta = inv_logit(alpha_cult[1] + alpha_cult[2] * pfu_hat[n]);
          ll_pfu = bernoulli_lpmf(to_int(pfu[n]) | theta);
          pfu_rep[n] = bernoulli_rng(theta);
        }
        log_lik[n] += ll_pfu;
      } else {
        pfu_rep[n] = lod_pfu[source[n]];
      }

      // ---- LFD log-likelihood + posterior predictive ----
      lfd_hat[n] = inv_logit(tau0_lfd_val + tau_lfd[1] * rna_hat[n] + tau_lfd[2] * pfu_hat[n] + tau_lfd[3] * post_peak_n + tau_lfd[4] * post_peak_n * rna_hat[n]);
      if (source_lfd) lfd_hat[n] = inv_logit(logit(lfd_hat[n]) + lfd_k[source[n]]);
      if (adj_lfd)    lfd_hat[n] = inv_logit(logit(lfd_hat[n]) + x[id[n], ] * beta_lfd);

      if (lfd_exist[n] == 1) {
        log_lik[n] += bernoulli_lpmf(lfd[n] | lfd_hat[n]);
        lfd_rep[n] = bernoulli_rng(lfd_hat[n]);
      } else {
        lfd_rep[n] = 0;
      }

      // ---- Symptom log-likelihood + posterior predictive ----
      if (sym_exist[n] == 1 && sym_at_risk[n] == 1) {
        real u_sym_gq = sigma_sym * z_sym[id[n]];
        real zeta_lin_gq = zeta_sym_intercept
                         + zeta_sym_pfu * (pfu_hat[n] / scale_vl)
                         + zeta_sym_rna * (rna_hat[n] / scale_vl)
                         + zeta_sym_postpeak * post_peak_n
                         + zeta_sym_postpeak_rna * post_peak_n * (rna_hat[n] / scale_vl)
                         + u_sym_gq;
        if (source_sym) zeta_lin_gq += to_k_sym[source[n]];
        if (adj_sym)    zeta_lin_gq += x[id[n], ] * beta_sym;
        zeta_lin_gq = fmin(zeta_lin_gq, 10.0);

        if (sym[n] == 1) {
          log_lik[n] += log1m_exp(-exp(zeta_lin_gq));
        } else {
          log_lik[n] += -exp(zeta_lin_gq);
        }
        real p_sym = 1 - exp(-exp(zeta_lin_gq));
        sym_rep[n] = bernoulli_rng(p_sym);
      } else {
        sym_rep[n] = 0;
      }
    }
  }
}
