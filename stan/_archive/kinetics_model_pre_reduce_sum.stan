functions {
  // This function counts the number of times a specific element appears in a given array.
  // vec: The array in which to search for the element.
  // elem: The element to count in the array.
  int count_elem(array[] int vec, int elem) {
    int count;
    count = 0;
    for (i in 1 : num_elements(vec)) {
      if (vec[i] == elem) {
        count = count + 1;
      }
    }
    return count;
  }
  
  array[] int first_instance(array[] int id) {
    array[num_elements(id)] int is_first_instance;
    int prev;
    prev = 0;
    for (i in 1 : num_elements(id)) {
      if (id[i] == prev) {
        is_first_instance[i] = 0;
      } else {
        is_first_instance[i] = 1;
      }
      prev = id[i];
    }
    return is_first_instance;
  }
  
  array[] int add(array[] int x, array[] int y) {
        int x_size = size(x);
        array[x_size] int z;
        for (i in 1:x_size){
          z[i] = x[i] + y[i];
        }
        return z;
    }
    
  // This function calculates a piece-wise exponential function
  // with an optional flat-top segment of width wf at the peak.
  // t: The time variable.
  // tp: The peak time (start of flat-top).
  // wp: The proliferation width (rising arm).
  // wr: The clearance width (falling arm).
  // dp: The peak height (log-scale).
  // wf: The flat-top duration at peak (0 = sharp peak).
  real piecewise(real t, real tp, real wp, real wr, real dp, real wf) {
    if (t <= tp) {
      return dp / wp * (t - (tp - wp)); // rising arm
    } else if (t <= tp + wf) {
      return dp;                         // flat top
    } else {
      return dp - dp / wr * (t - tp - wf); // falling arm
    }
  }

  // Smooth log-sum-exp approximation to the piecewise trajectory.
  // Same signature as piecewise(). Produces continuously differentiable
  // trajectories that round the peak transition. The falling arm is
  // shifted right by wf (natural scale, days), creating an approximate
  // plateau of width ~wf.  When wf = 0 reduces to the standard
  // smooth two-arm envelope.
  //
  // A scaled soft-min cap at dp prevents the log-sum-exp envelope from
  // overshooting the peak when wf > 0 (an inherent artifact of merging
  // two separated exponential arms).  For wf = 0 the raw value equals
  // dp exactly at the peak, so the cap introduces negligible rounding
  // (~0.07 log-units at k_cap = 10).  Gradients remain smooth.
  real smooth(real t, real tp, real wp, real wr, real dp, real wf) {
    real a = dp / wp; // proliferation rate
    real b = dp / wr; // clearance rate
    // Guard against overflow: cap exponents at ±50
    real arg1 = fmin(-a * (t - tp), 50.0);
    real arg2 = fmin( b * (t - (tp + wf)), 50.0);
    real raw = dp + log(
      (a + b) /
      (b * exp(arg1) + a * exp(arg2))
    );
    // Soft-cap at dp: soft_min(raw, dp) with sharpness k_cap = 10.
    // Prevents overshoot when wf > 0; transparent on the arms (raw << dp).
    // Stan's log1p_exp() is numerically stable for all input magnitudes.
    return raw - log1p_exp(10.0 * (raw - dp)) * 0.1;
  }

  // Numerically stable log-affine transformation:
  //   exp(a) * pow(x, b) = exp(a + b*log(x))
  // with overflow protection (cap result at exp(30) ~ 1e13).
  real log_affine(real intercept, real elasticity, real x) {
    return exp(fmin(intercept + elasticity * log(x), 30.0));
  }

  // Clamp a positive kinetic parameter to a safe range
  // (prevents overflow/underflow in piecewise divisions).
  real safe_pos(real x) {
    return fmax(fmin(x, 1e4), 1e-4);
  }

  // Clamp predicted viral load to a physically plausible range on the
  // log scale.  exp(50) ~ 5e21 copies/ml is impossible; exp(-50) ~ 0.
  real safe_vl(real x) {
    return fmax(fmin(x, 50.0), -50.0);
  }
}

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
  real<lower=0> prior_k_sd;
  real<lower=0> prior_beta_sd;

  real<lower=0> prior_sigma_sd;
  real<lower=0, upper=1> prior_lfd_mean;  
  
  real prior_fn; // prior mean on false negative probability
  real prior_fp; // prior mean on false positive probability
  real<lower=0> prior_fp_kappa; // precision (concentration) for fp/fn Beta priors

  real<lower=0> prior_wf_mean;  // prior mean for flat-top duration (days)
  real<lower=0> prior_wf_cv;    // prior CV for flat-top duration (NCP scale)
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
  real eta_sym_intercept; // baseline log-hazard
  real<lower=0> eta_sym_pfu; // log V_t coefficient
  real<lower=0> eta_sym_rna; // log R_t coefficient
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

  // individual effects: PFU
  array[sum(M)] real tp_i_pfu; // onset
  array[sum(M) && ind_effects ? sum(M) : 0] real dp_i_pfu; // peak
  array[sum(M) && ind_effects ? sum(M) : 0] real wp_i_pfu; // proliferation time
  array[sum(M) && ind_effects ? sum(M) : 0] real wr_i_pfu; // clearance time
  
  // source effects: PFU
  array[K && source_pfu ? K : 0] real tp_k_pfu; // onset
  array[K && source_pfu ? K : 0] real dp_k_pfu; // peak
  array[K && source_pfu ? K : 0] real wp_k_pfu; // proliferation time
  array[K && source_pfu ? K : 0] real wr_k_pfu; // clearance time
  
  // coefficient vectors for predictors of PFU peak, proliferation, and clearance 
  vector[P && adj_pfu ? P : 0] beta_dp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wp_pfu;
  vector[P && adj_pfu ? P : 0] beta_wr_pfu;
  
  // coefficient for transforming time to cell culture to PFUs
  vector[4] alpha_tcid50;
  
  // coefficient for binary culture result to PFUs 
  vector[2] alpha_cult;
  
  // non-centered parameterization for logistic intercept for LFD
  real tau0_lfd_raw;
  
  // coefficient vector for transformation of [RNA] and PFUs to LFD positivity
  vector[2] tau_lfd;
  
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

  // -----------------------------------------------------------------------
  // Total log-likelihood accumulator.
  // Replaces per-observation log_lik arrays — avoids writing ~100 K
  // columns per MCMC iteration.  Per-observation log_lik and posterior
  // predictive draws are computed in a post-hoc generated-quantities pass.
  // -----------------------------------------------------------------------
  real total_ll = 0;

  {
    // combined effects: RNA (non-centered)
    real tp_rna; // onset
    real dp_rna; // peak
    real wp_rna; // proliferation time
    real wr_rna; // clearance time

    // combined effects: PFU (non-centered)
    real tp_pfu; // onset
    real dp_pfu; // peak
    real wp_pfu; // proliferation time
    real wr_pfu; // clearance time

    // flat-top duration (shared between RNA and PFU)
    real wf; // 0 when use_wf=0

    // combined effects: symptoms (cloglog hazard)
    real u_sym; // individual random effect
    real eta_lin; // linear predictor for symptom hazard

    // alternative viral culture transformations
    real lambda; // time to TCID50 if positive
    real theta;  // probability positive / negative

    // per-observation fitted values (local — NOT saved to CSV)
    real pfu_hat_n;
    real rna_hat_n;
    real lfd_hat_n;

    // Effective RNA individual effects (from correlated NCP or independent params)
    array[sum(M)] real eff_tp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_dp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_wp_i_rna = rep_array(0.0, sum(M));
    array[sum(M)] real eff_wr_i_rna = rep_array(0.0, sum(M));

    if (ind_corr) {
      // NCP: eta = diag(sigma_ind_rna) * L_Omega_rna * z_ind_rna
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

    // main observation loop
    for (n in 1 : sum(N)) {

      // always include a random effect for peak time
      tp_rna = eff_tp_i_rna[id[n]];

      // RNA: start from population means (hoisted out of loop)
      dp_rna = dp_mean_rna;
      wp_rna = wp_mean_rna;
      wr_rna = wr_mean_rna;

      // RNA: update with individual random effects
      if (ind_effects) {
        dp_rna = dp_rna * exp(eff_dp_i_rna[id[n]]);
        wp_rna = wp_rna * exp(eff_wp_i_rna[id[n]]);
        wr_rna = wr_rna * exp(eff_wr_i_rna[id[n]]);
      }

      // RNA: update with source random effects
      if (source_rna) {
        tp_rna = tp_rna + tp_k_rna[source[n]];
        dp_rna = dp_rna * exp(dp_k_rna[source[n]]);
        wp_rna = wp_rna * exp(wp_k_rna[source[n]]);
        wr_rna = wr_rna * exp(wr_k_rna[source[n]]);
      }

      // update with covariate effects if requested
      if (adj_rna) {
        dp_rna = dp_rna * exp(x[id[n], ] * beta_dp_rna);
        wp_rna = wp_rna * exp(x[id[n], ] * beta_wp_rna);
        wr_rna = wr_rna * exp(x[id[n], ] * beta_wr_rna);
      }

      // Clamp RNA kinetics to safe range before PFU transform & piecewise
      dp_rna = safe_pos(dp_rna);
      wp_rna = safe_pos(wp_rna);
      wr_rna = safe_pos(wr_rna);

      // Flat-top duration (shared by RNA and PFU trajectories)
      if (use_wf) {
        wf = prior_wf_mean * exp(prior_wf_cv * wf_raw[1]);
        if (ind_effects) {
          wf = wf * exp(wf_i[id[n]]);
        }
        if (source_rna) {
          wf = wf * exp(wf_k[source[n]]);
        }
      } else {
        wf = 0.0;
      }

      // PFU: convert RNA to PFU parameters (log-affine transformation, overflow-safe)
      tp_pfu = tau_tp[1] + tau_tp[2] * tp_rna + tp_i_pfu[id[n]];
      dp_pfu = log_affine(tau_dp[1], tau_dp[2], dp_rna);
      wp_pfu = log_affine(tau_wp[1], tau_wp[2], wp_rna);
      wr_pfu = log_affine(tau_wr[1], tau_wr[2], wr_rna);

      // PFU: update with individual random effects
      if (ind_effects) {
        dp_pfu = dp_pfu * exp(dp_i_pfu[id[n]]);
        wp_pfu = wp_pfu * exp(wp_i_pfu[id[n]]);
        wr_pfu = wr_pfu * exp(wr_i_pfu[id[n]]);
      }

      // PFU: update with source random effects
      if (source_pfu) {
        tp_pfu = tp_pfu + tp_k_pfu[source[n]];
        dp_pfu = dp_pfu * exp(dp_k_pfu[source[n]]);
        wp_pfu = wp_pfu * exp(wp_k_pfu[source[n]]);
        wr_pfu = wr_pfu * exp(wr_k_pfu[source[n]]);
      }

      // PFU: update with covariate effects
      if (adj_pfu) {
        dp_pfu = dp_pfu * exp(x[id[n], ] * beta_dp_pfu);
        wp_pfu = wp_pfu * exp(x[id[n], ] * beta_wp_pfu);
        wr_pfu = wr_pfu * exp(x[id[n], ] * beta_wr_pfu);
      }

      // Clamp PFU kinetics to safe range
      dp_pfu = safe_pos(dp_pfu);
      wp_pfu = safe_pos(wp_pfu);
      wr_pfu = safe_pos(wr_pfu);

      // ---- Predicted PFU trajectory (local scalar) ----
      if (use_smooth) {
        pfu_hat_n = safe_vl(smooth(time[n], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      } else {
        pfu_hat_n = safe_vl(piecewise(time[n], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      }

      // if PFUs of culturable virus are not missing
      if (pfu_exist[n] == 1) {
        real ll_pfu;

        // if traditional PFU measurement
        if (pfu_type[n] == 1) {
          if (pfu[n] <= lod_pfu[source[n]]) { // if below limit of detection

            // log-likelihood for PFU: censored below limit
            ll_pfu = normal_lcdf(lod_pfu[source[n]] | pfu_hat_n, sigma_pfu);

            // test error mixture
            if (test_error) {
              ll_pfu += log1m(fn);
            }
          } else { // if above limit of detection

            // log-likelihood for PFU
            ll_pfu = normal_lpdf(pfu[n] | pfu_hat_n, sigma_pfu);

            // test error mixture
            if (test_error) {
              ll_pfu = log_sum_exp(
                log(fp) + exponential_lpdf(pfu[n] - lod_pfu[source[n]] | fp_mean[source[n]]),
                log1m(fp) + log1m(fn) + ll_pfu
              );
            }
          }

        // if days to 50% cell culture
        } else if (pfu_type[n] == 2) {
          theta = inv_logit(alpha_tcid50[1] + alpha_tcid50[2] * pfu_hat_n);

          // simple logistic model
          if (pfu[n] < 5.9) {
            ll_pfu = bernoulli_lpmf(1 | theta);
          } else {
            ll_pfu = bernoulli_lpmf(0 | theta);
          }

        // if simple positive/negative culture result
        } else if (pfu_type[n] == 3) {
          theta = inv_logit(alpha_cult[1] + alpha_cult[2] * pfu_hat_n);
          ll_pfu = bernoulli_lpmf(to_int(pfu[n]) | theta);
        }

        total_ll += ll_pfu;
      }

      // ---- Predicted RNA trajectory (local scalar) ----
      if (use_smooth) {
        rna_hat_n = safe_vl(smooth(time[n], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      } else {
        rna_hat_n = safe_vl(piecewise(time[n], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      }

      // if [RNA] is not missing
      if (rna_exist[n] == 1) {
        real ll_rna;
        if (rna[n] <= lod_rna[source[n]]) { // if below limit of detection

          // log-likelihood for [RNA]: censored below limit
          ll_rna = normal_lcdf(lod_rna[source[n]] | rna_hat_n, sigma_rna);

          // test error mixture
          if (test_error) {
            ll_rna += log1m(fn);
          }
        } else { // if above limit

          // log-likelihood for [RNA]
          ll_rna = normal_lpdf(rna[n] | rna_hat_n, sigma_rna);

          // test error mixture
          if (test_error) {
            ll_rna = log_sum_exp(
              log(fp) + exponential_lpdf(rna[n] - lod_rna[source[n]] | fp_mean[source[n]]),
              log1m(fp) + log1m(fn) + ll_rna
            );
          }
        }
        total_ll += ll_rna;
      }

      // ---- Predicted LFD probability (local scalar) ----
      lfd_hat_n = inv_logit(tau0_lfd + tau_lfd[1] * rna_hat_n + tau_lfd[2] * pfu_hat_n);

      // update with source effects if requested
      if (source_lfd) {
        lfd_hat_n = inv_logit(logit(lfd_hat_n) + lfd_k[source[n]]);
      }

      // update with covariate effects if requested
      if (adj_lfd) {
        lfd_hat_n = inv_logit(logit(lfd_hat_n) + x[id[n], ] * beta_lfd);
      }

      // if LFD is not missing
      if (lfd_exist[n] == 1) {
        total_ll += bernoulli_lpmf(lfd[n] | lfd_hat_n);
      }

      // ---- Symptom onset: discrete-time hazard with cloglog link ----
      // Individual random effect (non-centered)
      u_sym = sigma_sym * z_sym[id[n]];

      // Linear predictor on log-hazard scale (normalised by scale_vl)
      eta_lin = eta_sym_intercept
               + eta_sym_pfu * (pfu_hat_n / scale_vl)
               + eta_sym_rna * (rna_hat_n / scale_vl)
               + u_sym;

      if (source_sym) {
        eta_lin = eta_lin + to_k_sym[source[n]];
      }

      // update with covariate effects if requested
      if (adj_sym) {
        eta_lin = eta_lin + x[id[n], ] * beta_sym;
      }

      // Cap eta_lin to prevent exp(exp(eta)) overflow in cloglog link.
      eta_lin = fmin(eta_lin, 10.0);

      // Symptom likelihood (cloglog link: P(Y=1) = 1 - exp(-exp(eta)))
      if (sym_exist[n] == 1 && sym_at_risk[n] == 1) {
        if (sym[n] == 1) {
          total_ll += log1m_exp(-exp(eta_lin));  // log(1 - exp(-exp(eta)))
        } else {
          total_ll += -exp(eta_lin);  // log(exp(-exp(eta)))
        }
      }
    }
  }
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

  tp_i_pfu ~ normal(0, prior_i_sd);
  z_sym ~ std_normal();  // NCP for symptom random effects

  // RNA individual effects: correlated (Cholesky NCP) or independent
  if (ind_corr) {
    L_Omega_rna ~ lkj_corr_cholesky(2);  // LKJ(2): mild shrinkage toward identity
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
    dp_i_pfu ~ normal(0, prior_i_sd);
    wp_i_pfu ~ normal(0, prior_i_sd);
    wr_i_pfu ~ normal(0, prior_i_sd);
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
  eta_sym_intercept ~ normal(-3, 1);  // low baseline daily hazard ~exp(-3)≈0.05
  eta_sym_pfu ~ normal(0, 0.5) T[0, ];  // more virus -> more symptoms
  eta_sym_rna ~ normal(0, 0.5) T[0, ];  // more RNA -> more symptoms
  sigma_sym ~ normal(0, 1) T[0, ];  // individual heterogeneity in susceptibility

  tau_dp[1] ~ normal(-1, 1);  // log-affine intercept
  tau_wp[1] ~ normal(-1, 1);  // log-affine intercept
  tau_wr[1] ~ normal(-1, 1);  // log-affine intercept

  tau_dp[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)
  tau_wp[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)
  tau_wr[2] ~ normal(1, 0.5) T[0, ];  // log-affine elasticity (positive)

  alpha_tcid50 ~ normal(0, prior_beta_sd);
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
