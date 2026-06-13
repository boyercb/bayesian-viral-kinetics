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

  // Time derivative of the piecewise trajectory.
  // Returns dp/wp on the rising arm, 0 on flat-top, -dp/wr on the falling arm.
  real piecewise_deriv(real t, real tp, real wp, real wr, real dp, real wf) {
    if (t <= tp) {
      return dp / wp;
    } else if (t <= tp + wf) {
      return 0.0;
    } else {
      return -dp / wr;
    }
  }

  // Time derivative of the smooth trajectory, including the chain
  // rule through the soft-cap at dp.  Uses the same intermediates
  // as smooth() for numerical consistency.
  //
  // Raw derivative (before soft-cap):
  //   d(raw)/dt = a*b*(exp(arg1) - exp(arg2)) / (b*exp(arg1) + a*exp(arg2))
  //
  // Chain rule through soft-cap  f(raw) = raw - log1p_exp(k*(raw-dp))/k:
  //   f'(raw) = 1 - sigmoid(k*(raw - dp))   [≈ 1 on arms, → 0 at peak]
  //
  // Overall: d(smooth)/dt = f'(raw) * d(raw)/dt
  real smooth_deriv(real t, real tp, real wp, real wr, real dp, real wf) {
    real a = dp / wp;
    real b = dp / wr;
    real arg1 = fmin(-a * (t - tp), 50.0);
    real arg2 = fmin( b * (t - (tp + wf)), 50.0);
    real ea1 = exp(arg1);
    real ea2 = exp(arg2);
    real denom = b * ea1 + a * ea2;
    // Raw derivative of the log-sum-exp envelope
    real draw_dt = a * b * (ea1 - ea2) / denom;
    // Soft-cap correction: chain rule through log1p_exp(k*(raw-dp))/k
    real raw = dp + log((a + b) / denom);
    real cap_factor = 1.0 - inv_logit(10.0 * (raw - dp));
    return cap_factor * draw_dt;
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

  // ─────────────────────────────────────────────────────────────────────────
  // partial_sum_ll — reduce_sum worker for within-chain parallelism.
  //
  // Computes the partial log-likelihood for the slice id_slice (observations
  // [start, end] of the full id array).  All other per-observation arrays are
  // passed as additional arguments and indexed via start + n - 1.
  // ─────────────────────────────────────────────────────────────────────────
  real partial_sum_ll(
    array[] int    id_slice,        // slice of id[start:end]
    int start, int end,
    // ── per-observation data (full arrays, indexed by start+n-1) ───────────
    array[] real   time_obs,
    vector         rna_obs,
    data vector    pfu_obs,
    array[] int    lfd_obs,
    array[] int    sym_obs,
    array[] int    sym_at_risk_obs,
    array[] int    sym_exist_obs,
    array[] int    rna_exist_obs,
    array[] int    pfu_exist_obs,
    array[] int    lfd_exist_obs,
    array[] int    pfu_type_obs,
    array[] int    source_obs,
    // ── individual-level (indexed by id) ────────────────────────────────────
    matrix         eff_rna,         // [4 x N_ind]: rows = tp, dp, wp, wr
    array[] real   tp_i_pfu_arg,    // N_pfu_ind (PFU-informed individuals only)
    array[] real   dp_i_pfu_arg,    // N_pfu_ind or 0 (ind_effects)
    array[] real   wp_i_pfu_arg,
    array[] real   wr_i_pfu_arg,
    array[] int    pfu_ind_idx_arg, // [N_ind]: 0=no PFU RE, >0=index into PFU RE arrays
    array[] real   z_sym_arg,
    array[] real   wf_i_arg,        // N_ind or 0 (use_wf && ind_effects)
    array[] row_vector x_arg,
    // ── source-level ────────────────────────────────────────────────────────
    array[] real   lod_rna_arg,   array[] real lod_pfu_arg,
    array[] real   fp_mean_arg,
    array[] real   tp_k_rna_arg,  array[] real dp_k_rna_arg,
    array[] real   wp_k_rna_arg,  array[] real wr_k_rna_arg,
    array[] real   wf_k_arg,
    array[] real   tp_k_pfu_arg,  array[] real dp_k_pfu_arg,
    array[] real   wp_k_pfu_arg,  array[] real wr_k_pfu_arg,
    array[] real   lfd_k_arg,
    array[] real   to_k_sym_arg,
    // ── covariate betas ─────────────────────────────────────────────────────
    vector beta_dp_rna,  vector beta_wp_rna,  vector beta_wr_rna,
    vector beta_dp_pfu,  vector beta_wp_pfu,  vector beta_wr_pfu,
    vector beta_lfd_arg, vector beta_sym_arg,
    // ── population / transformation parameters ───────────────────────────────
    real dp_mean_rna,    real wp_mean_rna,    real wr_mean_rna,
    vector tau_tp_arg,   vector tau_dp_arg,
    vector tau_wp_arg,   vector tau_wr_arg,
    real tau0_lfd_arg,   vector tau_lfd_arg,
    real sigma_rna_arg,  real sigma_pfu_arg,
    real fp_arg,         real fn_arg,
    real zeta_sym_intercept_arg,
    real zeta_sym_pfu_arg, real zeta_sym_rna_arg,
    real zeta_sym_postpeak_arg, real zeta_sym_postpeak_rna_arg,
    real sigma_sym_arg,
    vector alpha_tcid50_arg, vector alpha_cult_arg,
    real wf_pop,            // pre-computed pop flat-top (0 when use_wf=0)
    // ── model flags ─────────────────────────────────────────────────────────
    int ind_effects,
    int source_rna,  int source_pfu,  int source_lfd,  int source_sym,
    int adj_rna,     int adj_pfu,     int adj_lfd,     int adj_sym,
    int test_error,  int use_wf,      int use_smooth,
    real scale_vl,
    real kappa_pp     // steepness for smooth post-peak sigmoid
  ) {
    real ll = 0.0;
    for (n in 1:size(id_slice)) {
      int nn = start + n - 1;   // global observation index
      int i  = id_slice[n];     // individual index
      int k  = source_obs[nn];  // source index

      // ── RNA trajectory ──────────────────────────────────────────────────
      real tp_rna = eff_rna[1, i];
      real dp_rna = dp_mean_rna;
      real wp_rna = wp_mean_rna;
      real wr_rna = wr_mean_rna;

      if (ind_effects) {
        dp_rna = dp_rna * exp(eff_rna[2, i]);
        wp_rna = wp_rna * exp(eff_rna[3, i]);
        wr_rna = wr_rna * exp(eff_rna[4, i]);
      }
      if (source_rna) {
        tp_rna = tp_rna + tp_k_rna_arg[k];
        dp_rna = dp_rna * exp(dp_k_rna_arg[k]);
        wp_rna = wp_rna * exp(wp_k_rna_arg[k]);
        wr_rna = wr_rna * exp(wr_k_rna_arg[k]);
      }
      if (adj_rna) {
        dp_rna = dp_rna * exp(x_arg[i] * beta_dp_rna);
        wp_rna = wp_rna * exp(x_arg[i] * beta_wp_rna);
        wr_rna = wr_rna * exp(x_arg[i] * beta_wr_rna);
      }
      dp_rna = safe_pos(dp_rna);
      wp_rna = safe_pos(wp_rna);
      wr_rna = safe_pos(wr_rna);

      // ── Flat-top duration ────────────────────────────────────────────────
      real wf = 0.0;
      if (use_wf) {
        wf = wf_pop;
        if (ind_effects) wf = wf * exp(wf_i_arg[i]);
        if (source_rna)  wf = wf * exp(wf_k_arg[k]);
      }

      // ── RNA hat (computed first — needed for PFU ≤ RNA constraint) ───────
      real rna_hat_n;
      if (use_smooth) {
        rna_hat_n = safe_vl(smooth(time_obs[nn], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      } else {
        rna_hat_n = safe_vl(piecewise(time_obs[nn], tp_rna, wp_rna, wr_rna, dp_rna, wf));
      }

      // ── Post-peak indicator (smooth sigmoid for HMC differentiability) ────
      real post_peak_n = inv_logit(kappa_pp * (time_obs[nn] - tp_rna));

      // ── PFU trajectory ───────────────────────────────────────────────────
      // PFU REs only exist for PFU-informed individuals (pidx > 0).
      // Non-informed individuals get population-level RNA→PFU transform only.
      int pidx = pfu_ind_idx_arg[i];  // 0 = no PFU RE
      real tp_pfu = tau_tp_arg[1] + tau_tp_arg[2] * tp_rna;
      if (pidx > 0) tp_pfu += tp_i_pfu_arg[pidx];
      real dp_pfu = log_affine(tau_dp_arg[1], tau_dp_arg[2], dp_rna);
      real wp_pfu = log_affine(tau_wp_arg[1], tau_wp_arg[2], wp_rna);
      real wr_pfu = log_affine(tau_wr_arg[1], tau_wr_arg[2], wr_rna);

      if (ind_effects && pidx > 0) {
        dp_pfu = dp_pfu * exp(dp_i_pfu_arg[pidx]);
        wp_pfu = wp_pfu * exp(wp_i_pfu_arg[pidx]);
        wr_pfu = wr_pfu * exp(wr_i_pfu_arg[pidx]);
      }
      if (source_pfu) {
        tp_pfu = tp_pfu + tp_k_pfu_arg[k];
        dp_pfu = dp_pfu * exp(dp_k_pfu_arg[k]);
        wp_pfu = wp_pfu * exp(wp_k_pfu_arg[k]);
        wr_pfu = wr_pfu * exp(wr_k_pfu_arg[k]);
      }
      if (adj_pfu) {
        dp_pfu = dp_pfu * exp(x_arg[i] * beta_dp_pfu);
        wp_pfu = wp_pfu * exp(x_arg[i] * beta_wp_pfu);
        wr_pfu = wr_pfu * exp(x_arg[i] * beta_wr_pfu);
      }
      dp_pfu = safe_pos(dp_pfu);
      wp_pfu = safe_pos(wp_pfu);
      wr_pfu = safe_pos(wr_pfu);

      // ── PFU hat ──────────────────────────────────────────────────────────
      real pfu_hat_n;
      if (use_smooth) {
        pfu_hat_n = safe_vl(smooth(time_obs[nn], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      } else {
        pfu_hat_n = safe_vl(piecewise(time_obs[nn], tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf));
      }
      // Biological constraint: infectious virus ≤ total viral RNA
      pfu_hat_n = fmin(pfu_hat_n, rna_hat_n);

      // ── PFU likelihood ───────────────────────────────────────────────────
      if (pfu_exist_obs[nn] == 1) {
        real ll_pfu;
        if (pfu_type_obs[nn] == 1) {
          if (pfu_obs[nn] <= lod_pfu_arg[k]) {
            ll_pfu = normal_lcdf(lod_pfu_arg[k] | pfu_hat_n, sigma_pfu_arg);
            if (test_error) ll_pfu += log1m(fn_arg);
          } else {
            ll_pfu = normal_lpdf(pfu_obs[nn] | pfu_hat_n, sigma_pfu_arg);
            if (test_error) {
              ll_pfu = log_sum_exp(
                log(fp_arg) + exponential_lpdf(pfu_obs[nn] - lod_pfu_arg[k] | fp_mean_arg[k]),
                log1m(fp_arg) + log1m(fn_arg) + ll_pfu
              );
            }
          }
        } else if (pfu_type_obs[nn] == 2) {
          // Interval-censored normal regression for TCID50 days-to-positivity.
          // Mechanistic: d* = a - b*log(V_t) + eps, eps ~ N(0, sigma)
          //   where b = exp(alpha[2]) > 0 is inverse culture growth rate,
          //   sigma = exp(alpha[3]) > 0 is scale of stochastic variation.
          // pfu_obs: 2,3,4,5 = days to CPE; 6 = negative (no CPE by day 5).
          real mu_tcid = alpha_tcid50_arg[1] - exp(alpha_tcid50_arg[2]) * pfu_hat_n;
          real sig_tcid = exp(alpha_tcid50_arg[3]);
          real d_obs = pfu_obs[nn];
          if (d_obs >= 5.5) {
            // Negative result: right-censored beyond day 5
            ll_pfu = normal_lccdf(5.0 | mu_tcid, sig_tcid);
          } else if (d_obs <= 2.5) {
            // Day 2: left-censored at day 2
            ll_pfu = normal_lcdf(2.0 | mu_tcid, sig_tcid);
          } else {
            // Days 3,4,5: interval (d-1, d]
            real d_val = round(d_obs);
            ll_pfu = log_diff_exp(
              normal_lcdf(d_val | mu_tcid, sig_tcid),
              normal_lcdf(d_val - 1.0 | mu_tcid, sig_tcid)
            );
          }
        } else if (pfu_type_obs[nn] == 3) {
          real theta = inv_logit(alpha_cult_arg[1] + alpha_cult_arg[2] * pfu_hat_n);
          ll_pfu = bernoulli_lpmf(to_int(pfu_obs[nn]) | theta);
        } else {
          ll_pfu = 0.0;
        }
        ll += ll_pfu;
      }

      // ── RNA likelihood ───────────────────────────────────────────────────
      if (rna_exist_obs[nn] == 1) {
        real ll_rna;
        if (rna_obs[nn] <= lod_rna_arg[k]) {
          ll_rna = normal_lcdf(lod_rna_arg[k] | rna_hat_n, sigma_rna_arg);
          if (test_error) ll_rna += log1m(fn_arg);
        } else {
          ll_rna = normal_lpdf(rna_obs[nn] | rna_hat_n, sigma_rna_arg);
          if (test_error) {
            ll_rna = log_sum_exp(
              log(fp_arg) + exponential_lpdf(rna_obs[nn] - lod_rna_arg[k] | fp_mean_arg[k]),
              log1m(fp_arg) + log1m(fn_arg) + ll_rna
            );
          }
        }
        ll += ll_rna;
      }

      // ── LFD hat & likelihood ─────────────────────────────────────────────
      real lfd_hat_n = inv_logit(  tau0_lfd_arg
                                 + tau_lfd_arg[1] * rna_hat_n
                                 + tau_lfd_arg[2] * pfu_hat_n
                                 + tau_lfd_arg[3] * post_peak_n
                                 + tau_lfd_arg[4] * post_peak_n * rna_hat_n);
      if (source_lfd) lfd_hat_n = inv_logit(logit(lfd_hat_n) + lfd_k_arg[k]);
      if (adj_lfd)    lfd_hat_n = inv_logit(logit(lfd_hat_n) + x_arg[i] * beta_lfd_arg);
      if (lfd_exist_obs[nn] == 1) ll += bernoulli_lpmf(lfd_obs[nn] | lfd_hat_n);

      // ── Symptom onset (discrete-time cloglog hazard) ─────────────────────
      real u_sym   = sigma_sym_arg * z_sym_arg[i];
      real zeta_lin = zeta_sym_intercept_arg
                    + zeta_sym_pfu_arg * (pfu_hat_n / scale_vl)
                    + zeta_sym_rna_arg * (rna_hat_n / scale_vl)
                    + zeta_sym_postpeak_arg * post_peak_n
                    + zeta_sym_postpeak_rna_arg * post_peak_n * (rna_hat_n / scale_vl)
                    + u_sym;
      if (source_sym) zeta_lin = zeta_lin + to_k_sym_arg[k];
      if (adj_sym)    zeta_lin = zeta_lin + x_arg[i] * beta_sym_arg;
      zeta_lin = fmin(zeta_lin, 10.0);
      if (sym_exist_obs[nn] == 1 && sym_at_risk_obs[nn] == 1) {
        if (sym_obs[nn] == 1) {
          ll += log1m_exp(-exp(zeta_lin));
        } else {
          ll += -exp(zeta_lin);
        }
      }
    }
    return ll;
  }
}
