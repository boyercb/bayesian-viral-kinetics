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
    
  // This function calculates a piece-wise exponential function.
  // t: The time variable.
  // tp: The peak time.
  // wp: The width of the peak.
  // wr: The width of the right side of the function.
  // dp: The depth of the peak.
  real piecewise(real t, real tp, real wp, real wr, real dp) {
    if (t <= tp) {
      return dp / wp * (t - (tp - wp)); // Viral load rises before peak
    } else {
      return dp - dp / wr * (t - tp); // Viral load falls after peak
    }
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
  // array[sum(N)] int<lower=0, upper=1> sym; // time-varying symptom indicator
  array[sum(N)] real sym_onset; // symptom onset time
  array[sum(N)] int<lower=0, upper=1> sym_ever; // symptomatic

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
  // int<lower=0, upper=1> ind_corr; // should individual effects be modeled with correlation?
  
  int<lower=0, upper=1> adj_pfu; // should PFU be adjusted?
  int<lower=0, upper=1> adj_rna; // should RNA be adjusted?
  int<lower=0, upper=1> adj_lfd; // should antigen test positivity be adjusted?
  // int<lower=0, upper=1> adj_sym; // should symptom onset be adjusted

  int<lower=0, upper=1> source_pfu; // should source effects be modeled?
  int<lower=0, upper=1> source_rna; // should source effects be modeled?
  int<lower=0, upper=1> source_lfd; // should source effects be modeled?
  int<lower=0, upper=1> source_sym; // should source effects be modeled?

  
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
  
  real prior_fn; // prior on false negative probability
  real prior_fp; // prior on false positive probability
}

transformed data {
  array[sum(N)] int<lower=0, upper=1> first_obs_per_id;
  first_obs_per_id = first_instance(id);
}

parameters {
  // observation noise for [RNA] and PFUs of culturable virus
  real<lower=0> sigma_rna;
  real<lower=0> sigma_pfu;
  
  // random effect sd for symptom onset
  real<lower=0> sigma_to_sym;

  // test characteristics 
  real<lower=0,upper=1> fp; // false positive probability
  real<lower=0,upper=1> fn; // false negative probability
  
  // cholesky_factor_corr[ind_corr ? K : 0] Omega;
  
  // population effects: RNA (non-ceneterd)
  real dp_raw; // peak
  real wp_raw; // proliferation time
  real wr_raw; // clearance time
  
  // individual effects: RNA
  array[sum(M)] real tp_i_rna; // onset
  array[sum(M) && ind_effects ? sum(M) : 0] real dp_i_rna; // peak
  array[sum(M) && ind_effects ? sum(M) : 0] real wp_i_rna; // proliferation time
  array[sum(M) && ind_effects ? sum(M) : 0] real wr_i_rna; // clearance time
  
  // source effects: RNA
  array[K && source_rna ? K : 0] real tp_k_rna; // onset
  array[K && source_rna ? K : 0] real dp_k_rna; // peak
  array[K && source_rna ? K : 0] real wp_k_rna; // proliferation time
  array[K && source_rna ? K : 0] real wr_k_rna; // clearance time
  
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
  
  // individual effects: symptom onset time
  array[sum(M)] real to_i_sym;
  
  // source effects: symptom onset time
  array[K && source_sym ? K : 0] real to_k_sym;

  vector[2] tau_sym;
  
  // real tau0_sym_raw;
  // vector[2] tau_sym;
  // vector[P && adj_sym ? P : 0] beta_sym;
  
}

transformed parameters {
  // predictions for RNA copies, viral culture, antigen positivity, 
  // and symptom onset
  array[sum(N)] real rna_hat, pfu_hat, lfd_hat; 
  
  // population effects: RNA 
  real<lower=0> dp_mean_rna; // peak
  real<lower=0> wp_mean_rna; // proliferation time
  real<lower=0> wr_mean_rna; // clearance time
  
  // population effects: PFU 
  real<lower=0> dp_mean_pfu; // peak
  real<lower=0> wp_mean_pfu; // proliferation time
  real<lower=0> wr_mean_pfu; // clearance time
  
  // arrays for the log-likelihoods 
  array[count_elem(rna_exist, 1)] real log_lik_rna;
  array[count_elem(pfu_exist, 1)] real log_lik_pfu;
  array[count_elem(lfd_exist, 1)] real log_lik_lfd;
  // array[count_elem(sym_exist, 1)] real log_lik_sym;
  
  array[count_elem(add(add(first_obs_per_id, sym_exist), sym_ever), 1)] real log_lik_to_sym;


  // calculate intercept using noncentered parameterization
  real tau0_lfd = logit(prior_lfd_mean) + 1 * tau0_lfd_raw;
  
  
   // matrix[P, K] eta;
   //  
   //  if (ind_corr) {
   //    // Cholesky factor of the covariance matrix
   //    matrix[K, K] L;
   //    Sigma = diag_pre_multiply(sigma, Omega);
   //    
   //    // Calculate per-infection correlated effects
   //    eta = (Sigma * ind_eta)';
   //  } else if (ind_var_m) {
   //    // All effects are independent
   //    for (i in 1 : K) {
   //      eta[1 : P, i] = to_vector(ind_eta[i, 1 : P]) * ind_var[i];
   //    }
   //  } else {
   //    // No infection level differences
   //    eta = rep_matrix(0, P, K);
   //  }
    
  
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
    
    // combined effects: symptoms
    real to_sym; // onset
    
    // alternative viral culture transformations
    real lambda; // time to TCID50 if positive
    real theta;  // probability positive / negative
    
    // counters for log-likelihoods
    int ii_rna = 1, ii_pfu = 1, ii_lfd = 1, ii_to_sym = 1;

    // sample statement
    for (n in 1 : sum(N)) {
      
      // always include a random effect for peak time
      tp_rna = tp_i_rna[id[n]]; 
      
      // RNA: calculate population means
      dp_mean_rna = prior_dp_mean * exp(prior_dp_cv * dp_raw);
      wp_mean_rna = prior_wp_mean * exp(prior_wp_cv * wp_raw);
      wr_mean_rna = prior_wr_mean * exp(prior_wr_cv * wr_raw);
      
      dp_rna = dp_mean_rna;
      wp_rna = wp_mean_rna;
      wr_rna = wr_mean_rna;
      
      // RNA: update with individual random effects 
      if (ind_effects) {
        dp_rna = dp_rna * exp(dp_i_rna[id[n]]);
        wp_rna = wp_rna * exp(wp_i_rna[id[n]]);
        wr_rna = wr_rna * exp(wr_i_rna[id[n]]);
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
      
      // PFU: calculate population means
      dp_mean_pfu = exp(tau_dp[1]) + exp(tau_dp[2]) * dp_mean_rna;
      wp_mean_pfu = exp(tau_wp[1]) + exp(tau_wp[2]) * wp_mean_rna;
      wr_mean_pfu = exp(tau_wr[1]) + exp(tau_wr[2]) * wr_mean_rna;
      
      // PFU: convert RNA to PFU parameters 
      tp_pfu = tau_tp[1] + tau_tp[2] * tp_rna + tp_i_pfu[id[n]];
      dp_pfu = exp(tau_dp[1]) + exp(tau_dp[2]) * dp_rna;
      wp_pfu = exp(tau_wp[1]) + exp(tau_wp[2]) * wp_rna;
      wr_pfu = exp(tau_wr[1]) + exp(tau_wr[2]) * wr_rna;
      
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
      
      // calculate piecewise exponential for predicted PFUs
      pfu_hat[n] = piecewise(time[n], tp_pfu, wp_pfu, wr_pfu, dp_pfu);
                   
      // if PFUs of culturable virus are not missing
      if (pfu_exist[n] == 1) {
        
        // if tradional PFU measurement
        if (pfu_type[n] == 1) {
          if (pfu[n] <= lod_pfu[source[n]]) { // if below limit of detection
          
            // update log-likelihood for PFU: censored below limit
            log_lik_pfu[ii_pfu] = normal_lcdf(lod_pfu[source[n]] | pfu_hat[n], sigma_pfu);
            
            // if test error, assume mixture based on test error characteristics
            if (test_error) {
              log_lik_pfu[ii_pfu] += log1m(fn);
            }
          } else { // if above limit of detection
            
            // update log-likelihood for PFU
            log_lik_pfu[ii_pfu] = normal_lpdf(pfu[n] | pfu_hat[n], sigma_pfu);
            
            // if test error, assume mixture based on test error characteristics
            if (test_error) {
              log_lik_pfu[ii_pfu] = log_sum_exp(
                log(fp) + exponential_lpdf(pfu[n] - lod_pfu[source[n]] | fp_mean[source[n]]),
                log1m(fp) + log1m(fn) + log_lik_pfu[ii_pfu]
              );
            }
          }
          
        // if days to 50% cell culture
        } else if (pfu_type[n] == 2) {
          theta = inv_logit(alpha_tcid50[1] + alpha_tcid50[2] * pfu_hat[n]);
          // lambda = 1 / exp(alpha_tcid50[3] + alpha_tcid50[4] * pfu_hat[n]);
          // 
          // // update log-likelihood for PFU using poisson hurdle model
          // if (pfu[n] < 5.9) {
          //   log_lik_pfu[ii_pfu] = bernoulli_lpmf(0 | theta) + poisson_lpmf(to_int(pfu[n]) | lambda) - log1m_exp(-lambda);
          // } else {
          //   log_lik_pfu[ii_pfu] = bernoulli_lpmf(1 | theta) + poisson_lccdf(6 | lambda);
          // }

          // update log-likelihood for PFU using simple logistic model
          if (pfu[n] < 5.9) {
            log_lik_pfu[ii_pfu] = bernoulli_lpmf(1 | theta);
          } else {
            log_lik_pfu[ii_pfu] = bernoulli_lpmf(0 | theta);
          }

          // // update log-likelihood for PFU using latent specification
          // if (pfu[n] < 5.9) {
          //   log_lik_pfu[ii_pfu] = normal_lccdf(alpha_tcid50[1] | pfu_hat[n], sigma_pfu);
          // } else {
          //   log_lik_pfu[ii_pfu] = normal_lcdf(alpha_tcid50[1] | pfu_hat[n], sigma_pfu);
          // }
        // if simple positive/negative culture result
        } else if (pfu_type[n] == 3) {
          theta = inv_logit(alpha_cult[1] + alpha_cult[2] * pfu_hat[n]);

          log_lik_pfu[ii_pfu] = bernoulli_lpmf(to_int(pfu[n]) | theta);
        }
        
        // if (is_inf(log_lik_pfu[ii_pfu])) {
        //   print("PFU:");
        //   print(time[n]);
        // }
        
        // increment counter
        ii_pfu += 1; 
      } 
  
      // calculate piecewise exponential for predicted [RNA]
      rna_hat[n] = piecewise(time[n], tp_rna, wp_rna, wr_rna, dp_rna);

      // if [RNA] is not missing
      if (rna_exist[n] == 1) {
        if (rna[n] <= lod_rna[source[n]]) { // if below limit of detection
          
          // update log-likelihood for [RNA]: censored below limit
          log_lik_rna[ii_rna] = normal_lcdf(lod_rna[source[n]] | rna_hat[n], sigma_rna);
          
          // if test error, assume mixture based on test error characteristics
          if (test_error) {
            log_lik_rna[ii_rna] += log1m(fn);
          }
        } else { // if above limit

          // update log-likelihood for [RNA]: perfect test
          log_lik_rna[ii_rna] = normal_lpdf(rna[n] | rna_hat[n], sigma_rna);
          
          // if test error, assume mixture based on test error characteristics
          if (test_error) {
            log_lik_rna[ii_rna] = log_sum_exp(
              log(fp) + exponential_lpdf(rna[n] - lod_rna[source[n]] | fp_mean[source[n]]),
              log1m(fp) + log1m(fn) + log_lik_rna[ii_rna]
            );
          }
        }
        if (is_inf(log_lik_rna[ii_rna])) {
          print("RNA:");
          print(rna_hat[n]);
        }

        // increment counter
        ii_rna += 1; 
      }
      
      // calculate predicted probability of LFD + based on logistic model
      lfd_hat[n] = inv_logit(tau0_lfd + tau_lfd[1] * rna_hat[n] + tau_lfd[2] * pfu_hat[n]);
      
      // update with covariate effects if requested
      if (adj_lfd) {
        lfd_hat[n] = inv_logit(logit(lfd_hat[n]) + x[id[n], ] * beta_lfd);
      }
      
      // if LFD is not missing
      if (lfd_exist[n] == 1) {
        // update log-likelihood for LFD
        log_lik_lfd[ii_lfd] = bernoulli_lpmf(lfd[n] | lfd_hat[n]);
        
        if (is_inf(log_lik_lfd[ii_lfd])) {
          // print("LFD:");
          // print(lfd_hat[n]);
          // print("ID:");
          // print(id[n]);
          // print("TIME:");
          // print(time[n]);
          // print("RNA");
          // print(rna_hat[n]);
          // print("PFU");
          // print(pfu_hat[n]);
          // print("INTERCEPT:");
          // print(tau0_lfd);
          // print("slope RNA:");
          // print(tau_lfd[1]);
          // print("slope PFU:");
          // print(tau_lfd[2]);
        }

        // increment counter
        ii_lfd += 1; 
      }
      
      // calculate symptom onset using random effect from peak
      to_sym = tau_sym[1] + tau_sym[2] * tp_rna + to_i_sym[id[n]];
      
      if (source_sym) {
        to_sym = to_sym + to_k_sym[source[n]];
      }

      if (sym_exist[n] == 1 && sym_ever[n] == 1 && first_obs_per_id[n] == 1) {
        log_lik_to_sym[ii_to_sym] = normal_lpdf(sym_onset[n] | to_sym, sigma_to_sym);
        
        if (is_inf(log_lik_to_sym[ii_to_sym])) {
          print("SYM:");
          print(to_sym);
        }
        // increment counter
        ii_to_sym += 1;
      }
      // // calculate predicted probability of symptom onset based on trajectory
      // sym_hat[n] = inv_logit();
      // 
      // // if symptom onset is not missing
      // if (sym_exist[n] == 1) {
      //   // update log-likelihood for sym
      //   log_lik_sym[ii_sym] = bernoulli_lpmf(sym[n] | sym_hat[n]);
      //   
      //   // increment counter
      //   ii_sym += 1; 
      // }
    }
  }
}

model {
  // PRIORS
  fp ~ beta_proportion(prior_fp, 50);
  fn ~ beta_proportion(prior_fn, 50);
  
  sigma_rna ~ normal(0, prior_sigma_sd) T[0, ]; 
  sigma_pfu ~ normal(0, prior_sigma_sd) T[0, ]; 
  
  dp_raw ~ std_normal();
  wp_raw ~ std_normal();
  wr_raw ~ std_normal();
   
  tp_i_pfu ~ normal(0, prior_i_sd);
  tp_i_rna ~ normal(0, prior_i_sd);
  to_i_sym ~ normal(0, prior_i_sd);
  
  if (ind_effects) {
    dp_i_pfu ~ normal(0, prior_i_sd);
    wp_i_pfu ~ normal(0, prior_i_sd);
    wr_i_pfu ~ normal(0, prior_i_sd);
    
    dp_i_rna ~ normal(0, prior_i_sd);
    wp_i_rna ~ normal(0, prior_i_sd);
    wr_i_rna ~ normal(0, prior_i_sd);
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
  
  tau_tp ~ std_normal();
  tau_sym ~ std_normal();
  
  tau_dp[1] ~ std_normal();  
  tau_wp[1] ~ std_normal();  
  tau_wr[1] ~ std_normal(); 

  tau_dp[2] ~ std_normal() T[, 0];  
  tau_wp[2] ~ std_normal() T[, 0];  
  tau_wr[2] ~ std_normal() T[, 0];  
  
  alpha_tcid50 ~ normal(0, prior_beta_sd);
  alpha_cult ~ normal(0, prior_beta_sd);
  
  tau0_lfd_raw ~ std_normal();
  tau_lfd ~ std_normal();
  
  // counters for log-likelihoods
  int ii_rna = 1, ii_pfu = 1, ii_lfd = 1, ii_to_sym = 1;
  
  // LIKELIHOOD
  for (n in 1 : sum(N)) {
    if (rna_exist[n] == 1) {
      target += log_lik_rna[ii_rna]; // update [RNA] likelihood
      ii_rna += 1; // increment counter
    }
    if (pfu_exist[n] == 1) {
      target += log_lik_pfu[ii_pfu]; // update PFU likelihood
      ii_pfu += 1; // increment counter
    }
    if (lfd_exist[n] == 1) {
      target += log_lik_lfd[ii_lfd]; // update LFD likelihood
      ii_lfd += 1; // increment counter
    }
    if (sym_exist[n] == 1 && sym_ever[n] == 1 && first_obs_per_id[n] == 1) {
      target += log_lik_to_sym[ii_to_sym]; // update symptom onset likelihood
      ii_to_sym += 1; // increment counter
    }
    
    // if (sym_exist[n] == 1) {
    //   target += log_lik_sym[ii_sym]; // update symptom likelihood
    //   ii_sym += 1; // increment counter
    // }
  }
}


