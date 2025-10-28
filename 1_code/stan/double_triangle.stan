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
  real prior_fn;
  real prior_fp;
  array[K] real fp_mean;
  
  // observed data -----------------------------------------------------------
  array[sum(N)] int id; // participant id
  array[sum(N)] real time; // time of observation
  vector[sum(N)] rna; // RNA copies per ml from PCR
  vector[sum(N)] pfu; // PFU of culturable virus 
  array[sum(N)] int<lower=0> source; // source indicator
  array[sum(N)] int<lower=0> pfu_type; // source indicator

  // missingness indicators
  array[sum(N)] int<lower=0, upper=1> rna_exist; 
  array[sum(N)] int<lower=0, upper=1> pfu_exist; 

  // model specification -----------------------------------------------------
  int<lower=0, upper=1> test_error; // should test error be included?
  // int<lower=0, upper=1> ind_effects; // should individual effects be modeled?
  // int<lower=0, upper=1> ind_corr; // should individual effects be modeled with correlation?
  // int<lower=0, upper=1> adj_pfu; // should PFU be adjusted
  // int<lower=0, upper=1> adj_rna; // should RNA be adjusted
  // int<lower=0, upper=1> adj_lfd; // should antigen test positivity be adjusted
  // int<lower=0, upper=1> adj_sym; // should symptom onset be adjusted
  
  // priors ------------------------------------------------------------------
  // array[K] real<lower=0, upper=1> sens_rna; // 
  // array[K] real<lower=0, upper=1> sens_pfu; // 
  // array[K] real<lower=0, upper=1> spec_rna; // 
  // array[K] real<lower=0, upper=1> spec_pfu; // 
  
  // real<lower=0> prior_dp; // prior on 
  // real<lower=0> prior_dp_sd; 
  // real<lower=0> prior_dpi_sd; 
  // 
  // real<lower=0> prior_wp; 
  // real<lower=0> prior_wp_sd; 
  // real<lower=0> prior_wpi_sd; 
  // 
  // real<lower=0> prior_wr; 
  // real<lower=0> prior_wr_sd; 
  // real<lower=0> prior_wri_sd; 
  // 
  // real<lower=0> prior_tp_sd; 
  // real<lower=0> prior_tpi_sd; 
  // 
  // real<lower=0> prior_sigma_sd; 
  
}

parameters {
  // observation noise for [RNA] and PFUs of culturable virus
  real<lower=0> sigma_rna;
  real<lower=0> sigma_pfu;
  
  real<lower=0,upper=1> fp;
  real<lower=0,upper=1> fn;
  
  // population effects: PFU 
  real<lower=0> dp_mean_pfu; // peak
  real<lower=0> wp_mean_pfu; // proliferation time
  real<lower=0> wr_mean_pfu; // clearance time
  
  // individual effects: PFU
  array[sum(M)] real tp_i_pfu; // onset
  array[sum(M)] real dp_i_pfu; // peak
  array[sum(M)] real wp_i_pfu; // proliferation time
  array[sum(M)] real wr_i_pfu; // clearance time
  
  // transformation parameters
  real tau_tp; // onset 
  real tau_dp; // peak 
  real tau_wp; // proliferation time 
  real tau_wr; // clearance time 
  
  real tau0_tp; // onset 
  real tau0_dp; // peak 
  real tau0_wp; // proliferation time 
  real tau0_wr; // clearance time 
  
  // individual effects: RNA
  array[sum(M)] real tp_i_rna;
  array[sum(M)] real dp_i_rna;
  array[sum(M)] real wp_i_rna;
  array[sum(M)] real wr_i_rna;
  
  // coefficient for transforming time to cell culture to PFUs
  vector[2] alpha;
}

transformed parameters {
  // predictions for RNA copies, viral culture, antigen positivity, 
  // and symptom onset
  array[sum(N)] real rna_hat, pfu_hat; 
  
  // // population effects: RNA 
  // real<lower=0> dp_mean_rna; // peak
  // real<lower=0> wp_mean_rna; // proliferation time
  // real<lower=0> wr_mean_rna; // clearance time
  
  // arrays for the log-likelihoods 
  array[count_elem(rna_exist, 1)] real log_lik_rna;
  array[count_elem(pfu_exist, 1)] real log_lik_pfu;
  
  {
    // combined effects: PFU (non-centered)
    real tp_pfu; // onset
    real dp_pfu; // peak
    real wp_pfu; // proliferation time
    real wr_pfu; // clearance time
    
    // combined effects: RNA (non-centered)
    real tp_rna; // onset
    real dp_rna; // peak
    real wp_rna; // proliferation time
    real wr_rna; // clearance time
  
    // counters for log-likelihoods
    int ii_rna = 1, ii_pfu = 1;

    // sample statement
    for (n in 1 : sum(N)) {
      
      // calculate individual random effects (centered)
      tp_pfu = tp_i_pfu[id[n]];
      dp_pfu = dp_mean_pfu * exp(dp_i_pfu[id[n]]);
      wp_pfu = wp_mean_pfu * exp(wp_i_pfu[id[n]]);
      wr_pfu = wr_mean_pfu * exp(wr_i_pfu[id[n]]);
      
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
          
          if (pfu[n] < 5.9) {
            log_lik_pfu[ii_pfu] = neg_binomial_lpmf(to_int(pfu[n]) | 1, inv_logit(alpha[1] + alpha[2] * pfu_hat[n])); 
          } else {
            log_lik_pfu[ii_pfu] = neg_binomial_lccdf(6 | 1, inv_logit(alpha[1] + alpha[2] * pfu_hat[n]));
          }
          // update log-likelihood for PFU using ordered logistic model
          
        }
        
        // increment counter
        ii_pfu += 1; 
      } 
  
      // convert [RNA] to PFUs parameters 
      tp_rna = tau0_tp + tp_pfu * tau_tp + tp_i_rna[id[n]];
      dp_rna = (exp(tau0_dp) + exp(tau_dp) * dp_pfu) * exp(dp_i_rna[id[n]]);
      wp_rna = (exp(tau0_wp) + exp(tau_wp) * wp_pfu) * exp(wp_i_rna[id[n]]);
      wr_rna = (exp(tau0_wr) + exp(tau_wr) * wr_pfu) * exp(wr_i_rna[id[n]]);
      
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

        // increment counter
        ii_rna += 1; 
      }
    }
  }
}

model {
  // PRIORS
  fp ~ beta_proportion(prior_fp, 50);
  fn ~ beta_proportion(prior_fn, 50);
   
  tp_i_pfu ~ normal(0, 1);
  dp_i_pfu ~ normal(0, 1);
  wp_i_pfu ~ normal(0, 1);
  wr_i_pfu ~ normal(0, 1);
  
  tp_i_rna ~ normal(0, 1);
  dp_i_rna ~ normal(0, 1);
  wp_i_rna ~ normal(0, 1);
  wr_i_rna ~ normal(0, 1);
  
  // counters for log-likelihoods
  int ii_rna = 1, ii_pfu = 1;
  
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
  }
}


