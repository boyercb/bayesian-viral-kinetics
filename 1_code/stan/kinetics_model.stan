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
  real pefun(real t, real tp, real wp, real wr, real dp) {
    if (t <= tp) {
      return dp / wp * (t - (tp - wp)); // Viral load rises before peak
    } else {
      return dp - dp / wr * (t - tp); // Viral load falls after peak
    }
  }
}

data {
  // fixed parameters --------------------------------------------------------
  // The number of data sources. 
  int<lower=0> N_sources;
  
  // The number of observations for each data source.
  array[N_sources] int<lower=0> N;
  
  // The number of unique ids for each data source.
  array[N_sources] int<lower=0> N_ids;
  
  // The number of covariates.
  int<lower=0> P;
  
  // The limit of detection for [RNA].
  array[N_sources] real lod_rna;
  
  // The limit of detection for PFU of culturable virus.
  array[N_sources] real lod_pfu;
  
  // The PCR sensitivity for each data source.
  array[N_sources] real<lower=0, upper=1> sens;
  
  // The mean [log(RNA)] for false positives for each data source.
  array[N_sources] real<lower=0> fpmean;

  // observations ------------------------------------------------------------
  // The id vector for all observations.
  array[sum(N)] int id;
  
  // The time vector for all observations.
  array[sum(N)] real time;
  
  // The RNA copies per ml for all observations.
  vector<lower=0>[sum(N)] rna;
  
  // The culturable virus (PFUs) for all observations.
  vector<lower=0>[sum(N)] pfu;
  
  // The lateral flow test results for all observations.
  array[sum(N)] int<lower=0, upper=1> lfd;
  
  // The individual characteristics for all unique ids.
  array[sum(N_ids)] row_vector[P] x;
  
  // The data source indicator for all observations.
  array[sum(N)] int<lower=0, upper=N_sources> source;
  
  // Missingness indicators
  array[sum(N)] int<lower=0, upper=1> rna_exist; 
  array[sum(N)] int<lower=0, upper=1> pfu_exist; 
  array[sum(N)] int<lower=0, upper=1> lfd_exist; 
  
  // priors ------------------------------------------------------------------
  // Prior mean and sd peak log(rna) copies
  real<lower=0> dp_prior_mean;   
  real<lower=0> dp_prior_sd;     
  
  // Prior mean and sd proliferation duration
  real<lower=0> wp_prior_mean;   
  real<lower=0> wp_prior_sd;     
  
  // Prior mean and sd clearance duration
  real<lower=0> wr_prior_mean;   
  real<lower=0> wr_prior_sd;     
  
  // Prior observation noise for [RNA] and PFUs
  real<lower=0> sigma_prior_sd;  
  
  // Prior sd for individual onset times (days)
  real<lower=0> tp_prior_sd;  
  
  // Prior for LFD logistic intercept
  real<lower=0, upper=1> lfd_prior;  

  // priors for random effects for peak, proliferation, and clearance
  real<lower=0> dpi_prior_sd;    
  real<lower=0> wpi_prior_sd;
  real<lower=0> wri_prior_sd;
}

parameters {
  // observation noise for [RNA] and PFUs of culturable virus
  real<lower=0> sigma_rna;
  real<lower=0> sigma_pfu;
  
  // non-centered parameterization for peak, proliferation, and clearance
  real log_dp_raw;
  real log_wp_raw;
  real log_wr_raw;

  // coefficient vectors for predictors of peak, proliferation, and clearance 
  vector[P] beta_dp;
  vector[P] beta_wp;
  vector[P] beta_wr;

  // population-level offset for peak timing 
  real theta_tp; 

  // coefficients for [RNA] to PFU transformations
  real rho_tp;
  real rho_dp;
  real rho_wp;
  real rho_wr;
   
  // non-centered parameterization for logistic intercept for LFD
  real gamma0_raw;
  
  // coefficient vector for transformation of [RNA] and PFUs to LFD positivity
  vector[2] gamma;

  // individual random effects for onset, peak, proliferation, and clearance
  array[sum(N_ids)] real tp_rna;
  array[sum(N_ids)] real log_dpi;
  array[sum(N_ids)] real log_wpi;
  array[sum(N_ids)] real log_wri;
  
  // cut-point vector for transforming time to cell culture to PFUs
  ordered[4] c;
  
  // coefficient for transforming time to cell culture to PFUs
  real alpha;
  
}

transformed parameters {
  // arrays for the predicted values for [RNA], PFUs of culturable virus, and probability LFD positive
  array[sum(N)] real rna_hat;
  array[sum(N)] real pfu_hat;
  array[sum(N)] real lfd_hat; 

  // arrays for the log-likelihoods 
  array[count_elem(rna_exist, 1)] vector[2] log_lik_rna;
  array[count_elem(pfu_exist, 1)] vector[1] log_lik_pfu;
  array[count_elem(lfd_exist, 1)] vector[1] log_lik_lfd;
  
  // mean peak, proliferation, and clearance duration
  real dp_mean; 
  real wp_mean;
  real wr_mean; 
  
  real gamma0;

  // random effects for peak, proliferation, and clearance duration for [RNA]
  array[sum(N_ids)] real<lower=0> dpi; 
  array[sum(N_ids)] real<lower=0> wpi; 
  array[sum(N_ids)] real<lower=0> wri; 
  
  // individual peak, proliferation, and clearance duration for [RNA]
  array[sum(N_ids)] real<lower=0> dp_rna; 
  array[sum(N_ids)] real<lower=0> wp_rna; 
  array[sum(N_ids)] real<lower=0> wr_rna; 
  
  // individual onset time, peak, proliferation, and clearance duration for PFUs
  array[sum(N_ids)] real tp_pfu;
  array[sum(N_ids)] real<lower=0> dp_pfu; 
  array[sum(N_ids)] real<lower=0> wp_pfu; 
  array[sum(N_ids)] real<lower=0> wr_pfu; 
  
  
  // calculate intercept using noncentered parameterization
  gamma0 = logit(lfd_prior) + 1 * gamma0_raw;

  // calculate mean peak, proliferation, and clearance duration using noncentered parameterization
  dp_mean = dp_prior_mean * exp(sqrt(log(dp_prior_sd / dp_prior_mean^2 + 1)) * log_dp_raw);
  wp_mean = wp_prior_mean * exp(sqrt(log(wp_prior_sd / wp_prior_mean^2 + 1)) * log_wp_raw);
  wr_mean = wr_prior_mean * exp(sqrt(log(wr_prior_sd / wr_prior_mean^2 + 1)) * log_wr_raw);
  
  {
    // counters for log-likelihoods
    int ii_rna;
    int ii_pfu;
    int ii_lfd;
    
    // initialize counters
    ii_rna = 1;
    ii_pfu = 1;
    ii_lfd = 1;
    
    // sample statement
    for (n in 1 : sum(N)) {
      
      // calculate individual random effects (centered)
      dpi[id[n]] = dp_mean * exp(log_dpi[id[n]]);
      wpi[id[n]] = wp_mean * exp(log_wpi[id[n]]);
      wri[id[n]] = wr_mean * exp(log_wri[id[n]]);
      
      // calculate individual peak, proliferation, and clearance from 
      // individual random effects and covariate values 
      dp_rna[id[n]] = dpi[id[n]] * exp(x[id[n], ] * beta_dp);
      wp_rna[id[n]] = wpi[id[n]] * exp(x[id[n], ] * beta_wp);
      wr_rna[id[n]] = wri[id[n]] * exp(x[id[n], ] * beta_wr);
      
      // calculate piecewise exponential for predicted [RNA]
      rna_hat[n] = pefun(time[n], tp_rna[id[n]], wp_rna[id[n]], wr_rna[id[n]], dp_rna[id[n]]);
      
      // if [RNA] is not missing
      if (rna_exist[n] == 1) {
        if (rna[n] <= lod_rna[source[n]]) { // if below limit of detection
          
          // update log-likelihood for [RNA]: mixture based on sensitivity
          log_lik_rna[ii_rna][1] = log1m(sens[source[n]]) + exponential_lpdf(0 | fpmean[source[n]]);
          log_lik_rna[ii_rna][2] = log(sens[source[n]]) + normal_lcdf(lod_rna[source[n]] | rna_hat[n], sigma_rna);
          
        } else { // if above limit

          // update log-likelihood for [RNA]: mixture based on sensitivity
          log_lik_rna[ii_rna][1] = log1m(sens[source[n]]) + exponential_lpdf(rna[n] - lod_rna[source[n]] | fpmean[source[n]]);
          log_lik_rna[ii_rna][2] = log(sens[source[n]]) + normal_lpdf(rna[n] | rna_hat[n], sigma_rna);
        }

        // increment counter
        ii_rna += 1; 
      }
      
      // convert [RNA] to PFUs parameters 
      tp_pfu[id[n]] = rho_tp * tp_rna[id[n]] + theta_tp;
      dp_pfu[id[n]] = dp_rna[id[n]] * exp(rho_dp);
      wp_pfu[id[n]] = wp_rna[id[n]] * exp(rho_wp);
      wr_pfu[id[n]] = wr_rna[id[n]] * exp(rho_wr);
  
      // calculate piecewise exponential for predicted PFUs
      pfu_hat[n] = pefun(time[n], tp_pfu[id[n]], wp_pfu[id[n]], wr_pfu[id[n]], dp_pfu[id[n]]);
  
      // if PFUs of culturable virus are not missing
      if (pfu_exist[n] == 1) {
        
        // if tradional PFU measurement
        if (source[n] == 2 || source[n] == 4) {
          if (pfu[n] <= lod_pfu[source[n]]) { // if below limit of detection
            
            // update log-likelihood for PFU (truncated)
            log_lik_pfu[ii_pfu][1] = normal_lcdf(lod_pfu[source[n]] | pfu_hat[n], sigma_pfu);
          } else { // if above limit of detection
            
            // update log-likelihood for PFU
            log_lik_pfu[ii_pfu][1] = normal_lpdf(pfu[n] | pfu_hat[n], sigma_pfu);
          }
          
        // if days to 50% cell culture
        } else if (source[n] == 3) {
          
          // update log-likelihood for PFU using ordered logistic model
          log_lik_pfu[ii_pfu][1] = ordered_logistic_lpmf(to_int(pfu[n]) | pfu_hat[n] * alpha, c);
        }
        
        // increment counter
        ii_pfu += 1; 
      }
      
      // calculate predicted probability of LFD + based on logistic model
      lfd_hat[n] = inv_logit(gamma0 + gamma[1] * rna_hat[n] + gamma[2] * pfu_hat[n]);
      
      // if LFD is not missing
      if (lfd_exist[n] == 1) {
        // update log-likelihood for LFD
        log_lik_lfd[ii_lfd][1] = bernoulli_lpmf(lfd[n] | lfd_hat[n]);
        
        // increment counter
        ii_lfd += 1; 
      }
        
    }
  }
}
model {
  
  // hierarchical priors:
  log_dp_raw ~ std_normal();
  log_wp_raw ~ std_normal();
  log_wr_raw ~ std_normal();

  // individual random effect priors
  log_dpi ~ normal(0, dpi_prior_sd);
  log_wpi ~ normal(0, wpi_prior_sd);
  log_wri ~ normal(0, wri_prior_sd);
  tp_rna ~ normal(0, tp_prior_sd);
  
  // process noise prior:
  sigma_rna ~ normal(0, sigma_prior_sd) T[0, ]; 
  sigma_pfu ~ normal(0, sigma_prior_sd) T[0, ]; 

  // priors for coefficients (ridge prior)
  beta_dp ~ std_normal(); 
  beta_wp ~ std_normal(); 
  beta_wr ~ std_normal(); 
  
  // priors for transformation coefficients
  rho_tp ~ std_normal();  
  rho_dp ~ std_normal() T[, 0];  
  rho_wp ~ std_normal() T[, 0];  
  rho_wr ~ std_normal() T[, 0];  

  theta_tp ~ std_normal();  

  // priors for LFD logistic regression
  gamma0_raw ~ std_normal();
  gamma ~ std_normal();
  
  // counters for log-likelihoods
  int ii_rna;
  int ii_pfu;
  int ii_lfd;
  
  // initialize counters
  ii_rna = 1;
  ii_pfu = 1;
  ii_lfd = 1;
  
  // sample statement
  for (n in 1 : sum(N)) {
    if (rna_exist[n] == 1) {
      target += log_sum_exp(log_lik_rna[ii_rna][1:2]); // update [RNA] likelihood
      ii_rna += 1; // increment counter
    }
    if (pfu_exist[n] == 1) {
      target += log_lik_pfu[ii_pfu][1]; // update PFU likelihood
      ii_pfu += 1; // increment counter
    }
    if (lfd_exist[n] == 1) {
      target += log_lik_lfd[ii_lfd][1]; // update LFD likelihood
      ii_lfd += 1; // increment counter
    }
  }
}

