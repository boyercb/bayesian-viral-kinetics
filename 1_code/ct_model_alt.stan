functions {
  // Ct calibration functions
  // convert Ct value to log rna copies
  real ct_to_rna(real ct, real kappa, real alpha) {
    return kappa + (40.0 - ct) / alpha;
  }
  
  // convert log rna copies to Ct value
  real rna_to_ct(real rna, real kappa, real alpha) {
    return 40.0 - (rna - kappa) * alpha;
  }
  
  // count number times elem appears in test set
  int count_elem(vector vec, int elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(vec))
      if(vec[i] == elem)
        count = count + 1;
    return(count);
  }
  
  // count number times elem appears in test set
  int count_unique(int[] vec) {
    int count = 1;
    int res[num_elements(vec)] = sort_desc(vec);

    for(i in 1:(num_elements(res) - 1))
      if(res[i] != res[i+1])
        count = count + 1;
    return(count);
  }
  
  // find elements in test which are equal to elem
  int[] which_elem(vector vec, int elem) {
    int res[count_elem(vec, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(vec))
      if(vec[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    return(res);
  }
  // count number of unique items in vector
  // int[] unique(int[] x) {
  //   // error occurs here with declaring res as int[]
  //   int res[num_elements(x)]  = sort_desc(x);
  //   for (i in 1:num_elements(res)) {
  //     res[i] = ...; // to be filled in later
  //   }
  // }
  
  // mean function
  real mufun(real t, real tp, real wp, real wr, real dp) {
    // Viral load rises before peak: 
    if (t <= tp)
      return dp / wp * (t - (tp - wp));
    // Viral load falls after peak: 
    else 
      return dp - dp / wr * (t - tp);
  }
  
}

data {
  // fixed parameters
  int<lower=0> N;                 // number of observations
  int<lower=0> P;                 // number of covariates
  real<lower=0> lod;              // limit of detection
  real kappa;                     // PCR Ct calibration intercept
  real alpha;                     // PCR Ct calibration slope
  real<lower=0, upper=1> sens;    // PCR sensitivity 
  real fpmean;                    // mean Ct value for false positives
  
  // observations
  int id[N];                        // id vector
  real time[N];                     // time vector
  vector<lower=0, upper=lod>[N] ct; // Ct value
  matrix[N, P] x;                   // individual characteristics
  
  // prior arguments
  real<lower=0> priorsd;    // Standard deviation for the lognormal priors
  real sigma_prior[2];      // Prior observation noise Cauchy scale
  real tp_prior[2];         // Prior sd for the onset time (days)
  
  real dp_midpoint;         // Prior mean peak Ct (delta from lod)
  real wp_midpoint;         // Prior mean proliferation duration
  real wr_midpoint;         // Prior mean clearance duration 

}

transformed data {
  // vector[N] rna;                     // RNA copies implied by Ct value
  vector<lower=0, upper=lod>[N] ctdrop; // Ct deviation from LOD
  int<lower=0> n_id;
  
  // rna = ct_to_rna(ct, kappa, alpha);  // calibration equation
  ctdrop = lod - ct;
  n_id = count_unique(id);
  
}

parameters{

  real dp_raw; // peak Ct value (noncentered trans.)
  real wp_raw; // proliferation phase duration (noncentered trans.)
  real wr_raw; // clearance phase duration (noncentered trans.)
 
  real<lower = 0> sigma; // observation noise
  
  real log_dp_mean; // population mean peak Ct value
  real log_wp_mean; // population mean proliferation phase duration
  real log_wr_mean; // population mean clearance phase duration
  
  real<lower = 0> log_dp_sd; // individual heterogeneity in peak Ct
  real<lower = 0> log_wp_sd; // individual heterogeneity proliferation phase
  real<lower = 0> log_wr_sd; // individual heterogeneity clearance phase
  
  vector[P] beta_dp; // coefficients for predictors of peak Ct 
  vector[P] beta_wp; // coefficients for predictors of proliferation
  vector[P] beta_wr; // coefficients for predictors of clearance
  
  real tp[n_id]; // onset time

}

transformed parameters {
  real dp_ref; // peak Ct value
  real<lower=0> wp_ref; // proliferation phase duration
  real<lower=0> wr_ref; // clearance phase duration
 
  real dp[n_id]; // mean peak Ct
  real wp[n_id]; // mean proliferation phase duration
  real wr[n_id]; // mean clearance phase duration
  real mu[N];    // mean Ct value for time and covariate pattern
  vector[2] log_lik[N]; // log likelihood
    
  // non-centered parameterization  
  dp_ref = dp_midpoint * exp(log_dp_mean + log_dp_sd * dp_raw);
  wp_ref = wp_midpoint * exp(log_wp_mean + log_wp_sd * wp_raw);
  wr_ref = wr_midpoint * exp(log_wr_mean + log_wr_sd * wr_raw);
  
   # sample statement
  for (n in 1:N) {
    
    dp[id[n]] = dp_ref * exp(x[id[n], ] * beta_dp);
    wp[id[n]] = wp_ref * exp(x[id[n], ] * beta_wp);
    wr[id[n]] = wr_ref * exp(x[id[n], ] * beta_wr);
    
    // calculate piecewise exponential for predicted Ct
    mu[n] = mufun(time[n], tp[id[n]], wp[id[n]], wr[id[n]], dp[id[n]]);

    if (ct[n] >= lod) { // if above limit of detection
      // target += normal_lcdf(0 | mu[i], sigma);
      log_lik[n][1] = log1m(sens) + exponential_lpdf(ctdrop[n] | 1 / fpmean);
      log_lik[n][2] = log(sens) + normal_lcdf(0 | mu[n], sigma);
    } else { // if below limit
      log_lik[n][1] = log1m(sens) + exponential_lpdf(ctdrop[n] | 1 / fpmean);
      log_lik[n][2] = log(sens) + normal_lpdf(ctdrop[n] | mu[n], sigma);
    }
     
  }
}

model {
  // define priors
  tp ~ normal(tp_prior[1], tp_prior[2]);

  log_dp_mean ~ normal(0, priorsd); // hierarchical mean
  log_dp_sd ~ normal(0, priorsd) T[0,]; // sd of the individual-level draws
  
  log_wp_mean ~ normal(0, priorsd); // hierarchical mean
  log_wp_sd ~ normal(0, priorsd) T[0,]; // sd of the individual-level draws
 
  log_wr_mean ~ normal(0, priorsd); // hierarchical mean
  log_wr_sd ~ normal(0, priorsd) T[0,]; // sd of the individual-level draws
  
  sigma ~ normal(sigma_prior[1], sigma_prior[2]) T[0,]; // process noise sd 
  
  dp_raw ~ std_normal();
  wp_raw ~ std_normal();
  wr_raw ~ std_normal();
    
  // priors for coefficients
  beta_dp ~ std_normal();//normal(0, priorsd);
  beta_wp ~ std_normal();//normal(0, priorsd);
  beta_wr ~ std_normal();//normal(0, priorsd);
  
  # sample statement
  for (n in 1:N) {
    target += log_sum_exp(log_lik[n]);
  }
 
}


