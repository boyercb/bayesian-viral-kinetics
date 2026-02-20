functions {
  // count number times elem appears in test set
  int count_unique(array[] int vec) {
    int count = 1;
    array[num_elements(vec)] int res = sort_desc(vec);

    for(i in 1:(num_elements(res) - 1))
      if(res[i] != res[i+1])
        count = count + 1;
    return(count);
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
  int<lower=0> N;   // number of observations
  array[N] int id;  // vector of ids
  vector[N] t;      // sample times
  vector[N] v_obs;  // observed outcome
  real lod;         // level of detection of instrument
  real mu;          // mean of fp error distribution
}

parameters {
  array[count_unique(id)] real<lower=5> dp;
  array[count_unique(id)] real<lower=0, upper=10> wp;
  array[count_unique(id)] real<lower=0, upper=15> wr;
  array[count_unique(id)] real tp;
  real<lower=0> sigma;
  real<lower=0,upper=1> fp;
  real<lower=0,upper=1> fn;
}

model {
   // PRIORS
   dp ~ normal(13, 1);
   wp ~ normal(4, 1);
   wr ~ normal(7, 1);
   tp ~ std_normal();
   
   sigma ~ normal(0, 2) T[0, ];
   fp ~ beta_proportion(0.1, 50);
   fn ~ beta_proportion(0.1, 50);

   // LIKELIHOOD 
   for (n in 1:N) {
     real v = pefun(t[n], tp[id[n]], wp[id[n]], wr[id[n]], dp[id[n]]);
     if (v_obs[n] == lod) {
       target += log1m(fn) + normal_lcdf(lod | v, sigma);
     } else {
       target += log_sum_exp(
         log(fp) + exponential_lpdf(v_obs[n] | mu),
         log1m(fp) + log1m(fn) + normal_lpdf(v_obs[n] | v, sigma)
       );
     }
   }
}

generated quantities {
  array[N] real v;
  
  for (n in 1:N) {
     v[n] = pefun(t[n], tp[id[n]], wp[id[n]], wr[id[n]], dp[id[n]]);
  }
}