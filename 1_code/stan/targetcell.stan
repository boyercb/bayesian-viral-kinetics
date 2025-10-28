functions {
  // This function specifies the ode system for the target cell model.
  // t: the time variable.
  // y: the state vector.
  // b: the transmission rate between virons and target cells.
  // d: the rate of death or clearance of infected cells.
  // p: the rate of virons produced by infected cells per time step.
  // c: the rate of death or clearance of virons.
  vector targetcell(real t,
                    vector y,     
                    real b,
                    real d,
                    real p,
                    real c) {  
    vector[3] dydt;
    dydt[1] = -b * y[3] * y[1];
    dydt[2] = b * y[3] * y[1] - d * y[2];
    dydt[3] = p * y[2] - c * y[3];
    return dydt;
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
  int<lower=1> T;         // number samples
  array[T] vector[2] y;   // observed values
  array[T] real ts;       // observed times
  vector[4] y0;           // initial state values
  array[2] real lod;      // limit of detection of the test
  int<lower=1> Tp;        // number of prediction time points
  array[Tp] real tsp;     // prediction times
}

parameters {
  // target cell parameters
  real<lower=0> sigma;
  real t0;
  real<lower=0> b;
  real<lower=0> d;
  real<lower=0> p;
  real<lower=0> c;
  
  // piece-wise exponential parameters
  real<lower=0> sigma_pe;
  real tp;
  real<lower=0> dp;
  real<lower=0> wp;
  real<lower=0> wr;

}

model {
  array[T] vector[4] mu = log(ode_rk45_tol(targetcell, y0, t0, ts,
                                 1e-9,
                                 1e-9,
                                 1000000,
                                 b, d, p, c));
  for (t in 1:T) {
    real pe = pefun(ts[t], tp, wp, wr, dp);
    real pe_rna = pefun(ts[t], tp_rna, wp_rna, wr_rna, dp_rna);

    if (y[t][1] > lod[1]) {
      target += normal_lpdf(y[t][1] | mu[t][3], sigma);
      target += normal_lpdf(y[t][1] | pe, sigma_pe);

    } else {
      target += normal_lcdf(lod[1] | mu[t][3], sigma);
      target += normal_lcdf(lod[1] | pe, sigma_pe);

    }
      
    if (y[t][2] > lod[2]) {
      target += normal_lpdf(y[t][2] | mu[t][4], sigma);
      target += normal_lpdf(y[t][2] | pe_rna, sigma_pe);

    } else {
      target += normal_lcdf(lod[2] | mu[t][4], sigma);
      target += normal_lcdf(lod[2] | pe_rna, sigma_pe);

    }
  }
}

generated quantities {
  array[Tp] vector[2] y_sim_pe;
  array[Tp] vector[2] y_sim_ode;
  
  array[Tp] vector[4] y_sim = ode_rk45(targetcell, y0, t0, tsp,
                                            b, d, p, c, q, e);

  for (t in 1:Tp) {
    y_sim_ode[t][1] = log(y_sim[t][3]);
    y_sim_ode[t][2] = log(y_sim[t][4]);
    y_sim_pe[t][1] = pefun(tsp[t], tp, wp, wr, dp);
    y_sim_pe[t][2] = pefun(tsp[t], tp_rna, wp_rna, wr_rna, dp_rna);
  }
}

