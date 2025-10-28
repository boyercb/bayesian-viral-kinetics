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
  // wf: The width of the flat period at peak
  real piecewise(real t, real tp, real wp, real wr, real dp, real wf) {
    if (t <= tp) { // Viral load rises before peak
      return dp / wp * (t - (tp - wp)); 
    } else if (t > tp && t <= (tp + wf)) { // Viral load flat for a period at peak
      return dp; 
    } else { // Viral load falls after peak
      return dp - dp / wr * (t - tp - wf); 
    }
  }
  
  // This function calculates a smooth mixed exponential function.
  // t: The time variable.
  // tp: The peak time.
  // wp: The width of the peak.
  // wr: The width of the right side of the function.
  // dp: The depth of the peak.
  // wf: The width of the flat period at peak
  real smooth(real t, real tp, real wp, real wr, real dp, real wf) {
    real a = dp / wp; // calculate proliferation rate
    real b = dp / wr; // calculate clearance rate
    
    return dp + log(
      (a + b) / 
      (b * exp(-a * (t - tp)) + a * exp(b * (t - tp)) + exp(wf))
    );
  }
}

data {
  int<lower=1> N;     // number of individuals
  int<lower=1> T;     // number of time points
  array[N] int id;    // participant id
  array[N] int start; // effective start time
  array[T] real time; // time point
  array[N*T] real y;  // observed values
  vector[3] y0;       // initial state values
  real lod;           // limit of detection of the test
  int<lower=0, upper=1> breakpoint; // add breakpoint?
}

parameters {
  // target cell parameters
  real<lower=0> sigma_tc;
  array[N] real t0;
  array[N] real<lower=0> b;
  array[N] real<lower=0> d;
  array[N] real<lower=0> p;
  array[N] real<lower=0> c;
  
  // piece-wise exponential parameters
  real<lower=0> sigma_pe;
  array[N] real tp_pe;
  array[N] real<lower=0> dp_pe;
  array[N] real<lower=0> wp_pe;
  array[N] real<lower=0> wr_pe;
  array[breakpoint ? N : 0] real<lower=0> wf_pe;

  
  // smooth mixed exponential parameters
  real<lower=0> sigma_sm;
  array[N] real tp_sm;
  array[N] real<lower=0> dp_sm;
  array[N] real<lower=0> wp_sm;
  array[N] real<lower=0> wr_sm;
  array[breakpoint ? N : 0] real<lower=0> wf_sm;

}

model {
  
  for (n in 1:N) {
    array[T-start[n]] vector[3] mu = log(ode_bdf(targetcell, y0, t0[n],
                                 time[(start[n] + 1):T],
                                 b[n], d[n], p[n], c[n]));
                                 
    for (t in 1:T) {
      int i = t + (n - 1) * T;
      int ts = t - start[n];
      
      real pe;
      real sm;
      
      if (breakpoint) {
        pe = piecewise(time[t], tp_pe[n], wp_pe[n], wr_pe[n], dp_pe[n], wf_pe[n]);
        sm = smooth(time[t], tp_sm[n], wp_sm[n], wr_sm[n], dp_sm[n], wf_sm[n]);
      } else {
        pe = piecewise(time[t], tp_pe[n], wp_pe[n], wr_pe[n], dp_pe[n], 0);
        sm = smooth(time[t], tp_sm[n], wp_sm[n], wr_sm[n], dp_sm[n], negative_infinity());
      }
      
      if (y[i] > lod) {
        if (ts >= 1) {
          target += normal_lpdf(y[i] | mu[ts][3], sigma_tc);
        }
        target += normal_lpdf(y[i] | pe, sigma_pe);
        target += normal_lpdf(y[i] | sm, sigma_sm);
      } else {
        if (ts >= 1) {
          target += normal_lcdf(lod | mu[ts][3], sigma_tc);
        }
        target += normal_lcdf(lod | pe, sigma_pe);
        target += normal_lcdf(lod | sm, sigma_sm);
      }
    }
                                 
  }
  
  
}

