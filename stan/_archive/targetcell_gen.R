library(deSolve)
library(tidyverse)

# target cell limited model ode system
targetcell <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dT <- -b * V * T
    dI <- b * V * T - d * I
    dV <- p * I - c * V
    #dR <- c * V - e * R
    
    list(c(dT, dI, dV))
  })
}

# observation process
observe <- function(df, sigma = 1) {
  n <- length(df$V)
  df$V_obs <- rnorm(n, log(df$V), 1)

  return(df)
}

fit <- function(df) {
  
}

# set true parameter values
params <- c(
  b = 0.000005788,
  d = 2.36,
  c = 974,
  p = 2.15
  # e = 0.5
)

# initialize stat variables
states <- c(
  T = 8e8,
  I = 0,
  V = 1
  # R = 0
)

# create vectors for plotting true value
times <- seq(0, 13, by = 0.01)

# samples for DGP
samples <- seq(0, 13, by = 1)

#
true <- ode(y = states, times = times, func = targetcell, parms = params)
true <- as.data.frame(true)

samples <- ode(y = states, times = samples, func = targetcell, parms = params)
samples <- as.data.frame(samples)

# observations process
samples <- observe(samples)


ggplot(samples) +
  geom_point(data = samples, aes(x = time, y = V_obs)) +
  geom_line(data = true, aes(x = time, y = log(V))) +
  # geom_line( aes(x = time, y = log(V) + log(R)), linetype = "dashed") +
  theme_bw()
