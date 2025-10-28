library(tidyverse)
library(cmdstanr)
library(progressr)
library(deSolve)

handlers("progress")

source("1_code/simstudy/bin.R")

# compile model
mod <- cmdstan_model("1_code/simstudy/sim.stan")

# setup parameter grid
param_grid <- expand.grid(
  sigma = c(1, 2, 3),
  b = 0.000005788,
  d = 2.36,
  c = 974,
  p = 2.15
)

results <- pmap(as.list(param_grid), function(sigma, b, d, c, p, model) {
  # set true parameter values
  params <- c(
    b = b,
    d = d,
    c = c,
    p = p
  )

  # initialize stat variables
  states <- c(
    T = 8e8,
    I = 0,
    V = 1
  )
  
  N_SAMPLES <- 500
  TIMES <- seq(0, 13, 1)
  
  df <- gendata(
    model = targetcell, 
    params = params, 
    states = states, 
    times = TIMES, 
    samples = N_SAMPLES, 
    sigma = sigma, 
    lod = 2.3
  )
  
  with_progress({
    pb <- progressor(steps = N_SAMPLES)
    k <- 0
    
    fits <- map(1:N_SAMPLES, function(x, df, pb, k) {
      dfx <- filter(df, id == x)
      
      fit <- suppressMessages(model$optimize(
        data = list(
          N = 1,
          T = length(TIMES),
          id = x,
          start = dfx |>
            group_by(id) |>
            summarise(
              start = min(time * I(V_obs != 2.3) + 99 * I(V_obs == 2.3))
            ) |>
            ungroup() |> 
            pull(start),
          time = TIMES,
          y = dfx$V_obs,
          y0 = c(8e8, 0, 1),
          lod = 2.3,
          breakpoint = 0
        ),
        init = list(
          list(
            t0 = array(rep(-0.5, 1), dim = 1),
            b = array(rep(b, 1), dim = 1),
            d = array(rep(d, 1), dim = 1),
            c = array(rep(c, 1), dim = 1),
            p = array(rep(p, 1), dim = 1),
            dp_pe = array(rep(13, 1), dim = 1),
            tp_pe = array(rep(0, 1), dim = 1),
            wr_pe = array(rep(6, 1), dim = 1),
            wp_pe = array(rep(2, 1), dim = 1),
            # wf_pe = NULL,
            dp_sm = array(rep(13, 1), dim = 1),
            tp_sm = array(rep(0, 1), dim = 1),
            wr_sm = array(rep(6, 1), dim = 1),
            wp_sm = array(rep(2, 1), dim = 1),
            # wf_sm = NULL,
            sigma_tc = sigma,
            sigma_pe = sigma,
            sigma_sm = sigma
          )
        ),
        iter = 5000,
        refresh = 0
      ))
      
      if (fit$return_codes() != 0) {
        k <- k + 1
        pb(sprintf("Failures %d", k), class = "sticky", amount = 0)
      }
      
      pb(message = sprintf("Runtime %0.3fs", fit$time()))
      
      return(fit)
    }, df = df, pb = pb, k = k)
  })
  
  
  list(
    df = df,
    fits = fits
  )
}, model = mod)
