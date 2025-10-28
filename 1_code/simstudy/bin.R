# target cell limited model ode system
targetcell <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dT <- -b * V * T
    dI <- b * V * T - d * I
    dV <- p * I - c * V

    list(c(dT, dI, dV))
  })
}

# observation process
observe <- function(y, sigma = 1, lod = 2.3) {
  y_obs <- rnorm(length(y), log(y), sigma)
  y_obs <- ifelse(y_obs <= lod, lod, y_obs)
  
  return(y_obs)
}

# data generation function
gendata <- function(model, params, states, times, samples, sigma, lod) {
  # if (x) {
  #   b <- rnorm(samples, log(params$b), cv)
  #   d <- rnorm(samples, log(params$d), cv)
  #   c <- rnorm(samples, log(params$c), cv)
  #   p <- rnorm(samples, log(params$p), cv)
  #   
  # }
  df <- ode(y = states, times = times, func = model, parms = params)
  df <- df[rep(1:nrow(df), samples),  ]
  df <- as.data.frame(df)
  
  df$V_obs <- observe(df$V, sigma, lod)
  df$id <- rep(1:samples, each = length(times))
  return(df)
}

# calc stats
calc_corners <- function(tp, wp, wr, dp, lod) {
  to <- lod * wp / dp + (tp - wp)
  tr <- (dp - lod) * wr / dp + tp 
  
  return(data.frame(to = to, tp = tp, tr = tr))
}


piecewise <- function(t, tp, wp, wr, dp, wf) {
  I(t <= tp) * (dp / wp * (t - (tp - wp))) +
    I(t > tp & t <= tp + wf) * dp +
    I(t > tp + wf) * (dp - dp / wr * (t - tp - wf))
}

smixed <- function(t, tp, wp, wr, dp, wf) {
  a <- dp / wp
  b <- dp / wr
  wf <- replace(wf, wf == 0, -Inf)
  wf <- wf
  dp + log((a + b) / (b * exp(-a * (t - tp)) + a * exp(b * (t - tp)) + exp(wf)))
}


get_truth <- function() {
  params <- c(
    b = 0.000005788,
    d = 2.36,
    c = 974,
    p = 2.15
  )
  
  states <- c(
    T = 8e8,
    I = 0,
    V = 1
  )
  
  p <- ode(
    y = states,
    times = seq(0, 13, 0.01),
    func = targetcell,
    parms = params
  )
  
  df <- as.data.frame(p)
  df$V <- ifelse(df$V < exp(2.3), exp(2.3), df$V)
  
  tibble(
    dp_truth = max(log(df$V)),
    tp_truth = max(I(log(df$V) == max(log(df$V))) * df$time),
    wp_truth = tp_truth - min(df$time * I(log(df$V) > 2.3) + I(log(df$V) <= 2.3) * 99),
    wr_truth = max(df$time * I(log(df$V) > 2.3)) - tp_truth,
  )
}
# plot


