specd <- function(x, k) trimws(format(round(x, k), nsmall = k))

pefun <- function(t, tp, wp, wr, dp) {
  # Viral load rises before peak: 
  ifelse(
    t <= tp, 
    dp / wp * (t - (tp - wp)), 
    dp - dp / wr * (t - tp) 
  )
}

rep.row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

rep.col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

# convert Ct to rna
ct_to_rna <- function(x, type) {
  
  if (type == "nba") {
    # Kissler, S. M. et al. (2021). Viral Dynamics of SARS-CoV-2 Variants in
    # Vaccinated and Unvaccinated Persons. New England Journal of Medicine,
    # 385(26), 2489–2491. https://doi.org/10.1056/NEJMc2102507
    
    # first convert to CDC equivalent Ct value
    # cdc <- -6.25 + 1.34 * x
    
    # then use CDC calibration equation to convert to log10 viral load
    rna <- (x - 40.93733) / (-3.60971) + log10(250)
    
    # convert to natural log scale?
    rna <- rna * log(10)
    
  } else if (type == "ata") {
    # Singanayagam, A. et al. (2022). Community transmission and viral load
    # kinetics of the SARS-CoV-2 delta (B.1.617.2) variant in vaccinated and
    # unvaccinated individuals in the UK: A prospective, longitudinal, cohort
    # study. The Lancet Infectious Diseases, 22(2), 183–195.
    # https://doi.org/10.1016/S1473-3099(21)00648-4
    
    
    # use PHE calibration equation 
    rna <- (40 - x) / 1.418 + 3.435
    
  } else if (type == "uiuc-cn") {
    # Ke, R., Martinez, P.P., Smith, R.L. et al. Daily longitudinal sampling of
    # SARS-CoV-2 infection reveals substantial heterogeneity in infectiousness.
    # Nat Microbiol 7, 640–652 (2022).
    # https://doi.org/10.1038/s41564-022-01105-z
    
    # use calibration equation for nasal swabs
    rna <- (11.35 - 0.25 * x)
      
    # convert to natural log scale?
    rna <- rna * log(10)
    
  } else if (type == "uiuc-ct") {
    # Ke, R., Martinez, P.P., Smith, R.L. et al. Daily longitudinal sampling of
    # SARS-CoV-2 infection reveals substantial heterogeneity in infectiousness.
    # Nat Microbiol 7, 640–652 (2022).
    # https://doi.org/10.1038/s41564-022-01105-z
    
    # use calibration equation for saliva swabs
    rna <- (14.24 - 0.28 * x)
    
    # convert to natural log scale?
    rna <- rna * log(10)
    
  } else if (type == "hct-cn") {
    # Killingley, B., Mann, A. J., Kalinova, M., Boyers, A., Goonawardane, N.,
    # Zhou, J., Lindsell, K., Hare, S. S., Brown, J., Frise, R., Smith, E.,
    # Hopkins, C., Noulin, N., Löndt, B., Wilkinson, T., Harden, S., McShane,
    # H., Baillet, M., Gilbert, A., … Chiu, C. (2022). Safety, tolerability and
    # viral kinetics during SARS-CoV-2 human challenge in young adults. Nature
    # Medicine, 28(5), Article 5. https://doi.org/10.1038/s41591-022-01780-9
    
  } else if (type == "hct-ct") {
    # Killingley, B., Mann, A. J., Kalinova, M., Boyers, A., Goonawardane, N.,
    # Zhou, J., Lindsell, K., Hare, S. S., Brown, J., Frise, R., Smith, E.,
    # Hopkins, C., Noulin, N., Löndt, B., Wilkinson, T., Harden, S., McShane,
    # H., Baillet, M., Gilbert, A., … Chiu, C. (2022). Safety, tolerability and
    # viral kinetics during SARS-CoV-2 human challenge in young adults. Nature
    # Medicine, 28(5), Article 5. https://doi.org/10.1038/s41591-022-01780-9
    
  }
  
    
  rna 
}

calc_corners <- function(tp, wp, wr, dp, lod) {
  to <- lod * wp / dp + (tp - wp)
  tr <- (dp - lod) * wr / dp + tp 
  
  return(data.frame(to = to, tp = tp, tr = tr))
}

add_corners <- function(id, tp, wp, wr, dp, lod, name = "rna") {
  corner_times <- calc_corners(tp, wp, wr, dp, lod)
  
  df <- data.frame(
    id = rep(id, 3),
    time = unlist(corner_times)
  )
  
  df[[name]] <- c(
    lod, 
    dp,
    lod
  )
  
  rownames(df) <- NULL
  return(df)
}


predict_kinetics <- function(fit, newdata) {
  variable_list <- c(
    "tp_rna",
    "tp_pfu",
    "dp_mean",
    "wp_mean",
    "wr_mean",
    "beta_dp",
    "beta_wp",
    "beta_wr",
    "rho_dp",
    "rho_wp",
    "rho_wr",
    "rho_tp",
    "theta_tp",
    "log_dpi",
    "log_wpi",
    "log_wri",
    "gamma0",
    "gamma"
  )
  
  x_vars <- c(
    "age_[30,50)",
    "age_[50,100)",
    "recurrence",
    # "alpha",
    "delta",
    "omicron",
    "ba4ba5",
    "other",
    "vaccinated_boosted",
    "vaccinated_unboosted",
    "vaccinated_unreported",
    "unvaccinated_unreported",
    "unreported_primary"
  )
  
  id <- newdata[["id"]]
  t <- newdata[["time"]]
  
  lod <- case_when(
    newdata[["source"]] == 1 ~ ct_to_rna(40, type = "nba") + 0.01,
    newdata[["source"]] == 2 ~ ct_to_rna(40, type = "ata") + 0.01,
    newdata[["source"]] == 3 ~ ct_to_rna(47, type = "uiuc-ct") + 0.01,
    newdata[["source"]] == 4 ~ log(1000) + 0.01
  )
  
  lod_pfu <- case_when(
    newdata[["source"]] == 1 ~ 2.3,
    newdata[["source"]] == 2 ~ 2.3,
    newdata[["source"]] == 3 ~ 2.3,
    newdata[["source"]] == 4 ~ log(5)
  )
  
  ids <- length(unique(id))
  count_ids <- table(id)
  
  x <- as.matrix(newdata[, x_vars])
  
  kinetics <- as_draws_rvars(fit$draws(variables = variable_list))
  
  dp_rna <- kinetics$dp_mean * exp(kinetics$log_dpi[id] + x %**% kinetics$beta_dp)
  wp_rna <- kinetics$wp_mean * exp(kinetics$log_wpi[id] + x %**% kinetics$beta_wp)
  wr_rna <- kinetics$wr_mean * exp(kinetics$log_wri[id] + x %**% kinetics$beta_wr)
  
  tp_rna <- kinetics$tp_rna[id]
  to_rna <- lod * wp_rna / dp_rna + (tp_rna - wp_rna)
  tr_rna <- (dp_rna - lod) * wr_rna / dp_rna + tp_rna
  
  rna_hat <- rvar_pefun(t, tp_rna, wp_rna, wr_rna, dp_rna)
  to_rna_hat <- rvar_pefun(to_rna, tp_rna, wp_rna, wr_rna, dp_rna)
  tp_rna_hat <- rvar_pefun(tp_rna, tp_rna, wp_rna, wr_rna, dp_rna)
  tr_rna_hat <- rvar_pefun(tr_rna, tp_rna, wp_rna, wr_rna, dp_rna)
  
  to_rna_hat_mean <- rvar_pefun(mean(to_rna), tp_rna, wp_rna, wr_rna, dp_rna)
  tp_rna_hat_mean <- rvar_pefun(mean(tp_rna), tp_rna, wp_rna, wr_rna, dp_rna)
  tr_rna_hat_mean <- rvar_pefun(mean(tr_rna), tp_rna, wp_rna, wr_rna, dp_rna)
  
  dp_pfu <- dp_rna * exp(kinetics$rho_dp)
  wp_pfu <- wp_rna * exp(kinetics$rho_wp)
  wr_pfu <- wr_rna * exp(kinetics$rho_wr)
  
  tp_pfu <- tp_rna * kinetics$rho_tp + kinetics$theta_tp
  to_pfu <- lod_pfu * wp_pfu / dp_pfu + (tp_pfu - wp_pfu)
  tr_pfu <- (dp_pfu - lod_pfu) * wr_pfu / dp_pfu + tp_pfu
  
  pfu_hat <- rvar_pefun(t, tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  to_pfu_hat <- rvar_pefun(to_pfu, tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  tp_pfu_hat <- rvar_pefun(tp_pfu, tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  tr_pfu_hat <- rvar_pefun(tr_pfu, tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  
  to_pfu_hat_mean <- rvar_pefun(mean(to_pfu), tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  tp_pfu_hat_mean <- rvar_pefun(mean(tp_pfu), tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  tr_pfu_hat_mean <- rvar_pefun(mean(tr_pfu), tp_pfu, wp_pfu, wr_pfu, dp_pfu)
  
  lfd_hat <- kinetics$gamma0 + kinetics$gamma[1] * rna_hat + kinetics$gamma[2] * pfu_hat
  lfd_hat <- exp(lfd_hat) / (1 + exp(lfd_hat))
  
  list(
    rna_hat = rna_hat,
    to_rna = to_rna,
    to_rna_hat = to_rna_hat,
    to_rna_hat_mean = to_rna_hat_mean,
    tp_rna = tp_rna,
    tp_rna_hat = tp_rna_hat,
    tp_rna_hat_mean = tp_rna_hat_mean,
    tr_rna = tr_rna,
    tr_rna_hat = tr_rna_hat,
    tr_rna_hat_mean = tr_rna_hat_mean,
    pfu_hat = pfu_hat,
    to_pfu = to_pfu,
    to_pfu_hat = to_pfu_hat,
    to_pfu_hat_mean = to_pfu_hat_mean,
    tp_pfu = tp_pfu,
    tp_pfu_hat = tp_pfu_hat,
    tp_pfu_hat_mean = tp_pfu_hat_mean,
    tr_pfu = tr_pfu,
    tr_pfu_hat = tr_pfu_hat,
    tr_pfu_hat_mean = tr_pfu_hat_mean,
    lfd_hat = lfd_hat
  )
  
}

expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

rvar_pefun <- rfun(pefun)
rvar_calc_corners <- rfun(calc_corners)


prior_predictive <- function(data, draws = 10) {
  
  id <- data$id
  source <- data$source
  
  src <- rep(1:data$K, data$M)
  
  # population parameter priors
  dp_raw <- rnorm(draws)
  wp_raw <- rnorm(draws)
  wr_raw <- rnorm(draws)
  
  dp_mean_rna <- data$prior_dp_mean * exp(data$prior_dp_cv * dp_raw)
  wp_mean_rna <- data$prior_wp_mean * exp(data$prior_wp_cv * wp_raw)
  wr_mean_rna <- data$prior_wr_mean * exp(data$prior_wr_cv * wr_raw)
  
  # error distribution
  sigma_rna <- truncnorm::rtruncnorm(draws, 0, Inf, 0, data$prior_sigma_sd) 
  sigma_pfu <- truncnorm::rtruncnorm(draws, 0, Inf, 0, data$prior_sigma_sd) 
  
  # transformation parameters 
  tau0_tp <- rnorm(draws)
  tau0_dp <- rnorm(draws)
  tau0_wp <- rnorm(draws)
  tau0_wr <- rnorm(draws)
  tau0_sym <- rnorm(draws)
  
  tau_tp <- rnorm(draws)
  tau_dp <- truncnorm::rtruncnorm(draws, -Inf, 0)
  tau_wp <- truncnorm::rtruncnorm(draws, -Inf, 0)
  tau_wr <- truncnorm::rtruncnorm(draws, -Inf, 0)
  tau_sym <- rnorm(draws)
  
  # individual onset or peak times
  tp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
  tp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
  to_i_sym <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
  
  # lfd coefficients
  tau0_lfd_raw <- rnorm(draws) 
  tau0_lfd <- logit(data$lfd_prior) + 1 * tau0_lfd_raw
  tau_lfd <- mvtnorm::rmvnorm(draws, rep(0, 2))  
  
  # test error rates
  if (data$test_error) {
    fp <- rbeta(draws, data$prior_fp * 50, (1 - data$prior_fp) * 50)
    fn <- rbeta(draws, data$prior_fn * 50, (1 - data$prior_fn) * 50)
  }
  
  # Individual effect priors
  if (data$ind_effects) {
    dp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
    wp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
    wr_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
    
    dp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
    wp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
    wr_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(data$prior_i_sd, draws))) 
  }
  
  # Source effects priors
  if (data$source_pfu) {
    tp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    dp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    wp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    wr_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
  }
  
  if (data$source_rna) {
    tp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    dp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    wp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
    wr_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
  }
  
  if (data$source_sym) {
    to_k_sym <- mvtnorm::rmvnorm(data$K, rep(0, draws), diag(rep(data$prior_k_sd, draws))) 
  }
  
  # Covariate effects priors
  if (data$adj_pfu) {
    beta_dp_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
    beta_wp_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
    beta_wr_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
  }
  
  if (data$adj_rna) {
    beta_dp_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
    beta_wp_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
    beta_wr_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P), diag(rep(data$prior_beta_sd, data$P)))
  }
  
  # initialize RNA trajectory parameters
  tp_rna <- tp_i_rna
  dp_rna <- dp_mean_rna
  wp_rna <- wp_mean_rna
  wr_rna <- wr_mean_rna
  
  # RNA: add individual effects
  if (data$ind_effects) {
    dp_rna <- dp_rna * exp(dp_i_rna)
    wp_rna <- wp_rna * exp(wp_i_rna)
    wr_rna <- wr_rna * exp(wr_i_rna)
  }
  
  # RNA: add source effects
  if (data$source_rna) {
    tp_rna <- tp_rna + tp_k_rna[src, ]
    dp_rna <- dp_rna * exp(dp_k_rna[src, ])
    wp_rna <- wp_rna * exp(wp_k_rna[src, ])
    wr_rna <- wr_rna * exp(wr_k_rna[src, ])
  }
  
  # RNA: add covariate effects
  if (data$adj_rna) {
    dp_rna <- dp_rna * as.vector(exp(as.matrix(data$x) %*% t(beta_dp_rna)))
    wp_rna <- wp_rna * as.vector(exp(as.matrix(data$x) %*% t(beta_wp_rna)))
    wr_rna <- wr_rna * as.vector(exp(as.matrix(data$x) %*% t(beta_wr_rna)))
  }
  
  # transform from RNA -> PFU
  tp_pfu <- rep.row(tau0_tp, sum(data$M)) + rep.row(tau_tp, sum(data$M)) * tp_rna + tp_i_pfu
  
  dp_mean_pfu <- exp(tau0_dp) + exp(tau_dp) * dp_mean_rna
  wp_mean_pfu <- exp(tau0_wp) + exp(tau_wp) * wp_mean_rna
  wr_mean_pfu <- exp(tau0_wr) + exp(tau_wr) * wr_mean_rna
  
  dp_pfu <- exp(rep.row(tau0_dp, sum(data$M))) + exp(rep.row(tau_dp, sum(data$M))) * dp_rna
  wp_pfu <- exp(rep.row(tau0_wp, sum(data$M))) + exp(rep.row(tau_wp, sum(data$M))) * wp_rna
  wr_pfu <- exp(rep.row(tau0_wr, sum(data$M))) + exp(rep.row(tau_wr, sum(data$M))) * wr_rna
  
  # PFU: add individual effects
  if (data$ind_effects) {
    dp_pfu <- dp_pfu * exp(dp_i_pfu) 
    wp_pfu <- wp_pfu * exp(wp_i_pfu) 
    wr_pfu <- wr_pfu * exp(wr_i_pfu) 
  }
  
  # PFU: add source effects
  if (data$source_pfu) {
    tp_pfu <- tp_pfu + tp_k_pfu[source]
    dp_pfu <- dp_pfu * exp(dp_k_pfu[source])
    wp_pfu <- wp_pfu * exp(wp_k_pfu[source])
    wr_pfu <- wr_pfu * exp(wr_k_pfu[source])
  }
  
  # PFU: add covariate effects
  if (data$adj_pfu) {
    dp_pfu <- dp_pfu * exp(data$x %*% t(beta_dp_pfu))
    wp_pfu <- wp_pfu * exp(data$x %*% t(beta_wp_pfu))
    wr_pfu <- wr_pfu * exp(data$x %*% t(beta_wr_pfu))
  }
  
  # transform symptom onset time
  to_sym <- rep.row(tau0_sym, sum(data$M)) + rep.row(tau_sym, sum(data$M)) * tp_rna + to_i_sym
  
  # SYM: add source effects
  if (data$source_sym) {
    to_sym <- to_sym + to_k_sym[src, ]
  }
  
  # initialize matrices 
  rna_hat <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  pfu_hat <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  lfd_hat <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  
  rna <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  pfu <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  lfd <- matrix(rep(NA, sum(data$N) * draws), nrow = sum(data$N))
  
  lod_rna <- data$lod_rna[data$source]
  lod_pfu <- data$lod_pfu[data$source]
  
  for (d in 1:draws) {
    # calculate RNA copies at each time step
    rna_hat[, d] <-  pefun(data$time, tp_rna[data$id, d], wp_rna[data$id, d], wr_rna[data$id, d], dp_rna[data$id, d])
    
    # calculate PFUs of culturable virus at each time step
    pfu_hat[, d] <-  pefun(data$time, tp_pfu[data$id, d], wp_pfu[data$id, d], wr_pfu[data$id, d], dp_pfu[data$id, d])
    
    # calculate LFD probabilities at each time step
    lfd_hat[, d] <- plogis(tau0_lfd[d] + tau_lfd[d, 1] * rna_hat[, d] + tau_lfd[d, 2] * pfu_hat[, d])
    
    # add sampling error
    rna[, d] <- rnorm(sum(data$N), rna_hat[, d], sigma_rna[d])
    pfu[, d] <- rnorm(sum(data$N), pfu_hat[, d], sigma_pfu[d])
    lfd[, d] <- rbinom(sum(data$N), 1, lfd_hat[, d])
    
    # add test error
    if (data$test_error) {
      fp_pfu <- rbinom(sum(data$N), 1, fp[d])
      fn_pfu <- rbinom(sum(data$N), 1, fn[d])
      
      fp_rna <- rbinom(sum(data$N), 1, fp[d])
      fn_rna <- rbinom(sum(data$N), 1, fn[d])
      
      pfu_error <- rexp(sum(data$N), data$fp_mean[data$source])
      rna_error <- rexp(sum(data$N), data$fp_mean[data$source])
      
      pfu[, d] <- (1 - fp_pfu) * (1 - fn_pfu) * pfu[, d] + 
        fp_pfu * (pfu_error + lod_pfu) + 
        fn_pfu * lod_pfu
      
      rna[, d] <- (1 - fp_rna) * (1 - fn_rna) * rna[, d] + 
        fp_rna * (rna_error + lod_rna) + 
        fn_rna * lod_rna
      
    }
    
    # add limits of detection
    rna[rna[, d] < lod_rna, d] <- 
      lod_rna[rna[, d] < lod_rna]
    
    pfu[pfu[, d] < lod_pfu, d] <- 
      lod_pfu[pfu[, d] < lod_pfu]
  }
  
  
  return(
    list(
      id = data$id,
      time = data$time,
      dp_mean_rna = dp_mean_rna,
      wp_mean_rna = wp_mean_rna,
      wr_mean_rna = wr_mean_rna,
      dp_mean_pfu = dp_mean_pfu,
      wp_mean_pfu = wp_mean_pfu,
      wr_mean_pfu = wr_mean_pfu,
      lod_rna = data$lod_rna[src],
      lod_pfu = data$lod_pfu[src],
      tau0_tp = tau0_tp,
      tau0_dp = tau0_dp,
      tau0_wp = tau0_wp,
      tau0_wr = tau0_wr,
      tau0_sym = tau0_sym,
      tau_tp = tau_tp,
      tau_dp = tau_dp,
      tau_wp = tau_wp,
      tau_wr = tau_wr,
      tau_sym = tau_sym,
      dp_rna = dp_rna,
      wp_rna = wp_rna,
      wr_rna = wr_rna,
      tp_rna = tp_rna,
      dp_pfu = dp_pfu,
      wp_pfu = wp_pfu,
      wr_pfu = wr_pfu,
      tp_pfu = tp_pfu,
      rna_hat = rna_hat,
      pfu_hat = pfu_hat,
      lfd_hat = lfd_hat,
      rna = rna,
      pfu = pfu,
      lfd = lfd
    )
  )
}
