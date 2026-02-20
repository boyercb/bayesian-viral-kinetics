# ──────────────────────────────────────────────────────────────────────────────
# plots.R — Prior predictive and ODE comparison plots
# ──────────────────────────────────────────────────────────────────────────────

#' Plot prior predictive check summaries
#'
#' @param pp        Output of prior_predictive()
#' @param stan_data Stan data list
#' @param out_dir   Output directory for figures
#' @return Character vector of saved file paths
plot_prior_predictive <- function(pp, stan_data, out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c()

  # population-level parameter draws
  pp_plot <- tibble::tibble(
    dp_mean_rna = pp$dp_mean_rna,
    wp_mean_rna = pp$wp_mean_rna,
    wr_mean_rna = pp$wr_mean_rna,
    dp_mean_pfu = pp$dp_mean_pfu,
    wp_mean_pfu = pp$wp_mean_pfu,
    wr_mean_pfu = pp$wr_mean_pfu,
    tau0_dp = pp$tau0_dp,
    tau0_wp = pp$tau0_wp,
    tau0_wr = pp$tau0_wr,
    tau_dp  = pp$tau_dp,
    tau_wp  = pp$tau_wp,
    tau_wr  = pp$tau_wr,
    tau0_tp = pp$tau0_tp,
    tau_tp  = pp$tau_tp,
    eta_sym_intercept = pp$eta_sym_intercept,
    eta_sym_pfu       = pp$eta_sym_pfu,
    eta_sym_rna       = pp$eta_sym_rna,
    sigma_sym         = pp$sigma_sym
  )

  # optionally add wf_mean prior draws
  if (!is.null(pp$wf_mean_pp)) {
    pp_plot$wf_mean <- pp$wf_mean_pp
  }

  # model type label for plot titles
  model_label <- paste0(
    if (isTRUE(pp$use_smooth)) "smooth" else "piecewise",
    if (isTRUE(pp$use_wf)) " + flat-top" else ""
  )

  # transformation parameters
  trans_vars <- c("tau0_dp", "tau_dp", "tau0_wp", "tau_wp",
                  "tau0_wr", "tau_wr")
  p1 <- pp_plot |>
    tidyr::pivot_longer(dplyr::all_of(trans_vars)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Prior: transformation parameters", model_label))

  f1 <- file.path(out_dir, "prior_trans_pe.pdf")
  ggplot2::ggsave(f1, p1, width = 10, height = 6)
  files <- c(files, f1)

  # population kinetics parameters
  pop_vars <- c("dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
                "dp_mean_pfu", "wp_mean_pfu", "wr_mean_pfu")
  if (!is.null(pp$wf_mean_pp)) pop_vars <- c(pop_vars, "wf_mean")
  p2 <- pp_plot |>
    tidyr::pivot_longer(dplyr::all_of(pop_vars)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Prior: population kinetics parameters", model_label))

  f2 <- file.path(out_dir, "prior_pe.pdf")
  ggplot2::ggsave(f2, p2, width = 10, height = 6)
  files <- c(files, f2)

  # symptom model parameters
  sym_vars <- c("eta_sym_intercept", "eta_sym_pfu",
                "eta_sym_rna", "sigma_sym")
  p3 <- pp_plot |>
    tidyr::pivot_longer(dplyr::all_of(sym_vars)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::facet_wrap(~ name, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Prior: symptom hazard parameters", model_label))

  f3 <- file.path(out_dir, "prior_sym.pdf")
  ggplot2::ggsave(f3, p3, width = 8, height = 5)
  files <- c(files, f3)

  # prior predictive trajectories (hex-binned)
  traj <- tibble::tibble(
    time = rep(stan_data$time, ncol(pp$rna)),
    rna  = as.vector(pp$rna),
    pfu  = as.vector(pp$pfu),
    id   = rep(stan_data$id, ncol(pp$rna))
  ) |>
    tidyr::pivot_longer(c(rna, pfu), names_to = "outcome", values_to = "value")

  p4 <- ggplot2::ggplot(traj, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_hex(bins = 50) +
    ggplot2::facet_wrap(~ outcome, scales = "free_y") +
    viridis::scale_fill_viridis(option = "C", trans = "log10") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Prior predictive trajectories", model_label),
                  x = "Days since peak", y = "Viral load (log scale)")

  f4 <- file.path(out_dir, "prior_trajectories.pdf")
  ggplot2::ggsave(f4, p4, width = 10, height = 5)
  files <- c(files, f4)

  files
}


#' Plot target-cell-limited ODE examples (for presentation/manuscript)
#'
#' @param out_dir  Output directory
#' @return Character vector of saved file paths
plot_ode_examples <- function(out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # single-compartment target cell model
  targetcell_ode <- function(t, states, params) {
    T <- states[1]; I <- states[2]; V <- states[3]
    with(as.list(params), {
      dT <- -beta * T * V
      dI <- beta * T * V - delta * I
      dV <- pi_param * I - c_param * V
      list(c(dT, dI, dV))
    })
  }

  params <- c(beta = 2.4e-5, delta = 2.0, pi_param = 1700, c_param = 10)
  states <- c(T = 4e8, I = 0, V = 0.4)
  times  <- seq(0, 20, by = 0.01)

  out <- deSolve::ode(y = states, times = times, func = targetcell_ode,
                       parms = params)
  df <- data.frame(out)
  df$logV <- log(df$V + 1)

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = logV)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time (days)", y = "log V",
                  title = "Target cell limited model")

  f1 <- file.path(out_dir, "tcl_example_1.pdf")
  ggplot2::ggsave(f1, p1, width = 6, height = 4)

  # with refractory compartment
  targetcell2 <- function(t, states, params) {
    T <- states[1]; R <- states[2]; I <- states[3]; V <- states[4]
    with(as.list(params), {
      dT <- -beta * T * V - phi * T + rho * R
      dR <- phi * T - rho * R
      dI <- beta * T * V - delta * I
      dV <- pi_param * I - c_param * V
      list(c(dT, dR, dI, dV))
    })
  }

  params2 <- c(beta = 2.4e-5, delta = 2.0, pi_param = 1700, c_param = 10,
               phi = 1, rho = 0.5)
  states2 <- c(T = 4e8, R = 0, I = 0, V = 0.4)

  out2 <- deSolve::ode(y = states2, times = times, func = targetcell2,
                        parms = params2)
  df2 <- data.frame(out2)
  df2$logV <- log(df2$V + 1)

  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = time, y = logV)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time (days)", y = "log V",
                  title = "Target cell limited model with refractory compartment")

  f2 <- file.path(out_dir, "tcl_example_2.pdf")
  ggplot2::ggsave(f2, p2, width = 6, height = 4)

  c(f1, f2)
}
