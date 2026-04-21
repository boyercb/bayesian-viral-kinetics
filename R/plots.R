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
    zeta_sym_intercept     = pp$zeta_sym_intercept,
    zeta_sym_pfu           = pp$zeta_sym_pfu,
    zeta_sym_rna           = pp$zeta_sym_rna,
    zeta_sym_postpeak      = pp$zeta_sym_postpeak,
    zeta_sym_postpeak_rna  = pp$zeta_sym_postpeak_rna,
    sigma_sym              = pp$sigma_sym
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
  cols <- journal_colors()
  p1 <- pp_plot |>
    tidyr::pivot_longer(dplyr::all_of(trans_vars)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30, fill = cols[["muted"]], color = "white") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    theme_journal() +
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
    ggplot2::geom_histogram(bins = 30, fill = cols[["muted"]], color = "white") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    theme_journal() +
    ggplot2::labs(title = paste("Prior: population kinetics parameters", model_label))

  f2 <- file.path(out_dir, "prior_pe.pdf")
  ggplot2::ggsave(f2, p2, width = 10, height = 6)
  files <- c(files, f2)

  # symptom model parameters
  sym_vars <- c("zeta_sym_intercept", "zeta_sym_pfu",
                "zeta_sym_rna", "zeta_sym_postpeak",
                "zeta_sym_postpeak_rna", "sigma_sym")
  p3 <- pp_plot |>
    tidyr::pivot_longer(dplyr::all_of(sym_vars)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30, fill = cols[["sym"]], color = "white") +
    ggplot2::facet_wrap(~ name, scales = "free") +
    theme_journal() +
    ggplot2::labs(title = paste("Prior: symptom hazard parameters", model_label))

  f3 <- file.path(out_dir, "prior_sym.pdf")
  ggplot2::ggsave(f3, p3, width = 8, height = 5)
  files <- c(files, f3)

  # prior predictive trajectories (2D kernel density)
  # Filter to visible range BEFORE KDE so density isn't diluted by extreme tails
  cols <- journal_colors()
  set.seed(42)

  # Use rna_hat/pfu_hat (latent trajectories) instead of rna/pfu
  # to avoid LOD censoring artefacts that create density spikes at LOD
  traj_full <- tibble::tibble(
    time = rep(stan_data$time, ncol(pp$rna_hat)),
    rna  = as.vector(pp$rna_hat),
    pfu  = as.vector(pp$pfu_hat)
  ) |>
    tidyr::pivot_longer(c(rna, pfu), names_to = "outcome", values_to = "value") |>
    dplyr::filter(value >= -1, value <= 22,
                  time >= -12, time <= 22) |>
    dplyr::mutate(
      outcome = factor(
        outcome,
        levels = c("rna", "pfu"),
        labels = c("Viral RNA", "Infectious virus (PFU)")
      )
    )

  # Subsample to ~500k points per outcome for tractable KDE
  n_per_outcome <- 500000L
  traj_rna <- traj_full |> dplyr::filter(outcome == "Viral RNA")
  traj_pfu <- traj_full |> dplyr::filter(outcome == "Infectious virus (PFU)")
  traj_rna <- traj_rna[sample.int(nrow(traj_rna), min(n_per_outcome, nrow(traj_rna))), ]
  traj_pfu <- traj_pfu[sample.int(nrow(traj_pfu), min(n_per_outcome, nrow(traj_pfu))), ]

  # Helper to build a contour panel with independent bin count
  make_panel <- function(dat, bins, title, show_ylab = TRUE) {
    n_fill <- bins  # number of fill bands = bins
    pal <- c("black", colorRampPalette(
      c("#0D0887", "#3B049A", "#7301A8", "#A62098",
        "#CC4678", "#E56B5D", "#F89441", "#FDC328", "#F0F921")
    )(n_fill))
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = time, y = value)) +
      ggplot2::stat_density_2d_filled(
        ggplot2::aes(fill = ggplot2::after_stat(level)),
        contour_var = "density", n = 200, bins = bins
      ) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::coord_cartesian(ylim = c(0, 20), xlim = c(-10, 20)) +
      theme_journal() +
      ggplot2::labs(
        title = title,
        x = "Days from peak",
        y = if (show_ylab) "log copies/mL" else NULL
      )
    if (!show_ylab) {
      p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank())
    }
    p
  }

  p4_rna <- make_panel(traj_rna, bins = 20, "Viral RNA", show_ylab = TRUE)
  p4_pfu <- make_panel(traj_pfu, bins = 16, "Infectious virus (PFU)", show_ylab = FALSE)
  p4 <- p4_rna + p4_pfu + patchwork::plot_layout(ncol = 2)

  f4 <- file.path(out_dir, "prior_trajectories.pdf")
  dims <- journal_dims("full")
  ggplot2::ggsave(f4, p4, width = dims$width, height = dims$height)
  files <- c(files, f4)

  files
}


#' Plot target-cell-limited ODE examples (for presentation/manuscript)
#'
#' Generates two figures comparing the target cell limited ODE solution
#' to a least-squares fitted piecewise exponential (smooth) approximation.
#' Example 1: standard TCL model. Example 2: TCL with refractory compartment.
#'
#' @param out_dir  Output directory
#' @return Character vector of saved file paths
plot_ode_examples <- function(out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cols <- journal_colors()

  # ── Helper: fit smfun to an ODE trajectory ──────────────────────────────
  fit_smfun <- function(time, logV) {
    # Subset to the detectable range for fitting
    detectable <- logV > 0.1
    t_fit <- time[detectable]
    y_fit <- logV[detectable]

    # Initial guesses from the ODE trajectory
    peak_idx <- which.max(y_fit)
    tp0  <- t_fit[peak_idx]
    dp0  <- y_fit[peak_idx]
    wp0  <- tp0 - t_fit[1]
    wr0  <- t_fit[length(t_fit)] - tp0

    # Optimize on log scale for positivity
    obj <- function(par) {
      tp <- par[1]
      dp <- exp(par[2])
      wp <- exp(par[3])
      wr <- exp(par[4])
      pred <- smfun(t_fit, tp, wp, wr, dp, wf = 0)
      sum((y_fit - pred)^2)
    }

    res <- optim(
      par = c(tp0, log(dp0), log(wp0), log(wr0)),
      fn  = obj,
      method = "Nelder-Mead",
      control = list(maxit = 5000)
    )

    list(tp = res$par[1], dp = exp(res$par[2]),
         wp = exp(res$par[3]), wr = exp(res$par[4]))
  }

  # ── Helper: build a comparison plot ─────────────────────────────────────
  make_tcl_plot <- function(df, title_text) {
    fit <- fit_smfun(df$time, df$logV)
    df$pe_fit <- smfun(df$time, fit$tp, fit$wp, fit$wr, fit$dp, wf = 0)
    # Zero out below-zero predictions
    df$pe_fit <- pmax(df$pe_fit, 0)

    ggplot2::ggplot(df, ggplot2::aes(x = time)) +
      ggplot2::geom_line(
        ggplot2::aes(y = logV, color = "ODE solution", linetype = "ODE solution"),
        linewidth = 0.8
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = pe_fit, color = "Piecewise exponential", linetype = "Piecewise exponential"),
        linewidth = 0.8
      ) +
      ggplot2::scale_color_manual(
        values = c("ODE solution" = cols[["rna"]], "Piecewise exponential" = cols[["pfu"]]),
        name = NULL
      ) +
      ggplot2::scale_linetype_manual(
        values = c("ODE solution" = "solid", "Piecewise exponential" = "dashed"),
        name = NULL
      ) +
      theme_journal() +
      ggplot2::theme(
        legend.position = "bottom",
        legend.key.width = ggplot2::unit(1.5, "cm")
      ) +
      ggplot2::labs(x = "Time (days)", y = "log V", title = title_text)
  }

  # ── Example 1: standard target cell limited model ───────────────────────
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

  p1 <- make_tcl_plot(df, "Target cell limited model")
  f1 <- file.path(out_dir, "tcl_example_1.pdf")
  ggplot2::ggsave(f1, p1, width = 6, height = 4)

  # ── Example 2: TCL with refractory compartment ─────────────────────────
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

  p2 <- make_tcl_plot(df2, "Target cell limited model with refractory compartment")
  f2 <- file.path(out_dir, "tcl_example_2.pdf")
  ggplot2::ggsave(f2, p2, width = 6, height = 4)

  c(f1, f2)
}

#' Plot antigen schematic: trajectory, derivative, convolution, and log-affine
#'
#' Generates a 4-panel supplementary figure illustrating the relationship
#' between the smooth trajectory, its derivative, the antigen convolution,
#' and the log-affine approximation. Uses representative parameters (no fit
#' dependency).
#'
#' @param style  Journal style for theming (default: current journal)
#' @return A patchwork composite ggplot object
plot_antigen_schematic <- function(style = NULL) {
  library(ggplot2)
  library(patchwork)

  cols <- journal_colors(style)

  # Representative parameters
  dp <- 10; wp <- 3; wr <- 7; tp <- 5; wf <- 0
  t_grid <- seq(0, 20, length.out = 500)

  # (a) Smooth trajectory
  g <- smfun(t_grid, tp, wp, wr, dp, wf)

  # (b) Derivative
  dg <- smfun_deriv(t_grid, tp, wp, wr, dp, wf)


  # (c) Antigen via numerical convolution vs V_t
  V <- exp(g)
  dt_step <- t_grid[2] - t_grid[1]
  alpha1 <- 1; alpha2 <- 0.5  # representative clearance rate
  A <- numeric(length(t_grid))
  for (j in seq_along(t_grid)) {
    # Discrete convolution: A_t = alpha1 * sum_{k<=j} exp(-alpha2*(t_j - t_k)) * V_k * dt
    A[j] <- alpha1 * sum(exp(-alpha2 * (t_grid[j] - t_grid[1:j])) * V[1:j]) * dt_step
  }
  # Normalize to unit peak for overlay
  V_norm <- V / max(V)
  A_norm <- A / max(A)

  # (d) Log-affine approximation vs convolution on log scale
  # The log-affine model produces a NEW piecewise exponential with transformed

  # parameters (different peak, prolif rate, clearance rate), not just a vertical
  # shift. Fit optimal log-affine parameters to the convolution.
  log_A_conv <- log(pmax(A, 1e-10))
  # Fit smfun() to log-convolution by least-squares optimization
  obj_fn <- function(par) {
    tp_a <- par[1]; dp_a <- par[2]; wp_a <- par[3]; wr_a <- par[4]
    pred <- smfun(t_grid, tp_a, wp_a, wr_a, dp_a, wf)
    sum((pred - log_A_conv)^2)
  }
  # Initialize: peak slightly delayed, clearance slower (= min(alpha2, dp/wr))
  fit_la <- optim(c(tp + 0.5, dp - 1, wp, wr * 2), obj_fn, method = "Nelder-Mead")
  log_A_logaffine <- smfun(t_grid, fit_la$par[1], fit_la$par[3],
                            fit_la$par[4], fit_la$par[2], wf)

  # Derivative model: gamma_1 + gamma_2 * g(t) + gamma_3 * dg/dt, fit to convolution
  # This matches the LFD observation model structure
  X_deriv <- cbind(1, g, dg)
  fit_deriv <- lm.fit(X_deriv, log_A_conv)
  log_A_deriv <- X_deriv %*% fit_deriv$coefficients

  # Indicator model: gamma_1 + gamma_2 * g(t) + gamma_3 * I(t >= tp) + gamma_4 * I(t >= tp) * g(t)
  # Allows different intercept AND slope pre vs post peak
  post_peak <- as.numeric(t_grid >= tp)
  X_ind <- cbind(1, g, post_peak, post_peak * g)
  fit_ind <- lm.fit(X_ind, log_A_conv)
  log_A_ind <- X_ind %*% fit_ind$coefficients

  # Sigmoid model: gamma_1 + gamma_2 * g(t) + gamma_3 * sigma_k(t - tp) + gamma_4 * sigma_k(t - tp) * g(t)
  # Smooth sigmoid approximation to the indicator (kappa = 5)
  kappa <- 5
  sig_pp <- soft_postpeak(t_grid, tp, kappa)
  X_sig <- cbind(1, g, sig_pp, sig_pp * g)
  fit_sig <- lm.fit(X_sig, log_A_conv)
  log_A_sig <- X_sig %*% fit_sig$coefficients

  # Build data frames
  df_a <- data.frame(t = t_grid, g = g)
  df_b <- data.frame(t = t_grid, dg = dg,
                      phase = ifelse(dg >= 0, "Proliferation", "Clearance"))
  df_c <- data.frame(
    t = rep(t_grid, 2),
    value = c(V_norm, A_norm),
    curve = rep(c("Viral load", "Antigen"), each = length(t_grid))
  )
  df_d <- data.frame(
    t = rep(t_grid, 5),
    value = c(log_A_conv, log_A_logaffine, log_A_deriv, log_A_ind, log_A_sig),
    curve = rep(c("Convolution", "Log-affine", "Derivative", "Indicator", "Sigmoid"),
                each = length(t_grid))
  )

  # Panel (a): trajectory
  pa <- ggplot(df_a, aes(t, g)) +
    geom_line(linewidth = 0.8, color = unname(cols["rna"])) +
    labs(x = "Days since infection", y = expression(g[s](t) ~ "(log copies/mL)"),
         title = "(a) Smooth trajectory") +
    theme_journal(style)

  # Panel (b): derivative with phase shading
  df_b_pos <- df_b[df_b$dg >= 0, ]
  df_b_neg <- df_b[df_b$dg < 0, ]
  pb <- ggplot(df_b, aes(t, dg)) +
    geom_area(data = df_b_pos, aes(t, dg), fill = unname(cols["rna"]),
              alpha = 0.25) +
    geom_area(data = df_b_neg, aes(t, dg), fill = unname(cols["pfu"]),
              alpha = 0.25) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(x = "Days since infection",
         y = expression(dot(g)[s](t) ~ "(log copies/mL/day)"),
         title = "(b) Trajectory derivative") +
    annotate("text", x = tp - wp/2, y = max(dg) * 0.7, label = "Proliferation",
             size = 3, hjust = 0.5) +
    annotate("text", x = tp + wr/2, y = min(dg) * 0.7, label = "Clearance",
             size = 3, hjust = 0.5) +
    theme_journal(style)

  # Panel (c): antigen vs virus (normalized)
  pc <- ggplot(df_c, aes(t, value, color = curve, linetype = curve)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Viral load" = unname(cols["pfu"]),
                                   "Antigen" = unname(cols["lfd"])),
                        name = NULL) +
    scale_linetype_manual(values = c("Viral load" = "dashed",
                                      "Antigen" = "solid"),
                           name = NULL) +
    labs(x = "Days since infection", y = "Normalized concentration",
         title = "(c) Antigen lag") +
    theme_journal(style) +
    theme(legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0))

  # Panel (d): all approximations vs convolution
  pd <- ggplot(df_d, aes(t, value, color = curve, linetype = curve)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Convolution" = "black",
                                   "Log-affine" = unname(cols["rna"]),
                                   "Derivative" = unname(cols["pfu"]),
                                   "Indicator" = unname(cols["lfd"]),
                                   "Sigmoid" = unname(cols["accent"])),
                        name = NULL) +
    scale_linetype_manual(values = c("Convolution" = "solid",
                                      "Log-affine" = "dashed",
                                      "Derivative" = "dotted",
                                      "Indicator" = "twodash",
                                      "Sigmoid" = "longdash"),
                           name = NULL) +
    labs(x = "Days since infection", y = "log(Antigen)",
         title = "(d) Approximations") +
    theme_journal(style) +
    theme(legend.position = "bottom",
          legend.margin = margin(0, 0, 0, 0))

  # Combine
  (pa | pb) / (pc | pd)
}


#' Save antigen schematic in one or more journal styles
#'
#' @param styles Character vector of journal styles
#' @param out_dir Output directory for figures
#' @return Character vector of saved file paths
save_antigen_schematic <- function(styles = c("pnas"),
                                   out_dir = "output/figures") {
  all_paths <- character(0)

  for (s in styles) {
    set_journal(s)
    p <- plot_antigen_schematic(style = s)
    paths <- save_journal_figure(
      p,
      name = "antigen_schematic",
      layout = "full",
      width = 12,
      height = 7,
      out_dir = out_dir,
      style = s
    )
    all_paths <- c(all_paths, paths)
  }

  unname(all_paths)
}
