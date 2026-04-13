# ──────────────────────────────────────────────────────────────────────────────
# predictions.R — Posterior prediction computation and trajectory plotting
# ──────────────────────────────────────────────────────────────────────────────

#' Sample observation IDs for trajectory plots
#'
#' @param stacked_dat   Stacked dataset
#' @param n_nba_sample  Number of NBA trajectories to sample (all others shown)
#' @param seed          Random seed for reproducibility
#' @return Filtered data frame for prediction
make_prediction_data <- function(stacked_dat, n_nba_sample = 42,
                                 seed = 7324650) {
  set.seed(seed)

  nba_smp  <- sample(
    unique(stacked_dat$id[stacked_dat$source == 1]),
    n_nba_sample
  )
  ata_smp  <- unique(stacked_dat$id[stacked_dat$source == 2])
  uiuc_smp <- unique(stacked_dat$id[stacked_dat$source == 3])
  hct_smp  <- unique(stacked_dat$id[stacked_dat$source == 4])
  legacy_smp <- unique(stacked_dat$id[stacked_dat$source == 5])

  dplyr::filter(
    stacked_dat,
    id %in% c(nba_smp, ata_smp, uiuc_smp, hct_smp, legacy_smp)
  )
}


#' Build a dense time grid per individual for smooth trajectory plotting
#'
#' Creates a fine grid of evaluation points for each individual in pred_dat.
#' Covariates are carried over from the first observation row per individual.
#'
#' @param pred_dat  Prediction data frame (from make_prediction_data)
#' @param dt        Time step for the dense grid (days)
#' @param pad       Padding beyond observed time range (days)
#' @return Data frame with columns: id, time, source, + covariate columns
make_dense_grid <- function(pred_dat, dt = 0.25, pad = 1.0) {
  covariates <- define_covariates()

  # One metadata row per individual (first observation carries covariates)
  meta <- pred_dat |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      source   = dplyr::first(source),
      t_min    = min(time) - pad,
      t_max    = max(time) + pad,
      dplyr::across(dplyr::all_of(covariates$x_vars), dplyr::first),
      .groups  = "drop"
    )

  # Build dense grid per individual
  grid_list <- lapply(seq_len(nrow(meta)), function(i) {
    row <- meta[i, ]
    tseq <- seq(row$t_min, row$t_max, by = dt)
    out <- data.frame(id = row$id, time = tseq, source = row$source)
    for (v in covariates$x_vars) out[[v]] <- row[[v]]
    out
  })

  dplyr::bind_rows(grid_list)
}


#' Compute posterior-predicted trajectories on a dense grid and observation points
#'
#' Runs predict_kinetics on a dense time grid (for smooth plotting lines) and
#' summarises the result into posterior mean/quantile columns. Observation-level
#' predictions (LFD, symptoms) are computed on the original pred_dat rows.
#'
#' @param fit       CmdStanMCMC fit object
#' @param pred_dat  Prediction data frame (from make_prediction_data)
#' @param stan_data Stan data list (from build_stan_data) — needed for flags
#' @param dt        Time resolution for the dense grid (days)
#' @return Named list with two elements:
#'   - obs: pred_dat augmented with observation-level predictions (lfd_hat, sym_hat, rna_hat, pfu_hat)
#'   - grid: dense grid with smooth RNA/PFU predictions and quantiles
compute_predictions <- function(fit, pred_dat, stan_data, dt = 0.25,
                                max_draws = 1000) {

  # ── 1. Observation-level predictions (LFD, symptoms, RNA/PFU) ───────────
  #    Uses rvar-based predict_kinetics on the (small) observation set
  message("compute_predictions: observation-level predictions...")
  p_obs <- predict_kinetics(fit, newdata = pred_dat, stan_data = stan_data,
                            max_draws = max_draws)

  pred_dat$rna_hat     <- posterior::E(p_obs$rna_hat)
  pred_dat$rna_hat_q1  <- posterior::quantile2(p_obs$rna_hat, 0.025)
  pred_dat$rna_hat_q99 <- posterior::quantile2(p_obs$rna_hat, 0.975)

  pred_dat$pfu_hat     <- posterior::E(p_obs$pfu_hat)
  pred_dat$pfu_hat_q1  <- posterior::quantile2(p_obs$pfu_hat, 0.025)
  pred_dat$pfu_hat_q99 <- posterior::quantile2(p_obs$pfu_hat, 0.975)

  pred_dat$lfd_hat <- posterior::E(p_obs$lfd_hat)
  pred_dat$sym_hat <- posterior::E(p_obs$sym_hat)
  rm(p_obs); gc()

  # ── 2. Dense grid predictions (smooth lines for RNA, PFU) ───────────────
  #    Computed per-individual to limit memory: extract draws once, then loop.
  message("compute_predictions: dense grid (dt=", dt, ")...")
  grid <- compute_grid_summaries(fit, pred_dat, stan_data, dt = dt,
                                  max_draws = max_draws)
  message("  grid rows: ", nrow(grid))

  list(obs = pred_dat, grid = grid)
}


#' Compute dense-grid trajectory summaries per individual (memory-efficient)
#'
#' Extracts posterior draws of kinetic parameters once, then loops over
#' individuals.  For each individual, evaluates the trajectory function on a
#' dense time grid for every posterior draw, summarises to mean / quantiles,
#' and discards the per-draw matrix immediately.  Peak memory is proportional
#' to (max_draws x grid_points_per_individual) rather than (max_draws x total_grid).
#'
#' @param fit       CmdStanMCMC fit object
#' @param pred_dat  Observation-level prediction data (to determine time ranges)
#' @param stan_data Stan data list (flags, priors)
#' @param dt        Grid spacing (days)
#' @param max_draws Max posterior draws to use
#' @param pad       Day padding beyond observed range
#' @return Data frame with id, time, source, rna_hat, rna_hat_q1, rna_hat_q99,
#'         pfu_hat, pfu_hat_q1, pfu_hat_q99
compute_grid_summaries <- function(fit, pred_dat, stan_data, dt = 0.25,
                                    max_draws = 1000, pad = 1.0) {

  covariates <- define_covariates()
  use_smooth <- as.logical(stan_data$use_smooth)

  # ── Extract and thin draws once ─────────────────────────────────────────
  var_list <- c(
    "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
    "tau_tp", "tau_dp", "tau_wp", "tau_wr",
    "tp_i_pfu"
  )
  if (isTRUE(stan_data$ind_corr == 1)) {
    var_list <- c(var_list, "z_ind_rna", "L_Omega_rna", "sigma_ind_rna")
  } else {
    var_list <- c(var_list, "tp_i_rna")
    if (stan_data$ind_effects) {
      var_list <- c(var_list, "dp_i_rna", "wp_i_rna", "wr_i_rna")
    }
  }
  if (stan_data$ind_effects) {
    var_list <- c(var_list, "dp_i_pfu", "wp_i_pfu", "wr_i_pfu")
  }
  if (stan_data$adj_rna) {
    var_list <- c(var_list, "beta_dp_rna", "beta_wp_rna", "beta_wr_rna")
  }
  if (stan_data$adj_pfu) {
    var_list <- c(var_list, "beta_dp_pfu", "beta_wp_pfu", "beta_wr_pfu")
  }
  if (stan_data$source_rna) {
    var_list <- c(var_list, "tp_k_rna", "dp_k_rna", "wp_k_rna", "wr_k_rna")
  }
  if (stan_data$source_pfu) {
    var_list <- c(var_list, "tp_k_pfu", "dp_k_pfu", "wp_k_pfu", "wr_k_pfu")
  }
  if (stan_data$use_wf) {
    var_list <- c(var_list, "wf_raw")
    if (stan_data$ind_effects) var_list <- c(var_list, "wf_i")
    if (stan_data$source_rna)  var_list <- c(var_list, "wf_k")
  }

  drws <- posterior::as_draws_matrix(fit$draws(variables = var_list))
  nd <- nrow(drws)
  if (!is.null(max_draws) && nd > max_draws) {
    idx <- round(seq(1, nd, length.out = max_draws))
    drws <- drws[idx, ]
    nd <- max_draws
  }

  # Helper to extract a parameter column as a numeric vector [nd]
  get_vec <- function(nm) as.numeric(drws[, nm])

  # Population means: [nd] vectors
  dp_mean <- get_vec("dp_mean_rna")
  wp_mean <- get_vec("wp_mean_rna")
  wr_mean <- get_vec("wr_mean_rna")

  # Log-affine transformation params
  tau_tp1 <- get_vec("tau_tp[1]"); tau_tp2 <- get_vec("tau_tp[2]")
  tau_dp1 <- get_vec("tau_dp[1]"); tau_dp2 <- get_vec("tau_dp[2]")
  tau_wp1 <- get_vec("tau_wp[1]"); tau_wp2 <- get_vec("tau_wp[2]")
  tau_wr1 <- get_vec("tau_wr[1]"); tau_wr2 <- get_vec("tau_wr[2]")

  # ── Pre-compute correlated RNA individual effects ───────────────────────
  N_ind <- sum(stan_data$M)
  if (isTRUE(stan_data$ind_corr == 1)) {
    k_rv <- posterior::as_draws_rvars(drws)
    sigma_d <- posterior::draws_of(k_rv$sigma_ind_rna)
    L_d     <- posterior::draws_of(k_rv$L_Omega_rna)
    z_d     <- posterior::draws_of(k_rv$z_ind_rna)
    eta_rna <- array(NA_real_, dim = c(nd, 4, N_ind))
    for (d in seq_len(nd)) {
      eta_rna[d, , ] <- (diag(sigma_d[d, ]) %*% L_d[d, , ]) %*% z_d[d, , ]
    }
    rm(k_rv, sigma_d, L_d, z_d)
  }

  # ── Flat-top base draws ─────────────────────────────────────────────────
  if (stan_data$use_wf) {
    wf_base <- stan_data$prior_wf_mean *
      exp(stan_data$prior_wf_cv * get_vec("wf_raw[1]"))
  }

  # Covariate betas
  P <- stan_data$P
  if (stan_data$adj_rna) {
    beta_dp_rna <- beta_wp_rna <- beta_wr_rna <- matrix(NA, nd, P)
    for (p in seq_len(P)) {
      beta_dp_rna[, p] <- get_vec(paste0("beta_dp_rna[", p, "]"))
      beta_wp_rna[, p] <- get_vec(paste0("beta_wp_rna[", p, "]"))
      beta_wr_rna[, p] <- get_vec(paste0("beta_wr_rna[", p, "]"))
    }
  }
  if (stan_data$adj_pfu) {
    beta_dp_pfu <- beta_wp_pfu <- beta_wr_pfu <- matrix(NA, nd, P)
    for (p in seq_len(P)) {
      beta_dp_pfu[, p] <- get_vec(paste0("beta_dp_pfu[", p, "]"))
      beta_wp_pfu[, p] <- get_vec(paste0("beta_wp_pfu[", p, "]"))
      beta_wr_pfu[, p] <- get_vec(paste0("beta_wr_pfu[", p, "]"))
    }
  }

  # ── Per-individual metadata ─────────────────────────────────────────────
  ind_meta <- pred_dat |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      source  = dplyr::first(source),
      t_min   = min(time) - pad,
      t_max   = max(time) + pad,
      dplyr::across(dplyr::all_of(covariates$x_vars), dplyr::first),
      .groups = "drop"
    )

  # Choose trajectory function
  tfun <- if (use_smooth) smfun else pefun

  # ── Loop over individuals ───────────────────────────────────────────────
  result_list <- vector("list", nrow(ind_meta))

  for (row_i in seq_len(nrow(ind_meta))) {
    ind    <- ind_meta$id[row_i]
    src    <- as.integer(ind_meta$source[row_i])
    tseq   <- seq(ind_meta$t_min[row_i], ind_meta$t_max[row_i], by = dt)
    n_t    <- length(tseq)
    x_ind  <- as.numeric(ind_meta[row_i, covariates$x_vars])

    # ── RNA parameters for this individual [nd vectors] ───────────────
    if (isTRUE(stan_data$ind_corr == 1)) {
      tp_rna <- eta_rna[, 1, ind]
      dp_rna <- dp_mean * exp(eta_rna[, 2, ind])
      wp_rna <- wp_mean * exp(eta_rna[, 3, ind])
      wr_rna <- wr_mean * exp(eta_rna[, 4, ind])
    } else {
      tp_rna <- get_vec(paste0("tp_i_rna[", ind, "]"))
      if (stan_data$ind_effects) {
        dp_rna <- dp_mean * exp(get_vec(paste0("dp_i_rna[", ind, "]")))
        wp_rna <- wp_mean * exp(get_vec(paste0("wp_i_rna[", ind, "]")))
        wr_rna <- wr_mean * exp(get_vec(paste0("wr_i_rna[", ind, "]")))
      } else {
        dp_rna <- dp_mean
        wp_rna <- wp_mean
        wr_rna <- wr_mean
      }
    }

    if (stan_data$source_rna) {
      tp_rna <- tp_rna + get_vec(paste0("tp_k_rna[", src, "]"))
      dp_rna <- dp_rna * exp(get_vec(paste0("dp_k_rna[", src, "]")))
      wp_rna <- wp_rna * exp(get_vec(paste0("wp_k_rna[", src, "]")))
      wr_rna <- wr_rna * exp(get_vec(paste0("wr_k_rna[", src, "]")))
    }

    if (stan_data$adj_rna) {
      xb_dp <- as.numeric(beta_dp_rna %*% x_ind)
      xb_wp <- as.numeric(beta_wp_rna %*% x_ind)
      xb_wr <- as.numeric(beta_wr_rna %*% x_ind)
      dp_rna <- dp_rna * exp(xb_dp)
      wp_rna <- wp_rna * exp(xb_wp)
      wr_rna <- wr_rna * exp(xb_wr)
    }

    if (stan_data$use_wf) {
      wf_d <- wf_base
      if (stan_data$ind_effects) {
        wf_d <- wf_d * exp(get_vec(paste0("wf_i[", ind, "]")))
      }
      if (stan_data$source_rna) {
        wf_d <- wf_d * exp(get_vec(paste0("wf_k[", src, "]")))
      }
    } else {
      wf_d <- rep(0, nd)
    }

    # ── Evaluate RNA trajectory: [nd x n_t] matrix ───────────────────
    rna_mat <- matrix(NA_real_, nd, n_t)
    for (d in seq_len(nd)) {
      rna_mat[d, ] <- tfun(tseq, tp_rna[d], wp_rna[d], wr_rna[d],
                            dp_rna[d], wf_d[d])
    }

    # ── PFU parameters (log-affine from RNA) ──────────────────────────
    dp_pfu <- exp(tau_dp1) * dp_rna^tau_dp2
    wp_pfu <- exp(tau_wp1) * wp_rna^tau_wp2
    wr_pfu <- exp(tau_wr1) * wr_rna^tau_wr2
    tp_pfu <- tau_tp1 + tau_tp2 * tp_rna

    # PFU individual effects: only for PFU-informed individuals
    pidx <- stan_data$pfu_ind_idx[ind]  # 0 = no RE, >0 = index
    if (pidx > 0) {
      tp_pfu <- tp_pfu + get_vec(paste0("tp_i_pfu[", pidx, "]"))
    }

    if (stan_data$ind_effects && pidx > 0) {
      dp_pfu <- dp_pfu * exp(get_vec(paste0("dp_i_pfu[", pidx, "]")))
      wp_pfu <- wp_pfu * exp(get_vec(paste0("wp_i_pfu[", pidx, "]")))
      wr_pfu <- wr_pfu * exp(get_vec(paste0("wr_i_pfu[", pidx, "]")))
    }

    if (stan_data$source_pfu) {
      tp_pfu <- tp_pfu + get_vec(paste0("tp_k_pfu[", src, "]"))
      dp_pfu <- dp_pfu * exp(get_vec(paste0("dp_k_pfu[", src, "]")))
      wp_pfu <- wp_pfu * exp(get_vec(paste0("wp_k_pfu[", src, "]")))
      wr_pfu <- wr_pfu * exp(get_vec(paste0("wr_k_pfu[", src, "]")))
    }

    if (stan_data$adj_pfu) {
      xb_dp <- as.numeric(beta_dp_pfu %*% x_ind)
      xb_wp <- as.numeric(beta_wp_pfu %*% x_ind)
      xb_wr <- as.numeric(beta_wr_pfu %*% x_ind)
      dp_pfu <- dp_pfu * exp(xb_dp)
      wp_pfu <- wp_pfu * exp(xb_wp)
      wr_pfu <- wr_pfu * exp(xb_wr)
    }

    # ── Evaluate PFU trajectory: [nd x n_t] matrix ───────────────────
    pfu_mat <- matrix(NA_real_, nd, n_t)
    for (d in seq_len(nd)) {
      pfu_mat[d, ] <- tfun(tseq, tp_pfu[d], wp_pfu[d], wr_pfu[d],
                            dp_pfu[d], wf_d[d])
    }

    # ── Summarise immediately ─────────────────────────────────────────
    result_list[[row_i]] <- data.frame(
      id          = ind,
      time        = tseq,
      source      = ind_meta$source[row_i],
      rna_hat     = colMeans(rna_mat),
      rna_hat_q1  = apply(rna_mat, 2, quantile, probs = 0.025),
      rna_hat_q99 = apply(rna_mat, 2, quantile, probs = 0.975),
      pfu_hat     = colMeans(pfu_mat),
      pfu_hat_q1  = apply(pfu_mat, 2, quantile, probs = 0.025),
      pfu_hat_q99 = apply(pfu_mat, 2, quantile, probs = 0.975)
    )

    if (row_i %% 50 == 0) {
      message(sprintf("  %d / %d individuals done", row_i, nrow(ind_meta)))
    }
  }

  dplyr::bind_rows(result_list)
}

#' Plot fitted trajectories for a single source
#'
#' Uses the dense grid for smooth RNA/PFU lines and ribbons, and the
#' observation-level data for data points, LFD tiles, and symptom markers.
#' Styling matches panel_individual() from pub_figures.R (Figure 2).
#'
#' @param predictions  List with elements `obs` and `grid` (from compute_predictions)
#' @param stan_data    Stan data list (for LODs and feature flags)
#' @param source_id    Numeric source ID (1=NBA, 2=ATACCC, 3=UIUC, 4=HCT, 5=Legacy)
#' @param style        Journal style for theme/colors (default "pnas")
#' @return ggplot object
plot_source_trajectories <- function(predictions, stan_data, source_id,
                                     style = "pnas") {

  obs  <- dplyr::filter(predictions$obs,  as.numeric(source) == source_id)
  grid <- dplyr::filter(predictions$grid, as.numeric(source) == source_id)
  if (nrow(obs) == 0) return(NULL)

  cols <- journal_colors(style)
  source_labels <- c("1" = "NBA", "2" = "ATACCC", "3" = "UIUC",
                      "4" = "HCT", "5" = "Legacy")

  lod_rna <- stan_data$lod_rna[source_id]
  lod_pfu <- stan_data$lod_pfu[source_id]

  has_pfu <- "pfu_exist" %in% names(obs) && any(obs$pfu_exist == 1)
  has_lfd <- "lfd_exist" %in% names(obs) && any(obs$lfd_exist == 1)
  has_sym <- "sym_exist" %in% names(obs) && any(obs$sym_exist == 1)

  # ---- clamp grid fitted values at LOD -----------------------------------
  clamp <- function(x, lod) replace(x, !is.na(x) & x < lod, lod)
  grid$rna_hat_c     <- clamp(grid$rna_hat,     lod_rna)
  grid$rna_hat_q1_c  <- clamp(grid$rna_hat_q1,  lod_rna)
  grid$rna_hat_q99_c <- clamp(grid$rna_hat_q99, lod_rna)

  if (has_pfu) {
    grid$pfu_hat_c     <- clamp(grid$pfu_hat,     lod_pfu)
    grid$pfu_hat_q1_c  <- clamp(grid$pfu_hat_q1,  lod_pfu)
    grid$pfu_hat_q99_c <- clamp(grid$pfu_hat_q99, lod_pfu)
  }

  # ---- y-axis limits for LFD strip placement -----------------------------
  y_vals <- c(obs$rna, obs$pfu, grid$rna_hat_q99)
  if (has_pfu) y_vals <- c(y_vals, grid$pfu_hat_q99)
  y_max <- max(y_vals, na.rm = TRUE) + 0.5

  # ---- base plot ---------------------------------------------------------
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = time, group = factor(id)))

  # ---- RNA ribbon + line + points ----------------------------------------
  rna_obs  <- dplyr::filter(obs, rna_exist == 1)
  rna_grid <- dplyr::filter(grid, !is.na(rna_hat))
  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rna_hat_q1_c, ymax = rna_hat_q99_c),
      data = rna_grid,
      alpha = 0.15, fill = cols[["rna"]]
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = rna_hat_c),
      data = rna_grid,
      color = cols[["rna"]], linewidth = 0.5
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = rna),
      data = rna_obs,
      color = cols[["rna"]], shape = 16, size = 1.5, alpha = 0.8
    )

  # ---- LOD line(s) for RNA -----------------------------------------------
  if (!is.na(lod_rna) && lod_rna > 0) {
    p <- p +
      ggplot2::geom_hline(yintercept = lod_rna, linetype = "dotted",
                           color = cols[["rna"]], alpha = 0.5,
                           linewidth = 0.3)
  }

  # ---- PFU ribbon + line + points ----------------------------------------
  if (has_pfu) {
    pfu_obs  <- dplyr::filter(obs, pfu_exist == 1)
    pfu_grid <- dplyr::filter(grid, !is.na(pfu_hat))
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = pfu_hat_q1_c, ymax = pfu_hat_q99_c),
        data = pfu_grid,
        alpha = 0.15, fill = cols[["pfu"]]
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = pfu_hat_c),
        data = pfu_grid,
        color = cols[["pfu"]], linewidth = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = pfu),
        data = pfu_obs,
        color = cols[["pfu"]], shape = 17, size = 1.5, alpha = 0.8
      )

    # LOD line for PFU
    if (!is.na(lod_pfu) && lod_pfu > 0) {
      p <- p +
        ggplot2::geom_hline(yintercept = lod_pfu, linetype = "dotted",
                             color = cols[["pfu"]], alpha = 0.5,
                             linewidth = 0.3)
    }
  }

  # ---- LFD tiles (matching Figure 2 style) -------------------------------
  if (has_lfd) {
    lfd_dat <- dplyr::filter(obs, lfd_exist == 1, !is.na(lfd_hat))
    if (nrow(lfd_dat) > 0) {
      lfd_dat$lfd_y <- y_max + 0.8
      p <- p +
        ggplot2::geom_tile(
          ggplot2::aes(y = lfd_y, fill = lfd_hat,
                       width = 0.8, height = 1.0),
          data = lfd_dat
        ) +
        ggplot2::geom_text(
          ggplot2::aes(y = lfd_y,
                       label = ifelse(lfd == 1, "+", "\u2013")),
          data = lfd_dat, size = 2, color = "white", fontface = "bold"
        ) +
        colorspace::scale_fill_continuous_sequential(
          name = "P(LFD+)", palette = "Greens 3",
          rev = FALSE, limits = c(0, 1), guide = "none"
        )
    }
  }

  # ---- Symptom onset (dashed vertical line, matching Figure 2) -----------
  if (has_sym) {
    sym_dat <- obs |>
      dplyr::filter(sym_exist == 1) |>
      dplyr::group_by(id) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::filter(!is.na(sym_onset), sym_ever == 1, sym_onset < 90)

    if (nrow(sym_dat) > 0) {
      p <- p +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = sym_onset),
          data = sym_dat,
          linetype = "dashed", color = cols[["sym"]], alpha = 0.6,
          linewidth = 0.3
        )
    }
  }

  # ---- facets & theme (matching Figure 2 journal style) ------------------
  p <- p +
    ggplot2::facet_wrap(~ id) +
    theme_journal(style, base_size = 7) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      strip.text       = ggplot2::element_text(size = ggplot2::rel(0.8))
    ) +
    ggplot2::labs(
      x = "Days from peak",
      y = "log copies/mL",
      title = source_labels[as.character(source_id)]
    ) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


#' Plot all trajectory figures and save to output/figures/
#'
#' One PDF per source, with RNA + PFU overlaid, LOD cutoffs, LFD heatmap,
#' and symptom-onset markers. Uses dense grid for smooth fitted lines.
#' Styling matches Figure 2 (panel_individual) from pub_figures.R.
#'
#' @param predictions  List with elements `obs` and `grid` (from compute_predictions)
#' @param stan_data    Stan data list (for LODs and feature flags)
#' @param out_dir      Output directory for figures
#' @param style        Journal style for theme/colors (default "pnas")
#' @return Character vector of saved file paths
plot_all_trajectories <- function(predictions, stan_data,
                                  out_dir = "output/figures",
                                  style = "pnas") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  source_names <- c("nba", "ataccc", "uiuc", "hct", "legacy")
  files <- c()

  for (i in seq_along(source_names)) {
    p <- plot_source_trajectories(predictions, stan_data, source_id = i,
                                  style = style)
    if (is.null(p)) next

    out_file <- file.path(out_dir, paste0(source_names[i], "_fit.pdf"))
    ggplot2::ggsave(out_file, p, width = 17, height = 10)
    files <- c(files, out_file)
  }

  files
}


#' Generate the inferred PFU figure (Figure 6)
#'
#' Produces a 3-row × 3-column figure comparing PFU inference quality
#' across data regimes: ATACCC (PFU observed), NBA (RNA only), Legacy
#' (RNA + symptoms). Requires fig_inferred_pfu() and panel_inferred_pfu()
#' from pub_figures.R.
#'
#' @param predictions  List with obs and grid (from compute_predictions)
#' @param stan_data    Stan data list (for LODs)
#' @param out_dir      Output directory
#' @param style        Journal style (default "pnas")
#' @return Character path to saved PDF (for targets file tracking)
plot_inferred_pfu <- function(predictions, stan_data,
                               out_dir = "output/figures",
                               style = "pnas") {
  dir.create(file.path(out_dir, style), recursive = TRUE,
             showWarnings = FALSE)

  p <- fig_inferred_pfu(predictions, stan_data, style = style,
                         ataccc_ids = c(1990, 1988, 1999),
                         nba_ids    = c(776, 253, 1650),
                         legacy_ids = c(2143, 2185, 2226))

  out_file <- file.path(out_dir, style, "fig6_inferred_pfu.pdf")
  ggplot2::ggsave(out_file, p, width = 14, height = 10, device = "pdf")

  out_file
}
