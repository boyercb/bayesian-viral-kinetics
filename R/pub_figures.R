# ──────────────────────────────────────────────────────────────────────────────
# pub_figures.R — Publication-ready main text figures
#
# Generates 5 main text figures in all journal styles (PNAS, PLoS, Annals).
# Run: Rscript R/pub_figures.R [style]
#   style = "pnas" (default), "plos", "annals", or "all"
#
# Depends on: targets cache (predictions, param_summary, kinetics_mcmc, etc.)
# Output: output/figures/{style}/fig{N}_{name}.{pdf,png}
# ──────────────────────────────────────────────────────────────────────────────

# ── Constants used by figure functions ────────────────────────────────────────
source_names <- c("1" = "NBA", "2" = "ATACCC", "3" = "UIUC",
                   "4" = "HCT", "5" = "Legacy")

# ── Helper: flatten 1-column matrix columns ──────────────────────────────────
flatten_mat_cols <- function(df) {
  for (col in names(df)) {
    if (is.matrix(df[[col]]) && ncol(df[[col]]) == 1) {
      df[[col]] <- as.numeric(df[[col]][, 1])
    }
  }
  df
}

# ── Helper: safe tar_read with qs2 fallback ──────────────────────────────────
safe_tar_read <- function(name) {
  tryCatch(
    targets::tar_read_raw(name),
    error = function(e) {
      message("  tar_read failed for '", name, "', reading directly with qs2...")
      qs2::qs_read(file.path("_targets", "objects", name))
    }
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: Example Individual Trajectory Fits
# ══════════════════════════════════════════════════════════════════════════════

#' Create a single individual trajectory panel
#'
#' Shows RNA trajectory (blue), PFU trajectory (red), LFD probability tiles,
#' and symptom onset marker for one individual.
#'
#' @param obs_i   Observation-level data for this individual
#' @param grid_i  Dense grid data for this individual
#' @param lods    Named list with lod_rna and lod_pfu for this source
#' @param cols    Color palette from journal_colors()
#' @param title   Panel title
#' @param source  Source identifier ("1"=NBA, "2"=ATACCC, "3"=UIUC, "4"=HCT)
#' @return ggplot object
panel_individual <- function(obs_i, grid_i, lods, cols, title = "", source = NULL) {

  has_pfu <- any(obs_i$pfu_exist == 1)
  has_lfd <- any(obs_i$lfd_exist == 1)
  has_sym <- any(obs_i$sym_exist == 1)

  # Clamp at LOD
  clamp <- function(x, lod) replace(x, !is.na(x) & x < lod, lod)
  grid_i$rna_hat_c     <- clamp(grid_i$rna_hat,     lods$rna)
  grid_i$rna_hat_q1_c  <- clamp(grid_i$rna_hat_q1,  lods$rna)
  grid_i$rna_hat_q99_c <- clamp(grid_i$rna_hat_q99, lods$rna)

  if (has_pfu) {
    grid_i$pfu_hat_c     <- clamp(grid_i$pfu_hat,     lods$pfu)
    grid_i$pfu_hat_q1_c  <- clamp(grid_i$pfu_hat_q1,  lods$pfu)
    grid_i$pfu_hat_q99_c <- clamp(grid_i$pfu_hat_q99, lods$pfu)
  }

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = time))

  # RNA ribbon + line
  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rna_hat_q1_c, ymax = rna_hat_q99_c),
      data = grid_i, alpha = 0.15, fill = cols["rna"]
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = rna_hat_c),
      data = grid_i, color = cols["rna"], linewidth = 0.7
    )

  # RNA observed points
  rna_obs <- dplyr::filter(obs_i, rna_exist == 1)
  if (nrow(rna_obs) > 0) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(y = rna), data = rna_obs,
      color = cols["rna"], shape = 16, size = 2.2, alpha = 0.8
    )
  }

  # PFU ribbon + line + points
  if (has_pfu) {
    pfu_grid <- dplyr::filter(grid_i, !is.na(pfu_hat_c))
    pfu_obs  <- dplyr::filter(obs_i, pfu_exist == 1)

    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = pfu_hat_q1_c, ymax = pfu_hat_q99_c),
        data = pfu_grid, alpha = 0.15, fill = cols["pfu"]
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = pfu_hat_c),
        data = pfu_grid, color = cols["pfu"], linewidth = 0.7
      )

    if (nrow(pfu_obs) > 0) {
      # UIUC (source "3") uses TCID50 assay: hollow triangles
      pfu_shape <- if (!is.null(source) && source == "3") 2 else 17
      p <- p + ggplot2::geom_point(
        ggplot2::aes(y = pfu), data = pfu_obs,
        color = cols["pfu"], shape = pfu_shape, size = 2.2, alpha = 0.8,
        stroke = if (pfu_shape == 2) 0.8 else 0.5
      )
    }
  }

  # LOD reference lines
  if (!is.na(lods$rna) && lods$rna > 0) {
    p <- p + ggplot2::geom_hline(
      yintercept = lods$rna, linetype = "dotted",
      color = cols["rna"], alpha = 0.5, linewidth = 0.3
    )
  }
  if (has_pfu && !is.na(lods$pfu) && lods$pfu > 0) {
    p <- p + ggplot2::geom_hline(
      yintercept = lods$pfu, linetype = "dotted",
      color = cols["pfu"], alpha = 0.5, linewidth = 0.3
    )
  }

  # LFD tiles at top
  if (has_lfd) {
    y_max <- max(c(obs_i$rna[obs_i$rna_exist == 1],
                   grid_i$rna_hat_q99_c), na.rm = TRUE) + 0.5
    lfd_dat <- dplyr::filter(obs_i, lfd_exist == 1, !is.na(lfd_hat))
    if (nrow(lfd_dat) > 0) {
      lfd_dat$lfd_y <- y_max + 0.8
      # Two-tone diverging scale: low prob = light/pale, high prob = dark green
      # Use discrete bins for clearer visual contrast
      lfd_dat$lfd_bin <- cut(lfd_dat$lfd_hat,
        breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
        labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8")
      )
      # Stepped green palette: light → dark
      tile_colors <- c(
        "<0.2"    = "#f0f0f0",
        "0.2-0.4" = "#a1d99b",
        "0.4-0.6" = "#41ab5d",
        "0.6-0.8" = "#238b45",
        ">0.8"    = "#005a32"
      )
      # Text color: dark on light tiles, white on dark tiles
      lfd_dat$text_col <- ifelse(
        lfd_dat$lfd_hat < 0.4, "grey30", "white"
      )
      p <- p +
        ggplot2::geom_tile(
          ggplot2::aes(y = lfd_y, fill = lfd_bin,
                       width = 0.8, height = 1.0),
          data = lfd_dat
        ) +
        ggplot2::geom_text(
          ggplot2::aes(y = lfd_y, label = ifelse(lfd == 1, "+", "\u2013"),
                       color = text_col),
          data = lfd_dat, size = 2.8, fontface = "bold",
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(
          values = tile_colors, name = "P(LFD+)", guide = "none",
          drop = FALSE
        ) +
        ggplot2::scale_color_identity()
    }
  }

  # Symptom onset marker
  if (has_sym) {
    sym_dat <- obs_i %>%
      dplyr::filter(sym_exist == 1, sym_ever == 1, sym_onset < 90) %>%
      dplyr::slice(1)
    if (nrow(sym_dat) > 0) {
      p <- p + ggplot2::annotate(
        "segment",
        x = sym_dat$sym_onset, xend = sym_dat$sym_onset,
        y = -Inf, yend = Inf,
        linetype = "dashed", color = cols["sym"], alpha = 1.0,
        linewidth = 0.7
      )
    }
  }

  p <- p +
    ggplot2::labs(x = "Days from peak", y = "log copies/mL", title = title) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


#' Generate Figure 2: Example trajectory fits
#'
#' Selects representative individuals from each multi-modal cohort
#' (ATACCC, UIUC, HCT) and a few from NBA to show the model fit.
#'
#' @param predictions  List with obs and grid
#' @param stan_data    Stan data list
#' @param style        Journal style
#' @return ggplot (patchwork composite)
fig_example_trajectories <- function(predictions, stan_data, style = "pnas") {

  cols <- unlist(journal_colors(style))
  obs  <- predictions$obs
  grid <- predictions$grid

  # Select representative individuals from each cohort with multi-modal data
  # ATACCC (source "2"): pick 3 with PFU + LFD + symptoms
  ataccc_ids <- obs %>%
    dplyr::filter(source == "2") %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      has_pfu = any(pfu_exist == 1),
      has_lfd = any(lfd_exist == 1),
      has_sym = any(sym_exist == 1),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(has_pfu, has_lfd, has_sym) %>%
    dplyr::arrange(desc(n_obs)) %>%
    dplyr::slice(c(1, 3, 5)) %>%
    dplyr::pull(id)

  # HCT (source "4"): pick 3 (all have multi-modal)
  hct_ids <- obs %>%
    dplyr::filter(source == "4") %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(n_obs = n(), .groups = "drop") %>%
    dplyr::arrange(desc(n_obs)) %>%
    dplyr::slice(c(1, 5, 10)) %>%
    dplyr::pull(id)

  # UIUC (source "3"): pick 2 with PFU + LFD + symptoms
  uiuc_ids <- obs %>%
    dplyr::filter(source == "3") %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      has_pfu = any(pfu_exist == 1),
      has_lfd = any(lfd_exist == 1),
      has_sym = any(sym_exist == 1),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(has_pfu, has_lfd, has_sym) %>%
    dplyr::arrange(desc(n_obs)) %>%
    dplyr::slice(c(1, 3)) %>%
    dplyr::pull(id)

  # NBA (source "1"): pick 2 (RNA only — shows contrast)
  nba_ids <- obs %>%
    dplyr::filter(source == "1") %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(n_obs = n(), .groups = "drop") %>%
    dplyr::arrange(desc(n_obs)) %>%
    dplyr::slice(c(2, 8)) %>%
    dplyr::pull(id)

  all_ids <- list(
    "2" = ataccc_ids,
    "4" = hct_ids,
    "3" = uiuc_ids,
    "1" = nba_ids
  )

  panels <- list()
  for (src in names(all_ids)) {
    src_num <- as.integer(src)
    lods <- list(
      rna = stan_data$lod_rna[src_num],
      pfu = stan_data$lod_pfu[src_num]
    )

    for (ind_id in all_ids[[src]]) {
      obs_i  <- dplyr::filter(obs, id == ind_id, source == src)
      grid_i <- dplyr::filter(grid, id == ind_id, source == src)

      title_str <- paste0(source_names[src], " #", obs_i$pid[1])

      panel <- panel_individual(obs_i, grid_i, lods, cols, title = title_str,
                                 source = src) +
        theme_journal(style, base_size = 11)

      panels <- c(panels, list(panel))
    }
  }

  # Compose with patchwork
  p <- patchwork::wrap_plots(panels, ncol = 5) +
    patchwork::plot_annotation(
      title = NULL,
      theme = theme_journal(style) +
        ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 2, 2, "pt"))
    )

  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6: Inferred PFU Trajectories from PCR-only Data
#
# Shows that the joint model infers infectious virus (PFU) trajectories even
# for individuals with only RNA data (NBA) or RNA + symptoms (Legacy), by
# comparing with ATACCC individuals whose PFU was directly measured.
#
# Three rows × 3 columns:
#   Row 1 — ATACCC: PFU directly observed (data points + narrow CI)
#   Row 2 — NBA: RNA-only (inferred PFU, moderate CI)
#   Row 3 — Legacy: RNA + symptoms (inferred PFU, wider CI)
# ══════════════════════════════════════════════════════════════════════════════

#' Create a single panel for the inferred PFU figure
#'
#' Like panel_individual() but always shows PFU trajectory (even without data)
#' and adds a row label annotation.
#'
#' @param obs_i   Observation-level data for this individual
#' @param grid_i  Dense grid data for this individual
#' @param lods    Named list with rna and pfu LODs
#' @param cols    Color palette from journal_colors()
#' @param title   Panel title (e.g. "ATACCC #9")
#' @return ggplot object
panel_inferred_pfu <- function(obs_i, grid_i, lods, cols, title = "") {

  has_pfu_data <- any(obs_i$pfu_exist == 1)
  has_sym <- any(obs_i$sym_exist == 1)

  # Clamp at LOD
  clamp <- function(x, lod) replace(x, !is.na(x) & x < lod, lod)
  grid_i$rna_hat_c     <- clamp(grid_i$rna_hat,     lods$rna)
  grid_i$rna_hat_q1_c  <- clamp(grid_i$rna_hat_q1,  lods$rna)
  grid_i$rna_hat_q99_c <- clamp(grid_i$rna_hat_q99, lods$rna)

  # Always show PFU — clamp at PFU LOD (use RNA LOD if PFU LOD is 0)
  pfu_lod <- if (!is.na(lods$pfu) && lods$pfu > 0) lods$pfu else 0
  grid_i$pfu_hat_c     <- clamp(grid_i$pfu_hat,     pfu_lod)
  grid_i$pfu_hat_q1_c  <- clamp(grid_i$pfu_hat_q1,  pfu_lod)
  grid_i$pfu_hat_q99_c <- clamp(grid_i$pfu_hat_q99, pfu_lod)

  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = time))

  # ---- RNA ribbon + line --------------------------------------------------
  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rna_hat_q1_c, ymax = rna_hat_q99_c),
      data = grid_i, alpha = 0.15, fill = cols[["rna"]]
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = rna_hat_c),
      data = grid_i, color = cols[["rna"]], linewidth = 0.7
    )

  # RNA observed points
  rna_obs <- dplyr::filter(obs_i, rna_exist == 1)
  if (nrow(rna_obs) > 0) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(y = rna), data = rna_obs,
      color = cols[["rna"]], shape = 16, size = 2.2, alpha = 0.8
    )
  }

  # ---- PFU ribbon + line (always shown) -----------------------------------
  pfu_grid <- dplyr::filter(grid_i, !is.na(pfu_hat_c))
  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pfu_hat_q1_c, ymax = pfu_hat_q99_c),
      data = pfu_grid, alpha = 0.12, fill = cols[["pfu"]]
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pfu_hat_c),
      data = pfu_grid, color = cols[["pfu"]], linewidth = 0.7
    )

  # PFU observed points (only for ATACCC/UIUC/HCT)
  if (has_pfu_data) {
    pfu_obs <- dplyr::filter(obs_i, pfu_exist == 1)
    if (nrow(pfu_obs) > 0) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(y = pfu), data = pfu_obs,
        color = cols[["pfu"]], shape = 17, size = 2.2, alpha = 0.8
      )
    }
  }

  # ---- LOD reference lines -----------------------------------------------
  if (!is.na(lods$rna) && lods$rna > 0) {
    p <- p + ggplot2::geom_hline(
      yintercept = lods$rna, linetype = "dotted",
      color = cols[["rna"]], alpha = 0.5, linewidth = 0.3
    )
  }
  if (!is.na(lods$pfu) && lods$pfu > 0) {
    p <- p + ggplot2::geom_hline(
      yintercept = lods$pfu, linetype = "dotted",
      color = cols[["pfu"]], alpha = 0.5, linewidth = 0.3
    )
  }

  # ---- Symptom onset marker -----------------------------------------------
  if (has_sym) {
    sym_dat <- obs_i %>%
      dplyr::filter(sym_exist == 1, sym_ever == 1, sym_onset < 90) %>%
      dplyr::slice(1)
    if (nrow(sym_dat) > 0) {
      p <- p + ggplot2::annotate(
        "segment",
        x = sym_dat$sym_onset, xend = sym_dat$sym_onset,
        y = -Inf, yend = Inf,
        linetype = "dashed", color = cols[["sym"]], alpha = 0.6,
        linewidth = 0.4
      )
    }
  }

  p <- p +
    ggplot2::labs(x = "Days from peak", y = "log copies/mL", title = title) +
    ggplot2::coord_cartesian(clip = "off")

  p
}


#' Generate Figure 6: Inferred PFU trajectories from RNA-only data
#'
#' Three-row figure comparing PFU inference quality across data regimes:
#' - Row 1: ATACCC (PFU directly observed — gold standard)
#' - Row 2: NBA (RNA only — inferred PFU)
#' - Row 3: Legacy (RNA + symptoms — inferred PFU)
#'
#' @param predictions  List with obs and grid
#' @param stan_data    Stan data list (for LODs)
#' @param style        Journal style
#' @param ataccc_ids   Individual IDs for ATACCC (default: auto-selected)
#' @param nba_ids      Individual IDs for NBA (default: auto-selected)
#' @param legacy_ids   Individual IDs for Legacy (default: auto-selected)
#' @return ggplot (patchwork composite)
fig_inferred_pfu <- function(predictions, stan_data, style = "pnas",
                              ataccc_ids = NULL, nba_ids = NULL,
                              legacy_ids = NULL) {

  cols <- unlist(journal_colors(style))
  obs  <- predictions$obs
  grid <- predictions$grid

  # ── Auto-select individuals if not specified ──────────────────────────────
  if (is.null(ataccc_ids)) {
    ataccc_ids <- obs %>%
      dplyr::filter(source == 2, pfu_exist == 1) %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(n_pfu = sum(pfu_exist), n_obs = n(), .groups = "drop") %>%
      dplyr::filter(n_pfu >= 10) %>%
      dplyr::arrange(desc(n_pfu)) %>%
      dplyr::slice(c(1, 4, 7)) %>%
      dplyr::pull(id)
  }

  if (is.null(nba_ids)) {
    nba_ids <- obs %>%
      dplyr::filter(source == 1, rna_exist == 1) %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        n_obs = n(),
        peak_rna = max(rna, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_obs >= 10, peak_rna >= 14, peak_rna <= 20) %>%
      dplyr::arrange(desc(n_obs)) %>%
      dplyr::slice(c(1, 3, 5)) %>%
      dplyr::pull(id)
  }

  if (is.null(legacy_ids)) {
    legacy_ids <- obs %>%
      dplyr::filter(source == 5, rna_exist == 1) %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        n_obs = n(),
        peak_rna = max(rna, na.rm = TRUE),
        has_sym = any(obs$sym_exist[obs$id == dplyr::first(id) & obs$source == 5] == 1),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_obs >= 7, peak_rna >= 16, has_sym) %>%
      dplyr::arrange(desc(n_obs)) %>%
      dplyr::slice(c(1, 3, 5)) %>%
      dplyr::pull(id)
  }

  # ── Build panels ──────────────────────────────────────────────────────────
  src_label <- c("1" = "NBA", "2" = "ATACCC", "3" = "UIUC",
                  "4" = "HCT", "5" = "Legacy")
  row_labels <- c(
    "2" = "A: Directly observed (ATACCC)",
    "1" = "B: Inferred from RNA (NBA)",
    "5" = "C: Inferred from RNA + symptoms (Legacy)"
  )

  all_ids <- list(
    "2" = ataccc_ids,
    "1" = nba_ids,
    "5" = legacy_ids
  )

  panels <- list()
  for (src in c("2", "1", "5")) {
    src_num <- as.integer(src)
    lods <- list(
      rna = stan_data$lod_rna[src_num],
      pfu = stan_data$lod_pfu[src_num]
    )

    ids <- all_ids[[src]]
    for (j in seq_along(ids)) {
      ind_id <- ids[j]
      obs_i  <- dplyr::filter(obs, id == ind_id, source == src_num)
      grid_i <- dplyr::filter(grid, id == ind_id, source == src_num)

      if (nrow(obs_i) == 0 || nrow(grid_i) == 0) next

      title_str <- paste0(src_label[src], " #", obs_i$pid[1])
      panel <- panel_inferred_pfu(obs_i, grid_i, lods, cols,
                                   title = title_str) +
        theme_journal(style, base_size = 9.5)

      panels <- c(panels, list(panel))
    }
  }

  # ── Compose with row labels ─────────────────────────────────────────────
  n_cols <- 3
  p <- patchwork::wrap_plots(panels, ncol = n_cols) +
    patchwork::plot_annotation(
      theme = theme_journal(style, base_size = 10) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(4, 4, 4, 4, "pt")
        )
    )

  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3: Population-Level Trajectory Summary  (two-panel design)
#
#   Panel A — Population-mean trajectory computed from median posterior
#             parameters, showing the smooth piecewise-exponential shape
#             for RNA + PFU, plus LFD positivity and symptom hazard on a
#             secondary y-axis.
#
#   Panel B — Spaghetti plot: latent trajectories drawn from the posterior
#             predictive by sampling new individual random effects for
#             each of 200 posterior draws.
# ══════════════════════════════════════════════════════════════════════════════

#' Extract thinned population parameter draws for trajectory plotting
#'
#' Pulls the subset of posterior draws needed by .trajectory_from_params and
#' .trajectory_with_re, thins to n_draws, and returns as a data.frame.
#'
#' @param fit      CmdStanMCMC fit object (must have accessible CSV output)
#' @param n_draws  Number of draws to retain (evenly thinned)
#' @param out_path If non-NULL, save as RDS to this path
#' @return data.frame with one row per draw, columns = population params
extract_pop_draws <- function(fit, n_draws = 200, out_path = NULL) {
  pop_vars <- c(
    "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
    paste0("tau_tp[", 1:2, "]"), paste0("tau_dp[", 1:2, "]"),
    paste0("tau_wp[", 1:2, "]"), paste0("tau_wr[", 1:2, "]"),
    "tau0_lfd", paste0("tau_lfd[", 1:4, "]"),
    "zeta_sym_intercept", "zeta_sym_pfu", "zeta_sym_rna",
    "zeta_sym_postpeak", "zeta_sym_postpeak_rna",
    "sigma_sym",
    paste0("sigma_ind_rna[", 1:4, "]"),
    paste0("sigma_ind_pfu[", 1:4, "]"),
    paste0("L_Omega_rna[", rep(1:4, each = 4), ",", rep(1:4, 4), "]")
  )
  drws <- posterior::as_draws_matrix(fit$draws(variables = pop_vars))
  n_total <- nrow(drws)
  if (n_total > n_draws) {
    idx <- round(seq(1, n_total, length.out = n_draws))
    drws <- drws[idx, ]
  }
  df <- as.data.frame(drws)
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(df, out_path)
    message(sprintf("Saved %d population draws to %s", nrow(df), out_path))
  }
  df
}

#' Compute a single latent trajectory + probabilities from a parameter vector
#'
#' @param pars  Named numeric — population params for one posterior draw
#' @param t     Time grid (days from RNA peak)
#' @param scale_vl  Normalization constant (prior_dp_mean)
#' @return data.frame with columns: time, rna, pfu, lfd_prob, sym_hazard
.trajectory_from_params <- function(pars, t, scale_vl = 17) {

  dp_rna <- pars["dp_mean_rna"]
  wp_rna <- pars["wp_mean_rna"]
  wr_rna <- pars["wr_mean_rna"]
  tp_rna <- 0                     # population mean peak at t = 0

  rna <- smfun(t, tp_rna, wp_rna, wr_rna, dp_rna, wf = 0)

  # PFU via log-affine transformation (no individual REs)
  dp_pfu <- exp(pars["tau_dp[1]"]) * dp_rna^pars["tau_dp[2]"]
  wp_pfu <- exp(pars["tau_wp[1]"]) * wp_rna^pars["tau_wp[2]"]
  wr_pfu <- exp(pars["tau_wr[1]"]) * wr_rna^pars["tau_wr[2]"]
  tp_pfu <- pars["tau_tp[1]"] + pars["tau_tp[2]"] * tp_rna

  pfu <- smfun(t, tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf = 0)
  pfu <- pmin(pfu, rna)           # PFU ≤ RNA constraint

  # LFD probability: inv_logit(tau0_lfd + tau_lfd[1]*rna + tau_lfd[2]*pfu + tau_lfd[3]*post_peak + tau_lfd[4]*post_peak*rna)
  post_peak <- soft_postpeak(t, tp_rna)
  lfd_logit <- pars["tau0_lfd"] +
    pars["tau_lfd[1]"] * rna +
    pars["tau_lfd[2]"] * pfu +
    pars["tau_lfd[3]"] * post_peak +
    pars["tau_lfd[4]"] * post_peak * rna
  lfd_prob <- expit(lfd_logit)

  # Symptom hazard (cloglog, no individual RE)
  zeta <- pars["zeta_sym_intercept"] +
    pars["zeta_sym_pfu"] * (pfu / scale_vl) +
    pars["zeta_sym_rna"] * (rna / scale_vl) +
    pars["zeta_sym_postpeak"] * post_peak +
    pars["zeta_sym_postpeak_rna"] * post_peak * (rna / scale_vl)
  sym_hz <- 1 - exp(-exp(pmin(zeta, 10)))

  data.frame(time = t, rna = rna, pfu = pfu,
             lfd_prob = lfd_prob, sym_hazard = sym_hz)
}


#' Sample one latent trajectory with new individual REs
#'
#' @param pars   Named numeric — one posterior draw of population params
#' @param t      Time grid
#' @param scale_vl  Normalization constant
#' @return data.frame with columns: time, rna, pfu, lfd_prob, sym_hazard
.trajectory_with_re <- function(pars, t, scale_vl = 17) {

  # RNA individual effects (correlated via Cholesky)
  si <- as.numeric(pars[paste0("sigma_ind_rna[", 1:4, "]")])
  L <- matrix(0, 4, 4)
  for (ii in 1:4) for (jj in 1:ii)
    L[ii, jj] <- pars[paste0("L_Omega_rna[", ii, ",", jj, "]")]
  LS <- diag(si) %*% L
  z4 <- rnorm(4)
  eta_re <- as.numeric(LS %*% z4)

  dp_rna <- pars["dp_mean_rna"] * exp(eta_re[2])
  wp_rna <- pars["wp_mean_rna"] * exp(eta_re[3])
  wr_rna <- pars["wr_mean_rna"] * exp(eta_re[4])
  tp_rna <- eta_re[1]

  rna <- smfun(t, tp_rna, wp_rna, wr_rna, dp_rna, wf = 0)

  # PFU via log-affine + independent PFU individual REs
  sd_pfu <- as.numeric(pars[paste0("sigma_ind_pfu[", 1:4, "]")])
  z_pfu  <- rnorm(4)

  dp_pfu <- exp(pars["tau_dp[1]"]) * dp_rna^pars["tau_dp[2]"] * exp(z_pfu[2] * sd_pfu[2])
  wp_pfu <- exp(pars["tau_wp[1]"]) * wp_rna^pars["tau_wp[2]"] * exp(z_pfu[3] * sd_pfu[3])
  wr_pfu <- exp(pars["tau_wr[1]"]) * wr_rna^pars["tau_wr[2]"] * exp(z_pfu[4] * sd_pfu[4])
  tp_pfu <- pars["tau_tp[1]"] + pars["tau_tp[2]"] * tp_rna + z_pfu[1] * sd_pfu[1]

  pfu <- smfun(t, tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf = 0)
  pfu <- pmin(pfu, rna)

  # LFD probability
  post_peak <- soft_postpeak(t, tp_rna)
  lfd_logit <- pars["tau0_lfd"] +
    pars["tau_lfd[1]"] * rna +
    pars["tau_lfd[2]"] * pfu +
    pars["tau_lfd[3]"] * post_peak +
    pars["tau_lfd[4]"] * post_peak * rna
  lfd_prob <- expit(lfd_logit)

  # Symptom hazard (with individual RE)
  z_sym <- rnorm(1)
  u_sym <- pars["sigma_sym"] * z_sym
  zeta <- pars["zeta_sym_intercept"] +
    pars["zeta_sym_pfu"] * (pfu / scale_vl) +
    pars["zeta_sym_rna"] * (rna / scale_vl) +
    pars["zeta_sym_postpeak"] * post_peak +
    pars["zeta_sym_postpeak_rna"] * post_peak * (rna / scale_vl) +
    u_sym
  sym_hz <- 1 - exp(-exp(pmin(zeta, 10)))

  data.frame(time = t, rna = rna, pfu = pfu,
             lfd_prob = lfd_prob, sym_hazard = sym_hz)
}


#' Generate Figure 3: Population trajectory summary (two panels)
#'
#' Panel A: population-mean smooth trajectory from median posterior params,
#'          showing RNA + PFU (left axis) and LFD/symptom probabilities
#'          (right axis).
#' Panel B: spaghetti plot of 200 latent trajectories sampled from the
#'          posterior predictive (new individual REs per draw).
#'
#' @param stan_data    Stan data list (for scale_vl and LODs)
#' @param draws_path   Path to cached posterior draws RDS (used if draws_df is NULL)
#' @param draws_df     Optional data.frame of population draws (overrides draws_path)
#' @param style        Journal style
#' @param n_spaghetti  Number of spaghetti trajectories
#' @return ggplot (patchwork composite)
fig_population_trajectories <- function(stan_data,
                                         draws_path = "output/pop_draws_200.rds",
                                         draws_df = NULL,
                                         style = "pnas",
                                         n_spaghetti = 200) {

  cols <- unlist(journal_colors(style))
  scale_vl <- stan_data$prior_dp_mean
  lod_rna <- min(stan_data$lod_rna, na.rm = TRUE)
  lod_pfu <- min(stan_data$lod_pfu[stan_data$lod_pfu > 0], na.rm = TRUE)

  if (is.null(draws_df)) draws_df <- readRDS(draws_path)
  t_grid <- seq(-12, 22, by = 0.25)
  clamp <- function(x, lod) pmax(x, lod)

  # ── Compute trajectories from each of the 200 posterior draws ───────────
  all_traj <- lapply(seq_len(nrow(draws_df)), function(d) {
    pars <- as.numeric(draws_df[d, ])
    names(pars) <- names(draws_df)
    tr <- .trajectory_from_params(pars, t_grid, scale_vl)
    tr$draw <- d
    tr
  })
  all_traj_df <- dplyr::bind_rows(all_traj)

  # Population summary: median + 50/80% CrI from population-parameter draws
  pop_summ <- all_traj_df %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      rna_med = median(rna),
      rna_q10 = quantile(rna, 0.10),
      rna_q25 = quantile(rna, 0.25),
      rna_q75 = quantile(rna, 0.75),
      rna_q90 = quantile(rna, 0.90),
      pfu_med = median(pfu),
      pfu_q10 = quantile(pfu, 0.10),
      pfu_q25 = quantile(pfu, 0.25),
      pfu_q75 = quantile(pfu, 0.75),
      pfu_q90 = quantile(pfu, 0.90),
      lfd_med = median(lfd_prob),
      sym_med = median(sym_hazard),
      .groups = "drop"
    )

  # Clamp RNA & PFU at LODs
  pop_summ <- pop_summ %>%
    dplyr::mutate(
      rna_med = clamp(rna_med, lod_rna),
      rna_q10 = clamp(rna_q10, lod_rna),
      rna_q25 = clamp(rna_q25, lod_rna),
      rna_q75 = clamp(rna_q75, lod_rna),
      rna_q90 = clamp(rna_q90, lod_rna),
      pfu_med = clamp(pfu_med, lod_pfu),
      pfu_q10 = clamp(pfu_q10, lod_pfu),
      pfu_q25 = clamp(pfu_q25, lod_pfu),
      pfu_q75 = clamp(pfu_q75, lod_pfu),
      pfu_q90 = clamp(pfu_q90, lod_pfu)
    )

  y_max <- max(pop_summ$rna_q90, na.rm = TRUE)

  # ── LFD probability tiles (at integer days, displayed above trajectory) ─
  lfd_days <- seq(-10, 20, by = 1)
  med_pars <- apply(draws_df, 2, median)
  lfd_traj <- .trajectory_from_params(med_pars, lfd_days, scale_vl)
  lfd_tiles <- data.frame(
    time     = lfd_days,
    lfd_prob = lfd_traj$lfd_prob,
    sym_hz   = lfd_traj$sym_hazard,
    y_lfd    = y_max + 1.2,
    y_sym    = y_max + 2.6
  )

  # ── Panel A: population-mean trajectory ─────────────────────────────────
  p_a <- ggplot2::ggplot(pop_summ, ggplot2::aes(x = time)) +
    # RNA 80% CI ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rna_q10, ymax = rna_q90),
      alpha = 0.12, fill = cols["rna"]
    ) +
    # RNA 50% CI ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rna_q25, ymax = rna_q75),
      alpha = 0.20, fill = cols["rna"]
    ) +
    # RNA median
    ggplot2::geom_line(
      ggplot2::aes(y = rna_med), color = cols["rna"],
      linewidth = 1.0
    ) +
    # PFU 80% CI ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pfu_q10, ymax = pfu_q90),
      alpha = 0.12, fill = cols["pfu"]
    ) +
    # PFU 50% CI ribbon
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pfu_q25, ymax = pfu_q75),
      alpha = 0.20, fill = cols["pfu"]
    ) +
    # PFU median
    ggplot2::geom_line(
      ggplot2::aes(y = pfu_med), color = cols["pfu"],
      linewidth = 1.0
    ) +
    # LOD reference lines
    ggplot2::geom_hline(
      yintercept = lod_rna, linetype = "dotted",
      color = cols["rna"], alpha = 0.4, linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = lod_pfu, linetype = "dotted",
      color = cols["pfu"], alpha = 0.4, linewidth = 0.3
    ) +
    # LFD probability tiles — 5-bin stepped green palette (matches Fig 2)
    ggplot2::geom_tile(
      data = lfd_tiles %>%
        dplyr::mutate(
          lfd_bin = cut(lfd_prob,
            breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
            labels = c("lfd_1", "lfd_2", "lfd_3", "lfd_4", "lfd_5")
          ),
          text_col = ifelse(lfd_prob < 0.4, "grey30", "white")
        ),
      ggplot2::aes(x = time, y = y_lfd, fill = lfd_bin,
                   width = 0.9, height = 1.0),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(lfd_tiles, lfd_prob >= 0.01) %>%
        dplyr::mutate(text_col = ifelse(lfd_prob < 0.4, "grey30", "white")),
      ggplot2::aes(x = time, y = y_lfd,
                   label = sprintf(".%02.0f", lfd_prob * 100),
                   color = text_col),
      inherit.aes = FALSE, size = 2.2, fontface = "bold",
      show.legend = FALSE
    ) +
    # Symptom hazard tiles — 5-bin stepped yellow palette
    ggplot2::geom_tile(
      data = lfd_tiles %>%
        dplyr::mutate(
          sym_bin = cut(sym_hz,
            breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
            labels = c("sym_1", "sym_2", "sym_3", "sym_4", "sym_5")
          ),
          text_col = ifelse(sym_hz < 0.4, "grey30", "white")
        ),
      ggplot2::aes(x = time, y = y_sym, fill = sym_bin,
                   width = 0.9, height = 1.0),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(lfd_tiles, sym_hz >= 0.01) %>%
        dplyr::mutate(text_col = ifelse(sym_hz < 0.4, "grey30", "white")),
      ggplot2::aes(x = time, y = y_sym,
                   label = sprintf(".%02.0f", sym_hz * 100),
                   color = text_col),
      inherit.aes = FALSE, size = 2.2, fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        # Green bins (LFD) — matches Fig 2
        "lfd_1" = "#f0f0f0", "lfd_2" = "#a1d99b",
        "lfd_3" = "#41ab5d", "lfd_4" = "#238b45", "lfd_5" = "#005a32",
        # Yellow/amber bins (symptoms)
        "sym_1" = "#f0f0f0", "sym_2" = "#fee391",
        "sym_3" = "#fec44f", "sym_4" = "#d95f0e", "sym_5" = "#993404"
      ),
      guide = "none"
    ) +
    ggplot2::scale_color_identity() +
    # Row labels
    ggplot2::annotate("text", x = -12.5, y = y_max + 1.2,
                      label = "P(LFD+)", size = 3.2, hjust = 1,
                      color = "grey30", fontface = "italic") +
    ggplot2::annotate("text", x = -12.5, y = y_max + 2.6,
                      label = "P(sym)", size = 3.2, hjust = 1,
                      color = "grey30", fontface = "italic") +
    # Legend annotations
    ggplot2::annotate("text", x = 18, y = y_max * 0.95,
                      label = "RNA", color = cols["rna"],
                      hjust = 0, size = 3.8, fontface = "bold") +
    ggplot2::annotate("text", x = 18, y = y_max * 0.85,
                      label = "PFU", color = cols["pfu"],
                      hjust = 0, size = 3.8, fontface = "bold") +
    ggplot2::labs(x = "Days from peak", y = "log copies/mL",
                  title = "A. Population mean") +
    ggplot2::coord_cartesian(clip = "off",
                              xlim = c(-12, 22),
                              ylim = c(lod_pfu - 0.5, y_max + 3.5)) +
    theme_journal(style, base_size = 11)

  # ── Panel B: spaghetti of sampled latent trajectories ───────────────────
  set.seed(42)
  n_use <- min(n_spaghetti, nrow(draws_df))

  spaghetti_list <- lapply(seq_len(n_use), function(d) {
    pars <- as.numeric(draws_df[d, ])
    names(pars) <- names(draws_df)
    traj <- .trajectory_with_re(pars, t_grid, scale_vl)
    traj$draw <- d
    traj
  })
  spaghetti <- dplyr::bind_rows(spaghetti_list)

  # Clamp spaghetti at LOD
  spaghetti$rna <- clamp(spaghetti$rna, lod_rna)
  spaghetti$pfu <- clamp(spaghetti$pfu, lod_pfu)

  p_b <- ggplot2::ggplot(spaghetti, ggplot2::aes(x = time, group = draw)) +
    # RNA spaghetti
    ggplot2::geom_line(
      ggplot2::aes(y = rna), alpha = 0.08, linewidth = 0.3,
      color = cols["rna"]
    ) +
    # PFU spaghetti
    ggplot2::geom_line(
      ggplot2::aes(y = pfu), alpha = 0.08, linewidth = 0.3,
      color = cols["pfu"]
    ) +
    # LODs
    ggplot2::geom_hline(
      yintercept = lod_rna, linetype = "dotted",
      color = cols["rna"], alpha = 0.4, linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = lod_pfu, linetype = "dotted",
      color = cols["pfu"], alpha = 0.4, linewidth = 0.3
    ) +
    # Population mean on top for reference
    ggplot2::geom_line(
      data = pop_summ, inherit.aes = FALSE,
      mapping = ggplot2::aes(x = time, y = rna_med),
      color = cols["rna"], linewidth = 0.9, alpha = 0.8
    ) +
    ggplot2::geom_line(
      data = pop_summ, inherit.aes = FALSE,
      mapping = ggplot2::aes(x = time, y = pfu_med),
      color = cols["pfu"], linewidth = 0.9, alpha = 0.8
    ) +
    ggplot2::annotate("text", x = 18, y = y_max * 0.95,
                      label = "RNA", color = cols["rna"],
                      hjust = 0, size = 3.8, fontface = "bold") +
    ggplot2::annotate("text", x = 18, y = y_max * 0.85,
                      label = "PFU", color = cols["pfu"],
                      hjust = 0, size = 3.8, fontface = "bold") +
    ggplot2::labs(x = "Days from peak",
                  y = "log copies/mL",
                  title = "B. Posterior predictive draws (N = 200)") +
    ggplot2::coord_cartesian(
      xlim = c(-12, 22),
      ylim = c(lod_pfu - 0.5, y_max + 3.5)
    ) +
    theme_journal(style, base_size = 11)

  p_a + p_b +
    patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(
      theme = theme_journal(style)
    )
}


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4: Covariate Effects Forest Plot
# ══════════════════════════════════════════════════════════════════════════════

#' Generate Figure 4: Forest plot of covariate effects on RNA kinetics
#'
#' Sorted within each group with visual separators between groups.
#'
#' @param param_summary  Output of summarize_parameters()
#' @param style          Journal style
#' @return ggplot
fig_forest_covariates <- function(param_summary, style = "pnas") {

  cols <- unlist(journal_colors(style))

  df <- param_summary$covariate_effects
  df$coef  <- as.numeric(df$coef)
  df$ci_lo <- as.numeric(df$ci_lo)
  df$ci_hi <- as.numeric(df$ci_hi)

  # Readable labels
  param_labels <- c(
    dp = "A. Peak viral load",
    wp = "B. Proliferation duration",
    wr = "C. Clearance duration"
  )
  df$param_label <- factor(
    param_labels[df$parameter],
    levels = c("A. Peak viral load", "B. Proliferation duration",
               "C. Clearance duration")
  )

  # Significance: CI excludes 1
  df$sig <- df$ci_lo > 1 | df$ci_hi < 1

  # Assign groups based on original label prefixes (before any cleaning)
  df$group <- dplyr::case_when(
    grepl("^Age:",        df$label) ~ "Age",
    grepl("^Recurrence:", df$label) ~ "Infection history",
    grepl("^Variant:",    df$label) ~ "Variant",
    grepl("^History:",    df$label) ~ "Vaccination",
    TRUE                            ~ "Other"
  )
  df$group <- factor(df$group, levels = c("Age", "Infection history",
                                           "Variant", "Vaccination"))

  # Clean display labels (strip prefixes for readability)
  df$display <- gsub("^Age: ",        "", df$label)
  df$display <- gsub("^Recurrence: ", "", df$display)
  df$display <- gsub("^Variant: ",    "", df$display)
  df$display <- gsub("^History: ",    "", df$display)

  # Explicit within-group sort order (top-to-bottom reading)
  # Use the *cleaned* display names here
  display_order <- c(
    # Age
    "[30,50)", "[50,100)",
    # Infection history
    "Yes",
    # Variant
    "Alpha", "Delta", "Omicron", "BA.4/BA.5", "Other",
    # Vaccination
    "Vaccinated boosted", "Vaccinated unboosted",
    "Vaccinated unreported", "Unreported",
    "Boosted unreported primary"
  )

  # Keep only display names present in data; factor levels bottom→top
  present <- display_order[display_order %in% df$display]

  # Build y-axis levels with group headers inserted between groups

  # Map each display name to its group
  display_to_group <- stats::setNames(
    as.character(df$group[match(present, df$display)]), present
  )

  # Insert group header labels before the first item in each group
  y_levels <- character()
  y_bold   <- character()   # track which levels are group headers
  prev_group <- ""
  for (d in present) {
    cur_group <- display_to_group[d]
    if (cur_group != prev_group) {
      y_levels <- c(y_levels, cur_group)
      y_bold   <- c(y_bold, cur_group)
      prev_group <- cur_group
    }
    y_levels <- c(y_levels, d)
  }

  df$display <- factor(df$display, levels = rev(y_levels))

  # Build y-axis label formatting: bold for group headers, plain for items
  y_label_faces <- ifelse(rev(y_levels) %in% y_bold, "bold", "plain")

  # Compute dashed separator positions: aligned with each group header
  y_num_map <- stats::setNames(seq_along(rev(y_levels)), rev(y_levels))
  sep_y <- numeric()
  for (hdr in y_bold[-1]) {  # skip first header (no group above it)
    hdr_y <- y_num_map[hdr]
    sep_y <- c(sep_y, hdr_y)
  }
  sep_df <- data.frame(y = sep_y)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = display)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "solid",
                         color = "grey80", linewidth = 0.3) +
    ggplot2::geom_hline(
      data = sep_df, ggplot2::aes(yintercept = y),
      linetype = "dashed", color = "grey70", linewidth = 0.3
    ) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi, color = sig),
      size = 0.5, linewidth = 0.7, fatten = 3.5
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = cols["sig"], "FALSE" = cols["nonsig"]),
      guide  = "none"
    ) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::facet_wrap(~ param_label, ncol = 3, scales = "free_x") +
    ggplot2::labs(
      x = "Multiplicative effect (95% CrI)",
      y = NULL
    ) +
    theme_journal(style, base_size = 11) +
    ggplot2::theme(
      strip.text  = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(
        size = 9,
        face = y_label_faces
      )
    )

  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5: Correlation Matrix Heatmap
# ══════════════════════════════════════════════════════════════════════════════

#' Generate Figure 5: RNA individual-effect correlation matrix
#'
#' @param param_summary  Output of summarize_parameters()
#' @param style          Journal style
#' @return ggplot
fig_correlation_matrix <- function(param_summary, style = "pnas") {

  cols <- unlist(journal_colors(style))
  corr_df <- param_summary$corr_params

  # Extract correlation rows (Omega_rna)
  corr_rows <- corr_df %>% dplyr::filter(grepl("Omega_rna", parameter))

  # Parse indices
  corr_rows <- corr_rows %>%
    dplyr::mutate(
      i = as.integer(gsub("Omega_rna\\[(\\d+),.*", "\\1", parameter)),
      j = as.integer(gsub("Omega_rna\\[\\d+,(\\d+)\\]", "\\1", parameter)),
      estimate = as.numeric(estimate),
      ci_lo    = as.numeric(ci_lo),
      ci_hi    = as.numeric(ci_hi)
    )

  re_labels <- c(
    expression(italic(t)[italic(p)]),
    expression(delta[italic(p)]),
    expression(omega[italic(p)]),
    expression(omega[italic(r)])
  )
  re_short <- c("tp", "dp", "wp", "wr")

  # Build full 4x4 matrix
  mat_df <- expand.grid(i = 1:4, j = 1:4) %>%
    dplyr::left_join(
      corr_rows %>% dplyr::select(i, j, estimate, ci_lo, ci_hi),
      by = c("i", "j")
    )
  # Fill diagonal
  mat_df$estimate[mat_df$i == mat_df$j] <- 1
  mat_df$ci_lo[mat_df$i == mat_df$j] <- NA
  mat_df$ci_hi[mat_df$i == mat_df$j] <- NA

  # Fill lower triangle by symmetry
  for (k in seq_len(nrow(mat_df))) {
    if (is.na(mat_df$estimate[k])) {
      mirror <- mat_df$i[k]
      orig   <- mat_df$j[k]
      match_row <- which(mat_df$i == orig & mat_df$j == mirror)
      if (length(match_row) == 1) {
        mat_df$estimate[k] <- mat_df$estimate[match_row]
        mat_df$ci_lo[k]    <- mat_df$ci_lo[match_row]
        mat_df$ci_hi[k]    <- mat_df$ci_hi[match_row]
      }
    }
  }

  mat_df$row_lab <- factor(re_short[mat_df$i], levels = re_short)
  mat_df$col_lab <- factor(re_short[mat_df$j], levels = re_short)

  # Annotation: show estimate [CI] for off-diagonal, "1" for diagonal
  mat_df$label <- ifelse(
    mat_df$i == mat_df$j, "",
    sprintf("%.2f", mat_df$estimate)
  )
  mat_df$ci_label <- ifelse(
    mat_df$i == mat_df$j, "",
    sprintf("[%.2f, %.2f]", mat_df$ci_lo, mat_df$ci_hi)
  )

  # Only show upper triangle + diagonal
  mat_df$show <- mat_df$i <= mat_df$j

  p <- ggplot2::ggplot(
    dplyr::filter(mat_df, show),
    ggplot2::aes(x = col_lab, y = row_lab)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = estimate), colour = "white", linewidth = 0.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 4.5, fontface = "bold", color = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = ci_label),
      size = 2.8, color = "black", nudge_y = -0.25
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation",
      guide = ggplot2::guide_colorbar(
        barwidth = 8, barheight = 0.5,
        title.position = "top", direction = "horizontal"
      )
    ) +
    ggplot2::scale_x_discrete(labels = re_labels) +
    ggplot2::scale_y_discrete(limits = rev(re_short), labels = rev(re_labels)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_journal(style, base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom"
    )

  # Add RE SD annotations alongside
  sigma_df <- corr_df %>%
    dplyr::filter(grepl("sigma_ind_rna", parameter)) %>%
    dplyr::mutate(
      idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", parameter)),
      param_lab = factor(re_short[idx], levels = rev(re_short)),
      estimate = as.numeric(estimate),
      ci_lo = as.numeric(ci_lo),
      ci_hi = as.numeric(ci_hi)
    )

  p_sigma <- ggplot2::ggplot(sigma_df, ggplot2::aes(x = estimate, y = param_lab)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
      color = cols["rna"], size = 0.4, fatten = 3
    ) +
    ggplot2::scale_y_discrete(labels = rev(re_labels)) +
    ggplot2::labs(x = "RE standard deviation", y = NULL,
                  title = "Individual-effect SDs") +
    theme_journal(style, base_size = 11) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(face = "bold", size = 13),
      plot.title   = ggplot2::element_text(size = 13),
      axis.title.x = ggplot2::element_text(size = 11)
    )

  p + p_sigma +
    patchwork::plot_layout(widths = c(3, 2)) +
    patchwork::plot_annotation(tag_levels = "A")
}


#' Generate main-text publication figures from existing targets objects
#'
#' Target-friendly wrapper that generates Figures 2--6 for the requested
#' journal styles without reading targets cache from disk.
#'
#' @param predictions   Output of \code{compute_predictions()}
#' @param param_summary Output of \code{summarize_parameters()}
#' @param stan_data     Stan data list
#' @param pop_draws_df  Data frame of population draws (e.g., from
#'   \code{extract_pop_draws()})
#' @param styles        Character vector of styles to generate
#' @param out_dir       Output directory
#' @return Character vector of saved file paths
generate_pub_figures <- function(predictions,
                                 param_summary,
                                 stan_data,
                                 pop_draws_df,
                                 styles = c("pnas"),
                                 out_dir = "output/figures") {

  predictions$obs <- flatten_mat_cols(predictions$obs)

  tmp_draws <- tempfile(pattern = "pop_draws_", fileext = ".rds")
  saveRDS(pop_draws_df, tmp_draws)
  on.exit(unlink(tmp_draws), add = TRUE)

  all_paths <- character(0)
  for (s in styles) {
    message(sprintf("Generating publication figures for style: %s", s))
    set_journal(s)

    p2 <- fig_example_trajectories(predictions, stan_data, style = s)
    paths2 <- save_journal_figure(
      p2, "fig2_trajectories", layout = "full",
      width = 14, height = 8, out_dir = out_dir, style = s
    )

    p3 <- fig_population_trajectories(
      stan_data,
      draws_path = tmp_draws,
      style = s,
      n_spaghetti = min(200, nrow(pop_draws_df))
    )
    paths3 <- save_journal_figure(
      p3, "fig3_population", layout = "full",
      width = 14, height = 6, out_dir = out_dir, style = s
    )

    p4 <- fig_forest_covariates(param_summary, style = s)
    paths4 <- save_journal_figure(
      p4, "fig4_forest", layout = "full",
      width = 10, height = 5, out_dir = out_dir, style = s
    )

    p5 <- fig_correlation_matrix(param_summary, style = s)
    paths5 <- save_journal_figure(
      p5, "fig5_correlations", layout = "full",
      width = 10, height = 5, out_dir = out_dir, style = s
    )

    p6 <- fig_inferred_pfu(predictions, stan_data, style = s)
    paths6 <- save_journal_figure(
      p6, "fig6_inferred_pfu", layout = "full",
      width = 14, height = 10, out_dir = out_dir, style = s
    )

    all_paths <- c(all_paths, paths2, paths3, paths4, paths5, paths6)
  }

  unname(all_paths)
}


# ══════════════════════════════════════════════════════════════════════════════
# SCRIPT EXECUTION — only when run directly (Rscript R/pub_figures.R [style])
# Guarded so tar_source("R/") does not trigger figure generation at load time.
# ══════════════════════════════════════════════════════════════════════════════

if (sys.nframe() == 0L) {

library(targets)
library(tidyverse)
library(patchwork)

# Load project functions (exclude this file to avoid recursion)
r_files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
r_files <- r_files[!grepl("pub_figures\\.R$", r_files)]
invisible(lapply(r_files, source))

# ── Parse command-line style argument ─────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
run_style <- if (length(args) > 0) args[1] else "all"

if (run_style == "all") {
  styles <- c("pnas", "plos", "annals")
} else {
  styles <- match.arg(run_style, c("pnas", "plos", "annals"))
}

# ── Load cached targets ──────────────────────────────────────────────────────
message("Loading targets cache...")
predictions   <- safe_tar_read("predictions")
param_summary <- safe_tar_read("param_summary")
stan_data     <- safe_tar_read("stan_data")
stacked_dat   <- safe_tar_read("stacked_dat")

# Try loading the fit (may fail if CSVs are remote)
fit <- tryCatch(targets::tar_read(kinetics_mcmc), error = function(e) NULL)

predictions$obs <- flatten_mat_cols(predictions$obs)

for (s in styles) {
  message(sprintf("\n=== Generating figures for style: %s ===", s))
  set_journal(s)

  # Figure 2: Example trajectories
  message("  Fig 2: Example trajectory fits...")
  p2 <- fig_example_trajectories(predictions, stan_data, style = s)
  save_journal_figure(p2, "fig2_trajectories", layout = "full",
                      width = 14, height = 8, style = s)

  # Figure 3: Population trajectories (from cached posterior draws)
  message("  Fig 3: Population trajectory summary...")
  p3 <- fig_population_trajectories(stan_data,
                                     draws_path = "output/pop_draws_200.rds",
                                     style = s, n_spaghetti = 200)
  save_journal_figure(p3, "fig3_population", layout = "full",
                      width = 14, height = 6, style = s)

  # Figure 4: Forest plot
  message("  Fig 4: Covariate forest plot...")
  p4 <- fig_forest_covariates(param_summary, style = s)
  save_journal_figure(p4, "fig4_forest", layout = "full",
                      width = 10, height = 5, style = s)

  # Figure 5: Correlation matrix
  message("  Fig 5: Correlation matrix...")
  p5 <- fig_correlation_matrix(param_summary, style = s)
  save_journal_figure(p5, "fig5_correlations", layout = "full",
                      width = 10, height = 5, style = s)

  # Figure 6: Inferred PFU trajectories
  message("  Fig 6: Inferred PFU trajectories...")
  p6 <- fig_inferred_pfu(predictions, stan_data, style = s)
  save_journal_figure(p6, "fig6_inferred_pfu", layout = "full",
                      width = 14, height = 10, style = s)
}

message("\nAll figures generated successfully.")

}  # end if (sys.nframe() == 0L)
