# ──────────────────────────────────────────────────────────────────────────────
# site_analysis.R — Illustrative analysis of anatomical site-specific viral
#                    kinetics in the HCT human challenge cohort
#
# Generates Supplementary Figure: figS_site_comparison.pdf
#   Panel A: Individual nose vs throat RNA trajectories for representative
#            HCT participants, with PFU (throat-only) overlay
#   Panel B: Population-level median (IQR) log(nose) − log(throat) difference
#            as a function of days since first positive
#
# Usage:
#   source("R/figure_themes.R")
#   source("R/site_analysis.R")
#   p <- plot_site_comparison()
#   save_journal_figure(p, "figS_site_comparison", layout = "full")
# ──────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Plot site-specific RNA trajectories and nose-throat differences (HCT)
#'
#' @param hct_path Path to the HCT data CSV
#' @param pids_panel_a Integer vector of participant IDs (pid) for Panel A.
#'   Defaults to c(1, 5, 7) — chosen to show a range of nose-throat offsets.
#' @param style Journal style for theming (NULL = current)
#' @return A patchwork ggplot (Panel A + Panel B stacked)
plot_site_comparison <- function(
    hct_path = "data/hct_dat.csv",
    pids_panel_a = c(1, 5, 7),
    style = NULL
) {

  # -- load journal theme --
  if (!is.null(style)) set_journal(style)
  cols <- journal_colors()

  # -- read and process HCT data --
  hct <- read.csv(hct_path, check.names = FALSE) |>
    mutate(
      pid = as.integer(pid),
      # Replace zeros with 1 before logging (consistent with clean_hct)
      rna_nose   = ifelse(qpcr_nose   == 0, NA_real_, log(qpcr_nose)),
      rna_throat = ifelse(qpcr_throat == 0, NA_real_, log(qpcr_throat)),
      pfu_throat = ifelse(throat == 0 | is.na(throat), NA_real_, log(throat)),
      # Current model average (log of arithmetic mean)
      rna_avg = case_when(
        !is.na(rna_nose) & !is.na(rna_throat) ~
          log(exp(rna_throat) / 2 + exp(rna_nose) / 2),
        !is.na(rna_throat) ~ rna_throat,
        !is.na(rna_nose)   ~ rna_nose,
        TRUE ~ NA_real_
      )
    )

  # Determine first positive day per individual (either site)
  first_pos <- hct |>
    filter(!is.na(rna_nose) | !is.na(rna_throat)) |>
    group_by(pid) |>
    summarise(day_first_pos = min(day), .groups = "drop")

  hct <- left_join(hct, first_pos, by = "pid") |>
    mutate(days_since_pos = day - day_first_pos)

  # ── Panel A: Individual trajectories ────────────────────────────────────────

  panel_a_dat <- hct |>
    filter(pid %in% pids_panel_a) |>
    mutate(pid_label = paste0("Participant ", pid))

  # Pivot longer for nose/throat/average
  rna_long <- panel_a_dat |>
    select(pid, pid_label, day, rna_nose, rna_throat, rna_avg) |>
    pivot_longer(
      cols = c(rna_nose, rna_throat, rna_avg),
      names_to = "site",
      values_to = "rna"
    ) |>
    filter(!is.na(rna)) |>
    mutate(
      site = recode(site,
        rna_nose   = "Nose",
        rna_throat = "Throat",
        rna_avg    = "Average (current model)"
      ),
      site = factor(site, levels = c("Nose", "Throat", "Average (current model)"))
    )

  # PFU data (throat only)
  pfu_dat <- panel_a_dat |>
    select(pid, pid_label, day, pfu_throat) |>
    filter(!is.na(pfu_throat))

  # LOD line
  lod_rna <- log(31.62278)  # HCT LOD = 10^1.5 copies/mL

  # -- Fit smoothed piecewise exponential to each pid × site ----
  # Uses smfun from R/utils.R
  fit_smfun <- function(t, y, lod = lod_rna) {
    # Only fit if enough above-LOD points to identify parameters
    above_lod <- y > lod
    if (sum(above_lod) < 4) return(NULL)
    # Initial parameter guesses from data
    tp0 <- t[which.max(y)]
    dp0 <- max(y)
    wp0 <- max(tp0 - min(t[above_lod]), 1)
    wr0 <- max(max(t[above_lod]) - tp0, 1)
    nll <- function(par) {
      tp <- par[1]; dp <- par[2]; wp <- exp(par[3]); wr <- exp(par[4])
      yhat <- smfun(t, tp, wp, wr, dp)
      sum((y - yhat)^2)
    }
    fit <- tryCatch(
      optim(c(tp0, dp0, log(wp0), log(wr0)), nll, method = "Nelder-Mead"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    list(tp = fit$par[1], dp = fit$par[2],
         wp = exp(fit$par[3]), wr = exp(fit$par[4]))
  }

  # Generate smooth curves for each pid × site (RNA)
  smooth_curves <- rna_long |>
    group_by(pid, pid_label, site) |>
    group_modify(function(.x, .key) {
      pars <- fit_smfun(.x$day, .x$rna)
      if (is.null(pars)) return(tibble())
      t_grid <- seq(min(.x$day) - 0.25, max(.x$day) + 0.25, by = 0.1)
      tibble(
        day = t_grid,
        rna = smfun(t_grid, pars$tp, pars$wp, pars$wr, pars$dp)
      )
    }) |>
    ungroup()

  # Fit smoothed piecewise exponential to PFU (throat-only) per pid
  lod_pfu <- log(5)  # HCT PFU LOD = 5 PFU/mL
  pfu_smooth <- pfu_dat |>
    group_by(pid, pid_label) |>
    group_modify(function(.x, .key) {
      pars <- fit_smfun(.x$day, .x$pfu_throat)
      if (is.null(pars)) return(tibble())
      t_grid <- seq(min(.x$day) - 0.25, max(.x$day) + 0.25, by = 0.1)
      tibble(
        day = t_grid,
        pfu = smfun(t_grid, pars$tp, pars$wp, pars$wr, pars$dp)
      )
    }) |>
    ungroup()

  pa <- ggplot(rna_long, aes(x = day, y = rna, color = site)) +
    geom_point(aes(shape = site), size = 1.4, alpha = 0.8) +
    geom_line(
      data = smooth_curves,
      aes(x = day, y = rna, color = site, linetype = site),
      linewidth = 0.9
    ) +
    geom_point(
      data = pfu_dat,
      aes(x = day, y = pfu_throat),
      color = cols$accent, shape = 17, size = 1.8,
      inherit.aes = FALSE
    ) +
    geom_line(
      data = pfu_smooth,
      aes(x = day, y = pfu),
      color = cols$accent, linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = lod_rna, linetype = "dotted", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = Inf, y = lod_rna, label = "RNA LOD",
             hjust = 1.1, vjust = -0.5, size = 2.2, color = "grey50") +
    facet_wrap(~pid_label, ncol = 3, scales = "free_x") +
    scale_color_manual(
      values = c(
        "Nose"   = cols$rna,
        "Throat" = cols$pfu,
        "Average (current model)" = "grey30"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Nose"   = 16,
        "Throat" = 15,
        "Average (current model)" = 4
      ),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(
        "Nose"   = "solid",
        "Throat" = "solid",
        "Average (current model)" = "dashed"
      ),
      name = NULL
    ) +
    labs(
      x = "Day since inoculation",
      y = "log(copies/mL) or log(PFU/mL)",
      tag = "A"
    ) +
    theme_journal() +
    theme(
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      plot.tag = element_text(face = "bold", size = 10)
    )

  # ── Panel B: Population-level nose-throat difference ────────────────────────

  diff_dat <- hct |>
    filter(!is.na(rna_nose) & !is.na(rna_throat)) |>
    mutate(diff = rna_nose - rna_throat) |>
    # Round to nearest 0.5 day for binning
    mutate(day_bin = round(days_since_pos * 2) / 2) |>
    group_by(day_bin) |>
    summarise(
      median_diff = median(diff),
      q25 = quantile(diff, 0.25),
      q75 = quantile(diff, 0.75),
      n = n(),
      .groups = "drop"
    ) |>
    filter(n >= 3)  # at least 3 observations per bin

  pb <- ggplot(diff_dat, aes(x = day_bin, y = median_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = cols$ci_fill, alpha = 0.5) +
    geom_line(color = cols$rna, linewidth = 0.7) +
    geom_point(color = cols$rna, size = 1.5) +
    labs(
      x = "Days since first positive",
      y = expression(log(RNA[nose]) - log(RNA[throat])),
      tag = "B"
    ) +
    annotate("text", x = Inf, y = 0.3, label = "Nose > Throat",
             hjust = 1.1, vjust = 0, size = 2.5, color = "grey40",
             fontface = "italic") +
    annotate("text", x = Inf, y = -0.3, label = "Throat > Nose",
             hjust = 1.1, vjust = 1, size = 2.5, color = "grey40",
             fontface = "italic") +
    theme_journal() +
    theme(
      plot.tag = element_text(face = "bold", size = 10)
    )

  # ── Combine ─────────────────────────────────────────────────────────────────

  combined <- pa / pb +
    plot_layout(heights = c(1.2, 1))

  combined
}

#' Run the full site analysis and save the figure
#'
#' @param out_dir Output directory for figures
#' @return File paths of saved figures (invisible)
run_site_analysis <- function(out_dir = "output/figures") {
  source("R/figure_themes.R", local = TRUE)
  p <- plot_site_comparison()
  paths <- save_journal_figure(p, "figS_site_comparison",
                               layout = "full", out_dir = out_dir)
  message("Saved: ", paste(paths, collapse = ", "))
  invisible(paths)
}


#' Generate site analysis figures in one or more journal styles
#'
#' @param styles Character vector of journal styles
#' @param out_dir Output directory for figures
#' @param hct_path Path to HCT data CSV
#' @param pids_panel_a Integer vector of participant IDs for panel A
#' @return Character vector of saved file paths
generate_site_analysis_figures <- function(styles = c("pnas"),
                                           out_dir = "output/figures",
                                           hct_path = "data/hct_dat.csv",
                                           pids_panel_a = c(1, 5, 7)) {
  all_paths <- character(0)

  for (s in styles) {
    p <- plot_site_comparison(
      hct_path = hct_path,
      pids_panel_a = pids_panel_a,
      style = s
    )
    paths <- save_journal_figure(
      p,
      "figS_site_comparison",
      layout = "full",
      out_dir = out_dir,
      style = s
    )
    all_paths <- c(all_paths, paths)
  }

  unname(all_paths)
}
