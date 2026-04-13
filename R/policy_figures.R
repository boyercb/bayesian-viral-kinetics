# ──────────────────────────────────────────────────────────────────────────────
# policy_figures.R — Publication figures for policy analysis results
#
# Figure 6: Probability of ongoing infectiousness over time
# Figure 7: Isolation and test-to-release analysis
# ──────────────────────────────────────────────────────────────────────────────


#' Figure 6: Probability of ongoing infectiousness
#'
#' Three-panel figure showing:
#'   A: P(culture+), P(LFD+), P(RNA detectable) vs days since first positive PCR
#'   B: Same aligned to symptom onset
#'   C: P(culture+ | LFD+) and P(culture+ | LFD-) vs days since first positive LFD
#'
#' @param policy_results  Output of \code{compute_all_policy()}
#' @param style           Journal style (NULL = use active)
#' @return patchwork ggplot object
fig_probability_curves <- function(policy_results, style = NULL) {

  cols <- journal_colors(style)

  # Marker and outcome colors
  col_culture <- cols[["pfu"]]
  col_lfd     <- cols[["lfd"]]
  col_rna     <- cols[["rna"]]
  col_neg     <- cols[["muted"]]

  dodge <- ggplot2::position_dodge(width = 0.5)

  # ── Panel A: Days since first positive PCR ─────────────────────────────
  pcr <- policy_results$curves$pcr_marginal

  pcr_long <- rbind(
    data.frame(
      day = pcr$day_since_landmark,
      prob = pcr$p_culture_med,
      lo = pcr$p_culture_lo, hi = pcr$p_culture_hi,
      marker = "Culture positive"
    ),
    data.frame(
      day = pcr$day_since_landmark,
      prob = pcr$p_lfd_med,
      lo = pcr$p_lfd_lo, hi = pcr$p_lfd_hi,
      marker = "LFD positive"
    ),
    data.frame(
      day = pcr$day_since_landmark,
      prob = pcr$p_rna_med,
      lo = pcr$p_rna_lo, hi = pcr$p_rna_hi,
      marker = "RNA detectable"
    )
  )
  pcr_long$marker <- factor(pcr_long$marker,
    levels = c("RNA detectable", "Culture positive", "LFD positive"))

  marker_colors <- c(
    "RNA detectable"  = col_rna,
    "Culture positive" = col_culture,
    "LFD positive"    = col_lfd
  )

  p_a <- ggplot2::ggplot(pcr_long,
    ggplot2::aes(x = day, y = prob, color = marker)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
      position = dodge, linewidth = 0.6, alpha = 0.35) +
    ggplot2::geom_point(position = dodge, size = 2.5) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.4) +
    ggplot2::annotate("text", x = 5, y = 1.02, label = "Day 5",
      size = 4.5, color = "grey40", hjust = -0.1) +
    ggplot2::annotate("text", x = 10, y = 1.02, label = "Day 10",
      size = 4.5, color = "grey40", hjust = -0.1) +
    ggplot2::scale_color_manual(values = marker_colors, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
      labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      x = "Days since first positive PCR",
      y = "Probability"
    ) +
    theme_journal(style, base_size = 13) +
    ggplot2::theme(legend.position = "bottom")

  # ── Panel B: Days since symptom onset ──────────────────────────────────
  sym <- policy_results$curves$sym_marginal

  sym_long <- rbind(
    data.frame(
      day = sym$day_since_landmark,
      prob = sym$p_culture_med,
      lo = sym$p_culture_lo, hi = sym$p_culture_hi,
      marker = "Culture positive"
    ),
    data.frame(
      day = sym$day_since_landmark,
      prob = sym$p_lfd_med,
      lo = sym$p_lfd_lo, hi = sym$p_lfd_hi,
      marker = "LFD positive"
    ),
    data.frame(
      day = sym$day_since_landmark,
      prob = sym$p_rna_med,
      lo = sym$p_rna_lo, hi = sym$p_rna_hi,
      marker = "RNA detectable"
    )
  )
  sym_long$marker <- factor(sym_long$marker,
    levels = c("RNA detectable", "Culture positive", "LFD positive"))

  p_b <- ggplot2::ggplot(sym_long,
    ggplot2::aes(x = day, y = prob, color = marker)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
      position = dodge, linewidth = 0.6, alpha = 0.35) +
    ggplot2::geom_point(position = dodge, size = 2.5) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.4) +
    ggplot2::annotate("text", x = 5, y = 1.02, label = "Day 5",
      size = 4.5, color = "grey40", hjust = -0.1) +
    ggplot2::annotate("text", x = 10, y = 1.02, label = "Day 10",
      size = 4.5, color = "grey40", hjust = -0.1) +
    ggplot2::scale_color_manual(values = marker_colors, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
      labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      x = "Days since symptom onset",
      y = "Probability"
    ) +
    theme_journal(style, base_size = 13) +
    ggplot2::theme(legend.position = "bottom")

  # (Panel C removed — conditional LFD analysis moved to Bayesian filtering)

  # ── Compose ────────────────────────────────────────────────────────────
  (p_a | p_b) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(legend.position = "bottom")
}


#' Figure 7: Isolation and test-to-release analysis
#'
#' Multi-panel figure showing:
#'   A: Expected residual infectious days vs isolation release day
#'      (marginal + selected profiles)
#'   B: Test-to-release: P(infectious | LFD-) at each day since PCR
#'   C: Isolation table heatmap: % still infectious at days 5, 7, 10
#'      by covariate profile
#'
#' @param policy_results  Output of \code{compute_all_policy()}
#' @param style           Journal style (NULL = use active)
#' @return patchwork ggplot object
fig_isolation_policy <- function(policy_results, style = NULL) {

  cols <- journal_colors(style)

  # ── Panel A: Residual infectiousness vs release day ────────────────────
  #   Show a few key profiles vs the marginal average

  # Marginal AUC (only has release_day and stats)
  auc_marginal <- policy_results$residual_auc$pcr_marginal
  auc_marginal$label <- "Population average"

  # Selected stratified profiles
  key_profiles <- c("Reference", "Boosted", "Delta", "Omicron",
                     "Boosted + Omicron (2022)", "Unboosted + Delta (2021)")
  auc_strat <- policy_results$residual_auc$pcr_by_label |>
    dplyr::filter(label %in% key_profiles)

  auc_combined <- rbind(
    auc_marginal[, c("label", "release_day",
                     "residual_days_med", "residual_days_lo",
                     "residual_days_hi")],
    auc_strat[, c("label", "release_day",
                  "residual_days_med", "residual_days_lo",
                  "residual_days_hi")]
  )

  # Order labels sensibly
  label_order <- c("Population average", key_profiles)
  auc_combined$label <- factor(auc_combined$label, levels = label_order)

  # Color palette: use a qualitative scale
  profile_colors <- c(
    "Population average"         = "black",
    "Reference"                  = cols[["muted"]],
    "Boosted"                    = cols[["lfd"]],
    "Delta"                      = cols[["pfu"]],
    "Omicron"                    = cols[["rna"]],
    "Boosted + Omicron (2022)"   = cols[["accent"]],
    "Unboosted + Delta (2021)"   = cols[["sym"]]
  )

  p_a <- ggplot2::ggplot(auc_combined,
    ggplot2::aes(x = release_day, y = residual_days_med,
                 color = label, fill = label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = residual_days_lo,
                                       ymax = residual_days_hi),
      alpha = 0.08, color = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = label == "Population average"),
      linewidth = 0.7) +
    ggplot2::scale_linetype_manual(values = c("TRUE" = "solid",
                                               "FALSE" = "solid"),
      guide = "none") +
    ggplot2::scale_color_manual(values = profile_colors, name = NULL) +
    ggplot2::scale_fill_manual(values = profile_colors, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(3, 14, 1)) +
    ggplot2::labs(
      x = "Isolation release day (days since first positive PCR)",
      y = "Expected remaining\ninfectious days"
    ) +
    theme_journal(style, base_size = 10) +
    ggplot2::theme(
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha("white", 0.85), color = NA),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.text = ggplot2::element_text(size = 7)
    )

  # ── Panel B: Test-to-release false reassurance rate ────────────────────
  ttr <- policy_results$test_release$pcr

  p_b <- ggplot2::ggplot(ttr,
    ggplot2::aes(x = day_since_landmark,
                 y = p_infectious_given_neg_lfd_med)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = p_infectious_given_neg_lfd_lo,
                   ymax = p_infectious_given_neg_lfd_hi),
      linewidth = 0.4, alpha = 0.35, color = cols[["lfd"]]) +
    ggplot2::geom_point(color = cols[["lfd"]], size = 1.5) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.3) +
    ggplot2::annotate("text", x = 5, y = max(ttr$p_infectious_given_neg_lfd_hi,
      na.rm = TRUE) * 1.05,
      label = "Day 5", size = 3.5, color = "grey40", hjust = -0.1) +
    ggplot2::annotate("text", x = 10, y = max(ttr$p_infectious_given_neg_lfd_hi,
      na.rm = TRUE) * 1.05,
      label = "Day 10", size = 3.5, color = "grey40", hjust = -0.1) +
    ggplot2::scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      x = "Days since first positive PCR",
      y = "P(culture+ | LFD negative)"
    ) +
    theme_journal(style, base_size = 10)

  # ── Panel C: Isolation table heatmap ───────────────────────────────────
  iso <- policy_results$isolation$pcr |>
    dplyr::filter(label %in% key_profiles)

  iso$label <- factor(iso$label, levels = rev(key_profiles))
  iso$day_label <- factor(paste0("Day ", iso$day),
                          levels = paste0("Day ", sort(unique(iso$day))))

  p_c <- ggplot2::ggplot(iso,
    ggplot2::aes(x = day_label, y = label,
                 fill = p_still_infectious_med)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%",
                                    p_still_infectious_med * 100)),
      size = 4, fontface = "bold") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("(%0.f\u2013%0.f%%)",
                                    p_still_infectious_lo * 100,
                                    p_still_infectious_hi * 100)),
      size = 2.5, vjust = 2.0, color = "grey40") +
    colorspace::scale_fill_continuous_sequential(
      palette = "Reds 3", rev = FALSE,
      limits = c(0, 1), name = "P(still infectious)",
      labels = scales::percent
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_journal(style, base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1)
    )

  # ── Compose ────────────────────────────────────────────────────────────
  top_row <- p_a + p_b + patchwork::plot_layout(widths = c(3, 2))
  top_row / p_c +
    patchwork::plot_layout(heights = c(3, 2)) +
    patchwork::plot_annotation(tag_levels = "A")
}


#' Generate and save all policy figures
#'
#' Called by the targets pipeline.  Generates Fig 6, Fig 7, and the
#' stratified supplement figure in all journal styles.
#'
#' @param policy_results  Output of \code{compute_all_policy()}
#' @param styles          Character vector of journal styles
#' @return Character vector of saved file paths
generate_policy_figures <- function(policy_results,
                                    styles = c("pnas", "plos", "annals")) {

  all_paths <- character(0)

  for (s in styles) {
    message(sprintf("  Policy figures for style: %s", s))
    set_journal(s)

    # Figure 6
    p6 <- fig_probability_curves(policy_results, style = s)
    paths6 <- save_journal_figure(p6, "fig6_probability_curves",
                                   layout = "full",
                                   width = 14, height = 7, style = s)

    # Figure 7
    p7 <- fig_isolation_policy(policy_results, style = s)
    paths7 <- save_journal_figure(p7, "fig7_isolation_policy",
                                   layout = "full",
                                   width = 14, height = 8, style = s)

    # Supplement: stratified probability curves
    p_strat <- fig_stratified_curves(policy_results, style = s)
    paths_s <- save_journal_figure(p_strat, "figS_stratified_curves",
                                    layout = "full",
                                    width = 14, height = 10, style = s)

    # Figure 8 (or supplement): conditional on LFD by symptom onset
    p8 <- fig_conditional_symptom(policy_results, style = s)
    paths8 <- save_journal_figure(p8, "fig8_conditional_symptom",
                                   layout = "full",
                                   width = 17, height = 5.5, style = s)

    all_paths <- c(all_paths, paths6, paths7, paths_s, paths8)
  }

  all_paths
}


#' Supplement Figure: Culture positivity stratified by covariate profile
#'
#' Shows P(culture+) curves for each of the 10 covariate profiles,
#' highlighting how variant, vaccination, and infection history
#' modulate the duration of infectiousness.
#'
#' @param policy_results  Output of \code{compute_all_policy()}
#' @param style           Journal style (NULL = use active)
#' @return ggplot object
fig_stratified_curves <- function(policy_results, style = NULL) {

  cols <- journal_colors(style)
  curves <- policy_results$curves$pcr_by_label
  if (is.null(curves)) {
    stop("pcr_by_label not found in policy_results$curves. ",
         "Re-run compute_all_policy().")
  }

  # Select profiles that maximise differences
  key_profiles <- c(
    "Reference",
    "Boosted",
    "Reinfection",
    "Delta",
    "Omicron",
    "Age 50+",
    "Boosted + Omicron (2022)",
    "Unboosted + Delta (2021)"
  )

  d <- curves |>
    dplyr::filter(label %in% key_profiles) |>
    dplyr::mutate(
      label = factor(label, levels = key_profiles)
    )

  # Qualitative palette — pick n distinct colors
  n_profiles <- length(key_profiles)
  pal <- colorspace::qualitative_hcl(n_profiles, palette = "Dark 3")
  names(pal) <- key_profiles

  dodge <- ggplot2::position_dodge(width = 0.5)

  p_cult <- ggplot2::ggplot(d,
    ggplot2::aes(x = day_since_landmark, y = p_culture_med,
                 color = label)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = p_culture_lo, ymax = p_culture_hi),
      position = dodge, linewidth = 0.35, alpha = 0.3) +
    ggplot2::geom_point(position = dodge, size = 1.3) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.3) +
    ggplot2::annotate("text", x = 5, y = 1.02, label = "Day 5",
      size = 3.5, color = "grey40", hjust = -0.1) +
    ggplot2::annotate("text", x = 10, y = 1.02, label = "Day 10",
      size = 3.5, color = "grey40", hjust = -0.1) +
    ggplot2::scale_color_manual(values = pal, name = "Covariate profile") +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
      labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      title = "P(culture positive) by covariate profile",
      x = "Days since first positive PCR",
      y = "P(culture positive)"
    ) +
    theme_journal(style, base_size = 10) +
    ggplot2::theme(legend.position = "right")

  p_lfd <- ggplot2::ggplot(d,
    ggplot2::aes(x = day_since_landmark, y = p_lfd_med,
                 color = label)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = p_lfd_lo, ymax = p_lfd_hi),
      position = dodge, linewidth = 0.35, alpha = 0.3) +
    ggplot2::geom_point(position = dodge, size = 1.3) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.3) +
    ggplot2::scale_color_manual(values = pal, name = "Covariate profile") +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
      labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      title = "P(LFD positive) by covariate profile",
      x = "Days since first positive PCR",
      y = "P(LFD positive)"
    ) +
    theme_journal(style, base_size = 10) +
    ggplot2::theme(legend.position = "right")

  (p_cult / p_lfd) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(legend.position = "right")
}


# ── Bayesian filtering figure ────────────────────────────────────────────────

#' Figure 9: Bayesian updating of infectiousness with test history
#'
#' Two-panel figure demonstrating personalized infectiousness estimation:
#' \describe{
#'   \item{Panel A}{Effect of initial viral load at diagnosis.  Shows how
#'     P(culture positive) differs for Ct = 20 (high), 25 (typical), and
#'     30 (low) diagnoses, compared to the population-averaged baseline.}
#'   \item{Panel B}{Sequential updating with negative LFD results.  Starting
#'     from a Ct = 25 diagnosis, shows how each additional negative LFD
#'     result narrows the posterior probability of ongoing infectiousness.}
#' }
#'
#' @param filter_results  Output of \code{compute_bayesian_filtering()}
#' @param style           Journal style (NULL = use active)
#' @return patchwork ggplot object
fig_bayesian_filtering <- function(filter_results, style = NULL) {

  cols <- journal_colors(style)
  pop  <- filter_results$population
  sc   <- filter_results$scenarios

  # ── Panel A: Initial viral load ──────────────────────────────────────────
  panel_a_labels <- c("Ct = 30", "Ct = 25", "Ct = 20")
  panel_a_data <- dplyr::bind_rows(
    sc$ct_20$curves |> dplyr::mutate(scenario = "Ct = 20"),
    sc$ct_25$curves |> dplyr::mutate(scenario = "Ct = 25"),
    sc$ct_30$curves |> dplyr::mutate(scenario = "Ct = 30")
  ) |>
    dplyr::mutate(scenario = factor(scenario, levels = panel_a_labels))

  pal_a <- c("Ct = 30" = "#4575b4", "Ct = 25" = "#fc8d59", "Ct = 20" = "#d73027")

  p_a <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pop,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi),
      fill = "grey80", alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = pop,
      ggplot2::aes(x = day, y = p_infectious_med),
      color = "grey50", linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::geom_ribbon(
      data = panel_a_data,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi,
                   fill = scenario),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = panel_a_data,
      ggplot2::aes(x = day, y = p_infectious_med, color = scenario),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(values = pal_a, name = "Diagnosis") +
    ggplot2::scale_fill_manual(values = pal_a, name = "Diagnosis") +
    ggplot2::scale_y_continuous(
      "P(culture positive)",
      limits = c(0, 1), expand = c(0.01, 0),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(
      "Days since diagnosis",
      breaks = seq(0, 20, 5)
    ) +
    ggplot2::labs(tag = "A") +
    theme_journal(style, base_size = 13) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2),
      fill  = ggplot2::guide_legend(nrow = 2)
    ) +
    ggplot2::theme(legend.position = "bottom")

  # ── Panel B: Serial LFD testing ─────────────────────────────────────────
  panel_b_labels <- c(
    "Diagnosis only (Ct 25)",
    "+ LFD\u2212 day 5",
    "+ LFD\u2212 days 5\u20136",
    "+ LFD\u2212 days 5\u20137"
  )
  # Full solid curves throughout; vertical lines mark conditioning events
  panel_b_data <- dplyr::bind_rows(
    sc$ct_25$curves       |> dplyr::mutate(scenario = panel_b_labels[1]),
    sc$lfd_neg_d5$curves  |> dplyr::mutate(scenario = panel_b_labels[2]),
    sc$lfd_neg_d56$curves |> dplyr::mutate(scenario = panel_b_labels[3]),
    sc$lfd_neg_d567$curves |> dplyr::mutate(scenario = panel_b_labels[4])
  ) |>
    dplyr::mutate(scenario = factor(scenario, levels = panel_b_labels))

  pal_b <- c(
    "Diagnosis only (Ct 25)"     = "#d73027",
    "+ LFD\u2212 day 5"          = "#fc8d59",
    "+ LFD\u2212 days 5\u20136"  = "#fee090",
    "+ LFD\u2212 days 5\u20137"  = "#4575b4"
  )

  # Vertical line annotations for conditioning events
  vline_data <- data.frame(
    day   = c(0, 5, 6, 7),
    label = c("PCR", "LFD\u2212", "LFD\u2212", "LFD\u2212"),
    stringsAsFactors = FALSE
  )

  p_b <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pop,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi),
      fill = "grey80", alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = pop,
      ggplot2::aes(x = day, y = p_infectious_med),
      color = "grey50", linetype = "dashed", linewidth = 0.7
    ) +
    # Vertical lines at conditioning events
    ggplot2::geom_vline(
      xintercept = c(0, 5, 6, 7),
      color = "grey40", linetype = "dotted", linewidth = 0.5
    ) +
    # Event labels at top
    ggplot2::annotate("text", x = 0, y = 0.98, label = "PCR",
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 5, y = 0.98,
                      label = paste0("LFD", "\u2212"),
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 6, y = 0.92,
                      label = paste0("LFD", "\u2212"),
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 7, y = 0.98,
                      label = paste0("LFD", "\u2212"),
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    # Ribbons and lines — solid throughout
    ggplot2::geom_ribbon(
      data = panel_b_data,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi,
                   fill = scenario),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = panel_b_data,
      ggplot2::aes(x = day, y = p_infectious_med, color = scenario),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(values = pal_b, name = "Test history") +
    ggplot2::scale_fill_manual(values = pal_b, name = "Test history") +
    ggplot2::scale_y_continuous(
      "P(culture positive)",
      limits = c(0, 1), expand = c(0.01, 0),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(
      "Days since diagnosis",
      breaks = seq(0, 20, 5)
    ) +
    ggplot2::labs(tag = "B") +
    theme_journal(style, base_size = 13) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2),
      fill  = ggplot2::guide_legend(nrow = 2)
    ) +
    ggplot2::theme(legend.position = "bottom")

  # ── Panel C: Serial POSITIVE LFD testing ────────────────────────────────
  panel_c_labels <- c(
    "Diagnosis only (Ct 25)",
    "+ LFD+ day 5",
    "+ LFD+ days 5\u20136",
    "+ LFD+ days 5\u20137"
  )
  panel_c_data <- dplyr::bind_rows(
    sc$ct_25$curves       |> dplyr::mutate(scenario = panel_c_labels[1]),
    sc$lfd_pos_d5$curves  |> dplyr::mutate(scenario = panel_c_labels[2]),
    sc$lfd_pos_d56$curves |> dplyr::mutate(scenario = panel_c_labels[3]),
    sc$lfd_pos_d567$curves |> dplyr::mutate(scenario = panel_c_labels[4])
  ) |>
    dplyr::mutate(scenario = factor(scenario, levels = panel_c_labels))

  pal_c <- c(
    "Diagnosis only (Ct 25)"     = "#d73027",
    "+ LFD+ day 5"              = "#fc8d59",
    "+ LFD+ days 5\u20136"      = "#fee090",
    "+ LFD+ days 5\u20137"      = "#4575b4"
  )

  p_c <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pop,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi),
      fill = "grey80", alpha = 0.4
    ) +
    ggplot2::geom_line(
      data = pop,
      ggplot2::aes(x = day, y = p_infectious_med),
      color = "grey50", linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = c(0, 5, 6, 7),
      color = "grey40", linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::annotate("text", x = 0, y = 0.98, label = "PCR",
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 5, y = 0.98,
                      label = "LFD+",
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 6, y = 0.92,
                      label = "LFD+",
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::annotate("text", x = 7, y = 0.98,
                      label = "LFD+",
                      size = 3.5, color = "grey30", hjust = 0.5,
                      fontface = "italic") +
    ggplot2::geom_ribbon(
      data = panel_c_data,
      ggplot2::aes(x = day, ymin = p_infectious_lo, ymax = p_infectious_hi,
                   fill = scenario),
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      data = panel_c_data,
      ggplot2::aes(x = day, y = p_infectious_med, color = scenario),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(values = pal_c, name = "Test history") +
    ggplot2::scale_fill_manual(values = pal_c, name = "Test history") +
    ggplot2::scale_y_continuous(
      "P(culture positive)",
      limits = c(0, 1), expand = c(0.01, 0),
      labels = scales::percent
    ) +
    ggplot2::scale_x_continuous(
      "Days since diagnosis",
      breaks = seq(0, 20, 5)
    ) +
    ggplot2::labs(tag = "C") +
    theme_journal(style, base_size = 13) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2),
      fill  = ggplot2::guide_legend(nrow = 2)
    ) +
    ggplot2::theme(legend.position = "bottom")

  # ── Combine ──────────────────────────────────────────────────────────────
  patchwork::wrap_plots(p_a, p_b, p_c, ncol = 3)
}


#' Generate Bayesian filtering figure (targets entry point)
#'
#' @param filter_results  Output of \code{compute_bayesian_filtering()}
#' @param styles          Journal styles to generate
#' @return Character vector of saved file paths
generate_filtering_figures <- function(filter_results,
                                       styles = c("pnas", "plos", "annals")) {
  all_paths <- character(0)
  for (s in styles) {
    set_journal(s)
    p <- fig_bayesian_filtering(filter_results, style = s)
    paths <- save_journal_figure(p, "fig9_bayesian_filtering",
                                  layout = "full",
                                  width = 18, height = 6.5, style = s)
    all_paths <- c(all_paths, paths)
  }
  all_paths
}


#' Figure 8: Conditional culture positivity by symptom onset,
#' stratified by covariate profile
#'
#' Three-column faceted figure (Reference, Boosted, Reinfection)
#' showing P(culture+ | LFD+), P(culture+ | LFD-), and
#' P(culture+ | 2 consecutive LFD-) as a function of days since
#' symptom onset.
#'
#' @param policy_results  Output of \code{compute_all_policy()}
#' @param style           Journal style (NULL = use active)
#' @return patchwork ggplot object
fig_conditional_symptom <- function(policy_results, style = NULL) {

  cols <- journal_colors(style)

  cond <- policy_results$conditional_sym_ext
  if (is.null(cond)) {
    stop("conditional_sym_ext not found in policy_results. ",
         "Re-run compute_all_policy().")
  }

  # Filter to three profiles of interest
  profiles <- c("Reinfection", "Reference", "Boosted + Omicron (2022)")
  d <- cond |>
    dplyr::filter(label %in% profiles) |>
    dplyr::mutate(label = factor(label, levels = profiles))

  # Pivot to long format for the three conditions
  d_long <- rbind(
    data.frame(
      day   = d$day_since_landmark,
      label = d$label,
      prob  = d$p_cult_given_lfd_pos_med,
      lo    = d$p_cult_given_lfd_pos_lo,
      hi    = d$p_cult_given_lfd_pos_hi,
      condition = "LFD positive"
    ),
    data.frame(
      day   = d$day_since_landmark,
      label = d$label,
      prob  = d$p_cult_given_lfd_neg_med,
      lo    = d$p_cult_given_lfd_neg_lo,
      hi    = d$p_cult_given_lfd_neg_hi,
      condition = "Single LFD negative"
    ),
    data.frame(
      day   = d$day_since_landmark,
      label = d$label,
      prob  = d$p_cult_given_2neg_med,
      lo    = d$p_cult_given_2neg_lo,
      hi    = d$p_cult_given_2neg_hi,
      condition = "Two consecutive LFD negatives"
    )
  )

  d_long$condition <- factor(d_long$condition,
    levels = c("LFD positive", "Single LFD negative",
               "Two consecutive LFD negatives"))

  # Remove NA rows (day 0 for 2-consecutive where lag is unavailable)
  d_long <- d_long[!is.na(d_long$prob), ]

  # Colors
  cond_colors <- c(
    "LFD positive"                 = cols[["pfu"]],
    "Single LFD negative"          = cols[["lfd"]],
    "Two consecutive LFD negatives" = cols[["rna"]]
  )

  dodge <- ggplot2::position_dodge(width = 0.5)

  ggplot2::ggplot(d_long,
    ggplot2::aes(x = day, y = prob, color = condition)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi),
      position = dodge, linewidth = 0.4, alpha = 0.35) +
    ggplot2::geom_point(position = dodge, size = 1.3) +
    ggplot2::geom_vline(xintercept = c(5, 10), linetype = "dashed",
      color = "grey60", linewidth = 0.3) +
    ggplot2::facet_wrap(~ label, nrow = 1) +
    ggplot2::scale_color_manual(values = cond_colors, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
      labels = scales::percent) +
    ggplot2::scale_x_continuous(breaks = seq(0, 25, 5)) +
    ggplot2::labs(
      x = "Days since symptom onset",
      y = "P(culture positive)"
    ) +
    theme_journal(style, base_size = 10) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", size = 10)
    )
}
