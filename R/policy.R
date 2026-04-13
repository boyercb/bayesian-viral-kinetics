# ──────────────────────────────────────────────────────────────────────────────
# policy.R — Policy-relevant derived quantities from the joint posterior
#
# All functions take the output of sample_trajectories() (with
# include_noise = TRUE) and compute aggregated probability curves,
# isolation duration tables, and related quantities.
# ──────────────────────────────────────────────────────────────────────────────


#' Ensure CmdStanMCMC fit can read its CSV files from the current machine
#'
#' If the fit was created on another machine (e.g., Dropbox path differs),
#' the stored output_files() paths won't resolve.  This function detects
#' that case and reconstructs the fit from local copies of the same CSVs.
#'
#' @param fit CmdStanMCMC object
#' @return CmdStanMCMC object with accessible CSV paths
#' @keywords internal
.ensure_local_fit <- function(fit) {
  csv_paths <- fit$output_files()
  if (all(file.exists(csv_paths))) return(fit)

  # Try to find the CSVs under the local output/stan_csv/ directory
  basenames <- basename(csv_paths)
  local_dir <- "output/stan_csv"
  local_paths <- file.path(local_dir, basenames)

  if (all(file.exists(local_paths))) {
    message("  Remapping fit CSV paths to local directory: ", local_dir)
    return(cmdstanr::as_cmdstan_fit(sort(local_paths)))
  }

  stop("Cannot locate Stan CSV files. Looked in:\n  ",
       paste(csv_paths[1], collapse = "\n  "), "\n  and\n  ",
       paste(local_paths[1], collapse = "\n  "))
}


#' Build the agent grid for policy trajectory sampling
#'
#' Constructs a set of covariate profiles for generating policy-relevant
#' trajectory samples.  Uses a one-at-a-time design (varying each covariate
#' from a reference profile) plus policy-relevant combined profiles.
#'
#' @return data.frame with one row per agent and columns:
#'   label, age_group, variant, vaccination, prior_infection
build_policy_agents <- function() {

  # Reference profile: pre-Alpha, unvaccinated, age <30, naive
  ref <- data.frame(
    label           = "Reference",
    age_group       = "0-29",
    variant         = "prealpha",
    vaccination     = "unvaccinated",
    prior_infection = FALSE,
    stringsAsFactors = FALSE
  )

  one_at_a_time <- data.frame(
    label = c(
      "Boosted",
      "Unboosted vaccinated",
      "Delta",
      "Omicron",
      "BA.4/BA.5",
      "Age 50+",
      "Reinfection"
    ),
    age_group = c(
      "0-29", "0-29", "0-29", "0-29", "0-29", "50+", "0-29"
    ),
    variant = c(
      "prealpha", "prealpha", "delta", "omicron", "ba4ba5",
      "prealpha", "prealpha"
    ),
    vaccination = c(
      "boosted", "unboosted", "unvaccinated", "unvaccinated", "unvaccinated",
      "unvaccinated", "unvaccinated"
    ),
    prior_infection = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE
    ),
    stringsAsFactors = FALSE
  )

  # Policy-relevant combined profiles
  combined <- data.frame(
    label = c(
      "Boosted + Omicron (2022)",
      "Unboosted + Delta (2021)"
    ),
    age_group = c("30-49", "0-29"),
    variant = c("omicron", "delta"),
    vaccination = c("boosted", "unboosted"),
    prior_infection = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  rbind(ref, one_at_a_time, combined)
}


#' Re-index trajectories to a landmark event
#'
#' Takes the output of \code{sample_trajectories(include_noise = TRUE)} and
#' re-aligns each trajectory so that day 0 corresponds to a clinically
#' meaningful landmark event.
#'
#' @param traj  Tibble from \code{sample_trajectories()}
#' @param landmark  One of \code{"pcr"} (first day RNA detectable),
#'   \code{"lfd"} (first day LFD positive), \code{"symptom"} (first day
#'   symptomatic)
#' @return Tibble with added column \code{day_since_landmark}, filtered to
#'   days >= 0.  An attribute \code{"excluded_frac"} records the fraction of
#'   trajectories that never reached the landmark and were dropped.
align_to_landmark <- function(traj, landmark = c("pcr", "lfd", "symptom")) {
  landmark <- match.arg(landmark)

  # Identify the event column
  event_col <- switch(landmark,
    pcr     = "rna_detectable",
    lfd     = "lfd",
    symptom = "symptomatic"
  )

  # For each trajectory (agent_id × draw), find the first day the event occurs
  landmark_days <- traj |>
    dplyr::group_by(agent_id, draw) |>
    dplyr::summarise(
      landmark_day = {
        hits <- day[.data[[event_col]] == 1 | .data[[event_col]] == TRUE]
        if (length(hits) == 0) NA_real_ else min(hits)
      },
      .groups = "drop"
    )

  n_total <- nrow(landmark_days)
  n_excluded <- sum(is.na(landmark_days$landmark_day))

  # Join and re-index

  out <- traj |>
    dplyr::left_join(landmark_days, by = c("agent_id", "draw")) |>
    dplyr::filter(!is.na(landmark_day)) |>
    dplyr::mutate(day_since_landmark = day - landmark_day) |>
    dplyr::filter(day_since_landmark >= 0)

  attr(out, "excluded_frac") <- n_excluded / n_total
  attr(out, "landmark") <- landmark
  out
}


#' Compute probability curves over time since a landmark
#'
#' For each day since the landmark, computes the probability of culture
#' positivity (PFU detectable), LFD positivity, and RNA detectability.
#' Uncertainty is propagated via within-draw aggregation then across-draw
#' quantiles.
#'
#' @param aligned  Output of \code{align_to_landmark()}
#' @param by       Character vector of column names to stratify by
#'   (e.g., \code{c("vaccination")}).  NULL for marginal curves.
#' @param max_day  Maximum day since landmark to include (default 25)
#' @return Tibble with columns: day_since_landmark, (by cols),
#'   p_culture_med, p_culture_lo, p_culture_hi,
#'   p_lfd_med, p_lfd_lo, p_lfd_hi,
#'   p_rna_med, p_rna_lo, p_rna_hi, n_traj
compute_probability_curves <- function(aligned, by = NULL, max_day = 25) {

  aligned <- dplyr::filter(aligned, day_since_landmark <= max_day)

  # Fix survivorship bias: agents whose detection window ends before max_day

  # would otherwise be absent at later days (inflating proportions).  Pad
  # the grid so every agent × draw has an entry at every day, with missing
  # days filled as undetectable (0).
  id_cols <- c("agent_id", "draw")
  if (!is.null(by)) id_cols <- c(id_cols, by)

  aligned <- aligned |>
    tidyr::complete(
      tidyr::nesting(!!!rlang::syms(id_cols)),
      day_since_landmark = seq(0L, max_day),
      fill = list(pfu_detectable = FALSE, lfd = 0L, rna_detectable = FALSE)
    )

  group_vars <- c("draw", "day_since_landmark")
  if (!is.null(by)) group_vars <- c(group_vars, by)

  # Step 1: within each draw, average across agents/trajectories
  draw_level <- aligned |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      p_culture = mean(as.numeric(pfu_detectable)),
      p_lfd     = mean(as.numeric(lfd)),
      p_rna     = mean(as.numeric(rna_detectable)),
      n_traj    = dplyr::n(),
      .groups   = "drop"
    )

  # Step 2: across draws, get median and 95% CrI
  summary_vars <- c("day_since_landmark")
  if (!is.null(by)) summary_vars <- c(summary_vars, by)

  draw_level |>
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
    dplyr::summarise(
      p_culture_med = stats::median(p_culture),
      p_culture_lo  = stats::quantile(p_culture, 0.025),
      p_culture_hi  = stats::quantile(p_culture, 0.975),
      p_lfd_med     = stats::median(p_lfd),
      p_lfd_lo      = stats::quantile(p_lfd, 0.025),
      p_lfd_hi      = stats::quantile(p_lfd, 0.975),
      p_rna_med     = stats::median(p_rna),
      p_rna_lo      = stats::quantile(p_rna, 0.025),
      p_rna_hi      = stats::quantile(p_rna, 0.975),
      n_traj        = stats::median(n_traj),
      .groups       = "drop"
    )
}


#' Compute conditional probability of culture positivity given LFD result
#'
#' At each day since the landmark, computes P(culture+ | LFD+) and
#' P(culture+ | LFD-) from the joint trajectory draws.
#'
#' @param aligned  Output of \code{align_to_landmark()}
#' @param by       Stratification columns (NULL for marginal)
#' @param max_day  Maximum day since landmark
#' @return Tibble with columns: day_since_landmark, (by cols),
#'   p_cult_given_lfd_pos_med/lo/hi, p_cult_given_lfd_neg_med/lo/hi,
#'   n_lfd_pos, n_lfd_neg
compute_conditional_curves <- function(aligned, by = NULL, max_day = 25) {

  aligned <- dplyr::filter(aligned, day_since_landmark <= max_day)

  group_vars <- c("draw", "day_since_landmark")
  if (!is.null(by)) group_vars <- c(group_vars, by)

  # Within each draw: direct proportions from joint samples
  draw_level <- aligned |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      # P(culture+ | LFD+)
      n_lfd_pos       = sum(lfd == 1),
      n_cult_lfd_pos  = sum(pfu_detectable & lfd == 1),
      p_cult_lfd_pos  = ifelse(n_lfd_pos > 0,
                               n_cult_lfd_pos / n_lfd_pos, NA_real_),
      # P(culture+ | LFD-)
      n_lfd_neg       = sum(lfd == 0),
      n_cult_lfd_neg  = sum(pfu_detectable & lfd == 0),
      p_cult_lfd_neg  = ifelse(n_lfd_neg > 0,
                               n_cult_lfd_neg / n_lfd_neg, NA_real_),
      .groups = "drop"
    )

  summary_vars <- c("day_since_landmark")
  if (!is.null(by)) summary_vars <- c(summary_vars, by)

  draw_level |>
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
    dplyr::summarise(
      p_cult_given_lfd_pos_med = stats::median(p_cult_lfd_pos, na.rm = TRUE),
      p_cult_given_lfd_pos_lo  = stats::quantile(p_cult_lfd_pos, 0.025,
                                                  na.rm = TRUE),
      p_cult_given_lfd_pos_hi  = stats::quantile(p_cult_lfd_pos, 0.975,
                                                  na.rm = TRUE),
      p_cult_given_lfd_neg_med = stats::median(p_cult_lfd_neg, na.rm = TRUE),
      p_cult_given_lfd_neg_lo  = stats::quantile(p_cult_lfd_neg, 0.025,
                                                  na.rm = TRUE),
      p_cult_given_lfd_neg_hi  = stats::quantile(p_cult_lfd_neg, 0.975,
                                                  na.rm = TRUE),
      n_lfd_pos                = stats::median(n_lfd_pos),
      n_lfd_neg                = stats::median(n_lfd_neg),
      .groups = "drop"
    )
}


#' Compute conditional culture-positivity curves with consecutive LFD negatives
#'
#' Extends \code{compute_conditional_curves} by also computing
#' P(culture+ | two consecutive LFD-), i.e. LFD negative on both the
#' current day and the previous day.
#'
#' @param aligned  Output of \code{align_to_landmark()}
#' @param by       Stratification columns (NULL for marginal)
#' @param max_day  Maximum day since landmark
#' @return Tibble with the same columns as \code{compute_conditional_curves}
#'   plus p_cult_given_2neg_med/lo/hi and n_2neg.
compute_conditional_curves_extended <- function(aligned, by = NULL,
                                                 max_day = 25) {

  aligned <- dplyr::filter(aligned, day_since_landmark <= max_day)

  # Add lagged LFD: previous-day LFD result within each trajectory

  aligned <- aligned |>
    dplyr::arrange(agent_id, draw, day_since_landmark) |>
    dplyr::group_by(agent_id, draw) |>
    dplyr::mutate(lfd_prev = dplyr::lag(lfd, 1, default = NA_integer_)) |>
    dplyr::ungroup()

  # Flag: two consecutive LFD negatives (today = 0 AND yesterday = 0)
  aligned <- aligned |>
    dplyr::mutate(
      two_neg = (!is.na(lfd_prev) & lfd == 0 & lfd_prev == 0)
    )

  group_vars <- c("draw", "day_since_landmark")
  if (!is.null(by)) group_vars <- c(group_vars, by)

  draw_level <- aligned |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      # P(culture+ | LFD+)
      n_lfd_pos       = sum(lfd == 1),
      n_cult_lfd_pos  = sum(pfu_detectable & lfd == 1),
      p_cult_lfd_pos  = ifelse(n_lfd_pos > 0,
                               n_cult_lfd_pos / n_lfd_pos, NA_real_),
      # P(culture+ | LFD-)
      n_lfd_neg       = sum(lfd == 0),
      n_cult_lfd_neg  = sum(pfu_detectable & lfd == 0),
      p_cult_lfd_neg  = ifelse(n_lfd_neg > 0,
                               n_cult_lfd_neg / n_lfd_neg, NA_real_),
      # P(culture+ | 2 consecutive LFD-)
      n_2neg          = sum(two_neg),
      n_cult_2neg     = sum(pfu_detectable & two_neg),
      p_cult_2neg     = ifelse(n_2neg > 0,
                               n_cult_2neg / n_2neg, NA_real_),
      .groups = "drop"
    )

  summary_vars <- c("day_since_landmark")
  if (!is.null(by)) summary_vars <- c(summary_vars, by)

  draw_level |>
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
    dplyr::summarise(
      p_cult_given_lfd_pos_med = stats::median(p_cult_lfd_pos, na.rm = TRUE),
      p_cult_given_lfd_pos_lo  = stats::quantile(p_cult_lfd_pos, 0.025,
                                                  na.rm = TRUE),
      p_cult_given_lfd_pos_hi  = stats::quantile(p_cult_lfd_pos, 0.975,
                                                  na.rm = TRUE),
      p_cult_given_lfd_neg_med = stats::median(p_cult_lfd_neg, na.rm = TRUE),
      p_cult_given_lfd_neg_lo  = stats::quantile(p_cult_lfd_neg, 0.025,
                                                  na.rm = TRUE),
      p_cult_given_lfd_neg_hi  = stats::quantile(p_cult_lfd_neg, 0.975,
                                                  na.rm = TRUE),
      p_cult_given_2neg_med    = stats::median(p_cult_2neg, na.rm = TRUE),
      p_cult_given_2neg_lo     = stats::quantile(p_cult_2neg, 0.025,
                                                  na.rm = TRUE),
      p_cult_given_2neg_hi     = stats::quantile(p_cult_2neg, 0.975,
                                                  na.rm = TRUE),
      n_lfd_pos                = stats::median(n_lfd_pos),
      n_lfd_neg                = stats::median(n_lfd_neg),
      n_2neg                   = stats::median(n_2neg),
      .groups = "drop"
    )
}


#' Compute isolation duration table
#'
#' For each covariate profile, returns the probability of remaining
#' culture-positive at fixed isolation durations (days since landmark).
#'
#' @param aligned   Output of \code{align_to_landmark()}
#' @param days      Integer vector of isolation durations to evaluate
#' @param by        Stratification columns (NULL uses "label" if present)
#' @return Tibble with columns: (by cols), day, p_still_infectious_med/lo/hi
compute_isolation_table <- function(aligned,
                                    days = c(5, 7, 10),
                                    by = NULL) {

  if (is.null(by) && "label" %in% names(aligned)) by <- "label"

  # Filter to the target days
  at_days <- aligned |>
    dplyr::filter(day_since_landmark %in% days)

  group_vars <- c("draw", "day_since_landmark")
  if (!is.null(by)) group_vars <- c(group_vars, by)

  # Within draw: fraction still culture positive
  draw_level <- at_days |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      p_infectious = mean(pfu_detectable),
      .groups = "drop"
    )

  summary_vars <- c("day_since_landmark")
  if (!is.null(by)) summary_vars <- c(summary_vars, by)

  draw_level |>
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
    dplyr::summarise(
      p_still_infectious_med = stats::median(p_infectious),
      p_still_infectious_lo  = stats::quantile(p_infectious, 0.025),
      p_still_infectious_hi  = stats::quantile(p_infectious, 0.975),
      .groups = "drop"
    ) |>
    dplyr::rename(day = day_since_landmark) |>
    dplyr::arrange(day)
}


#' Compute expected residual infectious days after isolation release
#'
#' For each release day d, computes the expected number of remaining days
#' that an individual would still be shedding infectious virus (PFU above
#' detection threshold).
#'
#' @param aligned      Output of \code{align_to_landmark()}
#' @param release_days Integer vector of candidate release days
#' @param by           Stratification columns
#' @return Tibble with columns: release_day, (by cols),
#'   residual_days_med/lo/hi
compute_residual_auc <- function(aligned,
                                  release_days = 3:14,
                                  by = NULL) {

  if (is.null(by) && "label" %in% names(aligned)) by <- "label"
  has_by <- length(by) > 0

  # For each trajectory, compute remaining culture-positive days after release
  results_list <- lapply(release_days, function(d) {

    residual <- aligned |>
      dplyr::filter(day_since_landmark > d) |>
      dplyr::group_by(agent_id, draw) |>
      dplyr::summarise(
        remaining_days = sum(pfu_detectable),
        .groups = "drop"
      )

    # Join back agent characteristics for stratification
    if (has_by) {
      agent_info <- aligned |>
        dplyr::select(agent_id, draw, dplyr::all_of(by)) |>
        dplyr::distinct()
      residual <- dplyr::left_join(residual, agent_info,
                                    by = c("agent_id", "draw"))
    }

    group_vars <- "draw"
    if (has_by) group_vars <- c(group_vars, by)

    # Within draw: average across agents
    draw_level <- residual |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        mean_remaining = mean(remaining_days),
        .groups = "drop"
      )

    summary_vars <- if (has_by) by else character(0)

    out <- draw_level |>
      dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
      dplyr::summarise(
        residual_days_med = stats::median(mean_remaining),
        residual_days_lo  = stats::quantile(mean_remaining, 0.025),
        residual_days_hi  = stats::quantile(mean_remaining, 0.975),
        .groups = "drop"
      )

    out$release_day <- d
    out
  })

  dplyr::bind_rows(results_list)
}


#' Compute test-to-release false reassurance rate
#'
#' At each day since the landmark, among trajectories where LFD is negative
#' on that day, computes the fraction that are still culture-positive.
#' This gives P(infectious | negative LFD at day d) — the probability
#' that a test-based release decision would be premature.
#'
#' @param aligned  Output of \code{align_to_landmark()}
#' @param by       Stratification columns
#' @param max_day  Maximum day since landmark
#' @return Tibble with columns: day_since_landmark, (by cols),
#'   p_infectious_given_neg_lfd_med/lo/hi, n_negative_tests
compute_test_to_release <- function(aligned, by = NULL, max_day = 25) {

  aligned <- dplyr::filter(aligned, day_since_landmark <= max_day)

  group_vars <- c("draw", "day_since_landmark")
  if (!is.null(by)) group_vars <- c(group_vars, by)

  # Within each draw: among LFD- at that day, fraction still culture+
  draw_level <- aligned |>
    dplyr::filter(lfd == 0) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      p_false_reassurance = mean(pfu_detectable),
      n_neg               = dplyr::n(),
      .groups = "drop"
    )

  summary_vars <- c("day_since_landmark")
  if (!is.null(by)) summary_vars <- c(summary_vars, by)

  draw_level |>
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_vars))) |>
    dplyr::summarise(
      p_infectious_given_neg_lfd_med = stats::median(p_false_reassurance,
                                                      na.rm = TRUE),
      p_infectious_given_neg_lfd_lo  = stats::quantile(p_false_reassurance,
                                                        0.025, na.rm = TRUE),
      p_infectious_given_neg_lfd_hi  = stats::quantile(p_false_reassurance,
                                                        0.975, na.rm = TRUE),
      n_negative_tests               = stats::median(n_neg),
      .groups = "drop"
    )
}


#' Run all policy analyses from sampled trajectories
#'
#' Wrapper function that generates trajectories, aligns them to all three
#' landmarks, and computes all derived policy quantities.  Designed to be
#' called as a single \code{targets} target.
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list
#' @param n_draws   Posterior draws per agent (default 500)
#' @param n_reps    Number of replicate individuals per profile per draw
#'                  (default 50).  Higher values give smoother within-draw
#'                  probability estimates but increase compute time.
#' @param seed      Random seed
#' @return Named list with all policy analysis results
compute_all_policy <- function(fit, stan_data, n_draws = 500,
                                n_reps = 50, seed = 2026) {

  # ── Remap CSV paths if running on a different machine ───────────────────
  fit <- .ensure_local_fit(fit)

  agents_unique <- build_policy_agents()

  # Replicate each profile n_reps times so that within each posterior draw,

  # we get n_reps independent trajectories per profile (each with fresh
  # individual REs).  This makes the within-draw probability estimate
  # meaningful — a single binary outcome per draw gives meaningless CrI.
  agents <- agents_unique[rep(seq_len(nrow(agents_unique)), each = n_reps), ]
  agents$rep_id  <- rep(seq_len(n_reps), times = nrow(agents_unique))
  # Assign a unique agent_id that encodes both profile and replicate
  agents$profile <- rep(agents_unique$label, each = n_reps)
  rownames(agents) <- NULL

  message(sprintf(
    "Policy analysis: %d profiles x %d reps x %d draws = %d trajectories",
    nrow(agents_unique), n_reps, n_draws,
    nrow(agents) * n_draws))

  # ── Generate trajectories ───────────────────────────────────────────────
  traj <- sample_trajectories(
    fit, stan_data, agents,
    n_draws = n_draws, dt = 1,
    include_noise = TRUE, seed = seed
  )
  message(sprintf("  Generated %d trajectory-days", nrow(traj)))

  # ── Align to landmarks ─────────────────────────────────────────────────
  aligned_pcr <- align_to_landmark(traj, "pcr")
  aligned_sym <- align_to_landmark(traj, "symptom")
  aligned_lfd <- align_to_landmark(traj, "lfd")

  message(sprintf("  Excluded fractions: PCR=%.1f%%, symptom=%.1f%%, LFD=%.1f%%",
                  attr(aligned_pcr, "excluded_frac") * 100,
                  attr(aligned_sym, "excluded_frac") * 100,
                  attr(aligned_lfd, "excluded_frac") * 100))

  # ── Probability curves (marginal + stratified) ─────────────────────────
  curves_pcr_marginal <- compute_probability_curves(aligned_pcr)
  curves_sym_marginal <- compute_probability_curves(aligned_sym)

  curves_pcr_by_label <- compute_probability_curves(
    aligned_pcr, by = "label"
  )
  curves_pcr_by_vax <- compute_probability_curves(
    aligned_pcr, by = "vaccination"
  )
  curves_pcr_by_var <- compute_probability_curves(
    aligned_pcr, by = "variant"
  )

  # ── Conditional on LFD result (aligned to first positive LFD) ──────────
  cond_lfd <- compute_conditional_curves(aligned_lfd)

  # ── Conditional on LFD result (aligned to symptom onset, by profile) ────
  #    Extended version includes P(culture+ | 2 consecutive LFD-)
  cond_sym_by_label <- compute_conditional_curves_extended(
    aligned_sym, by = "label"
  )

  # ── Isolation duration table ────────────────────────────────────────────
  isolation_pcr <- compute_isolation_table(aligned_pcr)
  isolation_sym <- compute_isolation_table(aligned_sym)

  # ── Residual infectiousness AUC ─────────────────────────────────────────
  # Stratified by label (covariate profile)
  auc_pcr <- compute_residual_auc(aligned_pcr, by = "label")
  auc_sym <- compute_residual_auc(aligned_sym, by = "label")

  # Marginal: pool across all agent profiles (use character(0) to prevent
  # the automatic by="label" default)
  auc_pcr_marginal <- compute_residual_auc(aligned_pcr, by = character(0))
  auc_sym_marginal <- compute_residual_auc(aligned_sym, by = character(0))

  # ── Test-to-release ─────────────────────────────────────────────────────
  test_release_pcr <- compute_test_to_release(aligned_pcr)
  test_release_sym <- compute_test_to_release(aligned_sym)

  # ── Return everything ──────────────────────────────────────────────────
  list(
    agents         = agents,
    trajectories   = traj,
    aligned_pcr    = aligned_pcr,
    aligned_sym    = aligned_sym,
    aligned_lfd    = aligned_lfd,
    excluded_fracs = c(
      pcr     = attr(aligned_pcr, "excluded_frac"),
      symptom = attr(aligned_sym, "excluded_frac"),
      lfd     = attr(aligned_lfd, "excluded_frac")
    ),
    curves = list(
      pcr_marginal = curves_pcr_marginal,
      sym_marginal = curves_sym_marginal,
      pcr_by_label = curves_pcr_by_label,
      pcr_by_vax   = curves_pcr_by_vax,
      pcr_by_var   = curves_pcr_by_var
    ),
    conditional_lfd     = cond_lfd,
    conditional_sym_ext = cond_sym_by_label,
    isolation = list(
      pcr = isolation_pcr,
      sym = isolation_sym
    ),
    residual_auc = list(
      pcr_by_label = auc_pcr,
      sym_by_label = auc_sym,
      pcr_marginal = auc_pcr_marginal,
      sym_marginal = auc_sym_marginal
    ),
    test_release = list(
      pcr = test_release_pcr,
      sym = test_release_sym
    )
  )
}
