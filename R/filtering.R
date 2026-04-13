# ──────────────────────────────────────────────────────────────────────────────
# filtering.R — Bayesian filtering: update infectiousness given test history
#
# Implements importance-sampling based Bayesian updating of viral trajectory
# probabilities conditioned on observed test results (PCR Ct values and/or
# LFD results).  Produces personalized infectiousness estimates that account
# for an individual's specific diagnostic history, in contrast to the
# population-averaged probability curves in policy.R.
#
# Approach:
#   1. Generate a large pool of latent trajectories from the population prior
#      (via sample_trajectories with include_noise = FALSE).
#   2. For each trajectory, marginalize over candidate "diagnosis times"
#      (all days when RNA > LOD, i.e., when a screening test would detect
#      the infection).
#   3. At each candidate diagnosis time, compute the likelihood of the
#      observed test history given the trajectory.
#   4. Self-normalize importance weights; compute weighted P(infectious)
#      at each clinical day since diagnosis.
#   5. Propagate posterior uncertainty by computing per-draw probability
#      estimates, then taking quantiles across draws.
# ──────────────────────────────────────────────────────────────────────────────


# ── Utility ──────────────────────────────────────────────────────────────────

#' Numerically stable log-sum-exp
#' @param x Numeric vector (may contain -Inf)
#' @return Scalar: log(sum(exp(x)))
#' @keywords internal
.log_sum_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}


# ── Default scenarios ────────────────────────────────────────────────────────

#' Define default Bayesian filtering scenarios
#'
#' Creates six scenarios in two groups for the manuscript figure:
#' \describe{
#'   \item{Panel A — Initial viral load}{Three scenarios varying the Ct value
#'     at diagnosis (Ct 20, 25, 30), each with a single RNA observation on
#'     day 0.}
#'   \item{Panel B — Serial LFD testing}{Three scenarios showing cumulative
#'     updating from a Ct 25 diagnosis followed by negative LFD results on
#'     days 5, 6, and 7.  The Ct-25–only scenario from Panel A serves as
#'     the baseline for comparison.}
#' }
#'
#' @return Named list of scenarios, each a list with \code{label} (character)
#'   and \code{test_history} (data.frame with columns \code{day}, \code{type},
#'   \code{value}).
default_filter_scenarios <- function() {
  rna_20 <- ct_to_rna(20, "nba")
  rna_25 <- ct_to_rna(25, "nba")
  rna_30 <- ct_to_rna(30, "nba")

  list(
    # ── Panel A: effect of initial viral load ───────────────────────────────
    ct_20 = list(
      label = "Ct = 20",
      test_history = data.frame(day = 0, type = "rna", value = rna_20,
                                stringsAsFactors = FALSE)
    ),
    ct_25 = list(
      label = "Ct = 25",
      test_history = data.frame(day = 0, type = "rna", value = rna_25,
                                stringsAsFactors = FALSE)
    ),
    ct_30 = list(
      label = "Ct = 30",
      test_history = data.frame(day = 0, type = "rna", value = rna_30,
                                stringsAsFactors = FALSE)
    ),
    # ── Panel B: serial LFD testing after Ct 25 diagnosis ───────────────────
    lfd_neg_d5 = list(
      label = "Ct 25 + LFD\u2212 day 5",
      test_history = data.frame(
        day   = c(0, 5),
        type  = c("rna", "lfd"),
        value = c(rna_25, 0),
        stringsAsFactors = FALSE
      )
    ),
    lfd_neg_d56 = list(
      label = "Ct 25 + LFD\u2212 days 5\u20136",
      test_history = data.frame(
        day   = c(0, 5, 6),
        type  = c("rna", "lfd", "lfd"),
        value = c(rna_25, 0, 0),
        stringsAsFactors = FALSE
      )
    ),
    lfd_neg_d567 = list(
      label = "Ct 25 + LFD\u2212 days 5\u20137",
      test_history = data.frame(
        day   = c(0, 5, 6, 7),
        type  = c("rna", "lfd", "lfd", "lfd"),
        value = c(rna_25, 0, 0, 0),
        stringsAsFactors = FALSE
      )
    ),
    # ── Panel C: serial POSITIVE LFD testing after Ct 25 diagnosis ─────────
    lfd_pos_d5 = list(
      label = "Ct 25 + LFD+ day 5",
      test_history = data.frame(
        day   = c(0, 5),
        type  = c("rna", "lfd"),
        value = c(rna_25, 1),
        stringsAsFactors = FALSE
      )
    ),
    lfd_pos_d56 = list(
      label = "Ct 25 + LFD+ days 5\u20136",
      test_history = data.frame(
        day   = c(0, 5, 6),
        type  = c("rna", "lfd", "lfd"),
        value = c(rna_25, 1, 1),
        stringsAsFactors = FALSE
      )
    ),
    lfd_pos_d567 = list(
      label = "Ct 25 + LFD+ days 5\u20137",
      test_history = data.frame(
        day   = c(0, 5, 6, 7),
        type  = c("rna", "lfd", "lfd", "lfd"),
        value = c(rna_25, 1, 1, 1),
        stringsAsFactors = FALSE
      )
    )
  )
}


# ── Core filtering ───────────────────────────────────────────────────────────

#' Filter trajectories by importance-weighting on a test history
#'
#' For each latent trajectory (agent \eqn{\times} draw), marginalizes over
#' candidate diagnosis times (all days when RNA exceeds the LOD) and weights
#' by the likelihood of the observed test results.  Returns per-draw
#' weighted probability curves with across-draw uncertainty.
#'
#' @param traj  Tibble from \code{sample_trajectories(include_noise = FALSE)}.
#'   Must have columns \code{agent_id}, \code{draw}, \code{day},
#'   \code{log_rna}, \code{log_pfu}, \code{lfd} (probability).
#' @param test_history  data.frame with columns:
#'   \describe{
#'     \item{day}{Integer day relative to diagnosis (0 = diagnosis)}
#'     \item{type}{"rna" or "lfd"}
#'     \item{value}{Log-RNA for "rna", 0/1 for "lfd"}
#'   }
#' @param sigma_rna_vec  Numeric vector of RNA observation SDs, one per
#'   posterior draw (indexed by draw number)
#' @param lod_rna  RNA limit of detection (natural-log scale)
#' @param lod_pfu  PFU limit of detection (natural-log scale)
#' @param max_day  Maximum clinical day to report (default 20)
#' @return List with:
#'   \describe{
#'     \item{curves}{Tibble with \code{day}, \code{p_infectious_med/lo/hi},
#'       \code{p_lfd_med/lo/hi}}
#'     \item{ess}{Per-draw effective sample size}
#'     \item{ess_median}{Median ESS across draws}
#'   }
filter_trajectories <- function(traj, test_history, sigma_rna_vec,
                                lod_rna, lod_pfu, max_day = 20) {

  draws   <- sort(unique(traj$draw))
  n_draws <- length(draws)
  clinical_days <- 0L:max_day
  n_days  <- length(clinical_days)
  has_obs <- nrow(test_history) > 0

  # Pre-split by draw for efficient subsetting
  traj_by_draw <- split(traj, traj$draw)

  # Per-draw weighted probability estimates
  draw_p_infectious <- matrix(NA_real_, n_draws, n_days)
  draw_p_lfd        <- matrix(NA_real_, n_draws, n_days)
  draw_ess          <- numeric(n_draws)

  for (di in seq_along(draws)) {
    d      <- draws[di]
    s_rna  <- sigma_rna_vec[di]
    traj_d <- traj_by_draw[[as.character(d)]]
    if (is.null(traj_d) || nrow(traj_d) == 0L) next

    agents_d <- unique(traj_d$agent_id)

    # Accumulators for weighted sums within this draw
    w_total      <- numeric(n_days)
    w_infectious <- numeric(n_days)
    w_lfd        <- numeric(n_days)
    agent_log_w  <- rep(-Inf, length(agents_d))

    for (ai in seq_along(agents_d)) {
      a  <- agents_d[ai]
      tr <- traj_d[traj_d$agent_id == a, ]
      tr <- tr[order(tr$day), ]
      n_t <- nrow(tr)
      if (n_t < 2L) next

      # Time step (days between consecutive trajectory points)
      dt_traj <- tr$day[2] - tr$day[1]

      # Candidate alignment indices: days where RNA is detectable
      align_idx <- which(tr$log_rna > lod_rna)
      n_align   <- length(align_idx)
      if (n_align < 1L) next

      # Index offsets for clinical days (in trajectory index space)
      day_offset <- round(clinical_days / dt_traj)

      # ── Compute log-likelihood for each candidate alignment ──────────────
      log_liks <- rep(-Inf, n_align)

      for (ki in seq_len(n_align)) {
        k <- align_idx[ki]

        if (!has_obs) {
          log_liks[ki] <- 0   # no observations → uniform weight
          next
        }

        ll    <- 0
        valid <- TRUE

        for (oi in seq_len(nrow(test_history))) {
          obs_idx <- k + round(test_history$day[oi] / dt_traj)
          if (obs_idx < 1L || obs_idx > n_t) { valid <- FALSE; break }

          ll <- ll + switch(test_history$type[oi],
            rna = stats::dnorm(test_history$value[oi], tr$log_rna[obs_idx],
                               s_rna, log = TRUE),
            lfd = {
              p <- min(max(tr$lfd[obs_idx], 1e-10), 1 - 1e-10)
              stats::dbinom(as.integer(test_history$value[oi]), 1L, p,
                            log = TRUE)
            },
            0  # unknown test types contribute nothing
          )
        }

        log_liks[ki] <- if (valid) ll else -Inf
      }

      # Marginalize over alignment: log_w = logSumExp(ll) - log(n_align)
      log_w <- .log_sum_exp(log_liks) - log(n_align)
      agent_log_w[ai] <- log_w
      if (!is.finite(log_w)) next

      # ── Accumulate weighted outcomes across alignment points ─────────────
      for (ki in seq_len(n_align)) {
        if (!is.finite(log_liks[ki])) next
        w_ali <- exp(log_liks[ki]) / n_align

        # Vectorized clinical-day lookup
        target_idx <- align_idx[ki] + day_offset
        in_range   <- target_idx >= 1L & target_idx <= n_t

        if (any(in_range)) {
          tidx <- target_idx[in_range]
          w_total[in_range]      <- w_total[in_range]      + w_ali
          w_infectious[in_range] <- w_infectious[in_range] +
            w_ali * (tr$log_pfu[tidx] > lod_pfu)
          w_lfd[in_range]        <- w_lfd[in_range]        +
            w_ali * tr$lfd[tidx]
        }
      }
    }

    # Within-draw probability
    nonzero <- w_total > 0
    if (any(nonzero)) {
      draw_p_infectious[di, nonzero] <- w_infectious[nonzero] / w_total[nonzero]
      draw_p_lfd[di, nonzero]        <- w_lfd[nonzero]        / w_total[nonzero]
    }

    # Effective sample size for this draw
    finite_lw <- agent_log_w[is.finite(agent_log_w)]
    if (length(finite_lw) > 1L) {
      nw <- exp(finite_lw - max(finite_lw))
      nw <- nw / sum(nw)
      draw_ess[di] <- 1 / sum(nw^2)
    } else {
      draw_ess[di] <- length(finite_lw)
    }
  }

  # ── Across-draw summaries ──────────────────────────────────────────────────
  .qfun <- function(mat, prob) {
    apply(mat, 2, stats::quantile, probs = prob, na.rm = TRUE)
  }

  curves <- tibble::tibble(
    day              = clinical_days,
    p_infectious_med = apply(draw_p_infectious, 2, stats::median, na.rm = TRUE),
    p_infectious_lo  = .qfun(draw_p_infectious, 0.025),
    p_infectious_hi  = .qfun(draw_p_infectious, 0.975),
    p_lfd_med        = apply(draw_p_lfd, 2, stats::median, na.rm = TRUE),
    p_lfd_lo         = .qfun(draw_p_lfd, 0.025),
    p_lfd_hi         = .qfun(draw_p_lfd, 0.975)
  )

  list(
    curves     = curves,
    ess        = draw_ess,
    ess_median = stats::median(draw_ess, na.rm = TRUE)
  )
}


# ── Population baseline ─────────────────────────────────────────────────────

#' Population-averaged P(infectious) aligned to first detection
#'
#' Computes the unfiltered population-averaged probability of culture
#' positivity at each day since first PCR detection.  Serves as the
#' baseline reference for comparing filtered (conditioned) curves.
#'
#' @param traj    Tibble from \code{sample_trajectories(include_noise = FALSE)}
#' @param lod_rna RNA LOD (natural-log scale)
#' @param lod_pfu PFU LOD (natural-log scale)
#' @param max_day Maximum day to report
#' @return Tibble with \code{day}, \code{p_infectious_med/lo/hi},
#'   \code{p_lfd_med/lo/hi}
compute_population_baseline <- function(traj, lod_rna, lod_pfu, max_day = 20) {

  traj <- dplyr::filter(traj, day_since_detectable <= max_day)

  # Pad missing days (short trajectories) with "not infectious"
  traj <- traj |>
    tidyr::complete(
      tidyr::nesting(agent_id, draw),
      day_since_detectable = 0:max_day,
      fill = list(log_pfu = -50, lfd = 0, log_rna = -50)
    )

  # Within each draw: average across agents
  draw_level <- traj |>
    dplyr::group_by(draw, day_since_detectable) |>
    dplyr::summarise(
      p_infectious = mean(log_pfu > lod_pfu),
      p_lfd        = mean(lfd),
      .groups      = "drop"
    )

  # Across draws: quantiles
  draw_level |>
    dplyr::group_by(day_since_detectable) |>
    dplyr::summarise(
      day              = day_since_detectable[1],
      p_infectious_med = stats::median(p_infectious),
      p_infectious_lo  = stats::quantile(p_infectious, 0.025),
      p_infectious_hi  = stats::quantile(p_infectious, 0.975),
      p_lfd_med        = stats::median(p_lfd),
      p_lfd_lo         = stats::quantile(p_lfd, 0.025),
      p_lfd_hi         = stats::quantile(p_lfd, 0.975),
      .groups          = "drop"
    ) |>
    dplyr::filter(day <= max_day) |>
    dplyr::select(-day_since_detectable)
}


# ── High-level wrapper ───────────────────────────────────────────────────────

#' Run Bayesian filtering analysis with default or custom scenarios
#'
#' Generates a large pool of latent trajectories from the population prior,
#' then runs importance-sampling filtering for each test-history scenario.
#' Returns results suitable for \code{fig_bayesian_filtering()}.
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list (from \code{build_stan_data()})
#' @param scenarios Named list of scenarios (NULL for defaults).  Each
#'   element is a list with \code{label} and \code{test_history}.
#' @param agent     Single-row data.frame with covariate profile (NULL for
#'   default: Omicron, boosted, age 30–49, first infection)
#' @param n_draws   Number of posterior draws (default 200)
#' @param n_reps    Replicate individuals per draw (default 500)
#' @param max_day   Maximum clinical day to report (default 20)
#' @param seed      Random seed for reproducibility
#' @return Named list with:
#'   \describe{
#'     \item{scenarios}{Per-scenario filtered results (curves, ESS)}
#'     \item{population}{Unfiltered population baseline curves}
#'     \item{agent}{Covariate profile used}
#'     \item{n_draws, n_reps, n_traj}{Sampling metadata}
#'   }
compute_bayesian_filtering <- function(fit, stan_data,
                                       scenarios = NULL,
                                       agent     = NULL,
                                       n_draws   = 200,
                                       n_reps    = 500,
                                       max_day   = 20,
                                       seed      = 2026) {

  # ── Remap CSV paths if needed ──────────────────────────────────────────────
  fit <- .ensure_local_fit(fit)

  # ── Default agent: Omicron, boosted, 30–49, first infection ────────────────
  if (is.null(agent)) {
    agent <- data.frame(
      age_group       = "30-49",
      variant         = "omicron",
      vaccination     = "boosted",
      prior_infection = FALSE,
      stringsAsFactors = FALSE
    )
  }

  # ── Replicate agent for n_reps independent individuals ─────────────────────
  agents <- agent[rep(1L, n_reps), ]
  rownames(agents) <- NULL

  # ── LODs ───────────────────────────────────────────────────────────────────
  lod_rna <- ct_to_rna(40, type = "nba") + 0.01
  lod_pfu <- 2.3

  # ── Default scenarios ──────────────────────────────────────────────────────
  if (is.null(scenarios)) scenarios <- default_filter_scenarios()

  # ── Extract sigma_rna per draw (same thinning as sample_trajectories) ──────
  #
  # sample_trajectories() calls set.seed(seed) then sort(sample(n_total, n_draws))

  # as its first RNG operation.  We replicate this to get the same draw indices.
  sigma_drws <- posterior::as_draws_matrix(
    fit$draws(variables = "sigma_rna")
  )
  n_total <- nrow(sigma_drws)
  if (n_total > n_draws) {
    set.seed(seed)
    sigma_idx  <- sort(sample(n_total, n_draws))
    sigma_drws <- sigma_drws[sigma_idx, ]
  }
  sigma_rna_vec <- as.numeric(sigma_drws[, "sigma_rna"])

  # ── Generate latent trajectories ───────────────────────────────────────────
  message(sprintf(
    "Bayesian filtering: %d agents \u00d7 %d draws = %s trajectories",
    n_reps, n_draws, format(n_reps * n_draws, big.mark = ",")
  ))

  traj <- sample_trajectories(
    fit, stan_data, agents,
    n_draws       = n_draws,
    dt            = 1,
    t_range       = c(-20, 50),
    include_noise = FALSE,
    seed          = seed
  )
  message(sprintf("  Generated %s trajectory-days",
                  format(nrow(traj), big.mark = ",")))

  # ── Population baseline (unfiltered) ───────────────────────────────────────
  pop <- compute_population_baseline(traj, lod_rna, lod_pfu, max_day)
  message(sprintf("  Population baseline computed (max day %d)", max_day))

  # ── Filter each scenario ───────────────────────────────────────────────────
  results <- lapply(names(scenarios), function(sc_name) {
    sc <- scenarios[[sc_name]]
    message(sprintf("  Filtering: %s (%d observations)",
                    sc$label, nrow(sc$test_history)))

    res <- filter_trajectories(
      traj          = traj,
      test_history  = sc$test_history,
      sigma_rna_vec = sigma_rna_vec,
      lod_rna       = lod_rna,
      lod_pfu       = lod_pfu,
      max_day       = max_day
    )
    res$label <- sc$label
    res
  })
  names(results) <- names(scenarios)

  # Report ESS
  for (nm in names(results)) {
    message(sprintf("    %s: median ESS = %.1f",
                    results[[nm]]$label, results[[nm]]$ess_median))
  }

  list(
    scenarios  = results,
    population = pop,
    agent      = agent,
    n_draws    = n_draws,
    n_reps     = n_reps,
    n_traj     = n_reps * n_draws
  )
}
