# ──────────────────────────────────────────────────────────────────────────────
# trajectories.R — Sample viral shedding trajectories for agent-based models
# ──────────────────────────────────────────────────────────────────────────────


#' Convert agent characteristics to model covariate matrix
#'
#' Maps human-readable agent attributes to the binary indicator matrix
#' expected by the kinetic model.  Reference categories (all zeros):
#' age \code{[0,30)}, first infection, pre-alpha variant,
#' unvaccinated/unboosted.
#'
#' @param agents  A \code{data.frame} with optional columns:
#'   \code{age_group}, \code{variant}, \code{vaccination},
#'   \code{prior_infection}.  Missing columns default to reference category.
#'   Extra columns (e.g., \code{label}) are silently ignored.
#' @return Numeric matrix (n_agents x P) of covariate indicators.
build_agent_covariates <- function(agents) {
  covariates <- define_covariates()
  n <- nrow(agents)
  P <- length(covariates$x_vars)
  x <- matrix(0L, nrow = n, ncol = P)
  colnames(x) <- covariates$x_vars

  # Age group (ref = [0,30))
  if ("age_group" %in% names(agents)) {
    age <- agents$age_group
    x[age %in% c("30-49", "[30,50)"), "age_[30,50)"]  <- 1L
    x[age %in% c("50+",   "[50,100)"), "age_[50,100)"] <- 1L
  }

  # Prior infection / recurrence (ref = first infection)
  if ("prior_infection" %in% names(agents)) {
    x[agents$prior_infection == TRUE, "recurrence"] <- 1L
  }

  # Variant (ref = prealpha)
  if ("variant" %in% names(agents)) {
    v <- tolower(agents$variant)
    for (vv in c("alpha", "delta", "omicron", "ba4ba5", "other")) {
      x[v == vv, vv] <- 1L
    }
  }

  # Vaccination history (ref = unvaccinated, unboosted)
  if ("vaccination" %in% names(agents)) {
    vax <- tolower(agents$vaccination)
    x[vax %in% c("boosted",  "vaccinated_boosted"),   "vaccinated_boosted"]    <- 1L
    x[vax %in% c("unboosted","vaccinated_unboosted"), "vaccinated_unboosted"]  <- 1L
    x[vax == "vaccinated_unreported",   "vaccinated_unreported"]   <- 1L
    x[vax == "unvaccinated_unreported", "unvaccinated_unreported"] <- 1L
    x[vax == "unreported_primary",      "unreported_primary"]      <- 1L
  }

  x
}


#' Sample viral shedding trajectories for new agents
#'
#' Generates realistic viral shedding trajectories by drawing from the
#' posterior predictive distribution.
#'
#' For each agent x posterior draw, **new** individual random effects are
#' sampled from the fitted population distribution—capturing true
#' inter-individual variability (not just posterior uncertainty for a
#' specific observed person).  This makes the output suitable for
#' initialising agents in an ABM with biologically realistic heterogeneity.
#'
#' @param fit        A \code{CmdStanMCMC} fit object (e.g., from
#'                   \code{targets::tar_read(kinetics_mcmc)}).
#' @param stan_data  Stan data list (from \code{build_stan_data()}).  Needed
#'                   for model flags, prior hyperparameters, and LOD values.
#' @param agents     A \code{data.frame} with one row per agent.  Optional
#'                   columns (missing columns default to the reference
#'                   category):
#'   \describe{
#'     \item{age_group}{\code{"0-29"}, \code{"30-49"}, or \code{"50+"}}
#'     \item{variant}{\code{"prealpha"}, \code{"alpha"}, \code{"delta"},
#'                    \code{"omicron"}, \code{"ba4ba5"}, or \code{"other"}}
#'     \item{vaccination}{\code{"unvaccinated"} (ref), \code{"boosted"},
#'                        \code{"unboosted"}, \code{"vaccinated_unreported"},
#'                        \code{"unvaccinated_unreported"},
#'                        \code{"unreported_primary"}}
#'     \item{prior_infection}{\code{TRUE} / \code{FALSE} (is this a
#'                            reinfection?)}
#'   }
#'   Additional columns (e.g., \code{label}) are carried through to the
#'   output for easy grouping/faceting.
#' @param n_draws    Number of posterior draws per agent (default 100).
#'                   Total trajectories generated =
#'                   \code{nrow(agents) * n_draws}.
#' @param dt         Time step in days (default 1).
#' @param t_range    Evaluation window in days relative to RNA peak.
#'                   Default \code{c(-20, 35)}.  Output is trimmed to the
#'                   detectable window (RNA or PFU above LOD).
#' @param lod_rna    RNA limit of detection (natural-log scale).  Default
#'                   uses Ct 40 on the NBA calibration curve
#'                   (~6.1 ln copies/mL).
#' @param lod_pfu    PFU limit of detection (natural-log scale).  Default
#'                   2.3 (~10 PFU/mL).
#' @param include_noise  If \code{TRUE} (default), RNA and PFU include
#'                   Gaussian measurement noise and LOD censoring; LFD and
#'                   symptoms are stochastic binary draws.  If \code{FALSE},
#'                   returns latent trajectories and exact probabilities.
#' @param seed       Random seed for reproducibility.
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{agent_id}{Row index into \code{agents}}
#'     \item{draw}{Posterior draw index}
#'     \item{day}{Days relative to RNA peak (day 0 = peak)}
#'     \item{day_since_detectable}{Days since RNA or PFU first exceeded LOD
#'           (starts at 0)}
#'     \item{log_rna}{Latent log RNA (ln copies/mL, true trajectory)}
#'     \item{log_pfu}{Latent log PFU/mL (true trajectory)}
#'     \item{rna_detectable}{Logical: is the true latent RNA above LOD?}
#'     \item{pfu_detectable}{Logical: is the true latent PFU above LOD?}
#'     \item{rna}{Measured RNA: with noise + LOD censoring when
#'           \code{include_noise = TRUE}; equals \code{log_rna} otherwise}
#'     \item{pfu}{Measured PFU: with noise + LOD censoring when
#'           \code{include_noise = TRUE}; equals \code{log_pfu} otherwise}
#'     \item{lfd}{LFD result: 0/1 when \code{include_noise = TRUE};
#'           P(positive) otherwise}
#'     \item{symptomatic}{1 from symptom onset onward when
#'           \code{include_noise = TRUE}; daily hazard otherwise}
#'   }
#'   Agent characteristic columns from \code{agents} are appended.
#'
#' @examples
#' \dontrun{
#' fit <- targets::tar_read(kinetics_mcmc)
#' sd  <- targets::tar_read(stan_data)
#'
#' agents <- data.frame(
#'   age_group       = c("30-49", "50+"),
#'   variant         = c("omicron", "delta"),
#'   vaccination     = c("boosted", "unvaccinated"),
#'   prior_infection = c(FALSE, TRUE)
#' )
#'
#' traj <- sample_trajectories(fit, sd, agents, n_draws = 50, seed = 1)
#' }
sample_trajectories <- function(
    fit,
    stan_data,
    agents,
    n_draws       = 100,
    dt            = 1,
    t_range       = c(-20, 35),
    lod_rna       = NULL,
    lod_pfu       = NULL,
    include_noise = TRUE,
    seed          = NULL) {

  if (!is.null(seed)) set.seed(seed)
  n_agents <- nrow(agents)

  # ── Defaults ──────────────────────────────────────────────────────────────
  if (is.null(lod_rna)) lod_rna <- ct_to_rna(40, type = "nba") + 0.01
  if (is.null(lod_pfu)) lod_pfu <- 2.3

  # ── Covariates ────────────────────────────────────────────────────────────
  x_mat <- build_agent_covariates(agents)

  # ── Time grid ─────────────────────────────────────────────────────────────
  t_grid <- seq(t_range[1], t_range[2], by = dt)

  # ── Model flags ───────────────────────────────────────────────────────────
  use_smooth  <- as.logical(stan_data$use_smooth)
  has_ind     <- as.logical(stan_data$ind_effects)
  has_corr    <- as.logical(stan_data$ind_corr)
  has_adj_rna <- as.logical(stan_data$adj_rna)
  has_adj_pfu <- as.logical(stan_data$adj_pfu)
  has_wf      <- as.logical(stan_data$use_wf)
  prior_i_sd  <- stan_data$prior_i_sd
  scale_vl    <- stan_data$prior_dp_mean

  # ── Variables to extract from posterior ────────────────────────────────────
  vars <- c(
    "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
    "tau_tp", "tau_dp", "tau_wp", "tau_wr",
    "tau0_lfd", "tau_lfd",
    "sigma_rna", "sigma_pfu",
    "eta_sym_intercept", "eta_sym_pfu", "eta_sym_rna", "sigma_sym"
  )
  if (has_corr)    vars <- c(vars, "sigma_ind_rna", "L_Omega_rna")
  if (has_adj_rna) vars <- c(vars, "beta_dp_rna", "beta_wp_rna", "beta_wr_rna")
  if (has_adj_pfu) vars <- c(vars, "beta_dp_pfu", "beta_wp_pfu", "beta_wr_pfu")
  if (has_wf)      vars <- c(vars, "wf_raw")

  # ── PFU individual-effect SDs ───────────────────────────────────────────────

  # Newer Stan model (mode 2):  sigma_resid_pfu is an estimated parameter that
  # captures residual PFU individual variation beyond what is propagated from
  # RNA via the log-affine transform.  If present, we use it directly.
  # Older fits:  PFU REs had a fixed prior_i_sd, so sigma_resid_pfu is missing.
  #   Fall back to extracting all individual PFU REs and computing empirical
  #   per-draw SDs (slower, but correct).
  all_param_names <- fit$metadata()$model_params
  has_sigma_resid_pfu <- "sigma_resid_pfu" %in% all_param_names

  # Pre-compute PFU-informed individual IDs for the fallback path
  pfu_informed <- NULL
  if (!has_sigma_resid_pfu) {
    pfu_informed <- which(tapply(
      stan_data$pfu_exist, stan_data$id, function(x) any(x == 1)
    ))
  }

  if (has_sigma_resid_pfu) {
    vars <- c(vars, "sigma_resid_pfu")
  } else {
    message("sigma_resid_pfu not found in posterior \u2014 using variance decomposition workaround")
    # PFU REs are extracted separately below (not added to vars)
  }

  # ── Extract and thin draws ────────────────────────────────────────────────
  drws <- posterior::as_draws_matrix(fit$draws(variables = vars))
  n_total <- nrow(drws)
  if (n_total > n_draws) {
    idx <- sort(sample(n_total, n_draws))
    drws <- drws[idx, ]
  }
  # Convert to plain matrix to avoid array-recycling deprecation warnings
  drws <- matrix(as.numeric(drws), nrow = nrow(drws),
                 dimnames = dimnames(drws))
  nd <- nrow(drws)  # actual draws used
  P  <- stan_data$P

  # ── Resolve PFU RE SDs per draw ───────────────────────────────────────────
  if (has_sigma_resid_pfu) {
    # Direct: sigma_resid_pfu[1..4] (or [1] when no ind_effects)
    sd_pfu <- list(
      tp = drws[, "sigma_resid_pfu[1]"],
      dp = if (has_ind) drws[, "sigma_resid_pfu[2]"] else rep(0, nd),
      wp = if (has_ind) drws[, "sigma_resid_pfu[3]"] else rep(0, nd),
      wr = if (has_ind) drws[, "sigma_resid_pfu[4]"] else rep(0, nd)
    )
  } else {
    # Fallback: estimate population SDs of PFU individual effects via
    # variance decomposition (needed until sigma_resid_pfu is fitted).
    #
    # Problem: with prior_i_sd=1 and sparse PFU data (~5-10 obs per person),
    # each individual's posterior RE is wide.  The per-draw empirical SD
    # across individuals conflates TRUE population variation with posterior
    # uncertainty:
    #     Var_between(draw d) = Var_pop + mean(Var_posterior_i)
    #
    # We correct this by subtracting the average within-individual posterior
    # variance (computed across MCMC draws) from the per-draw between-
    # individual variance:
    #     Var_pop ≈ max(0, Var_between - mean(Var_posterior_i))
    N_ind <- sum(stan_data$M)
    message(sprintf("  Using %d PFU-informed individuals (of %d total) for variance decomposition",
                    length(pfu_informed), N_ind))

    # We need ALL draws (not thinned) for accurate per-individual variance
    # estimation.  Extract the PFU RE columns from the full fit.
    pfu_re_nms <- list(
      tp = paste0("tp_i_pfu[", pfu_informed, "]"),
      dp = if (has_ind) paste0("dp_i_pfu[", pfu_informed, "]") else NULL,
      wp = if (has_ind) paste0("wp_i_pfu[", pfu_informed, "]") else NULL,
      wr = if (has_ind) paste0("wr_i_pfu[", pfu_informed, "]") else NULL
    )
    all_pfu_vars <- unlist(pfu_re_nms, use.names = FALSE)
    full_pfu_drws <- posterior::as_draws_matrix(fit$draws(variables = all_pfu_vars))
    full_pfu_drws <- matrix(as.numeric(full_pfu_drws), nrow = nrow(full_pfu_drws),
                            dimnames = dimnames(full_pfu_drws))

    .decompose_sd <- function(mat) {
      # Per-individual posterior variance (across draws)
      post_var_i <- apply(mat, 2, var)
      avg_post_var <- mean(post_var_i)
      # Per-draw between-individual variance
      between_var <- apply(mat, 1, var)
      # Population variance = between - average posterior uncertainty
      pop_var <- pmax(between_var - avg_post_var, 0)
      sqrt(pop_var)
    }

    sd_pfu <- list(
      tp = .decompose_sd(full_pfu_drws[, pfu_re_nms$tp, drop = FALSE]),
      dp = if (has_ind) .decompose_sd(full_pfu_drws[, pfu_re_nms$dp, drop = FALSE]) else rep(0, nrow(full_pfu_drws)),
      wp = if (has_ind) .decompose_sd(full_pfu_drws[, pfu_re_nms$wp, drop = FALSE]) else rep(0, nrow(full_pfu_drws)),
      wr = if (has_ind) .decompose_sd(full_pfu_drws[, pfu_re_nms$wr, drop = FALSE]) else rep(0, nrow(full_pfu_drws))
    )

    message(sprintf("  Variance-decomposed PFU RE SDs (mean): tp=%.3f dp=%.3f wp=%.3f wr=%.3f",
                    mean(sd_pfu$tp), mean(sd_pfu$dp), mean(sd_pfu$wp), mean(sd_pfu$wr)))

    # Thin sd_pfu to match the thinned draws
    if (n_total > n_draws) {
      sd_pfu <- lapply(sd_pfu, function(s) s[idx])
    }
    nd_check <- nrow(drws)
    stopifnot(length(sd_pfu$tp) == nd_check)
  }

  message(sprintf(
    "Sampling %d trajectories (%d agents \u00d7 %d draws)...",
    n_agents * nd, n_agents, nd
  ))

  # ── Main loop: draws x agents ─────────────────────────────────────────────
  results <- vector("list", nd * n_agents)
  ri <- 0L

  for (d in seq_len(nd)) {

    # ── Population parameters for draw d ────────────────────────────────────
    dp_pop <- drws[d, "dp_mean_rna"]
    wp_pop <- drws[d, "wp_mean_rna"]
    wr_pop <- drws[d, "wr_mean_rna"]

    ttp <- c(drws[d, "tau_tp[1]"], drws[d, "tau_tp[2]"])
    tdp <- c(drws[d, "tau_dp[1]"], drws[d, "tau_dp[2]"])
    twp <- c(drws[d, "tau_wp[1]"], drws[d, "tau_wp[2]"])
    twr <- c(drws[d, "tau_wr[1]"], drws[d, "tau_wr[2]"])

    lfd0  <- drws[d, "tau0_lfd"]
    lfd_c <- c(drws[d, "tau_lfd[1]"], drws[d, "tau_lfd[2]"])

    s_rna <- drws[d, "sigma_rna"]
    s_pfu <- drws[d, "sigma_pfu"]

    ei <- drws[d, "eta_sym_intercept"]
    ep <- drws[d, "eta_sym_pfu"]
    er <- drws[d, "eta_sym_rna"]
    ss <- drws[d, "sigma_sym"]

    # Covariate betas (RNA)
    if (has_adj_rna) {
      bdp_rna <- as.numeric(drws[d, paste0("beta_dp_rna[", 1:P, "]")])
      bwp_rna <- as.numeric(drws[d, paste0("beta_wp_rna[", 1:P, "]")])
      bwr_rna <- as.numeric(drws[d, paste0("beta_wr_rna[", 1:P, "]")])
    }

    # Covariate betas (PFU)
    if (has_adj_pfu) {
      bdp_pfu <- as.numeric(drws[d, paste0("beta_dp_pfu[", 1:P, "]")])
      bwp_pfu <- as.numeric(drws[d, paste0("beta_wp_pfu[", 1:P, "]")])
      bwr_pfu <- as.numeric(drws[d, paste0("beta_wr_pfu[", 1:P, "]")])
    }

    # Cholesky factor for correlated RNA individual effects
    if (has_corr) {
      si <- as.numeric(drws[d, paste0("sigma_ind_rna[", 1:4, "]")])
      L <- matrix(0, 4, 4)
      for (ii in 1:4) for (jj in 1:ii)
        L[ii, jj] <- drws[d, paste0("L_Omega_rna[", ii, ",", jj, "]")]
      LS <- diag(si) %*% L   # diag(sigma) * L_Omega
    }

    # Flat-top population value
    if (has_wf) {
      wf_pop <- stan_data$prior_wf_mean *
        exp(stan_data$prior_wf_cv * drws[d, "wf_raw[1]"])
    }

    for (a in seq_len(n_agents)) {

      # ── Sample NEW individual random effects ──────────────────────────────
      # RNA effects: from fitted population distribution
      if (has_corr) {
        z4 <- rnorm(4)
        eta_re <- LS %*% z4
        tp_re <- eta_re[1]; dp_re <- eta_re[2]
        wp_re <- eta_re[3]; wr_re <- eta_re[4]
      } else {
        tp_re <- rnorm(1, 0, prior_i_sd)
        dp_re <- if (has_ind) rnorm(1, 0, prior_i_sd) else 0
        wp_re <- if (has_ind) rnorm(1, 0, prior_i_sd) else 0
        wr_re <- if (has_ind) rnorm(1, 0, prior_i_sd) else 0
      }

      # PFU effects: independent, using posterior-learned population SDs
      tp_pfu_re <- rnorm(1, 0, sd_pfu$tp[d])
      dp_pfu_re <- if (has_ind) rnorm(1, 0, sd_pfu$dp[d]) else 0
      wp_pfu_re <- if (has_ind) rnorm(1, 0, sd_pfu$wp[d]) else 0
      wr_pfu_re <- if (has_ind) rnorm(1, 0, sd_pfu$wr[d]) else 0

      # Symptom RE
      z_sym <- rnorm(1)

      # Flat-top individual effect
      if (has_wf) {
        wf <- wf_pop * (if (has_ind) exp(rnorm(1, 0, prior_i_sd)) else 1)
      } else {
        wf <- 0
      }

      # ── RNA trajectory ────────────────────────────────────────────────────
      dp_rna <- dp_pop * exp(dp_re)
      wp_rna <- wp_pop * exp(wp_re)
      wr_rna <- wr_pop * exp(wr_re)
      tp_rna <- tp_re

      if (has_adj_rna) {
        xa <- x_mat[a, ]
        dp_rna <- dp_rna * exp(sum(xa * bdp_rna))
        wp_rna <- wp_rna * exp(sum(xa * bwp_rna))
        wr_rna <- wr_rna * exp(sum(xa * bwr_rna))
      }

      rna <- traj_fun(t_grid, tp_rna, wp_rna, wr_rna, dp_rna, wf, use_smooth)

      # ── PFU trajectory (log-affine from RNA) ──────────────────────────────
      dp_pfu <- exp(tdp[1]) * dp_rna^tdp[2]
      wp_pfu <- exp(twp[1]) * wp_rna^twp[2]
      wr_pfu <- exp(twr[1]) * wr_rna^twr[2]
      tp_pfu <- ttp[1] + ttp[2] * tp_rna + tp_pfu_re

      if (has_ind) {
        dp_pfu <- dp_pfu * exp(dp_pfu_re)
        wp_pfu <- wp_pfu * exp(wp_pfu_re)
        wr_pfu <- wr_pfu * exp(wr_pfu_re)
      }

      if (has_adj_pfu) {
        xa <- x_mat[a, ]
        dp_pfu <- dp_pfu * exp(sum(xa * bdp_pfu))
        wp_pfu <- wp_pfu * exp(sum(xa * bwp_pfu))
        wr_pfu <- wr_pfu * exp(sum(xa * bwr_pfu))
      }

      pfu <- traj_fun(t_grid, tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf, use_smooth)

      # Clamp to safe range (matches Stan's safe_vl)
      rna <- pmax(pmin(rna, 50), -50)
      # PFU clamped to RNA: infectious virus cannot exceed total viral RNA
      pfu <- pmax(pmin(pfu, rna), -50)

      # ── Detection window ──────────────────────────────────────────────────
      det <- rna > lod_rna | pfu > lod_pfu
      if (!any(det)) next
      i1 <- min(which(det));  i2 <- max(which(det))
      w  <- i1:i2;  nw <- length(w)

      # ── LFD probability ───────────────────────────────────────────────────
      lfd_p <- expit(lfd0 + lfd_c[1] * rna[w] + lfd_c[2] * pfu[w])

      # ── Symptom hazard (cloglog) ──────────────────────────────────────────
      u_sym  <- ss * z_sym
      eta_l  <- ei + ep * (pfu[w] / scale_vl) + er * (rna[w] / scale_vl) + u_sym
      sym_hz <- 1 - exp(-exp(pmin(eta_l, 10)))

      # ── Generate observations ─────────────────────────────────────────────
      if (include_noise) {
        rna_m <- pmax(rnorm(nw, rna[w], s_rna), lod_rna)
        pfu_m <- pmax(rnorm(nw, pfu[w], s_pfu), lod_pfu)
        lfd_o <- rbinom(nw, 1, lfd_p)

        # Symptom onset survival process: once symptomatic, stays symptomatic
        sym_o <- integer(nw)
        for (tt in seq_len(nw)) {
          if (tt > 1L && sym_o[tt - 1L] == 1L) {
            sym_o[tt] <- 1L
          } else {
            sym_o[tt] <- rbinom(1L, 1L, sym_hz[tt])
          }
        }
      } else {
        rna_m <- rna[w]
        pfu_m <- pfu[w]
        lfd_o <- lfd_p
        sym_o <- sym_hz
      }

      # ── Store result ──────────────────────────────────────────────────────
      ri <- ri + 1L
      tw <- t_grid[w]
      results[[ri]] <- tibble::tibble(
        agent_id             = a,
        draw                 = d,
        day                  = tw,
        day_since_detectable = tw - tw[1],
        log_rna              = rna[w],
        log_pfu              = pfu[w],        rna_detectable       = rna[w] > lod_rna,
        pfu_detectable       = pfu[w] > lod_pfu,        rna                  = rna_m,
        pfu                  = pfu_m,
        lfd                  = lfd_o,
        symptomatic          = sym_o
      )
    }
  }

  out <- dplyr::bind_rows(results[seq_len(ri)])

  # Append agent characteristics for easy grouping/faceting
  agent_chars <- dplyr::mutate(
    tibble::as_tibble(agents),
    agent_id = dplyr::row_number()
  )
  out <- dplyr::left_join(out, agent_chars, by = "agent_id")

  out
}
