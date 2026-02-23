# ──────────────────────────────────────────────────────────────────────────────
# diagnostics.R — Model convergence checks, LOO-CV, posterior predictive checks
# ──────────────────────────────────────────────────────────────────────────────

#' Check MCMC convergence diagnostics
#'
#' @param fit CmdStanMCMC fit object
#' @return List with divergences, rhat, ess_bulk, ess_tail summaries
check_convergence <- function(fit) {

  diag <- fit$diagnostic_summary()
  summ <- fit$summary()

  list(
    n_divergent   = sum(diag$num_divergent),
    n_max_treedepth = sum(diag$num_max_treedepth),
    rhat_max      = max(summ$rhat, na.rm = TRUE),
    rhat_warnings = summ |>
      dplyr::filter(rhat > 1.01) |>
      dplyr::select(variable, rhat),
    ess_bulk_min  = min(summ$ess_bulk, na.rm = TRUE),
    ess_tail_min  = min(summ$ess_tail, na.rm = TRUE),
    ess_warnings  = summ |>
      dplyr::filter(ess_bulk < 400 | ess_tail < 400) |>
      dplyr::select(variable, ess_bulk, ess_tail)
  )
}


#' Compute LOO-CV using PSIS-LOO (per-outcome and combined)
#'
#' Extracts per-outcome log-lik arrays (log_lik_rna, log_lik_pfu,
#' log_lik_lfd, log_lik_sym) from transformed parameters and the
#' combined log_lik from generated quantities. Returns a named list
#' of loo objects, one per outcome type plus "combined".
#'
#' Per-outcome LOO is more interpretable than the combined version
#' because the combined log_lik accumulates contributions from
#' multiple outcomes into a single value, inflating p_loo.
#'
#' @param fit CmdStanMCMC fit object
#' @return Named list of loo objects: rna, pfu, lfd, sym, combined
compute_loo <- function(fit) {
  outcomes <- c(rna = "log_lik_rna", pfu = "log_lik_pfu",
                lfd = "log_lik_lfd", sym = "log_lik_sym",
                combined = "log_lik")

  results <- list()

  for (nm in names(outcomes)) {
    varname <- outcomes[[nm]]
    ll <- tryCatch(
      fit$draws(varname, format = "draws_array"),
      error = function(e) NULL
    )
    if (is.null(ll) || prod(dim(ll)) == 0) {
      results[[nm]] <- NULL
      next
    }
    # Skip outcomes with very few observations (< 10)
    if (dim(ll)[3] < 10) {
      results[[nm]] <- NULL
      next
    }
    r_eff <- loo::relative_eff(exp(ll))
    results[[nm]] <- loo::loo(ll, r_eff = r_eff)
  }

  results
}


#' Summarize per-outcome LOO-CV results as a tibble
#'
#' Extracts key metrics from each loo object in the list returned by
#' compute_loo(): ELPD, p_loo, number of observations, and Pareto k
#' diagnostics (number and fraction with k > 0.7).
#'
#' @param loo_list Named list of loo objects (from compute_loo())
#' @return Tibble with one row per outcome
summarize_loo <- function(loo_list) {
  rows <- list()
  for (nm in names(loo_list)) {
    obj <- loo_list[[nm]]
    if (is.null(obj)) next
    est <- obj$estimates
    pk  <- obj$diagnostics$pareto_k
    rows[[nm]] <- tibble::tibble(
      outcome   = nm,
      n_obs     = length(pk),
      elpd      = est["elpd_loo", "Estimate"],
      elpd_se   = est["elpd_loo", "SE"],
      p_loo     = est["p_loo", "Estimate"],
      p_loo_se  = est["p_loo", "SE"],
      looic     = est["looic", "Estimate"],
      looic_se  = est["looic", "SE"],
      n_k_high  = sum(pk > 0.7, na.rm = TRUE),
      pct_k_high = 100 * mean(pk > 0.7, na.rm = TRUE)
    )
  }
  dplyr::bind_rows(rows)
}


#' Compute per-outcome WAIC as a complement to PSIS-LOO
#'
#' WAIC (Widely Applicable Information Criterion) is more stable than
#' PSIS-LOO when many Pareto k values exceed 0.7, which is expected
#' in models with dense individual random effects.
#'
#' @param fit CmdStanMCMC fit object
#' @return Named list of loo::waic objects: rna, pfu, lfd, sym, combined
compute_waic <- function(fit) {
  outcomes <- c(rna = "log_lik_rna", pfu = "log_lik_pfu",
                lfd = "log_lik_lfd", sym = "log_lik_sym",
                combined = "log_lik")

  results <- list()

  for (nm in names(outcomes)) {
    varname <- outcomes[[nm]]
    ll <- tryCatch(
      fit$draws(varname, format = "draws_matrix"),
      error = function(e) NULL
    )
    if (is.null(ll) || prod(dim(ll)) == 0) {
      results[[nm]] <- NULL
      next
    }
    if (ncol(ll) < 10) {
      results[[nm]] <- NULL
      next
    }
    results[[nm]] <- loo::waic(ll)
  }

  results
}


#' Summarize WAIC results as a tibble
#'
#' @param waic_list Named list of waic objects (from compute_waic())
#' @return Tibble with one row per outcome
summarize_waic <- function(waic_list) {
  rows <- list()
  for (nm in names(waic_list)) {
    obj <- waic_list[[nm]]
    if (is.null(obj)) next
    est <- obj$estimates
    rows[[nm]] <- tibble::tibble(
      outcome   = nm,
      n_obs     = nrow(obj$pointwise),
      elpd_waic = est["elpd_waic", "Estimate"],
      elpd_se   = est["elpd_waic", "SE"],
      p_waic    = est["p_waic", "Estimate"],
      p_waic_se = est["p_waic", "SE"],
      waic      = est["waic", "Estimate"],
      waic_se   = est["waic", "SE"]
    )
  }
  dplyr::bind_rows(rows)
}


#' Grouped k-fold CV (by individual)
#'
#' Partitions individuals into K folds, refits the model K times with
#' each fold held out, and computes out-of-fold log-likelihood.
#' This is computationally expensive (K × model fit time) but avoids
#' the Pareto k issues of PSIS-LOO with dense individual random effects.
#'
#' @param stan_file  Path to Stan model file
#' @param stan_data  Stan data list
#' @param stacked_dat Stacked analysis dataset (for individual IDs)
#' @param K_folds   Number of folds (default: 10)
#' @param seed      Random seed for fold assignment
#' @param ...       Additional arguments passed to fit_model()
#' @return List with: fold_logliks (per-obs held-out log-lik),
#'   elpd (total expected log predictive density), fold_assignments
kfold_cv <- function(stan_file, stan_data, stacked_dat, K_folds = 10,
                     seed = 42, ...) {
  set.seed(seed)

  # Assign individuals (not observations) to folds
  ind_ids <- unique(stacked_dat$id)
  fold_assign <- sample(rep(1:K_folds, length.out = length(ind_ids)))
  names(fold_assign) <- ind_ids

  # Map observations to folds via their individual
  obs_folds <- fold_assign[as.character(stacked_dat$id)]

  # Accumulate held-out log-likelihoods
  held_out_ll <- rep(NA_real_, nrow(stacked_dat))

  for (k in seq_len(K_folds)) {
    message(sprintf("k-fold CV: fold %d / %d", k, K_folds))

    # Build training data: exclude fold k individuals
    train_ids <- ind_ids[fold_assign != k]
    # Subsample stan_data to only include training individuals
    # This requires rebuilding the stan data — use subsample_stan_data
    # with a custom set of IDs
    train_data <- subset_stan_data_by_ids(stan_data, stacked_dat, train_ids)

    # Fit model on training data
    fit_k <- fit_model(stan_file, train_data, ...)

    # Extract log_lik for held-out observations
    # This requires computing log-lik on the held-out data using
    # the posterior draws from the training fit.
    # NOTE: Stan's generated quantities can be used for this via
    #   fit_k$generate_quantities(data = test_data)
    # but requires the test data to have the same structure.
    test_ids  <- ind_ids[fold_assign == k]
    test_data <- subset_stan_data_by_ids(stan_data, stacked_dat, test_ids)

    gq <- fit_k$generate_quantities(data = test_data,
                                     fitted_params = fit_k$draws())
    ll_k <- gq$draws("log_lik", format = "draws_matrix")
    held_out_ll[obs_folds == k] <- colMeans(ll_k)
  }

  list(
    elpd_kfold      = sum(held_out_ll, na.rm = TRUE),
    pointwise       = held_out_ll,
    fold_assignments = obs_folds,
    K_folds         = K_folds
  )
}


#' Subset Stan data to a specific set of individual IDs
#'
#' Rebuilds the Stan data list keeping only observations belonging to
#' the specified individuals. Updates M, N, id mappings, and all
#' observation-level arrays. Used internally by kfold_cv().
#'
#' @param stan_data    Full Stan data list
#' @param stacked_dat  Stacked analysis dataset
#' @param keep_ids     Integer vector of individual IDs to retain
#' @return Stan data list with only the specified individuals
subset_stan_data_by_ids <- function(stan_data, stacked_dat, keep_ids) {
  # Identify which observations belong to kept individuals
  keep_obs <- stacked_dat$id %in% keep_ids

  # Rebuild observation-level arrays
  obs_vars <- c("rna", "pfu", "lfd", "sym",
                "rna_exist", "pfu_exist", "lfd_exist", "sym_exist",
                "sym_at_risk", "pfu_type", "time", "source", "id")

  d <- stan_data
  for (v in obs_vars) {
    if (v %in% names(d) && length(d[[v]]) == length(keep_obs)) {
      d[[v]] <- d[[v]][keep_obs]
    }
  }

  # Re-index individual IDs to be sequential 1:N_new
  old_ids <- sort(unique(d$id))
  id_map  <- setNames(seq_along(old_ids), old_ids)
  d$id    <- as.integer(id_map[as.character(d$id)])

  # Update M (per-source counts) and N (per-source obs)
  src_tab <- table(factor(d$source, levels = 1:d$K))
  d$N     <- as.integer(src_tab)

  id_src  <- stacked_dat$source[match(old_ids, stacked_dat$id)]
  m_tab   <- table(factor(id_src, levels = 1:d$K))
  d$M     <- as.integer(m_tab)

  # Subset covariate matrix
  if ("x" %in% names(d) && is.matrix(d$x)) {
    d$x <- d$x[old_ids, , drop = FALSE]
  }

  d
}


#' Generate posterior predictive check data
#'
#' Extracts replicated outcomes (rna_rep, pfu_rep, lfd_rep) from the
#' generated quantities block and compares to observed data.
#'
#' @param fit         CmdStanMCMC fit object
#' @param stacked_dat Stacked analysis dataset
#' @return List with observed and replicated data for each outcome
posterior_predictive_check <- function(fit, stacked_dat) {

  rna_rep <- fit$draws("rna_rep", format = "draws_matrix")
  pfu_rep <- fit$draws("pfu_rep", format = "draws_matrix")
  lfd_rep <- fit$draws("lfd_rep", format = "draws_matrix")

  list(
    rna = list(
      y     = stacked_dat$rna[stacked_dat$rna_exist == 1],
      y_rep = rna_rep[, stacked_dat$rna_exist == 1],
      idx   = which(stacked_dat$rna_exist == 1)
    ),
    pfu = list(
      y     = stacked_dat$pfu[stacked_dat$pfu_exist == 1],
      y_rep = pfu_rep[, stacked_dat$pfu_exist == 1],
      idx   = which(stacked_dat$pfu_exist == 1)
    ),
    lfd = list(
      y     = stacked_dat$lfd[stacked_dat$lfd_exist == 1],
      y_rep = lfd_rep[, stacked_dat$lfd_exist == 1],
      idx   = which(stacked_dat$lfd_exist == 1)
    )
  )
}


#' Plot MCMC trace plots for key parameters
#'
#' @param fit       CmdStanMCMC fit object
#' @param params    Character vector of parameter names to plot
#' @param stan_data Stan data list (optional; used to include wf_raw when use_wf=1)
#' @param out_file  Path to save the PDF (or NULL for interactive)
#' @return File path (character) if out_file is specified; ggplot otherwise
plot_traces <- function(fit, params = NULL, stan_data = NULL, out_file = NULL) {
  if (is.null(params)) {
    params <- c("dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
                "sigma_rna", "sigma_pfu",
                "tau_dp", "tau_wp",
                "eta_sym_intercept")
    if (!is.null(stan_data) && isTRUE(stan_data$use_wf == 1)) {
      params <- c(params, "wf_raw")
    }
    if (!is.null(stan_data) && isTRUE(stan_data$ind_corr == 1)) {
      params <- c(params, "sigma_ind_rna")
    }
  }
  draws <- fit$draws(variables = params)
  p <- bayesplot::mcmc_trace(draws)

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(out_file, p, width = 12, height = 8)
    return(out_file)
  }

  invisible(p)
}


#' Check parameter recovery from a fit on simulated data
#'
#' Compares posterior summaries to the known true values returned by
#' \code{\link{simulate_data}}.  For each scalar population parameter the
#' function reports the true value, posterior median, 90 % CI, and whether
#' the truth falls inside that interval ("covered").
#'
#' @param fit       CmdStanMCMC fit object (fitted on \code{sim$sim_data})
#' @param truth     The \code{$truth} element returned by
#'                  \code{\link{simulate_data}}
#' @param stan_data Stan data list (optional; used to detect use_wf flag)
#' @param prob      Credible interval width (default 0.90)
#' @return A tibble with one row per parameter
check_recovery <- function(fit, truth, stan_data = NULL, prob = 0.90) {

  alpha <- (1 - prob) / 2

  # Scalar population parameters to check (name in truth → Stan variable)
  scalar_map <- c(
    dp_mean_rna       = "dp_mean_rna",
    wp_mean_rna       = "wp_mean_rna",
    wr_mean_rna       = "wr_mean_rna",
    sigma_rna         = "sigma_rna",
    sigma_pfu         = "sigma_pfu",
    eta_sym_intercept = "eta_sym_intercept",
    eta_sym_pfu       = "eta_sym_pfu",
    eta_sym_rna       = "eta_sym_rna",
    sigma_sym         = "sigma_sym",
    fp                = "fp",
    fn                = "fn"
  )

  # tau_* are vector[2] in Stan — map component-wise
  vec2_map <- list(
    list(truth_name = "tau0_tp", stan_var = "tau_tp", index = 1),
    list(truth_name = "tau_tp",  stan_var = "tau_tp", index = 2),
    list(truth_name = "tau0_dp", stan_var = "tau_dp", index = 1),
    list(truth_name = "tau_dp",  stan_var = "tau_dp", index = 2),
    list(truth_name = "tau0_wp", stan_var = "tau_wp", index = 1),
    list(truth_name = "tau_wp",  stan_var = "tau_wp", index = 2),
    list(truth_name = "tau0_wr", stan_var = "tau_wr", index = 1),
    list(truth_name = "tau_wr",  stan_var = "tau_wr", index = 2)
  )

  # alpha_tcid50[1:3] and alpha_cult[1:2] — viral culture model
  # (truth stores as vector; Stan also vector)
  if (!is.null(truth$alpha_tcid50)) {
    vec2_map <- c(vec2_map, list(
      list(truth_name = "alpha_tcid50_a", stan_var = "alpha_tcid50", index = 1,
           truth_val = truth$alpha_tcid50[1]),
      list(truth_name = "alpha_tcid50_log_b", stan_var = "alpha_tcid50", index = 2,
           truth_val = truth$alpha_tcid50[2]),
      list(truth_name = "alpha_tcid50_log_sig", stan_var = "alpha_tcid50", index = 3,
           truth_val = truth$alpha_tcid50[3])
    ))
  }
  if (!is.null(truth$alpha_cult)) {
    vec2_map <- c(vec2_map, list(
      list(truth_name = "alpha_cult_1", stan_var = "alpha_cult", index = 1,
           truth_val = truth$alpha_cult[1]),
      list(truth_name = "alpha_cult_2", stan_var = "alpha_cult", index = 2,
           truth_val = truth$alpha_cult[2])
    ))
  }

  summ <- fit$summary()

  # --- scalars ---
  scalar_rows <- purrr::map_dfr(names(scalar_map), function(nm) {
    stan_nm <- scalar_map[[nm]]
    row <- dplyr::filter(summ, variable == stan_nm)
    if (nrow(row) == 0) return(NULL)
    draws <- as.numeric(fit$draws(stan_nm, format = "matrix"))
    qi <- stats::quantile(draws, c(alpha, 1 - alpha))
    true_val <- truth[[nm]]
    tibble::tibble(
      parameter = nm,
      true_value = true_val,
      posterior_median = stats::median(draws),
      ci_lo = qi[[1]],
      ci_hi = qi[[2]],
      covered = true_val >= qi[[1]] && true_val <= qi[[2]],
      bias = stats::median(draws) - true_val,
      width = qi[[2]] - qi[[1]]
    )
  })

  # --- vector[2] components ---
  vec2_rows <- purrr::map_dfr(vec2_map, function(v) {
    stan_nm <- paste0(v$stan_var, "[", v$index, "]")
    draws <- tryCatch(
      as.numeric(fit$draws(stan_nm, format = "matrix")),
      error = function(e) NULL
    )
    if (is.null(draws)) return(NULL)
    qi <- stats::quantile(draws, c(alpha, 1 - alpha))
    # Use explicit truth_val if provided (for vector-valued truth entries),
    # otherwise look up by truth_name in the truth list
    true_val <- if (!is.null(v$truth_val)) v$truth_val else truth[[v$truth_name]]
    tibble::tibble(
      parameter = v$truth_name,
      true_value = true_val,
      posterior_median = stats::median(draws),
      ci_lo = qi[[1]],
      ci_hi = qi[[2]],
      covered = true_val >= qi[[1]] && true_val <= qi[[2]],
      bias = stats::median(draws) - true_val,
      width = qi[[2]] - qi[[1]]
    )
  })

  # --- derived wf_mean (only when use_wf is active) ---
  wf_rows <- NULL
  if (!is.null(stan_data) && isTRUE(stan_data$use_wf == 1) &&
      !is.null(truth$wf_mean)) {
    wf_raw_draws <- tryCatch(
      as.numeric(fit$draws("wf_raw[1]", format = "matrix")),
      error = function(e) NULL
    )
    if (!is.null(wf_raw_draws)) {
      wf_draws <- stan_data$prior_wf_mean *
        exp(stan_data$prior_wf_cv * wf_raw_draws)
      qi <- stats::quantile(wf_draws, c(alpha, 1 - alpha))
      true_val <- truth$wf_mean
      wf_rows <- tibble::tibble(
        parameter = "wf_mean",
        true_value = true_val,
        posterior_median = stats::median(wf_draws),
        ci_lo = qi[[1]],
        ci_hi = qi[[2]],
        covered = true_val >= qi[[1]] && true_val <= qi[[2]],
        bias = stats::median(wf_draws) - true_val,
        width = qi[[2]] - qi[[1]]
      )
    }
  }

  # --- correlated RE parameters (sigma_ind_rna and Omega_rna) ---
  corr_rows <- NULL
  if (!is.null(stan_data) && isTRUE(stan_data$ind_corr == 1)) {
    re_labels <- c("tp", "dp", "wp", "wr")

    # sigma_ind_rna[1..4]
    if (!is.null(truth$sigma_ind_rna)) {
      corr_rows <- purrr::map_dfr(1:4, function(i) {
        stan_nm <- paste0("sigma_ind_rna[", i, "]")
        draws <- tryCatch(
          as.numeric(fit$draws(stan_nm, format = "matrix")),
          error = function(e) NULL
        )
        if (is.null(draws)) return(NULL)
        qi <- stats::quantile(draws, c(alpha, 1 - alpha))
        true_val <- truth$sigma_ind_rna[i]
        tibble::tibble(
          parameter = paste0("sigma_ind_rna_", re_labels[i]),
          true_value = true_val,
          posterior_median = stats::median(draws),
          ci_lo = qi[[1]], ci_hi = qi[[2]],
          covered = true_val >= qi[[1]] && true_val <= qi[[2]],
          bias = stats::median(draws) - true_val,
          width = qi[[2]] - qi[[1]]
        )
      })
    }

    # Omega_rna off-diagonal correlations (6 unique pairs)
    if (!is.null(truth$Omega_rna)) {
      omega_rows <- purrr::map_dfr(
        list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)),
        function(ij) {
          i <- ij[1]; j <- ij[2]
          stan_nm <- paste0("Omega_rna[", i, ",", j, "]")
          draws <- tryCatch(
            as.numeric(fit$draws(stan_nm, format = "matrix")),
            error = function(e) NULL
          )
          if (is.null(draws)) return(NULL)
          qi <- stats::quantile(draws, c(alpha, 1 - alpha))
          true_val <- truth$Omega_rna[i, j]
          tibble::tibble(
            parameter = paste0("rho_rna_", re_labels[i], "_", re_labels[j]),
            true_value = true_val,
            posterior_median = stats::median(draws),
            ci_lo = qi[[1]], ci_hi = qi[[2]],
            covered = true_val >= qi[[1]] && true_val <= qi[[2]],
            bias = stats::median(draws) - true_val,
            width = qi[[2]] - qi[[1]]
          )
        }
      )
      corr_rows <- dplyr::bind_rows(corr_rows, omega_rows)
    }
  }

  dplyr::bind_rows(scalar_rows, vec2_rows, wf_rows, corr_rows)
}


#' Plot parameter recovery results
#'
#' @param recovery  Tibble from \code{\link{check_recovery}}
#' @param out_file  Path to save the PDF (or NULL for interactive)
#' @return ggplot object (invisible)
plot_recovery <- function(recovery, out_file = NULL) {

  recovery$parameter <- factor(recovery$parameter,
                                levels = rev(recovery$parameter))

  p <- ggplot2::ggplot(recovery, ggplot2::aes(y = parameter)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        x     = posterior_median,
        xmin  = ci_lo,
        xmax  = ci_hi,
        color = covered
      ),
      size = 0.6
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = true_value),
      shape = 4, size = 3, stroke = 1.2, color = "black"
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "steelblue", "FALSE" = "firebrick"),
      labels = c("TRUE" = "Covered", "FALSE" = "Missed"),
      name   = "90% CI"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      x = "Parameter value",
      y = NULL,
      title = "Parameter recovery (simulated data)",
      subtitle = "X = true value; point-range = posterior median + 90% CI"
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 8, height = 6)
    return(out_file)
  }

  invisible(p)
}


# ──────────────────────────────────────────────────────────────────────────────
# Publication-quality diagnostic figures and tables
# ──────────────────────────────────────────────────────────────────────────────

#' Forest plot of covariate effects on RNA kinetics
#'
#' Displays exponentiated beta coefficients (multiplicative effects) with
#' 95% credible intervals, grouped by kinetic parameter (peak, proliferation,
#' clearance). A vertical reference line at 1.0 marks no effect.
#'
#' @param param_summary Output of \code{\link{summarize_parameters}}
#' @param out_file      Path to save PDF/PNG (NULL for interactive)
#' @return ggplot object or file path
plot_forest <- function(param_summary, out_file = NULL) {

  df <- param_summary$covariate_effects
  df$coef  <- as.numeric(df$coef)
  df$ci_lo <- as.numeric(df$ci_lo)
  df$ci_hi <- as.numeric(df$ci_hi)

  # Readable facet labels
  param_labels <- c(
    dp = "Peak magnitude",
    wp = "Proliferation rate",
    wr = "Clearance rate"
  )
  df$param_label <- param_labels[df$parameter]
  df$param_label <- factor(df$param_label,
                           levels = c("Peak magnitude",
                                      "Proliferation rate",
                                      "Clearance rate"))

  # Significance flag: CI excludes 1
  df$sig <- df$ci_lo > 1 | df$ci_hi < 1

  # Reverse label order for top-down reading
  df$label <- factor(df$label, levels = rev(unique(df$label)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = coef, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lo, xmax = ci_hi, color = sig),
      size = 0.4
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
      guide  = "none"
    ) +
    ggplot2::facet_wrap(~ param_label, ncol = 3) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y      = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      x = "Multiplicative effect (95% CI)",
      y = NULL,
      title = "Covariate effects on RNA kinetics"
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 12, height = 5)
    return(out_file)
  }

  invisible(p)
}


#' PPC density overlay for RNA observations
#'
#' Plots the observed RNA distribution against posterior predictive
#' replications (semi-transparent density curves from draws).
#'
#' @param ppc      Output of \code{\link{posterior_predictive_check}}
#' @param n_draws  Number of posterior draws to overlay (default 100)
#' @param out_file Path to save PDF/PNG (NULL for interactive)
#' @return ggplot object or file path
plot_ppc_rna <- function(ppc, n_draws = 100, out_file = NULL) {

  y     <- ppc$rna$y
  y_rep <- as.matrix(ppc$rna$y_rep)

  # subsample draws for visual clarity
  draw_idx <- sample(nrow(y_rep), min(n_draws, nrow(y_rep)))

  # build long-format data
  obs_df <- data.frame(val = y, stringsAsFactors = FALSE)
  rep_list <- lapply(draw_idx, function(i) {
    data.frame(val  = as.numeric(y_rep[i, ]),
               draw = as.character(i),
               stringsAsFactors = FALSE)
  })
  rep_df <- do.call(rbind, rep_list)

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = rep_df,
      ggplot2::aes(x = val, group = draw),
      color = "steelblue", alpha = 0.08, linewidth = 0.3
    ) +
    ggplot2::geom_density(
      data = obs_df,
      ggplot2::aes(x = val),
      color = "black", linewidth = 0.8
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = expression(log[10]~RNA~copies/mL),
      y = "Density",
      title = "Posterior predictive check: RNA",
      subtitle = paste0("Black = observed; blue = ", length(draw_idx),
                        " replicated datasets")
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 7, height = 5)
    return(out_file)
  }

  invisible(p)
}


#' PPC density overlay for PFU observations
#'
#' @param ppc      Output of \code{\link{posterior_predictive_check}}
#' @param n_draws  Number of posterior draws to overlay (default 100)
#' @param out_file Path to save PDF/PNG (NULL for interactive)
#' @return ggplot object or file path
plot_ppc_pfu <- function(ppc, n_draws = 100, out_file = NULL) {

  y     <- ppc$pfu$y
  y_rep <- as.matrix(ppc$pfu$y_rep)

  draw_idx <- sample(nrow(y_rep), min(n_draws, nrow(y_rep)))

  obs_df <- data.frame(val = y, stringsAsFactors = FALSE)
  rep_list <- lapply(draw_idx, function(i) {
    data.frame(val  = as.numeric(y_rep[i, ]),
               draw = as.character(i),
               stringsAsFactors = FALSE)
  })
  rep_df <- do.call(rbind, rep_list)

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = rep_df,
      ggplot2::aes(x = val, group = draw),
      color = "firebrick", alpha = 0.08, linewidth = 0.3
    ) +
    ggplot2::geom_density(
      data = obs_df,
      ggplot2::aes(x = val),
      color = "black", linewidth = 0.8
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = expression(log[10]~PFU/mL),
      y = "Density",
      title = "Posterior predictive check: PFU",
      subtitle = paste0("Black = observed; red = ", length(draw_idx),
                        " replicated datasets")
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 7, height = 5)
    return(out_file)
  }

  invisible(p)
}


#' PPC calibration plot for LFD (binary outcome)
#'
#' Bins predicted probabilities and plots observed positive frequency
#' against mean predicted probability in each bin. Perfect calibration
#' falls on the diagonal.
#'
#' @param ppc       Output of \code{\link{posterior_predictive_check}}
#' @param n_bins    Number of probability bins (default 10)
#' @param out_file  Path to save PDF/PNG (NULL for interactive)
#' @return ggplot object or file path
plot_ppc_lfd <- function(ppc, n_bins = 10, out_file = NULL) {

  y     <- ppc$lfd$y
  y_rep <- as.matrix(ppc$lfd$y_rep)

  # posterior mean predicted probability per observation
  p_hat <- colMeans(y_rep)

  calib_df <- data.frame(
    y     = y,
    p_hat = p_hat,
    bin   = cut(p_hat, breaks = seq(0, 1, length.out = n_bins + 1),
                include.lowest = TRUE)
  )

  bin_summary <- calib_df |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      mean_pred = mean(p_hat),
      obs_freq  = mean(y),
      n         = dplyr::n(),
      se        = sqrt(obs_freq * (1 - obs_freq) / n),
      .groups   = "drop"
    )

  p <- ggplot2::ggplot(bin_summary,
                       ggplot2::aes(x = mean_pred, y = obs_freq)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                          linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = pmax(obs_freq - 1.96 * se, 0),
                    ymax = pmin(obs_freq + 1.96 * se, 1)),
      size = 0.5, color = "darkgreen"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = n),
      vjust = -1.2, size = 3, color = "grey40"
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x     = "Mean predicted P(LFD+)",
      y     = "Observed LFD+ frequency",
      title = "Calibration — LFD rapid test",
      subtitle = "Numbers = observations per bin"
    )

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(out_file, p, width = 6, height = 6)
    return(out_file)
  }

  invisible(p)
}


#' Convergence summary table (LaTeX)
#'
#' Produces a clean LaTeX table with key convergence metrics:
#' divergences, max treedepth, R-hat, ESS.
#'
#' @param convergence Output of \code{\link{check_convergence}}
#' @param out_file    Path to save .tex file
#' @return File path (invisibly)
save_convergence_table <- function(convergence, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

  tbl <- tibble::tibble(
    Metric = c(
      "Divergent transitions",
      "Max treedepth exceedances",
      "Max $\\hat{R}$",
      "Parameters with $\\hat{R} > 1.01$",
      "Min bulk ESS",
      "Min tail ESS",
      "Parameters with ESS $< 400$"
    ),
    Value = c(
      as.character(convergence$n_divergent),
      as.character(convergence$n_max_treedepth),
      sprintf("%.3f", convergence$rhat_max),
      as.character(nrow(convergence$rhat_warnings)),
      sprintf("%.0f", convergence$ess_bulk_min),
      sprintf("%.0f", convergence$ess_tail_min),
      as.character(nrow(convergence$ess_warnings))
    )
  )

  tex <- knitr::kable(
    tbl,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    col.names = c("Diagnostic", "Value"),
    caption   = "MCMC convergence diagnostics"
  )
  writeLines(as.character(tex), out_file)
  invisible(out_file)
}


#' Generate all publication-quality diagnostic figures and tables
#'
#' Orchestrator function that produces: forest plot, PPC density overlays,
#' LFD calibration, convergence table, and parameter recovery.
#'
#' @param param_summary  Output of \code{summarize_parameters()}
#' @param ppc            Output of \code{posterior_predictive_check()}
#' @param convergence    Output of \code{check_convergence()}
#' @param recovery_check Output of \code{check_recovery()} (optional)
#' @param out_dir        Output directory
#' @return Character vector of all generated file paths
plot_diagnostics <- function(param_summary, ppc, convergence,
                             recovery_check = NULL,
                             out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character()

  # 1. Forest plot
  f <- plot_forest(param_summary,
                   out_file = file.path(out_dir, "forest_covariates.pdf"))
  files <- c(files, f)

  # 2. PPC — RNA
  f <- plot_ppc_rna(ppc,
                    out_file = file.path(out_dir, "ppc_rna.pdf"))
  files <- c(files, f)

  # 3. PPC — PFU
  f <- plot_ppc_pfu(ppc,
                    out_file = file.path(out_dir, "ppc_pfu.pdf"))
  files <- c(files, f)

  # 4. PPC — LFD calibration
  f <- plot_ppc_lfd(ppc,
                    out_file = file.path(out_dir, "ppc_lfd_calibration.pdf"))
  files <- c(files, f)

  # 5. Convergence table
  f <- save_convergence_table(convergence,
                              out_file = file.path(dirname(out_dir),
                                                    "tables",
                                                    "convergence.tex"))
  files <- c(files, f)

  # 6. Parameter recovery (if available)
  if (!is.null(recovery_check)) {
    f <- plot_recovery(recovery_check,
                       out_file = file.path(out_dir, "param_recovery.pdf"))
    files <- c(files, f)
  }

  files
}


# ── Correlation matrix diagnostics ────────────────────────────────────────────

#' Summarize posterior distribution of RNA individual-effect correlations
#'
#' Extracts \code{Omega_rna} and \code{sigma_ind_rna} from the posterior
#' and returns a tidy tibble with medians and credible intervals.
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list (must have \code{ind_corr = 1})
#' @param prob      Width of the credible interval (default 0.90)
#' @return A tibble with columns: parameter, label, estimate, ci_lo, ci_hi
summarize_correlation <- function(fit, stan_data, prob = 0.90) {
  if (!isTRUE(stan_data$ind_corr == 1)) {
    return(tibble::tibble(
      parameter = character(), label = character(),
      estimate = numeric(), ci_lo = numeric(), ci_hi = numeric()
    ))
  }

  alpha <- (1 - prob) / 2
  re_labels <- c("tp (peak time)", "dp (peak height)",
                 "wp (proliferation)", "wr (clearance)")

  k <- posterior::as_draws_rvars(
    fit$draws(variables = c("sigma_ind_rna", "Omega_rna"))
  )

  # Helper to summarize an rvar
  .summ <- function(rv) {
    d <- posterior::draws_of(rv)
    tibble::tibble(
      estimate = round(stats::median(d), 3),
      ci_lo    = round(stats::quantile(d, alpha), 3),
      ci_hi    = round(stats::quantile(d, 1 - alpha), 3)
    )
  }

  # --- sigma_ind_rna (standard deviations) ---
  sigma_rows <- purrr::map_dfr(1:4, function(i) {
    dplyr::bind_cols(
      tibble::tibble(
        parameter = paste0("sigma_ind_rna[", i, "]"),
        label = paste0("SD of ", re_labels[i])
      ),
      .summ(k$sigma_ind_rna[i])
    )
  })

  # --- Omega_rna off-diagonal correlations ---
  pairs <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
  omega_rows <- purrr::map_dfr(pairs, function(ij) {
    i <- ij[1]; j <- ij[2]
    dplyr::bind_cols(
      tibble::tibble(
        parameter = paste0("Omega_rna[", i, ",", j, "]"),
        label = paste0("Corr(", re_labels[i], ", ", re_labels[j], ")")
      ),
      .summ(k$Omega_rna[i, j])
    )
  })

  dplyr::bind_rows(sigma_rows, omega_rows)
}


#' Plot posterior correlation matrix for RNA individual effects
#'
#' Produces a heatmap of the posterior median 4×4 correlation matrix
#' (\code{Omega_rna}) with point estimates annotated and an interval
#' summary in a companion panel.
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list (must have \code{ind_corr = 1})
#' @param prob      Credible interval width for annotations (default 0.90)
#' @param out_file  Path to save the PDF (or NULL for interactive)
#' @return ggplot object (invisible)
plot_correlation_matrix <- function(fit, stan_data, prob = 0.90,
                                    out_file = NULL) {
  if (!isTRUE(stan_data$ind_corr == 1)) {
    message("ind_corr = 0: no correlation matrix to plot")
    return(invisible(NULL))
  }

  alpha <- (1 - prob) / 2
  re_labels <- c("tp", "dp", "wp", "wr")

  k <- posterior::as_draws_rvars(
    fit$draws(variables = "Omega_rna")
  )

  # Build median correlation matrix
  med_mat <- matrix(NA, 4, 4)
  lo_mat  <- matrix(NA, 4, 4)
  hi_mat  <- matrix(NA, 4, 4)
  for (i in 1:4) {
    for (j in 1:4) {
      d <- posterior::draws_of(k$Omega_rna[i, j])
      med_mat[i, j] <- stats::median(d)
      lo_mat[i, j]  <- stats::quantile(d, alpha)
      hi_mat[i, j]  <- stats::quantile(d, 1 - alpha)
    }
  }

  # Convert to long format for ggplot
  df <- expand.grid(row = 1:4, col = 1:4)
  df$median <- as.vector(med_mat)
  df$lo     <- as.vector(lo_mat)
  df$hi     <- as.vector(hi_mat)
  df$row_lab <- factor(re_labels[df$row], levels = re_labels)
  df$col_lab <- factor(re_labels[df$col], levels = re_labels)
  df$label <- ifelse(
    df$row == df$col,
    "1",
    sprintf("%.2f\n[%.2f, %.2f]", df$median, df$lo, df$hi)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_lab, y = row_lab)) +
    ggplot2::geom_tile(ggplot2::aes(fill = median), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1),
      name = "Posterior\nmedian"
    ) +
    ggplot2::scale_y_discrete(limits = rev(re_labels)) +
    ggplot2::labs(
      title = expression(paste("Posterior correlation matrix ",
                                Omega[RNA])),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::coord_equal()

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(out_file, p, width = 6, height = 5)
    return(out_file)
  }
  invisible(p)
}


#' Plot posterior distributions of individual-effect SDs and correlations
#'
#' Shows ridge/density plots for the 4 \code{sigma_ind_rna} values and
#' the 6 off-diagonal correlations from \code{Omega_rna}.
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list (must have \code{ind_corr = 1})
#' @param out_file  Path to save the PDF (or NULL for interactive)
#' @return ggplot object (invisible)
plot_correlation_densities <- function(fit, stan_data, out_file = NULL) {
  if (!isTRUE(stan_data$ind_corr == 1)) {
    message("ind_corr = 0: no correlation parameters to plot")
    return(invisible(NULL))
  }

  re_labels <- c("tp", "dp", "wp", "wr")

  # --- sigma_ind_rna densities ---
  sigma_draws <- fit$draws(variables = "sigma_ind_rna", format = "draws_df")

  p_sigma <- bayesplot::mcmc_areas(
    fit$draws(variables = "sigma_ind_rna"),
    prob = 0.9
  ) +
    ggplot2::ggtitle("RNA individual-effect SDs") +
    ggplot2::theme_minimal(base_size = 11)

  # --- Omega_rna off-diagonal correlation densities ---
  pairs <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
  var_names <- sapply(pairs, function(ij) {
    paste0("Omega_rna[", ij[1], ",", ij[2], "]")
  })

  p_corr <- bayesplot::mcmc_areas(
    fit$draws(variables = var_names),
    prob = 0.9
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::ggtitle("RNA individual-effect correlations") +
    ggplot2::theme_minimal(base_size = 11)

  p <- patchwork::wrap_plots(p_sigma, p_corr, ncol = 1)

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(out_file, p, width = 8, height = 8)
    return(out_file)
  }
  invisible(p)
}
