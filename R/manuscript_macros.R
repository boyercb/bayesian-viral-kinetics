# ──────────────────────────────────────────────────────────────────────────────
# manuscript_macros.R — LaTeX macros for auto-updating manuscript numbers
# ──────────────────────────────────────────────────────────────────────────────

#' Save manuscript numeric results as LaTeX macros
#'
#' Builds a single \code{.tex} file with \code{\providecommand{...}{...}} entries
#' so main text numbers can be sourced from pipeline outputs.
#'
#' @param convergence   Output of \code{check_convergence()}
#' @param loo_result    Output of \code{compute_loo()}
#' @param waic_result   Output of \code{compute_waic()}
#' @param recovery_check Output of \code{check_recovery()}
#' @param recovery_coverage_summary Optional output of
#'   \code{summarize_recovery_coverage()} for replicate-aggregated metrics.
#' @param param_summary Output of \code{summarize_parameters()}
#' @param kinetics_mcmc Optional CmdStanMCMC fit object. If unavailable,
#'   results macros are derived from \code{param_summary} and convergence object.
#' @param stacked_dat   Stacked analysis dataset
#' @param credible_interval_level Credible interval mass used for results
#'   summaries in manuscript macros. Must be in (0, 1). Default is 0.95.
#' @param out_file      Output path for LaTeX macro file
#' @return Path to output file (invisibly)
save_manuscript_macros <- function(convergence,
                                   loo_result,
                                   waic_result,
                                   recovery_check,
                                   recovery_coverage_summary = NULL,
                                   param_summary,
                                   kinetics_mcmc = NULL,
                                   stacked_dat,
                                   credible_interval_level = 0.95,
                                   out_file = "output/tables/manuscript_numbers.tex") {

  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

  fmt_num <- function(x, digits = 2) {
    format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE, scientific = FALSE)
  }

  fmt_int <- function(x) {
    format(as.integer(round(as.numeric(x))), big.mark = ",", scientific = FALSE, trim = TRUE)
  }

  stan_var_to_latex <- function(x) {
    x <- as.character(x)
    if (grepl("^L_Omega_rna\\[[0-9]+,[0-9]+\\]$", x)) {
      ij <- sub("^L_Omega_rna\\[([0-9]+),([0-9]+)\\]$", "\\1,\\2", x)
      return(paste0("L_{\\Omega,\\text{rna}}[", ij, "]"))
    }
    x
  }

  if (!is.numeric(credible_interval_level) ||
      length(credible_interval_level) != 1 ||
      !is.finite(credible_interval_level) ||
      credible_interval_level <= 0 ||
      credible_interval_level >= 1) {
    stop("credible_interval_level must be a single number in (0, 1).")
  }

  q_prob <- (1 - credible_interval_level) / 2

  get_ci_from_param_summary <- function(var) {
    tables <- c(
      "pop_params", "transformation_params", "error_params",
      "symptom_params", "corr_params"
    )
    for (nm in tables) {
      if (!nm %in% names(param_summary)) next
      tab <- param_summary[[nm]]
      if (!all(c("parameter", "estimate", "ci_lo", "ci_hi") %in% names(tab))) next
      row <- tab[tab$parameter == var, , drop = FALSE]
      if (nrow(row) == 1) {
        return(c(
          med = as.numeric(row$estimate[1]),
          lo = as.numeric(row$ci_lo[1]),
          hi = as.numeric(row$ci_hi[1])
        ))
      }
    }
    NULL
  }

  q_ci <- function(var) {
    from_summary <- get_ci_from_param_summary(var)
    if (!is.null(from_summary) && all(is.finite(from_summary))) {
      return(from_summary)
    }

    if (is.null(kinetics_mcmc)) {
      stop(sprintf("Unable to derive summary for '%s': no fit and no param_summary row.", var))
    }

    d <- as.numeric(kinetics_mcmc$draws(var, format = "matrix"))
    if (length(d) == 0 || !all(is.finite(d))) {
      stop(sprintf("Unable to derive summary for '%s' from fit draws.", var))
    }
    c(
      med = stats::median(d),
      lo = as.numeric(stats::quantile(d, q_prob)),
      hi = as.numeric(stats::quantile(d, 1 - q_prob))
    )
  }

  # Convergence metrics
  rhat_max_var <- "NA"
  ess_bulk_var <- "NA"
  ess_tail_var <- "NA"

  if (!is.null(kinetics_mcmc)) {
    summ <- kinetics_mcmc$summary()
    i_rhat <- which.max(summ$rhat)
    i_essb <- which.min(summ$ess_bulk)
    i_esst <- which.min(summ$ess_tail)

    rhat_max_var <- stan_var_to_latex(summ$variable[i_rhat])
    ess_bulk_var <- stan_var_to_latex(summ$variable[i_essb])
    ess_tail_var <- stan_var_to_latex(summ$variable[i_esst])
  } else {
    if (!is.null(convergence$rhat_warnings) && nrow(convergence$rhat_warnings) > 0) {
      rhat_max_var <- stan_var_to_latex(convergence$rhat_warnings$variable[1])
    }
    if (!is.null(convergence$ess_warnings) && nrow(convergence$ess_warnings) > 0) {
      ess_bulk_var <- stan_var_to_latex(convergence$ess_warnings$variable[1])
      ess_tail_var <- stan_var_to_latex(convergence$ess_warnings$variable[1])
    }
  }

  n_rhat_warn <- if (is.null(convergence$rhat_warnings)) 0 else nrow(convergence$rhat_warnings)

  diag <- tryCatch(
    if (!is.null(kinetics_mcmc)) kinetics_mcmc$diagnostic_summary() else NULL,
    error = function(e) NULL
  )
  ebfmi_vals <- numeric(0)
  if (!is.null(diag)) {
    cand <- c("ebfmi", "e_bfmi")
    hit <- cand[cand %in% names(diag)]
    if (length(hit) > 0) ebfmi_vals <- as.numeric(diag[[hit[1]]])
  }

  # LOO/WAIC (combined)
  loo_comb <- loo_result$combined
  waic_comb <- waic_result$combined

  loo_est <- loo_comb$estimates
  waic_est <- waic_comb$estimates

  pk <- loo_comb$diagnostics$pareto_k
  n_obs <- length(pk)
  n_k_gt_07 <- sum(pk > 0.7, na.rm = TRUE)
  n_k_gt_1 <- sum(pk > 1, na.rm = TRUE)
  n_k_inf <- sum(is.infinite(pk))

  # Symptom-only enrichment among high-k observations
  symptom_only <- with(stacked_dat,
    sym_exist == 1 & rna_exist == 0 & pfu_exist == 0 & lfd_exist == 0
  )
  n_symptom_only <- sum(symptom_only)
  n_symptom_only_highk <- sum(symptom_only & pk > 1, na.rm = TRUE)
  n_symptom_only_inf <- sum(symptom_only & is.infinite(pk), na.rm = TRUE)

  prev_sym <- if (n_obs > 0) (n_symptom_only / n_obs) else NA_real_
  prev_sym_in_highk <- if (n_k_gt_1 > 0) (n_symptom_only_highk / n_k_gt_1) else NA_real_
  enrich_sym_highk <- prev_sym_in_highk / prev_sym

  # ATACCC fraction with k >= 1
  ataccc_mask <- stacked_dat$source == 2
  ataccc_n <- sum(ataccc_mask)
  ataccc_k_ge_1 <- sum(ataccc_mask & pk >= 1, na.rm = TRUE)

  # Recovery summary
  n_recovery <- nrow(recovery_check)
  n_recovery_covered <- sum(recovery_check$covered, na.rm = TRUE)

  recovery_empirical_pct <- 100 * n_recovery_covered / n_recovery
  recovery_empirical_lo <- NA_real_
  recovery_empirical_hi <- NA_real_
  recovery_n_reps <- 1L

  if (!is.null(recovery_coverage_summary) &&
      !is.null(recovery_coverage_summary$overall) &&
      !is.null(recovery_coverage_summary$by_parameter)) {
    recovery_empirical_pct <-
      100 * as.numeric(recovery_coverage_summary$overall$empirical_coverage[1])
    recovery_empirical_lo <-
      100 * as.numeric(recovery_coverage_summary$overall$ci_lo[1])
    recovery_empirical_hi <-
      100 * as.numeric(recovery_coverage_summary$overall$ci_hi[1])
    recovery_n_reps <- max(
      as.integer(recovery_coverage_summary$by_parameter$n_replicates),
      na.rm = TRUE
    )
  }

  # Population/results numbers (configurable CrI level)
  dp <- q_ci("dp_mean_rna")
  wp <- q_ci("wp_mean_rna")
  wr <- q_ci("wr_mean_rna")
  a1_dp <- q_ci("tau_dp[2]")
  a1_wr <- q_ci("tau_wr[2]")

  sigma_rna <- q_ci("sigma_rna")
  sigma_pfu <- q_ci("sigma_pfu")
  fp <- q_ci("fp")

  z0 <- q_ci("zeta_sym_intercept")
  z1 <- q_ci("zeta_sym_pfu")
  z2 <- q_ci("zeta_sym_rna")
  sig_sym <- q_ci("sigma_sym")

  rho_tp_wp <- q_ci("Omega_rna[1,3]")
  rho_tp_wr <- q_ci("Omega_rna[1,4]")
  rho_dp_wp <- q_ci("Omega_rna[2,3]")
  rho_dp_wr <- q_ci("Omega_rna[2,4]")
  rho_wp_wr <- q_ci("Omega_rna[3,4]")

  sd_tp <- q_ci("sigma_ind_rna[1]")
  sd_dp <- q_ci("sigma_ind_rna[2]")
  sd_wp <- q_ci("sigma_ind_rna[3]")
  sd_wr <- q_ci("sigma_ind_rna[4]")

  # Covariate effects (from summarize_parameters output)
  cov <- param_summary$covariate_effects
  covariates <- get("define_covariates", mode = "function")()
  get_cov <- function(parameter, var_name) {
    row <- cov[
      cov$parameter == parameter &
        cov$label == covariates$x_labels[[var_name]],
      ,
      drop = FALSE
    ]
    if (nrow(row) != 1) {
      stop(sprintf("Expected exactly one covariate row for %s/%s", parameter, var_name))
    }
    c(
      est = as.numeric(row$coef),
      lo = as.numeric(row$ci_lo),
      hi = as.numeric(row$ci_hi)
    )
  }

  pct_up <- function(x) as.integer(round(100 * (as.numeric(x) - 1)))
  pct_down <- function(x) as.integer(round(100 * (1 - as.numeric(x))))

  c_delta_dp <- get_cov("dp", "delta")
  c_delta_wp <- get_cov("wp", "delta")
  c_omic_dp <- get_cov("dp", "omicron")
  c_omic_wr <- get_cov("wr", "omicron")
  c_ba_dp <- get_cov("dp", "ba4ba5")
  c_ba_wr <- get_cov("wr", "ba4ba5")
  c_alpha_wp <- get_cov("wp", "alpha")
  c_boost_dp <- get_cov("dp", "vaccinated_boosted")
  c_boost_wp <- get_cov("wp", "vaccinated_boosted")
  c_boost_wr <- get_cov("wr", "vaccinated_boosted")
  c_vax_dp <- get_cov("dp", "vaccinated_unboosted")
  c_vax_wp <- get_cov("wp", "vaccinated_unboosted")
  c_vax_wr <- get_cov("wr", "vaccinated_unboosted")
  c_age50_wr <- get_cov("wr", "age_[50,100)")
  c_rec_dp <- get_cov("dp", "recurrence")
  c_rec_wp <- get_cov("wp", "recurrence")
  c_rec_wr <- get_cov("wr", "recurrence")

  macro <- function(name, value) {
    paste0("\\providecommand{\\", name, "}{", value, "}")
  }

  lines <- c(
    "% Auto-generated by save_manuscript_macros(); do not edit manually.",
    macro("MCNDivergent", fmt_int(convergence$n_divergent)),
    macro("MCNMaxTreeDepth", fmt_int(convergence$n_max_treedepth)),
    macro("MCRhatMax", fmt_num(convergence$rhat_max, 2)),
    macro("MCRhatMaxVar", rhat_max_var),
    macro("MCNRhatGTOneOhOne", fmt_int(n_rhat_warn)),
    macro("MCEssBulkMin", fmt_int(convergence$ess_bulk_min)),
    macro("MCEssBulkMinVar", ess_bulk_var),
    macro("MCEssTailMin", fmt_int(convergence$ess_tail_min)),
    macro("MCEssTailMinVar", ess_tail_var),
    macro("MCELPDLOO", fmt_int(loo_est["elpd_loo", "Estimate"])),
    macro("MCELPDLOOSE", fmt_int(loo_est["elpd_loo", "SE"])),
    macro("MCPLOO", fmt_int(loo_est["p_loo", "Estimate"])),
    macro("MCPLOOSE", fmt_int(loo_est["p_loo", "SE"])),
    macro("MCELPDWAIC", fmt_int(waic_est["elpd_waic", "Estimate"])),
    macro("MCELPDWAICSE", fmt_int(waic_est["elpd_waic", "SE"])),
    macro("MCPWAIC", fmt_int(waic_est["p_waic", "Estimate"])),
    macro("MCPWAICSE", fmt_int(waic_est["p_waic", "SE"])),
    macro("MCNObs", fmt_int(n_obs)),
    macro("MCNKGTSeven", fmt_int(n_k_gt_07)),
    macro("MCPctKGTSeven", fmt_num(100 * n_k_gt_07 / n_obs, 1)),
    macro("MCNKGTOne", fmt_int(n_k_gt_1)),
    macro("MCPctKGTOne", fmt_num(100 * n_k_gt_1 / n_obs, 1)),
    macro("MCNKInf", fmt_int(n_k_inf)),
    macro("MCNSymOnly", fmt_int(n_symptom_only)),
    macro("MCNSymOnlyHighK", fmt_int(n_symptom_only_highk)),
    macro("MCNSymOnlyInf", fmt_int(n_symptom_only_inf)),
    macro("MCPctSymOnlyInf", fmt_num(100 * n_symptom_only_inf / n_symptom_only, 1)),
    macro("MCEnrichSymHighK", fmt_num(enrich_sym_highk, 1)),
    macro("MCATACCCN", fmt_int(ataccc_n)),
    macro("MCATACCCKGEOne", fmt_int(ataccc_k_ge_1)),
    macro("MCPctATACCCKGEOne", fmt_num(100 * ataccc_k_ge_1 / ataccc_n, 1)),
    macro("MCRecoveryCovered", fmt_int(n_recovery_covered)),
    macro("MCRecoveryTotal", fmt_int(n_recovery)),
    macro("MCRecoveryPct", fmt_num(100 * n_recovery_covered / n_recovery, 1)),
    macro("MCRecoveryEmpiricalPct", fmt_num(recovery_empirical_pct, 1)),
    macro("MCRecoveryEmpiricalLo", fmt_num(recovery_empirical_lo, 1)),
    macro("MCRecoveryEmpiricalHi", fmt_num(recovery_empirical_hi, 1)),
    macro("MCRecoveryNReplicates", fmt_int(recovery_n_reps)),
    macro("ResDP", fmt_num(dp["med"], 2)),
    macro("ResDPLo", fmt_num(dp["lo"], 2)),
    macro("ResDPHi", fmt_num(dp["hi"], 2)),
    macro("ResWP", fmt_num(wp["med"], 2)),
    macro("ResWPLo", fmt_num(wp["lo"], 2)),
    macro("ResWPHi", fmt_num(wp["hi"], 2)),
    macro("ResWR", fmt_num(wr["med"], 2)),
    macro("ResWRLo", fmt_num(wr["lo"], 2)),
    macro("ResWRHi", fmt_num(wr["hi"], 2)),
    macro("ResClearToProlif", fmt_num(wr["med"] / wp["med"], 1)),
    macro("ResAOneDP", fmt_num(a1_dp["med"], 2)),
    macro("ResAOneDPLo", fmt_num(a1_dp["lo"], 2)),
    macro("ResAOneDPHi", fmt_num(a1_dp["hi"], 2)),
    macro("ResAOneWR", fmt_num(a1_wr["med"], 2)),
    macro("ResAOneWRLo", fmt_num(a1_wr["lo"], 2)),
    macro("ResAOneWRHi", fmt_num(a1_wr["hi"], 2)),
    macro("ResSigmaRNA", fmt_num(sigma_rna["med"], 2)),
    macro("ResSigmaRNALo", fmt_num(sigma_rna["lo"], 2)),
    macro("ResSigmaRNAHi", fmt_num(sigma_rna["hi"], 2)),
    macro("ResSigmaPFU", fmt_num(sigma_pfu["med"], 2)),
    macro("ResSigmaPFULo", fmt_num(sigma_pfu["lo"], 2)),
    macro("ResSigmaPFUHi", fmt_num(sigma_pfu["hi"], 2)),
    macro("ResFP", fmt_num(fp["med"], 3)),
    macro("ResFPLo", fmt_num(fp["lo"], 3)),
    macro("ResFPHi", fmt_num(fp["hi"], 3)),
    macro("ResZetaZero", fmt_num(z0["med"], 2)),
    macro("ResZetaZeroLo", fmt_num(z0["lo"], 2)),
    macro("ResZetaZeroHi", fmt_num(z0["hi"], 2)),
    macro("ResZetaOne", fmt_num(z1["med"], 2)),
    macro("ResZetaOneLo", fmt_num(z1["lo"], 2)),
    macro("ResZetaOneHi", fmt_num(z1["hi"], 2)),
    macro("ResZetaTwo", fmt_num(z2["med"], 2)),
    macro("ResZetaTwoLo", fmt_num(z2["lo"], 2)),
    macro("ResZetaTwoHi", fmt_num(z2["hi"], 2)),
    macro("ResSigmaSym", fmt_num(sig_sym["med"], 2)),
    macro("ResSigmaSymLo", fmt_num(sig_sym["lo"], 2)),
    macro("ResSigmaSymHi", fmt_num(sig_sym["hi"], 2)),
    macro("ResRhoTPWP", fmt_num(rho_tp_wp["med"], 2)),
    macro("ResRhoTPWPLo", fmt_num(rho_tp_wp["lo"], 2)),
    macro("ResRhoTPWPHi", fmt_num(rho_tp_wp["hi"], 2)),
    macro("ResRhoTPWR", fmt_num(rho_tp_wr["med"], 2)),
    macro("ResRhoTPWRLo", fmt_num(rho_tp_wr["lo"], 2)),
    macro("ResRhoTPWRHi", fmt_num(rho_tp_wr["hi"], 2)),
    macro("ResRhoDPWP", fmt_num(rho_dp_wp["med"], 2)),
    macro("ResRhoDPWPLo", fmt_num(rho_dp_wp["lo"], 2)),
    macro("ResRhoDPWPHi", fmt_num(rho_dp_wp["hi"], 2)),
    macro("ResRhoDPWR", fmt_num(rho_dp_wr["med"], 2)),
    macro("ResRhoDPWRLo", fmt_num(rho_dp_wr["lo"], 2)),
    macro("ResRhoDPWRHi", fmt_num(rho_dp_wr["hi"], 2)),
    macro("ResRhoWPWR", fmt_num(rho_wp_wr["med"], 2)),
    macro("ResRhoWPWRLo", fmt_num(rho_wp_wr["lo"], 2)),
    macro("ResRhoWPWRHi", fmt_num(rho_wp_wr["hi"], 2)),
    macro("ResSDTP", fmt_num(sd_tp["med"], 2)),
    macro("ResSDDP", fmt_num(sd_dp["med"], 2)),
    macro("ResSDWP", fmt_num(sd_wp["med"], 2)),
    macro("ResSDWR", fmt_num(sd_wr["med"], 2)),
    macro("ResCrILevelPct", fmt_int(100 * credible_interval_level)),
    macro("CovDeltaDP", fmt_num(c_delta_dp["est"], 2)),
    macro("CovDeltaDPLo", fmt_num(c_delta_dp["lo"], 2)),
    macro("CovDeltaDPHi", fmt_num(c_delta_dp["hi"], 2)),
    macro("CovDeltaDPPctUp", fmt_int(pct_up(c_delta_dp["est"]))),
    macro("CovDeltaWP", fmt_num(c_delta_wp["est"], 2)),
    macro("CovDeltaWPLo", fmt_num(c_delta_wp["lo"], 2)),
    macro("CovDeltaWPHi", fmt_num(c_delta_wp["hi"], 2)),
    macro("CovDeltaWPPctDown", fmt_int(pct_down(c_delta_wp["est"]))),
    macro("CovOmicronDP", fmt_num(c_omic_dp["est"], 2)),
    macro("CovOmicronDPLo", fmt_num(c_omic_dp["lo"], 2)),
    macro("CovOmicronDPHi", fmt_num(c_omic_dp["hi"], 2)),
    macro("CovOmicronDPPctUp", fmt_int(pct_up(c_omic_dp["est"]))),
    macro("CovOmicronWR", fmt_num(c_omic_wr["est"], 2)),
    macro("CovOmicronWRLo", fmt_num(c_omic_wr["lo"], 2)),
    macro("CovOmicronWRHi", fmt_num(c_omic_wr["hi"], 2)),
    macro("CovOmicronWRPctDown", fmt_int(pct_down(c_omic_wr["est"]))),
    macro("CovBAFourFiveDP", fmt_num(c_ba_dp["est"], 2)),
    macro("CovBAFourFiveDPLo", fmt_num(c_ba_dp["lo"], 2)),
    macro("CovBAFourFiveDPHi", fmt_num(c_ba_dp["hi"], 2)),
    macro("CovBAFourFiveDPPctUp", fmt_int(pct_up(c_ba_dp["est"]))),
    macro("CovBAFourFiveWR", fmt_num(c_ba_wr["est"], 2)),
    macro("CovBAFourFiveWRLo", fmt_num(c_ba_wr["lo"], 2)),
    macro("CovBAFourFiveWRHi", fmt_num(c_ba_wr["hi"], 2)),
    macro("CovBAFourFiveWRPctDown", fmt_int(pct_down(c_ba_wr["est"]))),
    macro("CovAlphaWP", fmt_num(c_alpha_wp["est"], 2)),
    macro("CovAlphaWPLo", fmt_num(c_alpha_wp["lo"], 2)),
    macro("CovAlphaWPHi", fmt_num(c_alpha_wp["hi"], 2)),
    macro("CovAlphaWPPctDown", fmt_int(pct_down(c_alpha_wp["est"]))),
    macro("CovBoostDP", fmt_num(c_boost_dp["est"], 2)),
    macro("CovBoostDPLo", fmt_num(c_boost_dp["lo"], 2)),
    macro("CovBoostDPHi", fmt_num(c_boost_dp["hi"], 2)),
    macro("CovBoostDPPctDown", fmt_int(pct_down(c_boost_dp["est"]))),
    macro("CovBoostWP", fmt_num(c_boost_wp["est"], 2)),
    macro("CovBoostWPLo", fmt_num(c_boost_wp["lo"], 2)),
    macro("CovBoostWPHi", fmt_num(c_boost_wp["hi"], 2)),
    macro("CovBoostWPPctUp", fmt_int(pct_up(c_boost_wp["est"]))),
    macro("CovBoostWR", fmt_num(c_boost_wr["est"], 2)),
    macro("CovBoostWRLo", fmt_num(c_boost_wr["lo"], 2)),
    macro("CovBoostWRHi", fmt_num(c_boost_wr["hi"], 2)),
    macro("CovBoostWRPctDown", fmt_int(pct_down(c_boost_wr["est"]))),
    macro("CovVaxDP", fmt_num(c_vax_dp["est"], 2)),
    macro("CovVaxDPLo", fmt_num(c_vax_dp["lo"], 2)),
    macro("CovVaxDPHi", fmt_num(c_vax_dp["hi"], 2)),
    macro("CovVaxDPPctDown", fmt_int(pct_down(c_vax_dp["est"]))),
    macro("CovVaxWP", fmt_num(c_vax_wp["est"], 2)),
    macro("CovVaxWPLo", fmt_num(c_vax_wp["lo"], 2)),
    macro("CovVaxWPHi", fmt_num(c_vax_wp["hi"], 2)),
    macro("CovVaxWPPctUp", fmt_int(pct_up(c_vax_wp["est"]))),
    macro("CovVaxWR", fmt_num(c_vax_wr["est"], 2)),
    macro("CovVaxWRLo", fmt_num(c_vax_wr["lo"], 2)),
    macro("CovVaxWRHi", fmt_num(c_vax_wr["hi"], 2)),
    macro("CovVaxWRPctDown", fmt_int(pct_down(c_vax_wr["est"]))),
    macro("CovAgeFiftyPlusWR", fmt_num(c_age50_wr["est"], 2)),
    macro("CovAgeFiftyPlusWRLo", fmt_num(c_age50_wr["lo"], 2)),
    macro("CovAgeFiftyPlusWRHi", fmt_num(c_age50_wr["hi"], 2)),
    macro("CovAgeFiftyPlusWRPctUp", fmt_int(pct_up(c_age50_wr["est"]))),
    macro("CovRecDP", fmt_num(c_rec_dp["est"], 2)),
    macro("CovRecDPLo", fmt_num(c_rec_dp["lo"], 2)),
    macro("CovRecDPHi", fmt_num(c_rec_dp["hi"], 2)),
    macro("CovRecDPPctDown", fmt_int(pct_down(c_rec_dp["est"]))),
    macro("CovRecWP", fmt_num(c_rec_wp["est"], 2)),
    macro("CovRecWPLo", fmt_num(c_rec_wp["lo"], 2)),
    macro("CovRecWPHi", fmt_num(c_rec_wp["hi"], 2)),
    macro("CovRecWPPctDown", fmt_int(pct_down(c_rec_wp["est"]))),
    macro("CovRecWR", fmt_num(c_rec_wr["est"], 2)),
    macro("CovRecWRLo", fmt_num(c_rec_wr["lo"], 2)),
    macro("CovRecWRHi", fmt_num(c_rec_wr["hi"], 2)),
    macro("CovRecWRPctDown", fmt_int(pct_down(c_rec_wr["est"]))),
    macro("MCEBFMIMin", "NA"),
    macro("MCEBFMIMax", "NA")
  )

  if (length(ebfmi_vals) > 0 && all(is.finite(ebfmi_vals))) {
    lines <- sub("\\providecommand{\\MCEBFMIMin}{NA}",
      macro("MCEBFMIMin", fmt_num(min(ebfmi_vals), 2)), lines, fixed = TRUE)
    lines <- sub("\\providecommand{\\MCEBFMIMax}{NA}",
      macro("MCEBFMIMax", fmt_num(max(ebfmi_vals), 2)), lines, fixed = TRUE)
  }

  writeLines(lines, out_file)
  invisible(out_file)
}
