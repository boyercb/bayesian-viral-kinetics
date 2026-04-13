# ──────────────────────────────────────────────────────────────────────────────
# summaries.R — Parameter and data summary tables
# ──────────────────────────────────────────────────────────────────────────────

#' Summarize posterior estimates for covariate effects
#'
#' @param fit       CmdStanMCMC fit object
#' @param stan_data Stan data list (optional; used to detect use_wf/use_smooth flags)
#' @return List of tibbles: covariate_effects (dp, wp, wr),
#'         transformation_params, lfd_params, etc.
summarize_parameters <- function(fit, stan_data = NULL) {
  covariates <- define_covariates()

  variable_list <- c(
    "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
    "beta_dp_rna", "beta_wp_rna", "beta_wr_rna",
    "tau_dp", "tau_wp", "tau_wr", "tau_tp",
    "tau0_lfd", "tau_lfd",
    "sigma_rna", "sigma_pfu",
    "zeta_sym_intercept", "zeta_sym_pfu", "zeta_sym_rna",
    "zeta_sym_postpeak", "zeta_sym_postpeak_rna", "sigma_sym",
    "fp", "fn"
  )
  if (!is.null(stan_data) && isTRUE(stan_data$use_wf == 1)) {
    variable_list <- c(variable_list, "wf_raw")
  }

  kinetics <- posterior::as_draws_rvars(fit$draws(variables = variable_list))

  # helper: summarize a scalar rvar to tibble row with median + 95% CI
  summarize_rv <- function(rv, digits = 2) {
    draws_vec <- posterior::draws_of(rv)
    tibble::tibble(
      estimate = specd(stats::median(draws_vec), digits),
      ci_lo    = specd(stats::quantile(draws_vec, 0.025), digits),
      ci_hi    = specd(stats::quantile(draws_vec, 0.975), digits)
    )
  }

  # --- population-level kinetics means ---
  pop_params <- dplyr::bind_rows(
    dplyr::bind_cols(
      tibble::tibble(parameter = "dp_mean_rna",
                     label = "Peak log-RNA (population mean)"),
      summarize_rv(kinetics$dp_mean_rna)
    ),
    dplyr::bind_cols(
      tibble::tibble(parameter = "wp_mean_rna",
                     label = "Proliferation rate (population mean)"),
      summarize_rv(kinetics$wp_mean_rna)
    ),
    dplyr::bind_cols(
      tibble::tibble(parameter = "wr_mean_rna",
                     label = "Clearance rate (population mean)"),
      summarize_rv(kinetics$wr_mean_rna)
    )
  )

  # --- covariate effects on RNA kinetics (exponentiated = multiplicative) ---
  covariate_effects <- purrr::map(
    c("dp", "wp", "wr"),
    function(x) {
      beta <- kinetics[[paste0("beta_", x, "_rna")]]
      tibble::tibble(
        parameter = x,
        label     = sapply(covariates$x_vars, \(v) covariates$x_labels[[v]]),
        coef      = specd(exp(posterior::E(beta)), 2),
        ci_lo     = specd(exp(posterior::quantile2(beta, 0.025)), 2),
        ci_hi     = specd(exp(posterior::quantile2(beta, 0.975)), 2)
      )
    }
  ) |> dplyr::bind_rows()

  # --- log-affine transformation parameters (RNA → PFU) ---
  tau_names <- c("tau_dp[1]", "tau_dp[2]",
                 "tau_wp[1]", "tau_wp[2]",
                 "tau_wr[1]", "tau_wr[2]",
                 "tau_tp[1]", "tau_tp[2]")
  tau_labels <- c("PFU peak intercept", "PFU peak elasticity",
                  "PFU prolif intercept", "PFU prolif elasticity",
                  "PFU clear intercept", "PFU clear elasticity",
                  "PFU peak time offset", "PFU peak time scaling")
  tau_rvars <- list(
    kinetics$tau_dp[1], kinetics$tau_dp[2],
    kinetics$tau_wp[1], kinetics$tau_wp[2],
    kinetics$tau_wr[1], kinetics$tau_wr[2],
    kinetics$tau_tp[1], kinetics$tau_tp[2]
  )
  transformation_params <- purrr::pmap_dfr(
    list(tau_names, tau_labels, tau_rvars),
    function(nm, lb, rv) {
      dplyr::bind_cols(
        tibble::tibble(parameter = nm, label = lb),
        summarize_rv(rv)
      )
    }
  )

  # --- LFD observation model parameters ---
  lfd_names  <- c("tau0_lfd", "tau_lfd[1]", "tau_lfd[2]",
                   "tau_lfd[3]", "tau_lfd[4]")
  lfd_labels <- c("LFD intercept (logit)", "LFD log-RNA coefficient",
                  "LFD log-PFU coefficient",
                  "LFD post-peak indicator",
                  "LFD post-peak $\\times$ RNA interaction")
  lfd_rvars  <- list(
    kinetics$tau0_lfd, kinetics$tau_lfd[1], kinetics$tau_lfd[2],
    kinetics$tau_lfd[3], kinetics$tau_lfd[4]
  )
  lfd_params <- purrr::pmap_dfr(
    list(lfd_names, lfd_labels, lfd_rvars),
    function(nm, lb, rv) {
      dplyr::bind_cols(
        tibble::tibble(parameter = nm, label = lb),
        summarize_rv(rv)
      )
    }
  )

  # --- symptom model parameters ---
  sym_names  <- c("zeta_sym_intercept", "zeta_sym_pfu",
                  "zeta_sym_rna", "zeta_sym_postpeak",
                  "zeta_sym_postpeak_rna", "sigma_sym")
  sym_labels <- c("Symptom intercept", "Symptom: log-PFU effect",
                  "Symptom: log-RNA effect",
                  "Symptom: post-peak indicator",
                  "Symptom: post-peak $\\times$ RNA interaction",
                  "Symptom RE SD")
  sym_rvars  <- list(
    kinetics$zeta_sym_intercept, kinetics$zeta_sym_pfu,
    kinetics$zeta_sym_rna, kinetics$zeta_sym_postpeak,
    kinetics$zeta_sym_postpeak_rna, kinetics$sigma_sym
  )
  symptom_params <- purrr::pmap_dfr(
    list(sym_names, sym_labels, sym_rvars),
    function(nm, lb, rv) {
      dplyr::bind_cols(
        tibble::tibble(parameter = nm, label = lb),
        summarize_rv(rv)
      )
    }
  )

  # --- measurement error ---
  err_names  <- c("sigma_rna", "sigma_pfu", "fp", "fn")
  err_labels <- c("RNA observation SD", "PFU observation SD",
                  "False positive rate", "False negative rate")
  err_rvars  <- list(
    kinetics$sigma_rna, kinetics$sigma_pfu,
    kinetics$fp, kinetics$fn
  )
  err_digits <- c(2, 2, 3, 3)
  error_params <- purrr::pmap_dfr(
    list(err_names, err_labels, err_rvars, err_digits),
    function(nm, lb, rv, dg) {
      dplyr::bind_cols(
        tibble::tibble(parameter = nm, label = lb),
        summarize_rv(rv, digits = dg)
      )
    }
  )

  # --- flat-top duration (if use_wf is on) ---
  wf_params <- NULL
  if (!is.null(stan_data) && isTRUE(stan_data$use_wf == 1)) {
    wf_mean_rv <- stan_data$prior_wf_mean *
      exp(stan_data$prior_wf_cv * kinetics$wf_raw[1])
    wf_params <- dplyr::bind_cols(
      tibble::tibble(parameter = "wf_mean",
                     label = "Flat-top duration (days)"),
      summarize_rv(wf_mean_rv)
    )
  }

  # --- correlated RNA individual-effect parameters (if ind_corr is on) ---
  corr_params <- NULL
  if (!is.null(stan_data) && isTRUE(stan_data$ind_corr == 1)) {
    corr_vars <- c("sigma_ind_rna", "Omega_rna")
    corr_draws <- posterior::as_draws_rvars(
      fit$draws(variables = corr_vars)
    )

    re_labels <- c("tp (peak time)", "dp (peak height)",
                   "wp (proliferation)", "wr (clearance)")

    # sigma_ind_rna
    sigma_rows <- purrr::map_dfr(1:4, function(i) {
      dplyr::bind_cols(
        tibble::tibble(
          parameter = paste0("sigma_ind_rna[", i, "]"),
          label     = paste0("RNA RE SD: ", re_labels[i])
        ),
        summarize_rv(corr_draws$sigma_ind_rna[i], digits = 3)
      )
    })

    # Omega_rna off-diagonal
    pairs <- list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4))
    omega_rows <- purrr::map_dfr(pairs, function(ij) {
      i <- ij[1]; j <- ij[2]
      dplyr::bind_cols(
        tibble::tibble(
          parameter = paste0("Omega_rna[", i, ",", j, "]"),
          label     = paste0("Corr(", re_labels[i], ", ", re_labels[j], ")")
        ),
        summarize_rv(corr_draws$Omega_rna[i, j], digits = 3)
      )
    })

    corr_params <- dplyr::bind_rows(sigma_rows, omega_rows)
  }

  # --- PFU individual-effect SDs (if ind_effects is on) ---
  pfu_re_params <- NULL
  if (!is.null(stan_data) && isTRUE(stan_data$ind_effects == 1)) {
    pfu_re_draws <- posterior::as_draws_rvars(
      fit$draws(variables = "sigma_ind_pfu")
    )

    re_labels_pfu <- c("tp (peak time)", "dp (peak height)",
                       "wp (proliferation)", "wr (clearance)")

    pfu_re_params <- purrr::map_dfr(1:4, function(i) {
      dplyr::bind_cols(
        tibble::tibble(
          parameter = paste0("sigma_ind_pfu[", i, "]"),
          label     = paste0("PFU RE SD: ", re_labels_pfu[i])
        ),
        summarize_rv(pfu_re_draws$sigma_ind_pfu[i], digits = 3)
      )
    })
  }

  # --- model type label ---
  model_type <- "piecewise"
  if (!is.null(stan_data)) {
    if (isTRUE(stan_data$use_smooth == 1)) model_type <- "smooth"
    if (isTRUE(stan_data$use_wf == 1)) {
      model_type <- paste0(model_type, " + flat-top")
    }
    if (isTRUE(stan_data$ind_corr == 1)) {
      model_type <- paste0(model_type, " + correlated REs")
    }
  }

  list(
    pop_params            = pop_params,
    covariate_effects     = covariate_effects,
    transformation_params = transformation_params,
    lfd_params            = lfd_params,
    symptom_params        = symptom_params,
    error_params          = error_params,
    wf_params             = wf_params,
    corr_params           = corr_params,
    pfu_re_params         = pfu_re_params,
    model_type            = model_type
  )
}


#' Create descriptive Table 1 (covariates by source)
#'
#' @param stacked_dat Stacked analysis dataset
#' @return kable object formatted for LaTeX
summarize_data <- function(stacked_dat) {
  covariates <- define_covariates()

  # covariate counts by source
  table1a <- stacked_dat |>
    dplyr::group_by(id) |>
    dplyr::slice(1) |>
    dplyr::group_by(source) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(covariates$x_vars_w_refs),
        sum,
        .names = "{.col}_n"
      ),
      dplyr::across(
        dplyr::all_of(covariates$x_vars_w_refs),
        function(x) sum(x) / dplyr::n(),
        .names = "{.col}_p"
      )
    ) |>
    tidyr::pivot_longer(-source) |>
    dplyr::mutate(
      type = stringr::str_extract(name, "_([np])$", group = 1),
      name = stringr::str_remove(name, "_[np]$")
    ) |>
    tidyr::pivot_wider(
      names_from = c(type, source),
      values_from = value
    ) |>
    dplyr::mutate(
      name = covariates$x_labels[name],
      dplyr::across(dplyr::starts_with("p_"), ~ .x * 100)
    )

  # observation counts by source
  table1b <- stacked_dat |>
    dplyr::group_by(source) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(covariates$miss_vars),
        sum,
        .names = "{.col}_n"
      ),
      n_obs = dplyr::n()
    )

  list(table1a = table1a, table1b = table1b)
}


#' Render a summary table as LaTeX and save to file
#'
#' @param summary_obj  Output of summarize_data() or summarize_parameters()
#' @param out_file     Path to write the .tex file
#' @return out_file path (invisibly, for targets file tracking)
save_table <- function(summary_obj, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

  if ("table1a" %in% names(summary_obj)) {
    # data summary
    tbl <- knitr::kable(
      summary_obj$table1a,
      format = "latex",
      booktabs = TRUE,
      digits = 1
    )
    writeLines(as.character(tbl), out_file)
  } else if ("covariate_effects" %in% names(summary_obj)) {
    # parameter summary — output all groups
    lines <- character()

    # population-level kinetics
    if ("pop_params" %in% names(summary_obj)) {
      tbl_pop <- knitr::kable(
        summary_obj$pop_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "Population-level kinetics parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Population-level parameters ---",
                 as.character(tbl_pop))
    }

    # covariate effects
    tbl_cov <- knitr::kable(
      summary_obj$covariate_effects,
      format = "latex", booktabs = TRUE,
      caption = "Covariate effects on RNA kinetics (multiplicative scale)"
    )
    lines <- c(lines, "", "% --- Covariate effects ---",
               as.character(tbl_cov))

    # transformation parameters
    if ("transformation_params" %in% names(summary_obj)) {
      tbl_tau <- knitr::kable(
        summary_obj$transformation_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "RNA-to-PFU transformation parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Transformation parameters ---",
                 as.character(tbl_tau))
    }

    # LFD observation model parameters
    if ("lfd_params" %in% names(summary_obj)) {
      tbl_lfd <- knitr::kable(
        summary_obj$lfd_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "LFD observation model parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- LFD parameters ---",
                 as.character(tbl_lfd))
    }

    # symptom parameters
    if ("symptom_params" %in% names(summary_obj)) {
      tbl_sym <- knitr::kable(
        summary_obj$symptom_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "Symptom onset model parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Symptom parameters ---",
                 as.character(tbl_sym))
    }

    # error parameters
    if ("error_params" %in% names(summary_obj)) {
      tbl_err <- knitr::kable(
        summary_obj$error_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "Measurement error parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Error parameters ---",
                 as.character(tbl_err))
    }

    # flat-top duration (if model includes wf)
    if (!is.null(summary_obj$wf_params)) {
      tbl_wf <- knitr::kable(
        summary_obj$wf_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "Flat-top duration parameters",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Flat-top parameters ---",
                 as.character(tbl_wf))
    }

    # correlated RE parameters (if ind_corr is on)
    if (!is.null(summary_obj$corr_params)) {
      tbl_corr <- knitr::kable(
        summary_obj$corr_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "RNA individual-effect correlation structure",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- Correlated RE parameters ---",
                 as.character(tbl_corr))
    }

    # PFU individual-effect SDs
    if (!is.null(summary_obj$pfu_re_params)) {
      tbl_pfu_re <- knitr::kable(
        summary_obj$pfu_re_params[, c("label", "estimate", "ci_lo", "ci_hi")],
        format = "latex", booktabs = TRUE,
        col.names = c("Parameter", "Median", "2.5\\%", "97.5\\%"),
        caption = "PFU individual-effect standard deviations",
        escape = FALSE
      )
      lines <- c(lines, "", "% --- PFU RE parameters ---",
                 as.character(tbl_pfu_re))
    }

    # model type annotation
    if (!is.null(summary_obj$model_type)) {
      lines <- c(lines, "",
                 paste0("% Model type: ", summary_obj$model_type))
    }

    writeLines(lines, out_file)
  }

  invisible(out_file)
}


# ──────────────────────────────────────────────────────────────────────────────
# Supplement table generation
# ──────────────────────────────────────────────────────────────────────────────

#' Generate LaTeX supplement tables from posterior summaries
#'
#' Produces individual .tex files for each supplement table (S3–S7)
#' that can be \input{} from supplement.tex. Each table is generated
#' from the param_summary object so tables stay in sync with the model fit.
#'
#' @param param_summary Output of summarize_parameters()
#' @param out_dir       Directory to write .tex files (default: "output/tables")
#' @return Character vector of output file paths (invisibly, for targets file tracking)
save_supplement_tables <- function(param_summary, out_dir = "output/tables") {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_files <- character()

  # helper: format number for LaTeX
  fmt <- function(x, digits = 2) {
    specd(as.numeric(x), digits)
  }

  # --- Extract population reference values ---
  pop   <- param_summary$pop_params
  dp    <- as.numeric(pop$estimate[1])
  wp    <- as.numeric(pop$estimate[2])
  wr    <- as.numeric(pop$estimate[3])
  dp_lo <- as.numeric(pop$ci_lo[1]); dp_hi <- as.numeric(pop$ci_hi[1])
  wp_lo <- as.numeric(pop$ci_lo[2]); wp_hi <- as.numeric(pop$ci_hi[2])
  wr_lo <- as.numeric(pop$ci_lo[3]); wr_hi <- as.numeric(pop$ci_hi[3])

  # LOD-based reference values (approximate; proper CrIs need posterior draws)
  lod_rna <- ct_to_rna(40, type = "nba") + 0.01
  frac    <- 1 - lod_rna / dp
  frac_lo <- 1 - lod_rna / dp_hi  # note: larger dp => larger frac

  frac_hi <- 1 - lod_rna / dp_lo

  # --- Covariate effects ---
  cov <- param_summary$covariate_effects

  # helper: build one covariate table
  build_cov_table <- function(par_code, title, ref_label, ref_val, ref_lo, ref_hi,
                               lod_label = NULL, lod_val = NULL, lod_lo = NULL, lod_hi = NULL) {
    sub <- cov[cov$parameter == par_code, ]
    rows <- paste0(
      "         ", sub$label, " & ", sub$coef, " & (", sub$ci_lo, ", ", sub$ci_hi, ")\\\\",
      collapse = "\n"
    )
    ref_row <- paste0("         ", ref_label, " & ", fmt(ref_val), " & (",
                       fmt(ref_lo), ", ", fmt(ref_hi), ")\\\\")
    lod_row <- ""
    if (!is.null(lod_label)) {
      lod_row <- paste0("\n         ", lod_label, " & ", fmt(lod_val),
                         " & (", fmt(lod_lo), ", ", fmt(lod_hi), ")\\\\")
    }
    paste0(rows, "\n         \\midrule\n", ref_row, lod_row)
  }

  # --- Individual RE SDs ---
  re_sd_section <- function(rna_corr, pfu_re, sym_params) {
    lines <- character()
    if (!is.null(rna_corr)) {
      # Extract just the SD rows (first 4)
      sd_rows <- rna_corr[grepl("^sigma_ind_rna", rna_corr$parameter), ]
      re_labels <- c("Peak time ($t_p$)", "Peak height ($\\delta$)",
                     "Proliferation ($\\omega_p$)", "Clearance ($\\omega_r$)")
      for (i in seq_len(nrow(sd_rows))) {
        lines <- c(lines, paste0(
          "         RNA: ", re_labels[i], " & ", sd_rows$estimate[i],
          " & (", sd_rows$ci_lo[i], ", ", sd_rows$ci_hi[i], ")\\\\"
        ))
      }
    }
    if (!is.null(pfu_re)) {
      re_labels <- c("Peak time ($t_p'$)", "Peak height ($\\delta'$)",
                     "Proliferation ($\\omega_p'$)", "Clearance ($\\omega_r'$)")
      for (i in seq_len(nrow(pfu_re))) {
        lines <- c(lines, paste0(
          "         PFU: ", re_labels[i], " & ", pfu_re$estimate[i],
          " & (", pfu_re$ci_lo[i], ", ", pfu_re$ci_hi[i], ")\\\\"
        ))
      }
    }
    # Symptom RE SD (scalar)
    if (!is.null(sym_params)) {
      sym_sd <- sym_params[sym_params$parameter == "sigma_sym", ]
      if (nrow(sym_sd) == 1) {
        lines <- c(lines, paste0(
          "         Symptom: Susceptibility ($u_i$) & ", sym_sd$estimate,
          " & (", sym_sd$ci_lo, ", ", sym_sd$ci_hi, ")\\\\"
        ))
      }
    }
    paste0(lines, collapse = "\n")
  }

  re_rows <- re_sd_section(param_summary$corr_params, param_summary$pfu_re_params,
                           param_summary$symptom_params)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S3: Covariate effects on peak RNA
  # ═══════════════════════════════════════════════════════════════════════════
  body_peak <- build_cov_table(
    "dp", "Peak value",
    "Reference value, log [RNA] per ml", dp, dp_lo, dp_hi
  )
  tab_peak <- paste0(
    "% Generated by save_supplement_tables() — do not edit manually\n",
    "\\begin{table}[p]\n",
    "    \\centering\n",
    "    \\caption{Posterior estimates of covariate effects on peak viral RNA concentration. ",
    "Each row reports the posterior median multiplicative effect $\\exp(\\hat\\beta)$ and 95\\% ",
    "credible interval (CrI) for the indicated covariate level relative to the reference category. ",
    "The reference value row reports the estimated population-mean peak concentration for the ",
    "reference category.}\n",
    "    \\label{tab:cov_peak}\n",
    "    \\begin{tabular}{lccc}\n",
    "    \\toprule\n",
    "     & \\multicolumn{2}{c}{Peak value} \\\\\n",
    "    Characteristic & $\\exp(\\beta)$ & 95\\% CrI\\\\\n",
    "    \\midrule\n",
    body_peak, "\n",
    "     \\bottomrule\n",
    "    \\end{tabular}\n",
    "\\end{table}\n"
  )
  f <- file.path(out_dir, "supp_cov_peak.tex")
  writeLines(tab_peak, f)
  out_files <- c(out_files, f)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S4: Covariate effects on proliferation
  # ═══════════════════════════════════════════════════════════════════════════
  body_prolif <- build_cov_table(
    "wp", "Proliferation duration",
    "Reference value, days", wp, wp_lo, wp_hi,
    "Reference value (lod), days", wp * frac, wp_lo * frac_hi, wp_hi * frac_lo
  )
  tab_prolif <- paste0(
    "% Generated by save_supplement_tables() — do not edit manually\n",
    "\\begin{table}[p]\n",
    "    \\caption{Posterior estimates of covariate effects on proliferation phase duration. ",
    "Format as in Table~\\ref{tab:cov_peak}. Values $>1$ indicate longer proliferation ",
    "(slower rise to peak) relative to the reference category; values $<1$ indicate faster ",
    "proliferation. Reference values are the estimated population-mean proliferation duration ",
    "(in days) for the reference category, reported both for the full trajectory (from infection ",
    "onset to peak) and from the limit of detection to peak.}\n",
    "    \\label{tab:cov_prolif}\n",
    "    \\centering\n",
    "    \\begin{tabular}{lcc}\n",
    "     \\toprule\n",
    "     & \\multicolumn{2}{c}{Proliferation duration } \\\\\n",
    "     Characteristic & $\\exp(\\beta)$ & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    body_prolif, "\n",
    "     \\bottomrule\n",
    "     \\end{tabular}\n",
    "\\end{table}\n"
  )
  f <- file.path(out_dir, "supp_cov_prolif.tex")
  writeLines(tab_prolif, f)
  out_files <- c(out_files, f)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S5: Covariate effects on clearance
  # ═══════════════════════════════════════════════════════════════════════════
  body_clear <- build_cov_table(
    "wr", "Clearance duration",
    "Reference value, days", wr, wr_lo, wr_hi,
    "Reference value (lod), days", wr * frac, wr_lo * frac_hi, wr_hi * frac_lo
  )
  tab_clear <- paste0(
    "% Generated by save_supplement_tables() — do not edit manually\n",
    "\\begin{table}[p]\n",
    "    \\centering\n",
    "    \\caption{Posterior estimates of covariate effects on clearance phase duration. ",
    "Format as in Table~\\ref{tab:cov_peak}. Values $>1$ indicate slower clearance ",
    "(longer time from peak to limit of detection) relative to the reference category; ",
    "values $<1$ indicate faster clearance.}\n",
    "    \\label{tab:cov_clearance}\n",
    "    \\begin{tabular}{lcc}\n",
    "     \\toprule\n",
    "     & \\multicolumn{2}{c}{Clearance duration} \\\\\n",
    "     Characteristic & $\\exp(\\beta)$ & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    body_clear, "\n",
    "     \\bottomrule\n",
    "     \\end{tabular}\n",
    "\\end{table}\n"
  )
  f <- file.path(out_dir, "supp_cov_clearance.tex")
  writeLines(tab_clear, f)
  out_files <- c(out_files, f)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S6: RNA-to-PFU transformation parameters
  # ═══════════════════════════════════════════════════════════════════════════
  tau <- param_summary$transformation_params
  # rows: dp intercept/elasticity, wp intercept/elasticity, wr intercept/elasticity,
  #        tp offset/scaling
  tab_trans <- paste0(
    "% Generated by save_supplement_tables() — do not edit manually\n",
    "\\begin{table}[p]\n",
    "    \\centering\n",
    "    \\caption{RNA-to-PFU transformation: posterior estimates of the log-affine ",
    "parameters linking infectious virus trajectory parameters to viral RNA trajectory ",
    "parameters. The transformation is $\\theta'_k = \\exp(a_{0,k} + a_{1,k} \\log \\theta_k)$ ",
    "for peak height, proliferation duration, and clearance duration, where $\\theta_k$ is the ",
    "RNA parameter and $\\theta'_k$ is the corresponding PFU parameter. The elasticity $a_1$ ",
    "represents the percent change in the PFU parameter for a 1\\% change in the RNA parameter; ",
    "values $<1$ indicate sub-proportional scaling. ",
    "The peak time transformation is affine: $t'_p = a_{0,tp} + a_{1,tp} \\cdot t_p$.}\n",
    "    \\label{tab:transformation}\n",
    "    \\begin{tabular}{lcccc}\n",
    "     \\toprule\n",
    "     Parameter & Intercept ($a_0$) & 95\\% CrI & Elasticity ($a_1$) & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    "     Peak ($\\delta'$) & $", tau$estimate[1], "$ & ($", tau$ci_lo[1], "$, $", tau$ci_hi[1], "$) & ",
      tau$estimate[2], " & (", tau$ci_lo[2], ", ", tau$ci_hi[2], ")\\\\\n",
    "     Proliferation ($\\omega_p'$) & $", tau$estimate[3], "$ & ($", tau$ci_lo[3], "$, $", tau$ci_hi[3], "$) & ",
      tau$estimate[4], " & (", tau$ci_lo[4], ", ", tau$ci_hi[4], ")\\\\\n",
    "     Clearance ($\\omega_r'$) & $", tau$estimate[5], "$ & ($", tau$ci_lo[5], "$, $", tau$ci_hi[5], "$) & ",
      tau$estimate[6], " & (", tau$ci_lo[6], ", ", tau$ci_hi[6], ")\\\\\n",
    "     \\midrule\n",
    "     & Offset ($a_0$) & 95\\% CrI & Scaling ($a_1$) & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    "     Peak time ($t_p'$) & $", tau$estimate[7], "$ & ($", tau$ci_lo[7], "$, $", tau$ci_hi[7], "$) & ",
      tau$estimate[8], " & (", tau$ci_lo[8], ", ", tau$ci_hi[8], ")\\\\\n",
    "     \\bottomrule\n",
    "     \\end{tabular}\n",
    "\\end{table}\n"
  )
  f <- file.path(out_dir, "supp_transformation.tex")
  writeLines(tab_trans, f)
  out_files <- c(out_files, f)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S7: LFD observation model parameters
  # ═══════════════════════════════════════════════════════════════════════════
  lfd <- param_summary$lfd_params
  # Row 1 = intercept (tau0_lfd), Row 2 = tau_lfd[1] (RNA), Row 3 = tau_lfd[2] (PFU)
  # Exponentiate the coefficients for odds ratios
  or_rna    <- fmt(exp(as.numeric(lfd$estimate[2])))
  or_rna_lo <- fmt(exp(as.numeric(lfd$ci_lo[2])))
  or_rna_hi <- fmt(exp(as.numeric(lfd$ci_hi[2])))
  or_pfu    <- fmt(exp(as.numeric(lfd$estimate[3])))
  or_pfu_lo <- fmt(exp(as.numeric(lfd$ci_lo[3])))
  or_pfu_hi <- fmt(exp(as.numeric(lfd$ci_hi[3])))

  # Post-peak indicator interaction terms (rows 4 and 5)
  or_ind    <- fmt(exp(as.numeric(lfd$estimate[4])))
  or_ind_lo <- fmt(exp(as.numeric(lfd$ci_lo[4])))
  or_ind_hi <- fmt(exp(as.numeric(lfd$ci_hi[4])))
  or_ix     <- fmt(exp(as.numeric(lfd$estimate[5])))
  or_ix_lo  <- fmt(exp(as.numeric(lfd$ci_lo[5])))
  or_ix_hi  <- fmt(exp(as.numeric(lfd$ci_hi[5])))

  tab_lfd <- paste0(
    "% Generated by save_supplement_tables() — do not edit manually\n",
    "\\begin{table}[p]\n",
    "    \\centering\n",
    "    \\caption{Posterior estimates of LFD logistic observation model parameters. ",
    "The exponentiated coefficients $\\exp(\\gamma)$ represent odds ratios. ",
    "The post-peak indicator $\\mathbbm{1}(t \\geq t_p)$ and its interaction with log RNA ",
    "capture the asymmetry between the proliferation and clearance phases of antigen ",
    "accumulation. The intercept is on the logit scale.}\n",
    "    \\label{tab:lfd}\n",
    "    \\begin{tabular}{lcc}\n",
    "     \\toprule\n",
    "     Predictor & $\\exp(\\gamma)$ & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    "     log RNA copies ($\\gamma_1$) & ", or_rna, " & (", or_rna_lo, ", ", or_rna_hi, ")\\\\\n",
    "     log PFU culturable virus ($\\gamma_2$) & ", or_pfu, " & (", or_pfu_lo, ", ", or_pfu_hi, ")\\\\\n",
    "     Post-peak indicator ($\\gamma_3$) & ", or_ind, " & (", or_ind_lo, ", ", or_ind_hi, ")\\\\\n",
    "     Post-peak $\\times$ log RNA ($\\gamma_4$) & ", or_ix, " & (", or_ix_lo, ", ", or_ix_hi, ")\\\\\n",
    "     \\midrule\n",
    "     & Coefficient & 95\\% CrI\\\\\n",
    "     \\midrule\n",
    "     Intercept ($\\gamma_0$) & ", lfd$estimate[1], " & (", lfd$ci_lo[1], ", ", lfd$ci_hi[1], ")\\\\\n",
    "     \\bottomrule\n",
    "     \\end{tabular}\n",
    "\\end{table}\n"
  )
  f <- file.path(out_dir, "supp_lfd.tex")
  writeLines(tab_lfd, f)
  out_files <- c(out_files, f)

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S8: Symptom onset model parameters
  # ═══════════════════════════════════════════════════════════════════════════
  if ("symptom_params" %in% names(param_summary)) {
    sym <- param_summary$symptom_params
    # Rows: 1=intercept, 2=log-PFU, 3=log-RNA, 4=post-peak, 5=post-peak×RNA, 6=RE SD
    # Exponentiate hazard coefficients for interpretability
    hr_pfu    <- fmt(exp(as.numeric(sym$estimate[2])))
    hr_pfu_lo <- fmt(exp(as.numeric(sym$ci_lo[2])))
    hr_pfu_hi <- fmt(exp(as.numeric(sym$ci_hi[2])))
    hr_rna    <- fmt(exp(as.numeric(sym$estimate[3])))
    hr_rna_lo <- fmt(exp(as.numeric(sym$ci_lo[3])))
    hr_rna_hi <- fmt(exp(as.numeric(sym$ci_hi[3])))
    hr_pp     <- fmt(exp(as.numeric(sym$estimate[4])))
    hr_pp_lo  <- fmt(exp(as.numeric(sym$ci_lo[4])))
    hr_pp_hi  <- fmt(exp(as.numeric(sym$ci_hi[4])))
    hr_ppx    <- fmt(exp(as.numeric(sym$estimate[5])))
    hr_ppx_lo <- fmt(exp(as.numeric(sym$ci_lo[5])))
    hr_ppx_hi <- fmt(exp(as.numeric(sym$ci_hi[5])))

    tab_sym <- paste0(
      "% Generated by save_supplement_tables() \u2014 do not edit manually\n",
      "\\begin{table}[p]\n",
      "    \\centering\n",
      "    \\caption{Posterior estimates of discrete-time complementary log-log symptom ",
      "onset model parameters. The hazard of symptom onset on day $t$ is ",
      "$h_t = 1 - \\exp\\bigl(-\\exp(\\zeta_0 + \\zeta_1 \\log V_t/s + \\zeta_2 \\log R_t/s ",
      "+ \\zeta_3 \\mathbbm{1}(t \\geq t_p) + \\zeta_4 \\mathbbm{1}(t \\geq t_p) \\log R_t/s + u_i)\\bigr)$, ",
      "where $V_t$ and $R_t$ are the latent infectious virus and viral RNA trajectories, ",
      "$\\mathbbm{1}(t \\geq t_p)$ is a post-peak indicator capturing the immune activation lag, ",
      "and $u_i \\sim N(0, \\sigma_{\\text{sym}}^2)$ is an individual random effect. ",
      "Exponentiated coefficients $\\exp(\\zeta)$ represent the multiplicative change in the ",
      "complementary log-log hazard per unit increase in the predictor.}\n",
      "    \\label{tab:symptom}\n",
      "    \\begin{tabular}{lcc}\n",
      "     \\toprule\n",
      "     Predictor & $\\exp(\\zeta)$ & 95\\% CrI\\\\\n",
      "     \\midrule\n",
      "     log PFU culturable virus ($\\zeta_1$) & ", hr_pfu, " & (", hr_pfu_lo, ", ", hr_pfu_hi, ")\\\\\n",
      "     log RNA copies ($\\zeta_2$) & ", hr_rna, " & (", hr_rna_lo, ", ", hr_rna_hi, ")\\\\\n",
      "     Post-peak indicator ($\\zeta_3$) & ", hr_pp, " & (", hr_pp_lo, ", ", hr_pp_hi, ")\\\\\n",
      "     Post-peak $\\times$ log RNA ($\\zeta_4$) & ", hr_ppx, " & (", hr_ppx_lo, ", ", hr_ppx_hi, ")\\\\\n",
      "     \\midrule\n",
      "     & Coefficient & 95\\% CrI\\\\\n",
      "     \\midrule\n",
      "     Intercept ($\\zeta_0$) & ", sym$estimate[1], " & (", sym$ci_lo[1], ", ", sym$ci_hi[1], ")\\\\\n",
      "     Individual RE SD ($\\sigma_{\\text{sym}}$) & ", sym$estimate[6], " & (", sym$ci_lo[6], ", ", sym$ci_hi[6], ")\\\\\n",
      "     \\bottomrule\n",
      "     \\end{tabular}\n",
      "\\end{table}\n"
    )
    f <- file.path(out_dir, "supp_symptom.tex")
    writeLines(tab_sym, f)
    out_files <- c(out_files, f)
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Table S9: Individual random effect standard deviations
  # ═══════════════════════════════════════════════════════════════════════════
  if (!is.null(param_summary$corr_params) || !is.null(param_summary$pfu_re_params)) {
    tab_re <- paste0(
      "% Generated by save_supplement_tables() — do not edit manually\n",
      "\\begin{table}[p]\n",
      "    \\centering\n",
      "    \\caption{Posterior estimates of individual random effect standard deviations ",
      "for the RNA and PFU trajectories. RNA random effects are correlated (Section~8.5 ",
      "of the main text); PFU random effects are independent with a half-normal(0, 0.3) ",
      "prior. Each row reports the posterior median and 95\\% credible interval for the ",
      "standard deviation of individual-level deviations from the population mean on the ",
      "log scale.}\n",
      "    \\label{tab:re_sd}\n",
      "    \\begin{tabular}{lcc}\n",
      "     \\toprule\n",
      "     Parameter & SD & 95\\% CrI\\\\\n",
      "     \\midrule\n",
      re_rows, "\n",
      "     \\bottomrule\n",
      "     \\end{tabular}\n",
      "\\end{table}\n"
    )
    f <- file.path(out_dir, "supp_re_sd.tex")
    writeLines(tab_re, f)
    out_files <- c(out_files, f)
  }

  invisible(out_files)
}
