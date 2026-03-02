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
    "tau_dp", "tau_wp", "tau_wr",
    "sigma_rna", "sigma_pfu",
    "eta_sym_intercept", "eta_sym_pfu", "eta_sym_rna", "sigma_sym",
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
                 "tau_wr[1]", "tau_wr[2]")
  tau_labels <- c("PFU peak intercept", "PFU peak elasticity",
                  "PFU prolif intercept", "PFU prolif elasticity",
                  "PFU clear intercept", "PFU clear elasticity")
  tau_rvars <- list(
    kinetics$tau_dp[1], kinetics$tau_dp[2],
    kinetics$tau_wp[1], kinetics$tau_wp[2],
    kinetics$tau_wr[1], kinetics$tau_wr[2]
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

  # --- symptom model parameters ---
  sym_names  <- c("eta_sym_intercept", "eta_sym_pfu",
                  "eta_sym_rna", "sigma_sym")
  sym_labels <- c("Symptom intercept", "Symptom: log-PFU effect",
                  "Symptom: log-RNA effect", "Symptom RE SD")
  sym_rvars  <- list(
    kinetics$eta_sym_intercept, kinetics$eta_sym_pfu,
    kinetics$eta_sym_rna, kinetics$sigma_sym
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
