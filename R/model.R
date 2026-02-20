# ──────────────────────────────────────────────────────────────────────────────
# model.R — Stan model fitting helpers, predictions, and prior-predictive
# ──────────────────────────────────────────────────────────────────────────────

#' Convert a flat CmdStan draw (bracket-notation names) into a nested init list
#'
#' Stan CSV / draws_matrix columns use names like \code{sigma_rna},
#' \code{z_sym[4]}, \code{L_Omega_rna[2,1]}, \code{z_ind_rna[3,100]}.
#' Stan's \code{init} argument expects a named R list where vectors are
#' numeric vectors and matrices are R matrices.  This function performs the
#' conversion.
#'
#' @param vals Numeric vector of parameter values.
#' @param nms  Character vector of bracket-notation parameter names (same
#'   length as \code{vals}).
#' @return Named list suitable for passing to CmdStan's \code{init}.
unflatten_stan_draw <- function(vals, nms) {
  # Parse each name into base + indices
  parsed <- regmatches(nms, regexec("^([^\\[]+)(\\[(.+)\\])?$", nms))

  bases <- vapply(parsed, `[`, character(1), 2L)
  idx_str <- vapply(parsed, function(p) ifelse(is.na(p[4]), "", p[4]), character(1))

  unique_bases <- unique(bases)
  out <- list()

  for (b in unique_bases) {
    sel <- which(bases == b)
    if (length(sel) == 1 && idx_str[sel] == "") {
      # Scalar parameter
      out[[b]] <- vals[sel]
    } else {
      # Array / vector / matrix parameter
      idxs <- strsplit(idx_str[sel], ",")
      ndim <- length(idxs[[1]])
      if (ndim == 1) {
        # 1-D: vector / array
        ix <- as.integer(unlist(idxs))
        v <- rep(NA_real_, max(ix))
        v[ix] <- vals[sel]
        out[[b]] <- v
      } else if (ndim == 2) {
        # 2-D: matrix (row, col)
        rows <- vapply(idxs, function(x) as.integer(x[1]), integer(1))
        cols <- vapply(idxs, function(x) as.integer(x[2]), integer(1))
        m <- matrix(NA_real_, nrow = max(rows), ncol = max(cols))
        for (k in seq_along(sel)) m[rows[k], cols[k]] <- vals[sel[k]]
        out[[b]] <- m
      }
    }
  }
  out
}


#' Build initialisation list for MCMC sampling
#'
#' @param stan_data  Stan data list (from build_stan_data)
#' @param chains     Number of chains
#' @return List of lists, one per chain
build_init <- function(stan_data, chains = 4) {
  N_ind <- sum(stan_data$M)
  P     <- stan_data$P
  K     <- stan_data$K
  lapply(seq_len(chains), function(x) {
    init <- list(
      # PFU individual effects (always independent)
      tp_i_pfu  = rep(0, N_ind),
      dp_i_pfu  = rep(0, N_ind),
      wp_i_pfu  = rep(0, N_ind),
      wr_i_pfu  = rep(0, N_ind),
      z_sym     = rep(0, N_ind),

      # population kinetics (NCP raw values → at 0 = prior mode)
      dp_raw = 0, wp_raw = 0, wr_raw = 0,

      # observation noise
      sigma_rna = 2,
      sigma_pfu = 3,

      # test error (near prior mode)
      fp = 0.02,
      fn = 0.01,

      # symptom hazard
      eta_sym_intercept = -3,
      eta_sym_pfu = 0.3,
      eta_sym_rna = 0.3,
      sigma_sym   = 0.5,

      # RNA→PFU transform (log-affine)
      tau_tp = c(0, 1),
      tau_dp = c(-1, 1),
      tau_wp = c(-1, 1),
      tau_wr = c(-1, 1),

      # LFD model
      tau0_lfd_raw = 0,
      tau_lfd  = c(0.1, 0.1),

      # viral culture models
      alpha_tcid50 = c(0, 0.1, 0, 0.1),
      alpha_cult   = c(0, 0.1)
    )

    # RNA individual effects: correlated (Cholesky NCP) vs independent
    if (isTRUE(stan_data$ind_corr == 1)) {
      init$z_ind_rna     <- matrix(0, 4, N_ind)
      init$L_Omega_rna   <- diag(4)
      init$sigma_ind_rna <- rep(1, 4)
    } else {
      init$tp_i_rna <- rep(0, N_ind)
      init$dp_i_rna <- rep(0, N_ind)
      init$wp_i_rna <- rep(0, N_ind)
      init$wr_i_rna <- rep(0, N_ind)
    }

    # covariate betas (only exist when adj_* flags are on)
    if (stan_data$adj_rna) {
      init$beta_dp_rna <- rep(0, P)
      init$beta_wp_rna <- rep(0, P)
      init$beta_wr_rna <- rep(0, P)
    }
    if (stan_data$adj_pfu) {
      init$beta_dp_pfu <- rep(0, P)
      init$beta_wp_pfu <- rep(0, P)
      init$beta_wr_pfu <- rep(0, P)
    }
    if (stan_data$adj_lfd) init$beta_lfd <- rep(0, P)
    if (stan_data$adj_sym) init$beta_sym <- rep(0, P)

    # source effects (only exist when source_* flags are on)
    if (stan_data$source_rna) {
      init$tp_k_rna <- rep(0, K)
      init$dp_k_rna <- rep(0, K)
      init$wp_k_rna <- rep(0, K)
      init$wr_k_rna <- rep(0, K)
    }
    if (stan_data$source_pfu) {
      init$tp_k_pfu <- rep(0, K)
      init$dp_k_pfu <- rep(0, K)
      init$wp_k_pfu <- rep(0, K)
      init$wr_k_pfu <- rep(0, K)
    }
    if (stan_data$source_lfd) init$lfd_k    <- rep(0, K)
    if (stan_data$source_sym) init$to_k_sym <- rep(0, K)

    # flat-top parameters (only when use_wf = 1)
    if (stan_data$use_wf) {
      init$wf_raw <- 0
      if (stan_data$ind_effects) init$wf_i <- rep(0, N_ind)
      if (stan_data$source_rna)  init$wf_k <- rep(0, K)
    }

    init
  })
}


#' Add random jitter to an init list so chains start at different points
#'
#' Perturbs each element of a single-chain init list by small random noise.
#' Respects constraints: bounded params are jittered on unconstrained scale,
#' Cholesky factors are left untouched (jitter goes into z_ind_rna instead).
#'
#' @param init   A single chain's init list (from build_init or MAP).
#' @param sd_scalar Jitter SD for scalar unconstrained params (default 0.1).
#' @param sd_vec    Jitter SD for vector/coefficient params (default 0.05).
#' @param sd_re     Jitter SD for individual random effects (default 0.01).
#' @param seed  Optional seed for reproducibility.
#' @return Jittered init list.
jitter_init <- function(init, sd_scalar = 0.1, sd_vec = 0.05,
                        sd_re = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Helper: jitter an array/vector by sd
  jit <- function(x, sd) {
    if (is.matrix(x)) {
      x + matrix(rnorm(length(x), 0, sd), nrow = nrow(x), ncol = ncol(x))
    } else {
      x + rnorm(length(x), 0, sd)
    }
  }

  # Unconstrained scalar population params
  for (nm in c("dp_raw", "wp_raw", "wr_raw", "wf_raw",
               "eta_sym_intercept", "tau0_lfd_raw")) {
    if (!is.null(init[[nm]])) init[[nm]] <- init[[nm]] + rnorm(length(init[[nm]]), 0, sd_scalar)
  }

  # Positive scalars: jitter on log scale
  for (nm in c("sigma_rna", "sigma_pfu", "sigma_sym")) {
    if (!is.null(init[[nm]])) {
      init[[nm]] <- init[[nm]] * exp(rnorm(1, 0, sd_scalar))
    }
  }

  # Bounded [0,1]: jitter on logit scale
  for (nm in c("fp", "fn")) {
    if (!is.null(init[[nm]])) {
      logit_val <- qlogis(init[[nm]])
      init[[nm]] <- plogis(logit_val + rnorm(1, 0, sd_scalar))
    }
  }

  # Positive coefficients (eta_sym_pfu, eta_sym_rna): jitter on log scale
  for (nm in c("eta_sym_pfu", "eta_sym_rna")) {
    if (!is.null(init[[nm]]) && init[[nm]] > 0) {
      init[[nm]] <- init[[nm]] * exp(rnorm(1, 0, sd_scalar))
    }
  }

  # sigma_ind_rna: jitter on log scale (positive vector)
  if (!is.null(init$sigma_ind_rna)) {
    init$sigma_ind_rna <- init$sigma_ind_rna *
      exp(rnorm(length(init$sigma_ind_rna), 0, sd_scalar))
  }

  # Vector/coefficient params
  for (nm in c("tau_tp", "tau_dp", "tau_wp", "tau_wr", "tau_lfd",
               "alpha_tcid50", "alpha_cult",
               "beta_dp_rna", "beta_wp_rna", "beta_wr_rna",
               "beta_dp_pfu", "beta_wp_pfu", "beta_wr_pfu",
               "beta_lfd", "beta_sym")) {
    if (!is.null(init[[nm]])) init[[nm]] <- jit(init[[nm]], sd_vec)
  }

  # Individual random effects (large arrays)
  for (nm in c("tp_i_pfu", "dp_i_pfu", "wp_i_pfu", "wr_i_pfu",
               "tp_i_rna", "dp_i_rna", "wp_i_rna", "wr_i_rna",
               "z_sym", "wf_i", "z_ind_rna")) {
    if (!is.null(init[[nm]])) init[[nm]] <- jit(init[[nm]], sd_re)
  }

  # Source effects
  for (nm in c("tp_k_rna", "dp_k_rna", "wp_k_rna", "wr_k_rna", "wf_k",
               "tp_k_pfu", "dp_k_pfu", "wp_k_pfu", "wr_k_pfu",
               "lfd_k", "to_k_sym")) {
    if (!is.null(init[[nm]])) init[[nm]] <- jit(init[[nm]], sd_vec)
  }

  # L_Omega_rna: leave untouched — valid Cholesky factor
  init
}


#' Find MAP estimate and return jittered init lists for MCMC
#'
#' Runs CmdStan's L-BFGS optimizer (\code{$optimize()}) starting from the
#' prior-mode init (\code{build_init()}), then creates \code{chains} jittered
#' copies of the MAP to give MCMC chains diverse starting points near the mode.
#'
#' @param mod        Compiled CmdStanModel object.
#' @param stan_data  Stan data list.
#' @param chains     Number of MCMC chains (number of jittered copies).
#' @param seed       RNG seed for the optimizer.
#' @return List of \code{chains} init lists (jittered copies of the MAP).
map_init <- function(mod, stan_data, chains = 4, seed = 42) {
  base <- build_init(stan_data, chains = 1)[[1]]
  message("Running L-BFGS optimizer to find MAP...")
  t0 <- proc.time()

  opt <- tryCatch(
    mod$optimize(
      data = stan_data,
      init = base,
      algorithm = "lbfgs",
      iter = 2000,
      seed = seed
    ),
    error = function(e) {
      warning("Optimizer failed (", conditionMessage(e),
              "); falling back to build_init()")
      NULL
    }
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  message(sprintf("Optimizer finished in %.1f s", elapsed))

  if (is.null(opt)) {
    message("Falling back to prior-mode init with jitter.")
    return(lapply(seq_len(chains), function(i)
      jitter_init(base, seed = seed + i)))
  }

  # Extract MAP draws and convert to init list
  draws_mat <- opt$draws(format = "draws_matrix")
  all_cols  <- colnames(draws_mat)

  # Only keep actual parameters — skip transformed params and GQ
  skip_cols <- c("lp__", "lp_approx__",
                 "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
                 "dp_mean_pfu", "wp_mean_pfu", "wr_mean_pfu",
                 "tau0_lfd", "total_ll")
  skip_pat  <- "^(Omega_rna)"
  keep      <- all_cols[!all_cols %in% skip_cols & !grepl(skip_pat, all_cols)]

  map_list <- unflatten_stan_draw(as.numeric(draws_mat[1, keep]), keep)

  map_lp <- as.numeric(draws_mat[1, "lp__"])
  message(sprintf("MAP lp__ = %.1f", map_lp))

  # Return jittered copies of the MAP for each chain
  lapply(seq_len(chains), function(i) jitter_init(map_list, seed = seed + i))
}


#' Fit the generative Stan model via MCMC
#'
#' @param stan_file   Path to the .stan file
#' @param stan_data   Stan data list
#' @param chains      Number of chains
#' @param iter_warmup Warmup iterations per chain
#' @param iter_sampling Sampling iterations per chain
#' @param adapt_delta Target acceptance rate
#' @param max_treedepth Maximum tree depth for NUTS
#' @param init_method Initialization strategy: \code{"map"} (default) runs
#'   L-BFGS optimizer then jitters, \code{"prior"} uses prior-mode init with
#'   jitter, \code{"pathfinder"} uses pathfinder (known to fail — see
#'   DIAGNOSIS_PATHFINDER.md).
#' @param threads_per_chain Number of within-chain threads for reduce_sum.
#'   Requires OpenMP (macOS: \code{brew install libomp}).
#'   When \code{> 1}, the model is recompiled with \code{stan_threads = TRUE}.
#'   Speedup is approximately linear up to the number of physical cores per
#'   chain.  With 4 parallel chains on a 16-core machine, use
#'   \code{threads_per_chain = 4} (4 chains × 4 threads = 16 cores).
#'   Default is 1 (no threading, sequential reduce_sum fallback).
#' @param ...         Additional arguments passed to $sample()
#' @return CmdStanMCMC fit object
fit_model <- function(stan_file, stan_data,
                      chains = 4, iter_warmup = 1000, iter_sampling = 2000,
                      adapt_delta = 0.95, max_treedepth = 12,
                      init_method = "map",
                      threads_per_chain = 1,
                      ...) {
  # Compile with threading support when threads_per_chain > 1.
  # reduce_sum falls back to sequential when STAN_THREADS is not defined,
  # so threads_per_chain = 1 (default) works without OpenMP.
  cpp_opts <- if (threads_per_chain > 1) list(stan_threads = TRUE) else list()
  mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = cpp_opts)

  # --- Initialization strategy ---
  if (init_method == "map") {
    init <- map_init(mod, stan_data, chains = chains)
  } else if (init_method == "pathfinder") {
    # WARNING: Pathfinder fails on this model — non-finite gradients,
    # converges to degenerate mode (dp ≈ 0 instead of ~17).
    message("Running pathfinder for MCMC initialisation...")
    init <- build_init(stan_data, chains)  # fallback default
    pf <- tryCatch(
      mod$pathfinder(
        data = stan_data,
        num_paths = 8,
        draws = 1000,
        psis_resample = FALSE,
        seed = 42
      ),
      error = function(e) {
        warning("Pathfinder failed (", conditionMessage(e),
                "); falling back to build_init()")
        NULL
      }
    )
    if (!is.null(pf)) {
      pf_init <- tryCatch({
        draws_mat <- pf$draws(format = "draws_matrix")
        lp_vals   <- as.vector(draws_mat[, "lp__"])
        best_idx  <- order(lp_vals, decreasing = TRUE)
        best_idx  <- unique(best_idx)[seq_len(min(chains, length(unique(best_idx))))]
        skip_cols <- c("lp__", "lp_approx__", "path__",
                       "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
                       "dp_mean_pfu", "wp_mean_pfu", "wr_mean_pfu",
                       "tau0_lfd", "total_ll")
        all_cols <- colnames(draws_mat)
        skip_pat <- "^(Omega_rna)"
        keep <- all_cols[!all_cols %in% skip_cols & !grepl(skip_pat, all_cols)]
        lapply(best_idx, function(i) {
          unflatten_stan_draw(as.numeric(draws_mat[i, keep]), keep)
        })
      }, error = function(e) { NULL })
      if (!is.null(pf_init) && length(pf_init) >= chains) {
        init <- pf_init[seq_len(chains)]
      }
    }
  } else {
    # "prior" — build_init at prior mode with jitter
    base <- build_init(stan_data, chains = 1)[[1]]
    init <- lapply(seq_len(chains), function(i)
      jitter_init(base, seed = 123 + i))
  }

  fit <- mod$sample(
    data = stan_data,
    seed = 123,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    init = init,
    refresh = 100,
    ...
  )
  # Persist CSV output files so the CmdStanMCMC object survives
  # serialization/deserialization by targets (qs format).
  csv_dir <- file.path("output", "stan_csv")
  if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)
  fit$save_output_files(dir = csv_dir)
  fit
}


#' Run generated quantities pass on a fitted model
#'
#' Compiles a separate GQ-only Stan model and runs it on the posterior draws
#' from a sampling fit.  This produces per-observation log-likelihoods and
#' posterior predictive replicates without inflating the MCMC output.
#'
#' @param fit A CmdStanMCMC object from \code{fit_model()}.
#' @param stan_data The same Stan data list used for sampling.
#' @param gq_stan_file Path to the GQ-only Stan file
#'   (default: \code{"stan/kinetics_model_gq.stan"}).
#' @return A CmdStanGQ object with draws for log_lik, *_hat, *_rep, etc.
run_gq <- function(fit, stan_data,
                   gq_stan_file = "stan/kinetics_model_gq.stan") {
  gq_mod <- cmdstanr::cmdstan_model(gq_stan_file)
  gq_fit <- gq_mod$generate_quantities(
    fitted_params = fit,
    data = stan_data,
    parallel_chains = fit$num_chains()
  )
  # Persist GQ CSV files
  gq_csv_dir <- file.path("output", "stan_csv_gq")
  if (!dir.exists(gq_csv_dir)) dir.create(gq_csv_dir, recursive = TRUE)
  gq_fit$save_output_files(dir = gq_csv_dir)
  gq_fit
}


#' Extract posterior predictions from a fitted model
#'
#' Uses parameter names from kinetics_model.stan.  Computes predicted RNA,
#' PFU, LFD and symptom-onset trajectories respecting all model flags.
#'
#' @param fit       CmdStanMCMC fit object
#' @param newdata   Data frame with columns: id, time, source, + covariate cols
#' @param stan_data Stan data list (from build_stan_data) — needed for flags
#' @param max_draws Maximum number of posterior draws to use (NULL = all).
#'                  Thinning to ~1000 draws is usually sufficient for posterior
#'                  means/quantiles and dramatically reduces memory and runtime.
#' @return Named list of rvar vectors (rna_hat, pfu_hat, lfd_hat, sym_hat,
#'         corner times, etc.)
predict_kinetics <- function(fit, newdata, stan_data, max_draws = 1000) {

  covariates <- define_covariates()

  # --- Build variable list matching Stan parameter names --------------------
  # Base variables (always present)
  variable_list <- c(
    "tp_i_pfu",
    "dp_mean_rna", "wp_mean_rna", "wr_mean_rna",
    # log-affine transformation: each is a vector[2] in Stan
    #   [1] = intercept,  [2] = elasticity / slope
    "tau_tp", "tau_dp", "tau_wp", "tau_wr",
    # LFD coefficients
    "tau0_lfd", "tau_lfd",
    # symptom onset hazard (cloglog)
    "eta_sym_intercept", "eta_sym_pfu", "eta_sym_rna",
    "sigma_sym", "z_sym"
  )

  # RNA individual effects: correlated (Cholesky NCP) or independent
  if (isTRUE(stan_data$ind_corr == 1)) {
    variable_list <- c(variable_list,
      "z_ind_rna", "L_Omega_rna", "sigma_ind_rna"
    )
  } else {
    variable_list <- c(variable_list, "tp_i_rna")
  }

  # Conditional parameters based on active flags
  if (stan_data$ind_effects) {
    if (!isTRUE(stan_data$ind_corr == 1)) {
      variable_list <- c(variable_list,
        "dp_i_rna", "wp_i_rna", "wr_i_rna"
      )
    }
    variable_list <- c(variable_list,
      "dp_i_pfu", "wp_i_pfu", "wr_i_pfu"
    )
  }
  if (stan_data$adj_rna) {
    variable_list <- c(variable_list,
      "beta_dp_rna", "beta_wp_rna", "beta_wr_rna"
    )
  }
  if (stan_data$adj_pfu) {
    variable_list <- c(variable_list,
      "beta_dp_pfu", "beta_wp_pfu", "beta_wr_pfu"
    )
  }
  if (stan_data$source_rna) {
    variable_list <- c(variable_list,
      "tp_k_rna", "dp_k_rna", "wp_k_rna", "wr_k_rna"
    )
  }
  if (stan_data$source_pfu) {
    variable_list <- c(variable_list,
      "tp_k_pfu", "dp_k_pfu", "wp_k_pfu", "wr_k_pfu"
    )
  }
  if (stan_data$source_lfd) {
    variable_list <- c(variable_list, "lfd_k")
  }
  if (stan_data$adj_lfd) {
    variable_list <- c(variable_list, "beta_lfd")
  }
  if (stan_data$source_sym) {
    variable_list <- c(variable_list, "to_k_sym")
  }
  if (stan_data$adj_sym) {
    variable_list <- c(variable_list, "beta_sym")
  }
  if (stan_data$use_wf) {
    variable_list <- c(variable_list, "wf_raw")
    if (stan_data$ind_effects) variable_list <- c(variable_list, "wf_i")
    if (stan_data$source_rna)  variable_list <- c(variable_list, "wf_k")
  }

  # --- Extract posterior draws ----------------------------------------------
  id  <- newdata[["id"]]
  t   <- newdata[["time"]]
  src <- as.integer(newdata[["source"]])

  lod <- dplyr::case_when(
    src == 1 ~ ct_to_rna(40, type = "nba")     + 0.01,
    src == 2 ~ ct_to_rna(40, type = "ata")     + 0.01,
    src == 3 ~ ct_to_rna(47, type = "uiuc-ct") + 0.01,
    src == 4 ~ log(1000) + 0.01,
    src == 5 ~ ct_to_rna(40, type = "nba")     + 0.01
  )

  lod_pfu <- dplyr::case_when(
    src == 1 ~ 0,
    src == 2 ~ 2.3,
    src == 3 ~ 2.3,
    src == 4 ~ log(5),
    src == 5 ~ 0
  )

  x <- as.matrix(newdata[, covariates$x_vars])

  # --- Extract and optionally thin draws -----------------------------------
  drws <- posterior::as_draws_matrix(fit$draws(variables = variable_list))
  n_total <- nrow(drws)
  if (!is.null(max_draws) && n_total > max_draws) {
    idx <- round(seq(1, n_total, length.out = max_draws))
    drws <- drws[idx, ]
    message(sprintf("predict_kinetics: thinned %d → %d draws", n_total, max_draws))
  }
  k <- posterior::as_draws_rvars(drws)
  rm(drws)  # free memory

  # --- RNA individual effects (correlated NCP → effective, or independent) --
  if (isTRUE(stan_data$ind_corr == 1)) {
    # NCP reconstruction: eta = diag(sigma) * L_Omega * z
    N_ind <- sum(stan_data$M)
    sigma_draws <- posterior::draws_of(k$sigma_ind_rna)  # [n_draws, 4]
    L_draws     <- posterior::draws_of(k$L_Omega_rna)    # [n_draws, 4, 4]
    z_draws     <- posterior::draws_of(k$z_ind_rna)      # [n_draws, 4, N_ind]
    n_draws     <- dim(sigma_draws)[1]

    eta <- array(NA_real_, dim = c(n_draws, 4, N_ind))
    for (d in seq_len(n_draws)) {
      L_Sigma <- diag(sigma_draws[d, ]) %*% L_draws[d, , ]
      eta[d, , ] <- L_Sigma %*% z_draws[d, , ]
    }
    eff_tp_i_rna <- posterior::rvar(eta[, 1, ])
    eff_dp_i_rna <- posterior::rvar(eta[, 2, ])
    eff_wp_i_rna <- posterior::rvar(eta[, 3, ])
    eff_wr_i_rna <- posterior::rvar(eta[, 4, ])
  } else {
    eff_tp_i_rna <- k$tp_i_rna
    if (stan_data$ind_effects) {
      eff_dp_i_rna <- k$dp_i_rna
      eff_wp_i_rna <- k$wp_i_rna
      eff_wr_i_rna <- k$wr_i_rna
    }
  }

  # --- RNA trajectory -------------------------------------------------------
  tp_rna <- eff_tp_i_rna[id]
  dp_rna <- k$dp_mean_rna
  wp_rna <- k$wp_mean_rna
  wr_rna <- k$wr_mean_rna

  if (stan_data$ind_effects) {
    dp_rna <- dp_rna * exp(eff_dp_i_rna[id])
    wp_rna <- wp_rna * exp(eff_wp_i_rna[id])
    wr_rna <- wr_rna * exp(eff_wr_i_rna[id])
  }

  if (stan_data$source_rna) {
    tp_rna <- tp_rna + k$tp_k_rna[src]
    dp_rna <- dp_rna * exp(k$dp_k_rna[src])
    wp_rna <- wp_rna * exp(k$wp_k_rna[src])
    wr_rna <- wr_rna * exp(k$wr_k_rna[src])
  }

  if (stan_data$adj_rna) {
    dp_rna <- dp_rna * exp(x %**% k$beta_dp_rna)
    wp_rna <- wp_rna * exp(x %**% k$beta_wp_rna)
    wr_rna <- wr_rna * exp(x %**% k$beta_wr_rna)
  }

  # --- Flat-top duration (shared RNA/PFU) -----------------------------------
  if (stan_data$use_wf) {
    wf <- stan_data$prior_wf_mean * exp(stan_data$prior_wf_cv * k$wf_raw[1])
    if (stan_data$ind_effects) wf <- wf * exp(k$wf_i[id])
    if (stan_data$source_rna)  wf <- wf * exp(k$wf_k[src])
  } else {
    wf <- posterior::rvar(0)
  }

  use_smooth_flag <- as.logical(stan_data$use_smooth)

  # RNA predictions
  rna_hat <- rvar_traj_fun(t, tp_rna, wp_rna, wr_rna, dp_rna, wf, use_smooth_flag)

  # --- PFU trajectory (log-affine from RNA) ---------------------------------
  # Stan: tau_dp = vector[2]; [1] = intercept, [2] = elasticity
  dp_pfu <- exp(k$tau_dp[1]) * dp_rna^k$tau_dp[2]
  wp_pfu <- exp(k$tau_wp[1]) * wp_rna^k$tau_wp[2]
  wr_pfu <- exp(k$tau_wr[1]) * wr_rna^k$tau_wr[2]

  # Stan: tp_pfu = tau_tp[1] + tau_tp[2] * tp_rna + tp_i_pfu[id]
  tp_pfu <- k$tau_tp[1] + k$tau_tp[2] * tp_rna + k$tp_i_pfu[id]

  if (stan_data$ind_effects) {
    dp_pfu <- dp_pfu * exp(k$dp_i_pfu[id])
    wp_pfu <- wp_pfu * exp(k$wp_i_pfu[id])
    wr_pfu <- wr_pfu * exp(k$wr_i_pfu[id])
  }

  if (stan_data$source_pfu) {
    tp_pfu <- tp_pfu + k$tp_k_pfu[src]
    dp_pfu <- dp_pfu * exp(k$dp_k_pfu[src])
    wp_pfu <- wp_pfu * exp(k$wp_k_pfu[src])
    wr_pfu <- wr_pfu * exp(k$wr_k_pfu[src])
  }

  if (stan_data$adj_pfu) {
    dp_pfu <- dp_pfu * exp(x %**% k$beta_dp_pfu)
    wp_pfu <- wp_pfu * exp(x %**% k$beta_wp_pfu)
    wr_pfu <- wr_pfu * exp(x %**% k$beta_wr_pfu)
  }

  pfu_hat <- rvar_traj_fun(t, tp_pfu, wp_pfu, wr_pfu, dp_pfu, wf, use_smooth_flag)

  # --- LFD probability -----------------------------------------------------
  # Stan: tau0_lfd (transformed param), tau_lfd (vector[2])
  lfd_hat <- expit(k$tau0_lfd +
                     k$tau_lfd[1] * rna_hat +
                     k$tau_lfd[2] * pfu_hat)

  if (stan_data$source_lfd) {
    lfd_hat <- expit(logit(lfd_hat) + k$lfd_k[src])
  }

  if (stan_data$adj_lfd) {
    lfd_hat <- expit(logit(lfd_hat) + x %**% k$beta_lfd)
  }

  # --- Symptom onset hazard (cloglog) ---------------------------------------
  # Normalise viral load by prior_dp_mean (matches Stan's scale_vl)
  scale_vl <- stan_data$prior_dp_mean
  u_sym   <- k$sigma_sym * k$z_sym[id]
  eta_lin <- k$eta_sym_intercept +
    k$eta_sym_pfu * (pfu_hat / scale_vl) +
    k$eta_sym_rna * (rna_hat / scale_vl) +
    u_sym

  if (stan_data$source_sym) {
    eta_lin <- eta_lin + k$to_k_sym[src]
  }

  if (stan_data$adj_sym) {
    eta_lin <- eta_lin + x %**% k$beta_sym
  }

  sym_hat <- 1 - exp(-exp(eta_lin))

  list(
    rna_hat = rna_hat,
    pfu_hat = pfu_hat,
    lfd_hat = lfd_hat,
    sym_hat = sym_hat
  )
}


#' Simulate from the prior predictive distribution
#'
#' Draws population parameters, individual effects, source/covariate effects,
#' and generates RNA/PFU/LFD/symptom predictions from the prior.
#'
#' @param data   Stan data list (from build_stan_data)
#' @param draws  Number of prior draws
#' @return Named list of matrices and vectors with prior predictive quantities
prior_predictive <- function(data, draws = 10) {

  id     <- data$id
  source <- data$source
  src    <- rep(1:data$K, data$M)

  # --- population parameter priors ---
  dp_raw <- rnorm(draws)
  wp_raw <- rnorm(draws)
  wr_raw <- rnorm(draws)

  dp_mean_rna <- data$prior_dp_mean * exp(data$prior_dp_cv * dp_raw)
  wp_mean_rna <- data$prior_wp_mean * exp(data$prior_wp_cv * wp_raw)
  wr_mean_rna <- data$prior_wr_mean * exp(data$prior_wr_cv * wr_raw)

  # error distributions
  sigma_rna <- truncnorm::rtruncnorm(draws, 0, Inf, 0, data$prior_sigma_sd)
  sigma_pfu <- truncnorm::rtruncnorm(draws, 0, Inf, 0, data$prior_sigma_sd)

  # --- transformation parameters (log-affine) ---
  tau0_tp <- rnorm(draws)
  tau_tp  <- rnorm(draws)
  tau0_dp <- rnorm(draws, -1, 1)
  tau0_wp <- rnorm(draws, -1, 1)
  tau0_wr <- rnorm(draws, -1, 1)
  tau_dp  <- truncnorm::rtruncnorm(draws, 0, Inf, 1, 0.5)
  tau_wp  <- truncnorm::rtruncnorm(draws, 0, Inf, 1, 0.5)
  tau_wr  <- truncnorm::rtruncnorm(draws, 0, Inf, 1, 0.5)

  # --- symptom onset hazard (cloglog) ---
  eta_sym_intercept <- rnorm(draws, -3, 1)
  eta_sym_pfu <- truncnorm::rtruncnorm(draws, 0, Inf, 0, 0.5)
  eta_sym_rna <- truncnorm::rtruncnorm(draws, 0, Inf, 0, 0.5)
  sigma_sym   <- truncnorm::rtruncnorm(draws, 0, Inf, 0, 1)
  z_sym <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws), diag(rep(1, draws)))

  # --- individual peak times (PFU only; RNA tp is handled above) ---
  tp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                diag(rep(data$prior_i_sd, draws)))
  if (!isTRUE(data$ind_corr == 1)) {
    tp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
  }

  # --- LFD coefficients ---
  tau0_lfd_raw <- rnorm(draws)
  tau0_lfd <- logit(data$prior_lfd_mean) + 1 * tau0_lfd_raw
  tau_lfd  <- mvtnorm::rmvnorm(draws, rep(0, 2))

  # --- test error rates ---
  if (data$test_error) {
    fp <- rbeta(draws, data$prior_fp * 50, (1 - data$prior_fp) * 50)
    fn <- rbeta(draws, data$prior_fn * 50, (1 - data$prior_fn) * 50)
  }

  # --- individual effects ---
  # When ind_corr = 1, draw RNA individual effects from MVN via Cholesky NCP.
  # For the prior predictive we approximate LKJ(2) by drawing from Wishart
  # and normalising to a correlation matrix, then sample MVN.
  if (isTRUE(data$ind_corr == 1)) {
    # sigma_ind_rna ~ normal(0, prior_i_sd) T[0,] (4 SDs per draw)
    sigma_ind_rna_pp <- matrix(
      truncnorm::rtruncnorm(4 * draws, 0, Inf, 0, data$prior_i_sd),
      nrow = 4, ncol = draws
    )
    # Approximate LKJ(2) correlation matrices via Wishart normalisation
    nu_lkj <- 4 + 2 * 2 - 2  # K + 2*eta - 2 = 6
    tp_i_rna <- matrix(NA, sum(data$M), draws)
    dp_i_rna <- matrix(NA, sum(data$M), draws)
    wp_i_rna <- matrix(NA, sum(data$M), draws)
    wr_i_rna <- matrix(NA, sum(data$M), draws)
    for (d in seq_len(draws)) {
      W <- stats::rWishart(1, nu_lkj, diag(4))[, , 1]
      D_inv <- diag(1 / sqrt(diag(W)))
      Omega <- D_inv %*% W %*% D_inv
      Sigma <- diag(sigma_ind_rna_pp[, d]) %*% Omega %*% diag(sigma_ind_rna_pp[, d])
      eff <- mvtnorm::rmvnorm(sum(data$M), rep(0, 4), Sigma)
      tp_i_rna[, d] <- eff[, 1]
      dp_i_rna[, d] <- eff[, 2]
      wp_i_rna[, d] <- eff[, 3]
      wr_i_rna[, d] <- eff[, 4]
    }
  }
  # When ind_corr = 0, tp_i_rna was already created in the peak-times section
  # above; dp/wp/wr are handled next.

  if (data$ind_effects && !isTRUE(data$ind_corr == 1)) {
    dp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
    wp_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
    wr_i_rna <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
  }

  if (data$ind_effects) {
    dp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
    wp_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
    wr_i_pfu <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                  diag(rep(data$prior_i_sd, draws)))
  }

  # --- source effects ---
  if (data$source_pfu) {
    tp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    dp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    wp_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    wr_k_pfu <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
  }

  if (data$source_rna) {
    tp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    dp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    wp_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
    wr_k_rna <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
  }

  if (data$source_sym) {
    to_k_sym <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                  diag(rep(data$prior_k_sd, draws)))
  }

  if (data$source_lfd) {
    lfd_k <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                               diag(rep(data$prior_k_sd, draws)))
  }

  # --- flat-top duration ---
  if (data$use_wf) {
    wf_raw <- rnorm(draws)
    wf_mean_pp <- data$prior_wf_mean * exp(data$prior_wf_cv * wf_raw)
    if (data$ind_effects) {
      wf_i <- mvtnorm::rmvnorm(sum(data$M), rep(0, draws),
                                 diag(rep(data$prior_i_sd, draws)))
    }
    if (data$source_rna) {
      wf_k <- mvtnorm::rmvnorm(data$K, rep(0, draws),
                                 diag(rep(data$prior_k_sd, draws)))
    }
  }

  # --- covariate effects ---
  if (data$adj_pfu) {
    beta_dp_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
    beta_wp_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
    beta_wr_pfu <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
  }

  if (data$adj_rna) {
    beta_dp_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
    beta_wp_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
    beta_wr_rna <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                     diag(rep(data$prior_beta_sd, data$P)))
  }

  if (data$adj_lfd) {
    beta_lfd <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                  diag(rep(data$prior_beta_sd, data$P)))
  }

  if (data$adj_sym) {
    beta_sym <- mvtnorm::rmvnorm(draws, rep(0, data$P),
                                  diag(rep(data$prior_beta_sd, data$P)))
  }

  # --- build RNA trajectories ---
  M_total <- sum(data$M)
  tp_rna <- tp_i_rna
  # Broadcast population means to M × draws matrices
  dp_rna <- rep_row(dp_mean_rna, M_total)
  wp_rna <- rep_row(wp_mean_rna, M_total)
  wr_rna <- rep_row(wr_mean_rna, M_total)

  if (data$ind_effects) {
    dp_rna <- dp_rna * exp(dp_i_rna)
    wp_rna <- wp_rna * exp(wp_i_rna)
    wr_rna <- wr_rna * exp(wr_i_rna)
  }

  if (data$source_rna) {
    tp_rna <- tp_rna + tp_k_rna[src, ]
    dp_rna <- dp_rna * exp(dp_k_rna[src, ])
    wp_rna <- wp_rna * exp(wp_k_rna[src, ])
    wr_rna <- wr_rna * exp(wr_k_rna[src, ])
  }

  if (data$adj_rna) {
    dp_rna <- dp_rna * exp(as.matrix(data$x) %*% t(beta_dp_rna))
    wp_rna <- wp_rna * exp(as.matrix(data$x) %*% t(beta_wp_rna))
    wr_rna <- wr_rna * exp(as.matrix(data$x) %*% t(beta_wr_rna))
  }

  # --- transform RNA → PFU (log-affine) ---
  tp_pfu <- rep_row(tau0_tp, sum(data$M)) +
    rep_row(tau_tp, sum(data$M)) * tp_rna + tp_i_pfu

  dp_mean_pfu <- exp(tau0_dp) * dp_mean_rna^tau_dp
  wp_mean_pfu <- exp(tau0_wp) * wp_mean_rna^tau_wp
  wr_mean_pfu <- exp(tau0_wr) * wr_mean_rna^tau_wr

  dp_pfu <- exp(rep_row(tau0_dp, sum(data$M))) *
    dp_rna^rep_row(tau_dp, sum(data$M))
  wp_pfu <- exp(rep_row(tau0_wp, sum(data$M))) *
    wp_rna^rep_row(tau_wp, sum(data$M))
  wr_pfu <- exp(rep_row(tau0_wr, sum(data$M))) *
    wr_rna^rep_row(tau_wr, sum(data$M))

  if (data$ind_effects) {
    dp_pfu <- dp_pfu * exp(dp_i_pfu)
    wp_pfu <- wp_pfu * exp(wp_i_pfu)
    wr_pfu <- wr_pfu * exp(wr_i_pfu)
  }

  if (data$source_pfu) {
    tp_pfu <- tp_pfu + tp_k_pfu[src, ]
    dp_pfu <- dp_pfu * exp(dp_k_pfu[src, ])
    wp_pfu <- wp_pfu * exp(wp_k_pfu[src, ])
    wr_pfu <- wr_pfu * exp(wr_k_pfu[src, ])
  }

  if (data$adj_pfu) {
    dp_pfu <- dp_pfu * exp(as.matrix(data$x) %*% t(beta_dp_pfu))
    wp_pfu <- wp_pfu * exp(as.matrix(data$x) %*% t(beta_wp_pfu))
    wr_pfu <- wr_pfu * exp(as.matrix(data$x) %*% t(beta_wr_pfu))
  }

  # --- build individual flat-top (N_ind x draws matrix) ---
  if (data$use_wf) {
    wf_pp <- rep_row(wf_mean_pp, sum(data$M))
    if (data$ind_effects) wf_pp <- wf_pp * exp(wf_i)
    if (data$source_rna)  wf_pp <- wf_pp * exp(wf_k[src, ])
  } else {
    wf_pp <- matrix(0, nrow = sum(data$M), ncol = draws)
  }

  # --- compute predicted values ---
  rna_hat <- matrix(NA, nrow = sum(data$N), ncol = draws)
  pfu_hat <- matrix(NA, nrow = sum(data$N), ncol = draws)
  lfd_hat <- matrix(NA, nrow = sum(data$N), ncol = draws)
  rna     <- matrix(NA, nrow = sum(data$N), ncol = draws)
  pfu     <- matrix(NA, nrow = sum(data$N), ncol = draws)
  lfd     <- matrix(NA, nrow = sum(data$N), ncol = draws)

  lod_rna_vec <- data$lod_rna[data$source]
  lod_pfu_vec <- data$lod_pfu[data$source]

  use_smooth_pp <- as.logical(data$use_smooth)

  for (d in seq_len(draws)) {
    rna_hat[, d] <- traj_fun(data$time, tp_rna[data$id, d], wp_rna[data$id, d],
                              wr_rna[data$id, d], dp_rna[data$id, d],
                              wf_pp[data$id, d], use_smooth = use_smooth_pp)
    pfu_hat[, d] <- traj_fun(data$time, tp_pfu[data$id, d], wp_pfu[data$id, d],
                              wr_pfu[data$id, d], dp_pfu[data$id, d],
                              wf_pp[data$id, d], use_smooth = use_smooth_pp)
    # Clamp non-finite values: NaN/Inf from extreme prior draws → 0 (below LOD)
    rna_hat[, d][!is.finite(rna_hat[, d])] <- 0
    pfu_hat[, d][!is.finite(pfu_hat[, d])] <- 0
    lfd_hat[, d] <- plogis(tau0_lfd[d] + tau_lfd[d, 1] * rna_hat[, d] +
                             tau_lfd[d, 2] * pfu_hat[, d])

    if (data$source_lfd) {
      lfd_hat[, d] <- plogis(qlogis(lfd_hat[, d]) + lfd_k[data$source, d])
    }

    if (data$adj_lfd) {
      lfd_hat[, d] <- plogis(qlogis(lfd_hat[, d]) +
                              as.matrix(data$x[data$id, ]) %*% beta_lfd[d, ])
    }

    rna[, d] <- rnorm(sum(data$N), rna_hat[, d], sigma_rna[d])
    pfu[, d] <- rnorm(sum(data$N), pfu_hat[, d], sigma_pfu[d])
    lfd[, d] <- rbinom(sum(data$N), 1, lfd_hat[, d])

    if (data$test_error) {
      fp_pfu <- rbinom(sum(data$N), 1, fp[d])
      fn_pfu <- rbinom(sum(data$N), 1, fn[d])
      fp_rna <- rbinom(sum(data$N), 1, fp[d])
      fn_rna <- rbinom(sum(data$N), 1, fn[d])
      pfu_error <- rexp(sum(data$N), data$fp_mean[data$source])
      rna_error <- rexp(sum(data$N), data$fp_mean[data$source])

      pfu[, d] <- (1 - fp_pfu) * (1 - fn_pfu) * pfu[, d] +
        fp_pfu * (pfu_error + lod_pfu_vec) + fn_pfu * lod_pfu_vec
      rna[, d] <- (1 - fp_rna) * (1 - fn_rna) * rna[, d] +
        fp_rna * (rna_error + lod_rna_vec) + fn_rna * lod_rna_vec
    }

    idx_rna <- which(rna[, d] < lod_rna_vec)
    if (length(idx_rna) > 0) rna[idx_rna, d] <- lod_rna_vec[idx_rna]
    idx_pfu <- which(pfu[, d] < lod_pfu_vec)
    if (length(idx_pfu) > 0) pfu[idx_pfu, d] <- lod_pfu_vec[idx_pfu]
  }

  # --- symptom onset (cloglog) ---
  sym_hat <- matrix(NA, nrow = sum(data$N), ncol = draws)
  sym     <- matrix(NA, nrow = sum(data$N), ncol = draws)
  u_sym   <- rep_row(sigma_sym, sum(data$M)) * z_sym

  # Normalise viral load by prior_dp_mean (matches Stan's scale_vl)
  scale_vl <- data$prior_dp_mean

  for (d in seq_len(draws)) {
    eta_lin <- eta_sym_intercept[d] +
      eta_sym_pfu[d] * (pfu_hat[, d] / scale_vl) +
      eta_sym_rna[d] * (rna_hat[, d] / scale_vl) +
      u_sym[data$id, d]

    if (data$source_sym) {
      eta_lin <- eta_lin + to_k_sym[data$source, d]
    }

    if (data$adj_sym) {
      eta_lin <- eta_lin +
        as.matrix(data$x[data$id, ]) %*% beta_sym[d, ]
    }

    sym_hat[, d] <- 1 - exp(-exp(eta_lin))
    sym[, d] <- rbinom(sum(data$N), 1, sym_hat[, d])
  }

  list(
    id = data$id, time = data$time,
    dp_mean_rna = dp_mean_rna, wp_mean_rna = wp_mean_rna,
    wr_mean_rna = wr_mean_rna,
    dp_mean_pfu = dp_mean_pfu, wp_mean_pfu = wp_mean_pfu,
    wr_mean_pfu = wr_mean_pfu,
    lod_rna = data$lod_rna[src], lod_pfu = data$lod_pfu[src],
    tau0_tp = tau0_tp, tau0_dp = tau0_dp,
    tau0_wp = tau0_wp, tau0_wr = tau0_wr,
    eta_sym_intercept = eta_sym_intercept,
    eta_sym_pfu = eta_sym_pfu, eta_sym_rna = eta_sym_rna,
    sigma_sym = sigma_sym,
    tau_tp = tau_tp, tau_dp = tau_dp, tau_wp = tau_wp, tau_wr = tau_wr,
    dp_rna = dp_rna, wp_rna = wp_rna, wr_rna = wr_rna, tp_rna = tp_rna,
    dp_pfu = dp_pfu, wp_pfu = wp_pfu, wr_pfu = wr_pfu, tp_pfu = tp_pfu,
    rna_hat = rna_hat, pfu_hat = pfu_hat, lfd_hat = lfd_hat,
    sym_hat = sym_hat, rna = rna, pfu = pfu, lfd = lfd, sym = sym,
    wf_mean_pp = if (data$use_wf) wf_mean_pp else NULL,
    use_wf = data$use_wf, use_smooth = data$use_smooth
  )
}


#' Build a default parameter list for simulate_data()
#'
#' Returns plausible population-level parameter values that can be passed
#' directly to \code{\link{simulate_data}}.  Edit any element to explore
#' sensitivity / identifiability.
#'
#' @param data  Stan data list (from build_stan_data) — used for dimensions
#' @return Named list of parameter values matching the Stan model
default_params <- function(data) {
  N_ind <- sum(data$M)  # total individuals
  K     <- data$K       # number of sources
  P     <- data$P       # number of covariates

  list(
    # --- RNA kinetics (population) ---
    # Values informed by posterior medians from production fit (Feb 2026)
    dp_mean_rna = 16.5,   # peak log-RNA (posterior ≈ 16.6)
    wp_mean_rna =  5.3,   # proliferation width in days (posterior ≈ 5.3)
    wr_mean_rna = 15.0,   # clearance width in days (posterior ≈ 15.1)

    # --- observation noise ---
    # Posterior: sigma_rna ≈ 2.2, sigma_pfu ≈ 2.2
    sigma_rna = 2.2,
    sigma_pfu = 2.2,

    # --- log-affine RNA → PFU transformation ---
    # tau_*[1] = log-scale intercept, tau_*[2] = elasticity
    # Posterior medians: tau_tp ≈ (-0.15, 1.24), tau_dp ≈ (-0.93, 1.11),
    #   tau_wp ≈ (-1.26, 1.19), tau_wr ≈ (-0.66, 0.84)
    tau0_tp = -0.15,  tau_tp = 1.2,
    tau0_dp = -0.9,   tau_dp = 1.1,
    tau0_wp = -1.3,   tau_wp = 1.2,
    tau0_wr = -0.7,   tau_wr = 0.8,

    # --- LFD model ---
    # Posterior: tau_lfd ≈ (0.35, 0.29)
    tau0_lfd = logit(data$prior_lfd_mean),
    tau_lfd  = c(0.35, 0.29),

    # --- symptom onset hazard (cloglog) ---
    # Posterior: intercept ≈ -1.76, pfu ≈ 0.40, rna ≈ 1.55, sigma ≈ 0.77
    eta_sym_intercept = -1.8,
    eta_sym_pfu       =  0.4,
    eta_sym_rna       =  1.5,
    sigma_sym         =  0.8,

    # --- test error ---
    # Posterior: fp ≈ 0.016, fn ≈ ~0
    fp = 0.015,
    fn = 0.001,

    # --- individual-level SDs (for generating individual effects) ---
    # prior_i_sd = 1 allows exp(±2σ) ≈ [0.14, 7.4]× pop mean.
    # Use 0.5 for simulation → exp(±2σ) ≈ [0.37, 2.7]× — more realistic.
    sd_tp_i    = 0.5,   # peak time RE SD
    sd_ind     = 0.5,   # dp/wp/wr individual RE SD
    sd_source  = data$prior_k_sd,   # source RE SD
    sd_beta    = 0.0,               # covariate effect magnitude (0 = no effect)

    # --- flat-top duration (shared RNA/PFU) ---
    wf_mean = 1.0,   # flat-top duration in days (only used when use_wf = 1)

    # --- correlated RNA individual effects (only used when ind_corr = 1) ---
    sigma_ind_rna = rep(0.5, 4),  # SD of (tp, dp, wp, wr) RNA effects
    Omega_rna = diag(4)           # 4×4 correlation matrix (identity = no corr)
  )
}


#' Simulate data under known/fixed parameter values
#'
#' Uses the same generative model as the Stan program but with parameters
#' supplied as a named list rather than drawn from priors.
#' Useful for:
#' \itemize{
#'   \item Simulation-based calibration (SBC)
#'   \item Parameter recovery / identifiability checks
#'   \item Power analysis
#'   \item Debugging the Stan model
#' }
#'
#' @param data    Stan data list (from build_stan_data) — provides structure,
#'                flags, covariates, and observation times
#' @param params  Named list of fixed parameter values (see
#'                \code{\link{default_params}} for the expected names).
#'                If NULL, uses \code{default_params(data)}.
#' @param seed    Random seed (for reproducibility of individual effects and
#'                observation noise)
#' @return A list with two elements:
#'   \describe{
#'     \item{sim_data}{A Stan-ready data list identical in structure to
#'                     \code{data} but with \code{rna}, \code{pfu}, \code{lfd},
#'                     \code{sym} replaced by simulated observations.}
#'     \item{truth}{A named list of all true parameter values (population,
#'                  individual, source, covariate effects) so recovery can be
#'                  checked.}
#'   }
simulate_data <- function(data, params = NULL, seed = 42) {

  set.seed(seed)
  if (is.null(params)) params <- default_params(data)

  N_obs  <- sum(data$N)
  N_ind  <- sum(data$M)
  K      <- data$K
  P      <- data$P
  src    <- rep(1:K, data$M)       # source index per individual
  id     <- data$id                # individual index per observation
  source <- data$source            # source index per observation

  # ─── Population parameters (fixed) ────────────────────────────────────────
  dp_mean_rna <- params$dp_mean_rna
  wp_mean_rna <- params$wp_mean_rna
  wr_mean_rna <- params$wr_mean_rna
  sigma_rna   <- params$sigma_rna
  sigma_pfu   <- params$sigma_pfu

  # ─── Individual peak times ────────────────────────────────────────────────
  tp_i_rna <- rnorm(N_ind, 0, params$sd_tp_i)
  tp_i_pfu <- rnorm(N_ind, 0, params$sd_tp_i)

  # ─── Symptom random effects ──────────────────────────────────────────────
  z_sym <- rnorm(N_ind)

  # ─── Individual effects (dp, wp, wr) ─────────────────────────────────────
  if (isTRUE(data$ind_corr == 1)) {
    # Correlated RNA individual effects: draw from MVN(0, Sigma)
    Sigma_rna <- diag(params$sigma_ind_rna) %*% params$Omega_rna %*%
                 diag(params$sigma_ind_rna)
    eff <- mvtnorm::rmvnorm(N_ind, rep(0, 4), Sigma_rna)
    tp_i_rna <- eff[, 1]
    dp_i_rna <- eff[, 2]
    wp_i_rna <- eff[, 3]
    wr_i_rna <- eff[, 4]
    # PFU individual effects remain independent
    if (data$ind_effects) {
      dp_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
      wp_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
      wr_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
    } else {
      dp_i_pfu <- wp_i_pfu <- wr_i_pfu <- rep(0, N_ind)
    }
  } else if (data$ind_effects) {
    dp_i_rna <- rnorm(N_ind, 0, params$sd_ind)
    wp_i_rna <- rnorm(N_ind, 0, params$sd_ind)
    wr_i_rna <- rnorm(N_ind, 0, params$sd_ind)
    dp_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
    wp_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
    wr_i_pfu <- rnorm(N_ind, 0, params$sd_ind)
  } else {
    dp_i_rna <- wp_i_rna <- wr_i_rna <- rep(0, N_ind)
    dp_i_pfu <- wp_i_pfu <- wr_i_pfu <- rep(0, N_ind)
  }

  # ─── Source effects ──────────────────────────────────────────────────────
  if (data$source_rna) {
    tp_k_rna <- rnorm(K, 0, params$sd_source)
    dp_k_rna <- rnorm(K, 0, params$sd_source)
    wp_k_rna <- rnorm(K, 0, params$sd_source)
    wr_k_rna <- rnorm(K, 0, params$sd_source)
  } else {
    tp_k_rna <- dp_k_rna <- wp_k_rna <- wr_k_rna <- rep(0, K)
  }

  if (data$source_pfu) {
    tp_k_pfu <- rnorm(K, 0, params$sd_source)
    dp_k_pfu <- rnorm(K, 0, params$sd_source)
    wp_k_pfu <- rnorm(K, 0, params$sd_source)
    wr_k_pfu <- rnorm(K, 0, params$sd_source)
  } else {
    tp_k_pfu <- dp_k_pfu <- wp_k_pfu <- wr_k_pfu <- rep(0, K)
  }

  if (data$source_lfd) {
    lfd_k <- rnorm(K, 0, params$sd_source)
  } else {
    lfd_k <- rep(0, K)
  }

  if (data$source_sym) {
    to_k_sym <- rnorm(K, 0, params$sd_source)
  } else {
    to_k_sym <- rep(0, K)
  }

  # ─── Flat-top duration ──────────────────────────────────────────────────
  wf_mean <- if (data$use_wf) params$wf_mean else 0
  if (data$use_wf) {
    if (data$ind_effects) {
      wf_i <- rnorm(N_ind, 0, params$sd_ind)
    } else {
      wf_i <- rep(0, N_ind)
    }
    if (data$source_rna) {
      wf_k <- rnorm(K, 0, params$sd_source)
    } else {
      wf_k <- rep(0, K)
    }
  } else {
    wf_i <- rep(0, N_ind)
    wf_k <- rep(0, K)
  }

  # ─── Covariate effects ──────────────────────────────────────────────────
  if (data$adj_rna && P > 0) {
    beta_dp_rna <- rnorm(P, 0, params$sd_beta)
    beta_wp_rna <- rnorm(P, 0, params$sd_beta)
    beta_wr_rna <- rnorm(P, 0, params$sd_beta)
  } else {
    beta_dp_rna <- beta_wp_rna <- beta_wr_rna <- rep(0, max(P, 1))
  }

  if (data$adj_pfu && P > 0) {
    beta_dp_pfu <- rnorm(P, 0, params$sd_beta)
    beta_wp_pfu <- rnorm(P, 0, params$sd_beta)
    beta_wr_pfu <- rnorm(P, 0, params$sd_beta)
  } else {
    beta_dp_pfu <- beta_wp_pfu <- beta_wr_pfu <- rep(0, max(P, 1))
  }

  if (data$adj_lfd && P > 0) {
    beta_lfd <- rnorm(P, 0, params$sd_beta)
  } else {
    beta_lfd <- rep(0, max(P, 1))
  }

  if (data$adj_sym && P > 0) {
    beta_sym <- rnorm(P, 0, params$sd_beta)
  } else {
    beta_sym <- rep(0, max(P, 1))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Build individual-level RNA kinetics
  # ═══════════════════════════════════════════════════════════════════════════
  tp_rna_ind <- tp_i_rna
  dp_rna_ind <- rep(dp_mean_rna, N_ind)
  wp_rna_ind <- rep(wp_mean_rna, N_ind)
  wr_rna_ind <- rep(wr_mean_rna, N_ind)

  if (data$ind_effects) {
    dp_rna_ind <- dp_rna_ind * exp(dp_i_rna)
    wp_rna_ind <- wp_rna_ind * exp(wp_i_rna)
    wr_rna_ind <- wr_rna_ind * exp(wr_i_rna)
  }

  if (data$source_rna) {
    tp_rna_ind <- tp_rna_ind + tp_k_rna[src]
    dp_rna_ind <- dp_rna_ind * exp(dp_k_rna[src])
    wp_rna_ind <- wp_rna_ind * exp(wp_k_rna[src])
    wr_rna_ind <- wr_rna_ind * exp(wr_k_rna[src])
  }

  if (data$adj_rna && P > 0) {
    x_mat <- as.matrix(data$x)  # N_ind x P
    dp_rna_ind <- dp_rna_ind * exp(as.numeric(x_mat %*% beta_dp_rna))
    wp_rna_ind <- wp_rna_ind * exp(as.numeric(x_mat %*% beta_wp_rna))
    wr_rna_ind <- wr_rna_ind * exp(as.numeric(x_mat %*% beta_wr_rna))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Build individual-level PFU kinetics (log-affine from RNA)
  # ═══════════════════════════════════════════════════════════════════════════
  tp_pfu_ind <- params$tau0_tp + params$tau_tp * tp_rna_ind + tp_i_pfu
  dp_pfu_ind <- exp(params$tau0_dp) * dp_rna_ind^params$tau_dp
  wp_pfu_ind <- exp(params$tau0_wp) * wp_rna_ind^params$tau_wp
  wr_pfu_ind <- exp(params$tau0_wr) * wr_rna_ind^params$tau_wr

  if (data$ind_effects) {
    dp_pfu_ind <- dp_pfu_ind * exp(dp_i_pfu)
    wp_pfu_ind <- wp_pfu_ind * exp(wp_i_pfu)
    wr_pfu_ind <- wr_pfu_ind * exp(wr_i_pfu)
  }

  if (data$source_pfu) {
    tp_pfu_ind <- tp_pfu_ind + tp_k_pfu[src]
    dp_pfu_ind <- dp_pfu_ind * exp(dp_k_pfu[src])
    wp_pfu_ind <- wp_pfu_ind * exp(wp_k_pfu[src])
    wr_pfu_ind <- wr_pfu_ind * exp(wr_k_pfu[src])
  }

  if (data$adj_pfu && P > 0) {
    if (!exists("x_mat")) x_mat <- as.matrix(data$x)
    dp_pfu_ind <- dp_pfu_ind * exp(as.numeric(x_mat %*% beta_dp_pfu))
    wp_pfu_ind <- wp_pfu_ind * exp(as.numeric(x_mat %*% beta_wp_pfu))
    wr_pfu_ind <- wr_pfu_ind * exp(as.numeric(x_mat %*% beta_wr_pfu))
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Build individual-level flat-top duration (shared RNA/PFU)
  # ═══════════════════════════════════════════════════════════════════════════
  wf_ind <- rep(wf_mean, N_ind)
  if (data$use_wf) {
    if (data$ind_effects) wf_ind <- wf_ind * exp(wf_i)
    if (data$source_rna)  wf_ind <- wf_ind * exp(wf_k[src])
  }

  # ═══════════════════════════════════════════════════════════════════════════
  # Generate observations (per time-point)
  # ═══════════════════════════════════════════════════════════════════════════
  use_smooth <- as.logical(data$use_smooth)
  lod_rna_vec <- data$lod_rna[source]
  lod_pfu_vec <- data$lod_pfu[source]

  rna_hat <- traj_fun(data$time, tp_rna_ind[id], wp_rna_ind[id],
                       wr_rna_ind[id], dp_rna_ind[id], wf_ind[id],
                       use_smooth = use_smooth)
  pfu_hat <- traj_fun(data$time, tp_pfu_ind[id], wp_pfu_ind[id],
                       wr_pfu_ind[id], dp_pfu_ind[id], wf_ind[id],
                       use_smooth = use_smooth)

  # Apply safe_vl clamping to match Stan model (stan/kinetics_model.stan L69-71)
  # Stan clamps predicted viral loads to [-50, 50] before evaluating likelihood
  rna_hat <- pmax(pmin(rna_hat, 50), -50)
  pfu_hat <- pmax(pmin(pfu_hat, 50), -50)

  # LFD probability
  lfd_logit <- params$tau0_lfd +
    params$tau_lfd[1] * rna_hat +
    params$tau_lfd[2] * pfu_hat

  if (data$source_lfd) {
    lfd_logit <- lfd_logit + lfd_k[source]
  }
  if (data$adj_lfd && P > 0) {
    if (!exists("x_mat")) x_mat <- as.matrix(data$x)
    lfd_logit <- lfd_logit + as.numeric(x_mat[id, , drop = FALSE] %*% beta_lfd)
  }
  lfd_hat <- plogis(lfd_logit)

  # Symptom hazard (cloglog) — normalise by prior_dp_mean (matches Stan)
  scale_vl <- data$prior_dp_mean
  u_sym <- params$sigma_sym * z_sym[id]
  eta_lin <- params$eta_sym_intercept +
    params$eta_sym_pfu * (pfu_hat / scale_vl) +
    params$eta_sym_rna * (rna_hat / scale_vl) +
    u_sym

  if (data$source_sym) {
    eta_lin <- eta_lin + to_k_sym[source]
  }
  if (data$adj_sym && P > 0) {
    if (!exists("x_mat")) x_mat <- as.matrix(data$x)
    eta_lin <- eta_lin + as.numeric(x_mat[id, , drop = FALSE] %*% beta_sym)
  }
  sym_hat <- 1 - exp(-exp(eta_lin))

  # ─── Draw noisy observations ─────────────────────────────────────────────
  rna_obs <- rnorm(N_obs, rna_hat, sigma_rna)
  pfu_obs <- rnorm(N_obs, pfu_hat, sigma_pfu)
  lfd_obs <- rbinom(N_obs, 1, lfd_hat)
  sym_obs <- rbinom(N_obs, 1, sym_hat)

  # Test error contamination
  # The Stan model treats fp/fn/true as mutually exclusive mixture components:
  #   P(obs | above LOD) = fp * Exp(fp_mean) + (1-fp)(1-fn) * Normal(rna_hat, sigma)
  #   P(obs | below LOD) = (1-fn) * Phi_normal(lod | rna_hat, sigma)
  # So we sample from the same three-component mixture:
  #   component 0 = true observation (weight (1-fp)*(1-fn))
  #   component 1 = false positive   (weight fp)
  #   component 2 = false negative   (weight fn*(1-fp))
  if (data$test_error) {
    fp_val <- params$fp
    fn_val <- params$fn
    mix_probs <- c((1 - fp_val) * (1 - fn_val), fp_val, fn_val * (1 - fp_val))
    component_rna <- sample(0:2, N_obs, replace = TRUE, prob = mix_probs)
    component_pfu <- sample(0:2, N_obs, replace = TRUE, prob = mix_probs)

    rna_error <- rexp(N_obs, data$fp_mean[source])
    pfu_error <- rexp(N_obs, data$fp_mean[source])

    rna_obs[component_rna == 1] <- rna_error[component_rna == 1] +
      lod_rna_vec[component_rna == 1]
    rna_obs[component_rna == 2] <- lod_rna_vec[component_rna == 2]

    pfu_obs[component_pfu == 1] <- pfu_error[component_pfu == 1] +
      lod_pfu_vec[component_pfu == 1]
    pfu_obs[component_pfu == 2] <- lod_pfu_vec[component_pfu == 2]
  }

  # Censor below LOD
  rna_obs <- pmax(rna_obs, lod_rna_vec)
  pfu_obs <- pmax(pfu_obs, lod_pfu_vec)

  # ═══════════════════════════════════════════════════════════════════════════
  # Package output
  # ═══════════════════════════════════════════════════════════════════════════

  # Build a Stan-ready data list by replacing observed outcomes
  sim_data       <- data
  sim_data$rna   <- rna_obs
  sim_data$pfu   <- pfu_obs
  sim_data$lfd   <- as.integer(lfd_obs)
  sim_data$sym   <- as.integer(sym_obs)

  # Ground truth for parameter recovery checks
  truth <- list(
    # population
    dp_mean_rna = dp_mean_rna,
    wp_mean_rna = wp_mean_rna,
    wr_mean_rna = wr_mean_rna,
    sigma_rna   = sigma_rna,
    sigma_pfu   = sigma_pfu,
    # log-affine transformation
    tau0_tp = params$tau0_tp, tau_tp = params$tau_tp,
    tau0_dp = params$tau0_dp, tau_dp = params$tau_dp,
    tau0_wp = params$tau0_wp, tau_wp = params$tau_wp,
    tau0_wr = params$tau0_wr, tau_wr = params$tau_wr,
    # LFD
    tau0_lfd = params$tau0_lfd, tau_lfd = params$tau_lfd,
    # symptom onset
    eta_sym_intercept = params$eta_sym_intercept,
    eta_sym_pfu = params$eta_sym_pfu,
    eta_sym_rna = params$eta_sym_rna,
    sigma_sym   = params$sigma_sym,
    # test error
    fp = params$fp, fn = params$fn,
    # individual effects
    tp_i_rna = tp_i_rna, tp_i_pfu = tp_i_pfu,
    dp_i_rna = dp_i_rna, wp_i_rna = wp_i_rna, wr_i_rna = wr_i_rna,
    dp_i_pfu = dp_i_pfu, wp_i_pfu = wp_i_pfu, wr_i_pfu = wr_i_pfu,
    z_sym    = z_sym,
    # source effects
    tp_k_rna = tp_k_rna, dp_k_rna = dp_k_rna,
    wp_k_rna = wp_k_rna, wr_k_rna = wr_k_rna,
    tp_k_pfu = tp_k_pfu, dp_k_pfu = dp_k_pfu,
    wp_k_pfu = wp_k_pfu, wr_k_pfu = wr_k_pfu,
    lfd_k    = lfd_k,    to_k_sym = to_k_sym,
    wf_i     = wf_i,    wf_k     = wf_k,
    # covariate effects
    beta_dp_rna = beta_dp_rna, beta_wp_rna = beta_wp_rna,
    beta_wr_rna = beta_wr_rna,
    beta_dp_pfu = beta_dp_pfu, beta_wp_pfu = beta_wp_pfu,
    beta_wr_pfu = beta_wr_pfu,
    beta_lfd    = beta_lfd,    beta_sym    = beta_sym,
    # correlated RNA individual effects (only when ind_corr = 1)
    sigma_ind_rna = if (isTRUE(data$ind_corr == 1)) params$sigma_ind_rna else NULL,
    Omega_rna     = if (isTRUE(data$ind_corr == 1)) params$Omega_rna else NULL,
    # derived individual kinetics
    tp_rna = tp_rna_ind, dp_rna = dp_rna_ind,
    wp_rna = wp_rna_ind, wr_rna = wr_rna_ind,
    tp_pfu = tp_pfu_ind, dp_pfu = dp_pfu_ind,
    wp_pfu = wp_pfu_ind, wr_pfu = wr_pfu_ind,
    wf = wf_ind,  wf_mean = wf_mean,
    # latent means
    rna_hat = rna_hat, pfu_hat = pfu_hat,
    lfd_hat = lfd_hat, sym_hat = sym_hat
  )

  list(sim_data = sim_data, truth = truth)
}
