# ──────────────────────────────────────────────────────────────────────────────
# simulation.R — Simulation study functions (piecewise vs target-cell)
# ──────────────────────────────────────────────────────────────────────────────

#' Target cell limited ODE system (for simulation study)
#'
#' @param t       Time
#' @param states  Named vector: T, I, V
#' @param params  Named vector: beta, delta, pi_param, c_param
#' @return List of derivatives
targetcell <- function(t, states, params) {
  T <- states[1]; I <- states[2]; V <- states[3]
  with(as.list(params), {
    dT <- -beta * T * V
    dI <- beta * T * V - delta * I
    dV <- pi_param * I - c_param * V
    list(c(dT, dI, dV))
  })
}

#' Add observation noise with LOD censoring
#'
#' @param y     True log-viral load
#' @param sigma Observation noise SD
#' @param lod   Limit of detection
#' @return Observed value (censored at LOD)
observe <- function(y, sigma = 1, lod = 2.3) {
  obs <- y + rnorm(length(y), 0, sigma)
  ifelse(obs < lod, lod, obs)
}

#' Generate simulated data from the target cell model
#'
#' @param model    ODE function
#' @param params   ODE parameters
#' @param states   Initial states
#' @param times    Observation times
#' @param samples  Number of simulated individuals
#' @param sigma    Observation noise SD
#' @param lod      Limit of detection
#' @return List: df (data frame), truth (ODE solution)
gendata <- function(model, params, states, times, samples, sigma, lod) {
  out <- deSolve::ode(y = states, times = times, func = model, parms = params)
  truth <- data.frame(out)
  truth$logV <- log(truth$V + 1)

  df <- purrr::map_dfr(seq_len(samples), function(i) {
    tibble::tibble(
      id   = i,
      time = times,
      y    = observe(truth$logV, sigma = sigma, lod = lod)
    )
  })

  list(df = df, truth = truth)
}

#' Calculate onset/peak/resolution from piecewise parameters (sim version)
#'
#' @param tp  Peak time
#' @param wp  Proliferation width
#' @param wr  Clearance width
#' @param dp  Peak height (log-scale)
#' @param lod Limit of detection
#' @return data.frame with to, tp, tr
sim_calc_corners <- function(tp, wp, wr, dp, lod) {
  to <- lod * wp / dp + (tp - wp)
  tr <- (dp - lod) * wr / dp + tp
  data.frame(to = to, tp = tp, tr = tr)
}

#' Piecewise-linear kinetics with optional flat top
piecewise_sim <- function(t, tp, wp, wr, dp, wf = 0) {
  ifelse(
    t <= tp - wf / 2,
    dp / wp * (t - (tp - wp)),
    ifelse(
      t <= tp + wf / 2,
      dp,
      dp - dp / wr * (t - tp - wf / 2)
    )
  )
}

#' Smooth-mixed approximation
smixed <- function(t, tp, wp, wr, dp, wf = 0) {
  a <- dp / wp
  b <- dp / wr
  k <- 5
  a * (t - (tp - wp)) * plogis(k * (tp - t)) +
    (dp - b * (t - tp)) * plogis(k * (t - tp))
}

#' Compute ground truth from ODE with fixed parameters
#'
#' @return data.frame with onset, peak, resolution from ODE trajectory
get_truth <- function() {
  params <- c(beta = 2.4e-5, delta = 2.0, pi_param = 1700, c_param = 10)
  states <- c(T = 4e8, I = 0, V = 0.4)
  times  <- seq(0, 20, by = 0.01)

  out <- deSolve::ode(y = states, times = times, func = targetcell,
                       parms = params)
  df <- data.frame(out)
  df$logV <- log(df$V + 1)
  df
}


#' Build simulation parameter grid
#'
#' @return tibble of simulation scenarios
build_sim_grid <- function() {
  tibble::tibble(
    sigma = c(0.5, 1.0, 2.0),
    beta  = 2.4e-5,
    delta = 2.0,
    pi_param = 1700,
    c_param  = 10
  )
}


#' Run simulation study across parameter grid
#'
#' @param stan_file Path to sim.stan
#' @param param_grid  Tibble of scenarios (from build_sim_grid)
#' @param n_sims      Number of simulated datasets per scenario
#' @return List of results (one per scenario)
run_simulation <- function(stan_file, param_grid, n_sims = 50) {
  mod <- cmdstanr::cmdstan_model(stan_file)

  states <- c(T = 4e8, I = 0, V = 0.4)
  times  <- seq(0, 20, by = 0.5)
  lod    <- 2.3
  n_times <- length(times)

  results <- purrr::pmap(param_grid, function(sigma, beta, delta,
                                               pi_param, c_param) {
    params <- c(beta = beta, delta = delta,
                pi_param = pi_param, c_param = c_param)
    simdata <- gendata(targetcell, params, states, times,
                        samples = n_sims, sigma = sigma, lod = lod)

    # Build the data structure that sim.stan expects:
    #   N individuals, T time points, y as N*T vector (row-major by individual)
    y_vec <- numeric(n_sims * n_times)
    start_vec <- integer(n_sims)
    for (i in seq_len(n_sims)) {
      d <- dplyr::filter(simdata$df, id == i)
      idx <- ((i - 1) * n_times + 1):(i * n_times)
      y_vec[idx] <- d$y
      # start = index of first observation above LOD (0-based for Stan)
      first_above <- which(d$y > lod)[1]
      start_vec[i] <- if (is.na(first_above)) 0L else max(0L, first_above - 2L)
    }

    stan_dat <- list(
      N = n_sims,
      T = n_times,
      id = seq_len(n_sims),
      start = start_vec,
      time = times,
      y = y_vec,
      y0 = as.array(states),
      lod = lod,
      breakpoint = 0L
    )

    fit <- tryCatch(
      mod$optimize(data = stan_dat, seed = 42, algorithm = "lbfgs",
                    iter = 2000, init = 0.1),
      error = function(e) {
        warning("sim optimization failed: ", conditionMessage(e))
        NULL
      }
    )

    # Package as a list of per-individual "fits" for process_sim_results
    fits <- list()
    if (!is.null(fit)) {
      for (i in seq_len(n_sims)) {
        tryCatch({
          mle_all <- fit$mle()
          fits[[i]] <- list(
            mle = function() {
              c(tp = mle_all[[paste0("tp_pe[", i, "]")]],
                dp = mle_all[[paste0("dp_pe[", i, "]")]],
                wp = mle_all[[paste0("wp_pe[", i, "]")]],
                wr = mle_all[[paste0("wr_pe[", i, "]")]])
            }
          )
        }, error = function(e) {
          fits[[i]] <<- NULL
        })
      }
    }

    list(df = simdata$df, truth = simdata$truth, fits = fits)
  })

  results
}


#' Process simulation results into summary statistics
#'
#' @param results  Output of run_simulation()
#' @return tibble with bias estimates for each scenario and model type
process_sim_results <- function(results) {
  truth <- get_truth()
  truth_corners <- sim_calc_corners(
    tp = truth$time[which.max(truth$logV)],
    wp = truth$time[which.max(truth$logV)] -
      truth$time[min(which(truth$logV > 2.3))],
    wr = truth$time[max(which(truth$logV > 2.3))] -
      truth$time[which.max(truth$logV)],
    dp = max(truth$logV),
    lod = 2.3
  )

  out <- purrr::map_dfr(seq_along(results), function(scenario) {
    res <- results[[scenario]]
    purrr::map_dfr(seq_along(res$fits), function(i) {
      fit <- res$fits[[i]]
      if (is.null(fit)) return(NULL)

      tryCatch({
        mle <- fit$mle()
        corners <- sim_calc_corners(mle[["tp"]], mle[["wp"]], mle[["wr"]],
                                     mle[["dp"]], lod = 2.3)
        tibble::tibble(
          scenario = scenario,
          sim      = i,
          tp_bias  = corners$tp - truth_corners$tp,
          dp_bias  = mle[["dp"]] - max(truth$logV),
          wp_bias  = mle[["wp"]] - (truth_corners$tp - truth_corners$to),
          wr_bias  = mle[["wr"]] - (truth_corners$tr - truth_corners$tp)
        )
      }, error = function(e) NULL)
    })
  })

  # Ensure we always return a tibble with the expected columns
  if (nrow(out) == 0) {
    out <- tibble::tibble(
      scenario = integer(), sim = integer(),
      tp_bias = double(), dp_bias = double(),
      wp_bias = double(), wr_bias = double()
    )
    warning("process_sim_results: all fits failed; returning empty tibble")
  }

  out
}


#' Plot simulation study results
#'
#' @param simsum   Output of process_sim_results()
#' @param out_dir  Output directory
#' @return Character vector of saved file paths
plot_sim_results <- function(simsum, out_dir = "output/figures") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  f <- file.path(out_dir, "simulation_results.pdf")

  # Guard: if no successful fits, write a placeholder and return early
  if (is.null(simsum) || nrow(simsum) == 0 ||
      !any(grepl("_bias$", names(simsum)))) {
    grDevices::pdf(f, width = 10, height = 6)
    plot.new()
    text(0.5, 0.5, "No simulation results to plot\n(all fits failed)",
         cex = 1.5)
    grDevices::dev.off()
    return(f)
  }

  bias_cols <- grep("_bias$", names(simsum), value = TRUE)

  plot_dat <- simsum |>
    tidyr::pivot_longer(
      tidyselect::all_of(bias_cols),
      names_to = "parameter",
      values_to = "bias"
    )

  p <- ggplot2::ggplot(plot_dat,
                        ggplot2::aes(x = factor(scenario), y = bias)) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Scenario (noise level)", y = "Bias",
                  title = "Simulation study: piecewise approximation bias")

  f <- file.path(out_dir, "simulation_results.pdf")
  ggplot2::ggsave(f, p, width = 10, height = 6)

  f
}
