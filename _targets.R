# ──────────────────────────────────────────────────────────────────────────────
# _targets.R — targets pipeline for Bayesian viral kinetics
# ──────────────────────────────────────────────────────────────────────────────
#
# Usage:
#   targets::tar_make()        — run the full pipeline
#   targets::tar_visnetwork()  — visualise the DAG
#   targets::tar_read(name)    — read a cached target
#   targets::tar_outdated()    — see what needs to be re-run
#
# ──────────────────────────────────────────────────────────────────────────────

library(targets)
library(tarchetypes)
library(stantargets)

tar_option_set(
  packages = c(
    "tidyverse",
    "cmdstanr",
    "posterior",
    "loo",
    "bayesplot",
    "knitr",
    "kableExtra",
    "truncnorm",
    "mvtnorm",
    "colorspace",
    "viridis",
    "deSolve",
    "purrr"
  ),
  # store CmdStanMCMC objects properly
  format = "qs"
)

# source all function definitions in R/
tar_source("R/")


# ══════════════════════════════════════════════════════════════════════════════
# Pipeline
# ══════════════════════════════════════════════════════════════════════════════

list(

  # ── Raw data ─────────────────────────────────────────────────────────────────

  tar_target(nba_raw,    readr::read_csv("data/nba_dat.csv",     show_col_types = FALSE)),
  tar_target(ataccc_raw, readr::read_csv("data/ataccc_dat.csv",  show_col_types = FALSE)),
  tar_target(ataccc_sym, readr::read_csv("data/ataccc_sx_dat.csv", show_col_types = FALSE)),
  tar_target(uiuc_raw,   readr::read_csv("data/uiuc_dat.csv",    show_col_types = FALSE)),
  tar_target(uiuc_sym,   readr::read_csv("data/uiuc_sx_dat.csv", show_col_types = FALSE)),
  tar_target(hct_raw,    readr::read_csv("data/hct_dat.csv",     show_col_types = FALSE)),
  tar_target(hct_sym,    readr::read_csv("data/hct_sx_dat.csv",  show_col_types = FALSE)),
  tar_target(legacy_raw, readr::read_csv("data/legacy_dat.csv",  show_col_types = FALSE)),


  # ── Clean each source ───────────────────────────────────────────────────────

  tar_target(nba_clean,    clean_nba(nba_raw)),
  tar_target(ataccc_clean, clean_ataccc(ataccc_raw, ataccc_sym)),
  tar_target(uiuc_clean,   clean_uiuc(uiuc_raw, uiuc_sym)),
  tar_target(hct_clean,    clean_hct(hct_raw, hct_sym)),
  tar_target(legacy_clean, clean_legacy(legacy_raw)),


  # ── Stack & build Stan data ─────────────────────────────────────────────────

  tar_target(
    stacked_dat,
    stack_sources(nba_clean, ataccc_clean, uiuc_clean, hct_clean, legacy_clean)
  ),

  tar_target(
    stan_data,
    build_stan_data(stacked_dat,
                    flags = list(
                      ind_effects = 1,
                      test_error  = 1,
                      adj_rna     = 1,
                      ind_corr    = 1,
                      use_smooth  = 1,
                      use_wf      = 0
                    ))
  ),


  # ── Prior predictive ────────────────────────────────────────────────────────

  tar_target(prior_pred, prior_predictive(stan_data, draws = 1000)),

  tar_target(
    fig_prior,
    plot_prior_predictive(prior_pred, stan_data),
    format = "file"
  ),


  # ── Model fitting ──────────────────────────────────────────────────────────
  #
  # Uses fit_model() which calls build_init() for proper initialisation.
  # Replaces tar_stan_mcmc() to allow init = build_init(stan_data).

  tar_target(
    kinetics_mcmc,
    fit_model("stan/kinetics_model.stan", stan_data,
              chains = 4, iter_warmup = 1000, iter_sampling = 4000,
              adapt_delta = 0.95, max_treedepth = 12,
              init_method = "map",
              threads_per_chain = 4)
  ),

  # ── Generated quantities (post-hoc) ────────────────────────────────────────
  #
  # Run a separate GQ-only Stan model on the posterior draws to compute

  # per-observation log-likelihoods and posterior predictive replicates.
  # This avoids inflating MCMC output during sampling.

  tar_target(
    gq_fit,
    run_gq(kinetics_mcmc, stan_data)
  ),


  # ── Diagnostics ─────────────────────────────────────────────────────────────

  tar_target(convergence, check_convergence(kinetics_mcmc)),
  tar_target(loo_result,  compute_loo(gq_fit)),
  tar_target(waic_result, compute_waic(gq_fit)),
  tar_target(ppc,         posterior_predictive_check(
    gq_fit, stacked_dat
  )),

  # k-fold CV grouped by individual — computationally expensive (~K × 43 hrs).
  # Uncomment to run; consider using subsample_stan_data() first.
  # tar_target(
  #   kfold_result,
  #   kfold_cv("stan/kinetics_model.stan", stan_data, stacked_dat,
  #            K_folds = 10, chains = 4, iter_warmup = 500,
  #            iter_sampling = 1000, adapt_delta = 0.95)
  # ),

  tar_target(
    fig_traces,
    plot_traces(kinetics_mcmc, stan_data = stan_data,
                out_file = "output/figures/trace_plots.pdf"),
    format = "file"
  ),

  tar_target(
    fig_diagnostics,
    plot_diagnostics(param_summary, ppc, convergence,
                     recovery_check = recovery_check),
    format = "file"
  ),


  # ── Correlation diagnostics (only meaningful when ind_corr = 1) ─────────────

  tar_target(
    corr_summary,
    summarize_correlation(kinetics_mcmc, stan_data)
  ),

  tar_target(
    fig_corr_matrix,
    plot_correlation_matrix(kinetics_mcmc, stan_data,
                            out_file = "output/figures/corr_matrix.pdf"),
    format = "file"
  ),

  tar_target(
    fig_corr_densities,
    plot_correlation_densities(kinetics_mcmc, stan_data,
                               out_file = "output/figures/corr_densities.pdf"),
    format = "file"
  ),


  # ── Summaries ───────────────────────────────────────────────────────────────

  tar_target(param_summary, summarize_parameters(
    kinetics_mcmc, stan_data
  )),
  tar_target(data_summary, summarize_data(stacked_dat)),

  tar_target(
    tex_table1,
    save_table(data_summary, "output/tables/table1.tex"),
    format = "file"
  ),
  tar_target(
    tex_params,
    save_table(param_summary, "output/tables/params.tex"),
    format = "file"
  ),
  tar_target(
    tex_supplement,
    save_supplement_tables(param_summary, "output/tables"),
    format = "file"
  ),


  # ── Predictions & trajectory plots ──────────────────────────────────────────

  tar_target(pred_dat,    make_prediction_data(stacked_dat)),
  tar_target(predictions, compute_predictions(
    kinetics_mcmc, pred_dat, stan_data
  )),
  tar_target(
    pop_draws,
    extract_pop_draws(kinetics_mcmc, n_draws = 200,
                      out_path = "output/pop_draws_200.rds")
  ),
  tar_target(
    fig_trajectories,
    plot_all_trajectories(predictions, stan_data),
    format = "file"
  ),
  tar_target(
    fig_inferred_pfu_plot,
    plot_inferred_pfu(predictions, stan_data),
    format = "file"
  ),


  # ── Policy analysis ────────────────────────────────────────────────────────
  #
  #  Generates trajectories for hypothetical agents under representative

  #  covariate profiles, aligns to clinically meaningful landmarks (first
  #  positive PCR, symptom onset, first positive LFD), and computes derived
  #  policy quantities: probability curves, conditional infectiousness,
  #  isolation duration tables, residual infectiousness AUC, and
  #  test-to-release false reassurance rates.

  tar_target(
    policy_results,
    compute_all_policy(kinetics_mcmc, stan_data,
                       n_draws = 200, n_reps = 50, seed = 2026)
  ),
  tar_target(
    fig_policy,
    generate_policy_figures(policy_results),
    format = "file"
  ),

  # ── Bayesian filtering (personalized infectiousness given test history) ─────
  #
  #  Generates latent trajectories from the population prior, then uses
  #  importance sampling to condition on hypothetical test histories (PCR Ct
  #  values and serial LFD results).  Produces Figure 9 showing how
  #  P(culture positive) updates with diagnostic information.

  tar_target(
    filter_results,
    compute_bayesian_filtering(kinetics_mcmc, stan_data,
                               n_draws = 200, n_reps = 500, seed = 2026)
  ),
  tar_target(
    fig_filter,
    generate_filtering_figures(filter_results),
    format = "file"
  ),


  # ── ODE comparison figures (for manuscript/presentation) ────────────────────

  tar_target(fig_ode, plot_ode_examples(), format = "file"),

  tar_target(
    fig_pub_main,
    generate_pub_figures(
      predictions = predictions,
      param_summary = param_summary,
      stan_data = stan_data,
      pop_draws_df = pop_draws,
      styles = c("pnas")
    ),
    format = "file"
  ),

  tar_target(
    fig_antigen_schematic,
    save_antigen_schematic(styles = c("pnas")),
    format = "file"
  ),

  tar_target(
    fig_site_analysis,
    generate_site_analysis_figures(styles = c("pnas")),
    format = "file"
  ),


  # ── Parameter recovery (simulated-data identifiability check) ──────────────
  #
  #  Uses 50% subsample of individuals for speed (~2× faster than full data).
  #  1. subsample the stan_data structure
  #  2. simulate_data() generates a synthetic dataset under known params
  #  3. fit the same Stan model on the simulated data
  #  4. check_recovery() compares posteriors to ground truth
  #  5. plot_recovery() visualises coverage

  tar_target(
    recovery_stan_data,
    subsample_stan_data(stan_data, frac = 0.5, seed = 99)
  ),

  tar_target(sim_params,   default_params(recovery_stan_data)),
  tar_target(sim_result,   simulate_data(recovery_stan_data, sim_params, seed = 42)),
  tar_target(sim_stan_data, sim_result$sim_data),
  tar_target(sim_truth,     sim_result$truth),

  tar_target(
    recovery_mcmc,
    fit_model("stan/kinetics_model.stan", sim_stan_data,
              chains = 4, iter_warmup = 1000, iter_sampling = 2000,
              adapt_delta = 0.95, max_treedepth = 12,
              init_method = "map",
              threads_per_chain = 4)
  ),

  tar_target(recovery_check, check_recovery(
    recovery_mcmc, sim_truth, stan_data = sim_stan_data
  )),

  tar_target(
    tex_manuscript_numbers,
    save_manuscript_macros(
      convergence = convergence,
      loo_result = loo_result,
      waic_result = waic_result,
      recovery_check = recovery_check,
      param_summary = param_summary,
      kinetics_mcmc = kinetics_mcmc,
      stacked_dat = stacked_dat,
      out_file = "output/tables/manuscript_numbers.tex"
    ),
    format = "file"
  ),

  tar_target(
    fig_recovery,
    plot_recovery(recovery_check,
                  out_file = "output/figures/param_recovery.pdf"),
    format = "file"
  )


  # ── Simulation study ────────────────────────────────────────────────────────

  # tar_target(sim_grid, build_sim_grid()),
  # tar_target(sim_results, run_simulation("stan/sim.stan", sim_grid, n_sims = 50)),
  # tar_target(sim_summary, process_sim_results(sim_results)),
  # tar_target(fig_sim, plot_sim_results(sim_summary), format = "file")
)
