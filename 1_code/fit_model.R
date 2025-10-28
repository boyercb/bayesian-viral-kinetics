# compile model
mod <- cmdstan_model("1_code/stan/generative_model.stan")

# run MCMC
fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  init = lapply(1:4, function(x) {
    list(
      tp_i_pfu = rep(0, sum(stan_data$M)),
      dp_i_pfu = rep(0, sum(stan_data$M)),
      wp_i_pfu = rep(0, sum(stan_data$M)),
      wr_i_pfu = rep(0, sum(stan_data$M)),
      tp_i_rna = rep(0, sum(stan_data$M)),
      dp_i_rna = rep(0, sum(stan_data$M)),
      wp_i_rna = rep(0, sum(stan_data$M)),
      wr_i_rna = rep(0, sum(stan_data$M)),
      # tp_k_pfu = rep(0, stan_data$K),
      # dp_k_pfu = rep(0, stan_data$K),
      # wp_k_pfu = rep(0, stan_data$K),
      # wr_k_pfu = rep(0, stan_data$K),
      # tp_k_rna = rep(0, stan_data$K),
      # dp_k_rna = rep(0, stan_data$K),
      # wp_k_rna = rep(0, stan_data$K),
      # wr_k_rna = rep(0, stan_data$K),
      sigma_rna = 2,
      sigma_pfu = 3,
      fp = 0.05,
      fn = 0.01,
      dp_raw = 0,
      wp_raw = 0,
      wr_raw = 0,
      beta_dp_rna = rep(0, stan_data$P),
      beta_wp_rna = rep(0, stan_data$P),
      beta_wr_rna = rep(0, stan_data$P),
      # beta_lfd = rep(0, stan_data$P),
      tau_tp = c(0, 1),
      tau_dp = c(0, -0.5),
      tau_wp = c(0, -0.5),
      tau_wr = c(0, -0.5),
      tau_lfd = rep(0, 2),
      tau0_lfd_raw = 0,
      
      alpha_tcid50 = rep(1, 4),
      alpha_cult = rep(1, 2),
      
      to_i_sym = rep(0, sum(stan_data$M)),
      # to_k_sym = rep(0, stan_data$K),
      tau_sym = c(0, 1),
      sigma_to_sym = 3
    )

  }),
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 2000
)

fit$save_object("1_code/stan/kinetics_model.rds")

# MLE for comparison
mle <- mod$optimize(
  data = stan_data,
  init = list(
    list(
      tp_i_pfu = rep(0, sum(stan_data$M)),
      dp_i_pfu = rep(0, sum(stan_data$M)),
      wp_i_pfu = rep(0, sum(stan_data$M)),
      wr_i_pfu = rep(0, sum(stan_data$M)),
      tp_i_rna = rep(0, sum(stan_data$M)),
      dp_i_rna = rep(0, sum(stan_data$M)),
      wp_i_rna = rep(0, sum(stan_data$M)),
      wr_i_rna = rep(0, sum(stan_data$M)),
      # tp_k_pfu = rep(0, stan_data$K),
      # dp_k_pfu = rep(0, stan_data$K),
      # wp_k_pfu = rep(0, stan_data$K),
      # wr_k_pfu = rep(0, stan_data$K),
      # tp_k_rna = rep(0, stan_data$K),
      # dp_k_rna = rep(0, stan_data$K),
      # wp_k_rna = rep(0, stan_data$K),
      # wr_k_rna = rep(0, stan_data$K),
      sigma_rna = 2,
      sigma_pfu = 3,
      fp = 0.05,
      fn = 0.01,
      dp_raw = 0,
      wp_raw = 0,
      wr_raw = 0,
      beta_dp_rna = rep(0, stan_data$P),
      beta_wp_rna = rep(0, stan_data$P),
      beta_wr_rna = rep(0, stan_data$P),
      # beta_lfd = rep(0, stan_data$P),
      tau_tp = c(0, 1),
      tau_dp = c(0, -0.5),
      tau_wp = c(0, -0.5),
      tau_wr = c(0, -0.5),
      tau_lfd = rep(0, 2),
      tau0_lfd_raw = 0,
      
      alpha_tcid50 = rep(1, 4),
      alpha_cult = rep(1, 2),
      
      to_i_sym = rep(0, sum(stan_data$M)),
      # to_k_sym = rep(0, stan_data$K),
      tau_sym = c(0, 1),
      sigma_to_sym = 3
    )
  ),
  iter = 10000
)


# d <- list()
# for (i in 1:15) {
#   # MLE for comparison
#   mle <- mod$optimize(
#     data = stan_data,
#     init = list(
#       list(
#         log_dpi = rep(0, stan_data$N_ids),
#         log_wpi = rep(0, stan_data$N_ids),
#         log_wri = rep(0, stan_data$N_ids),
#         tp_rna = rep(0, stan_data$N_ids),
#         beta_dp = rep(0, stan_data$P),
#         beta_wp = rep(0, stan_data$P),
#         beta_wr = rep(0, stan_data$P),
#         theta_tp = 0,
#         log_dp_raw = 0,
#         log_wp_raw = 0,
#         log_wr_raw = 0,
#         gamma = rep(0, 2),
#         gamma0_raw = 0,
#         rho_tp = 0,
#         rho_dp = -1,
#         rho_wp = -1,
#         rho_wr = -1,
#         sigma_rna = 3,
#         sigma_pfu = 2
#       )
#     ),
#     iter = 10000
#   )
# 
#   d[[i]] <- mle$summary(variables = c("dp_mean", "wp_mean", "wr_mean"))
# }
# 
# bind_rows(d, .id = "sim") |>
#   arrange(variable) |> View()
