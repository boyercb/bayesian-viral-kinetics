variable_list <- c(
  # "tp_rna",
  # "tp_pfu",
  "dp_mean_rna",
  "wp_mean_rna",
  "wr_mean_rna",
  "beta_dp_rna",
  "beta_wp_rna",
  "beta_wr_rna",
  "tau_dp",
  "tau_wp",
  "tau_wr",
  "tau_tp"
  # "log_dpi",
  # "log_wpi",
  # "log_wri",
  # "gamma0",
  # "gamma"
)


kinetics <- as_draws_rvars(fit$draws(variables = variable_list))


map(
  c("dp", "wp", "wr"),
  function(x) {
    df <- tibble(
      label = sapply(x_vars, \(x) x_labels[[x]]),
      coef = specd(exp(E(kinetics[[paste0("beta_", x,  "_rna")]])), 2),
      ci = paste0("(", 
                  specd(exp(quantile(kinetics[[paste0("beta_", x, "_rna")]], 0.025)), 2),
                  ", ",
                  specd(exp(quantile(kinetics[[paste0("beta_", x, "_rna")]], 0.975)), 2),
                  ")"
      )
    )
    
    df <- add_row(df,
            label = "Reference value", 
            coef = specd(E(kinetics[[paste0(x, "_mean",  "_rna")]]), 2),
            ci = paste0("(", 
                        specd(quantile(kinetics[[paste0(x, "_mean_rna")]], 0.025), 2),
                        ", ",
                        specd(quantile(kinetics[[paste0(x, "_mean_rna")]], 0.975), 2),
                        ")")
    )
    
    kable(
      x = df,
      format = "latex",
      col.names = c("Characteristic", "$\\exp(\\beta)$", "95\\% CrI"),
      align = "lccc",
      booktabs = TRUE,
      linesep = "",
      escape = FALSE
    ) 
  }
)




df <- tibble(
  label = c(
    "Slope"
  ),
  dp_coef = specd(E(exp(kinetics$rho_dp)), 2),
  dp_ci = paste0("(", 
              specd(exp(quantile(kinetics$rho_dp, 0.025)), 2),
              ", ",
              specd(exp(quantile(kinetics$rho_dp, 0.975)), 2),
              ")"
  ),
  wp_coef = specd(E(exp(kinetics$rho_wp)), 2),
  wp_ci = paste0("(", 
                 specd(exp(quantile(kinetics$rho_wp, 0.025)), 2),
                 ", ",
                 specd(exp(quantile(kinetics$rho_wp, 0.975)), 2),
                 ")"
  ),
  wr_coef = specd(E(exp(kinetics$rho_wr)), 2),
  wr_ci = paste0("(", 
                 specd(exp(quantile(kinetics$rho_wr, 0.025)), 2),
                 ", ",
                 specd(exp(quantile(kinetics$rho_wr, 0.975)), 2),
                 ")"
  ),
  tp_coef = specd(E(kinetics$rho_tp), 2),
  tp_ci = paste0("(", 
                 specd(quantile(kinetics$rho_tp, 0.025), 2),
                 ", ",
                 specd(quantile(kinetics$rho_tp, 0.975), 2),
                 ")"
  )
)

df <- add_row(df,
        label = c(
          "Intercept"
        ),
        dp_coef = "",
        dp_ci = "",
        wp_coef = "",
        wp_ci = "",
        wr_coef = "",
        wr_ci = "",
        tp_coef = specd(E(kinetics$theta_tp), 2),
        tp_ci = paste0("(", 
                       specd(quantile(kinetics$theta_tp, 0.025), 2),
                       ", ",
                       specd(quantile(kinetics$theta_tp, 0.975), 2),
                       ")"
        ))

pivot_longer(df, -label) |>
  separate(name, c("param", "type")) |>
  pivot_wider(names_from = c(label, type), values_from = value) |>
  kable(
    format = "latex",
    col.names = c("Parameter", "$\\rho$", "95\\% CrI", "$\\theta$", "95\\% CrI"),
    align = "lcccc",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  )


df <- 
  tibble(
    gamma1_coef = specd(E(exp(kinetics$gamma[1])), 2),
    gamma1_ci = paste0("(", 
                   specd(exp(quantile(kinetics$gamma[1], 0.025)), 2),
                   ", ",
                   specd(exp(quantile(kinetics$gamma[1], 0.975)), 2),
                   ")"
    ),
    gamma2_coef = specd(E(exp(kinetics$gamma[2])), 2),
    gamma2_ci = paste0("(", 
                   specd(exp(quantile(kinetics$gamma[2], 0.025)), 2),
                   ", ",
                   specd(exp(quantile(kinetics$gamma[2], 0.975)), 2),
                   ")"
    ),
    gamma0_coef = specd(E(kinetics$gamma0), 2),
    gamma0_ci = paste0("(", 
                   specd(quantile(kinetics$gamma0, 0.025), 2),
                   ", ",
                   specd(quantile(kinetics$gamma0, 0.975), 2),
                   ")"
    )
  )

df <- 
  pivot_longer(df, everything()) |>
  separate(name, c("param", "type")) |>
  mutate(
    param = case_when(
      param == "gamma1" ~ "log RNA copies",
      param == "gamma2" ~ "log PFU culturable virus",
      param == "gamma0" ~ "Intercept"
    )
  ) |>
  pivot_wider(names_from = c(type), values_from = value)

kable(
  x = df,
  format = "latex",
  col.names = c("Parameter", "$\\exp(\\gamma)$", "95\\% CrI"),
  align = "lcc",
  booktabs = TRUE,
  linesep = "",
  escape = FALSE
)
