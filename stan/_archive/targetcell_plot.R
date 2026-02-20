smp <- filter(ataccc_dat, id == 52 & time >= -5 & time <= 10)
mod <- cmdstan_model("1_code/stan/targetcell.stan")
times <- seq(-5, 10, 0.01)

fit <- mod$optimize(
  data = list(
    T = nrow(smp),
    y = smp[, c("pfu", "rna")],
    ts = smp$time,
    y0 = c(8e8, 0, 1, 0),
    lod = c(2.3, ct_to_rna(40, type = "ata") + 0.01),
    Tp = length(times),
    tsp = times
  ),
  init = list(
    list(
      b = 0.00000578,
      d = 2.36,
      c = 974,
      p = 2.15, 
      t0 = -5.5,
      q = 100,
      e = 5,
      sigma = 2,
      dp = 14,
      tp = -1,
      wp = 3,
      wr = 6,
      sigma_pe = 4,
      tau_tp = 0,
      tau_dp = 0.3,
      tau_wp = 0.7,
      tau_wr = 1,
      tp0 = 0,
      dp0 = 1,
      wp0 = 1,
      wr0 = 1
    )),
  iter = 10000
)

y_sim_ode <- fit$mle(variables = "y_sim_ode")
y_sim_pe <- fit$mle(variables = "y_sim_pe")

smp_pred <-
  data.frame(
    time = times, 
    pfu_ode = y_sim_ode[1:length(times)],
    rna_ode = y_sim_ode[(length(times) + 1):(length(times) * 2)],
    pfu_pe = y_sim_pe[1:length(times)],
    rna_pe = y_sim_pe[(length(times) + 1):(length(times) * 2)]
  )

smp_pred <-
  mutate(smp_pred, 
         pfu_ode = replace(pfu_ode, pfu_ode < 2.3, NA),
         pfu_pe = replace(pfu_pe, pfu_pe < 2.3, NA)
         )

ggplot(hct_dat, aes(x = time, y = rna)) +
  facet_wrap(~id) +
  geom_point(color = "blue") +
  # geom_point(aes(y = rna), color = "red") + 
  geom_point(aes(y = rna_ct), color = "green") + 
  geom_point(aes(y = rna_cn), color = "magenta") +
  theme_bw()
  
  geom_line(data = smp_pred, aes(y = pfu_ode)) +
  geom_line(data = smp_pred, aes(y = pfu_pe), linetype = "dashed") +
  geom_line(data = smp_pred, aes(y = rna_ode)) +
  geom_line(data = smp_pred, aes(y = rna_pe), linetype = "dashed") +
  labs(
    x = "Time",
    y = "count"
  ) +
  theme_bw() 
