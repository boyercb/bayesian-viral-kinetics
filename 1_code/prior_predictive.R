

pp <- prior_predictive(stan_data, draws = 1000)

pp_plot <- 
  tibble(
    draw = 1:1000,
    dp_mean_rna = pp$dp_mean_rna,
    wp_mean_rna = pp$wp_mean_rna,
    wr_mean_rna = pp$wr_mean_rna,
    dp_mean_pfu = pp$dp_mean_pfu,
    wp_mean_pfu = pp$wp_mean_pfu,
    wr_mean_pfu = pp$wr_mean_pfu,
    tau0_tp = pp$tau0_tp,
    tau0_dp = pp$tau0_dp,
    tau0_wp = pp$tau0_wp,
    tau0_wr = pp$tau0_wr,
    tau0_sym = pp$tau0_sym,
    tau_tp = pp$tau_tp,
    tau_dp = pp$tau_dp,
    tau_wp = pp$tau_wp,
    tau_wr = pp$tau_wr,
    tau_sym = pp$tau_sym,
  )


p1 <- pp_plot |>
  pivot_longer(-draw) |>
  filter(str_detect(name, "tau")) |>
  mutate(name = factor(
    name,
    levels = c(
      "tau0_tp",
      "tau0_dp",
      "tau0_wp",
      "tau0_wr",
      "tau0_sym",
      "tau_tp",
      "tau_dp",
      "tau_wp",
      "tau_wr",
      "tau_sym"
    ),
    labels = c(
      "a0: peak time",
      "a0: peak",
      "a0: proliferation",
      "a0: clearance",
      "a0: symptom onset",
      "a1: peak time",
      "a1: peak",
      "a1: proliferation",
      "a1: clearance",
      "a1: symptom onset"
    )
  ),
  value = exp(value)
  ) |>
  filter(
    !name %in% c(
      "a0: peak time",
      "a1: peak time",
      "a0: symptom onset",
      "a1: symptom onset"
    )
  ) |>
  ggplot(aes(x = value)) +
  facet_wrap(~name, scales = "free") +
  scale_y_continuous(n.breaks = 3) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) + 
  theme_minimal() 


p2 <- pp_plot |>
  pivot_longer(-draw) |>
  filter(str_detect(name, "tau")) |>
  mutate(name = factor(
    name,
    levels = c(
      "tau0_tp",
      "tau0_dp",
      "tau0_wp",
      "tau0_wr",
      "tau0_sym",
      "tau_tp",
      "tau_dp",
      "tau_wp",
      "tau_wr",
      "tau_sym"
    ),
    labels = c(
      "a0: peak time",
      "a0: peak",
      "a0: proliferation",
      "a0: clearance",
      "a0: symptom onset",
      "a1: peak time",
      "a1: peak",
      "a1: proliferation",
      "a1: clearance",
      "a1: symptom onset"
    )
  )
  ) |>
  filter(
    name %in% c(
      "a0: peak time",
      "a1: peak time",
      "a0: symptom onset",
      "a1: symptom onset"
    )
  ) |>
  ggplot(aes(x = value)) +
  facet_wrap(~name, scales = "free") +
  scale_y_continuous(n.breaks = 3) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) + 
  theme_minimal()


p3 <- pp_plot |>
  pivot_longer(-draw) |>
  filter(!str_detect(name, "tau")) |>
  mutate(name = factor(
    name,
    levels = c(
      "dp_mean_rna",
      "wp_mean_rna",
      "wr_mean_rna",
      "dp_mean_pfu",
      "wp_mean_pfu",
      "wr_mean_pfu"
    ),
    labels = c(
      "RNA: peak",
      "RNA: proliferation",
      "RNA: clearance",
      "PFU: peak",
      "PFU: proliferation",
      "PFU: clearance"
    )
  )
  ) |>
  ggplot(aes(x = value)) +
  facet_wrap(~name, scales = "free") +
  scale_y_continuous(n.breaks = 3) +
  geom_histogram(fill = "lightblue", color = "black", bins = 30) + 
  theme_minimal() 

ggsave(
  filename = "3_figures/prior_trans_pe.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 5
)

ggsave(
  filename = "3_figures/prior_trans_times.pdf",
  plot = p2,
  device = "pdf",
  width = 8,
  height = 5
)

ggsave(
  filename = "3_figures/prior_pe.pdf",
  plot = p3,
  device = "pdf",
  width = 8,
  height = 5
)



pp_plot <-
  tibble(
    id = pp$id,
    time = pp$time,
    rna = pp$rna,
    pfu = pp$pfu
  )

pp_plot <- as_tibble(do.call(data.frame, pp_plot))

p4 <- 
  pp_plot |>
  pivot_longer(-c(id, time),
               names_to = c("variable", "draw"),
               names_sep = "\\.") |>
  ggplot(aes(x = time, y = value)) +
  ylim(c(0, 30)) +
  facet_wrap(~variable) +
  geom_hex(binwidth = c(1, 1)) +
  scale_fill_viridis_c(trans = "log", option = "A") +
  #scale_fill_gradient(trans = "log") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = "3_figures/prior_trajectories.pdf",
  plot = p4,
  device = "pdf",
  width = 9,
  height = 5
)

