
set.seed(7324650)

# sample just a few NBA trajectories to plot
nba_plot_smp <- sample(unique(stacked_dat$id[stacked_dat$source == 1]), 42)

# ATACCC sample
ata_plot_smp <- unique(stacked_dat$id[stacked_dat$source == 2])

# UIUC sample
uiuc_plot_smp <- unique(stacked_dat$id[stacked_dat$source == 3])
  
# HCT sample
hct_plot_smp <- unique(stacked_dat$id[stacked_dat$source == 4])

# create prediction data frame
pred_dat <-
  filter(stacked_dat,
         id %in% c(ata_plot_smp, nba_plot_smp, uiuc_plot_smp, hct_plot_smp))

# sample from posterior
p <- predict_kinetics(fit, newdata = pred_dat)


# ATACCC plot -------------------------------------------------------------

# RNA copies
pred_dat$rna_hat <- mean(p$rna_hat)
pred_dat$rna_hat_q1 <- quantile2(p$rna_hat, 0.025)
pred_dat$rna_hat_q99 <- quantile2(p$rna_hat, 0.975)
pred_dat$rna_hat <- replace(pred_dat$rna_hat, pred_dat$rna_hat < stan_data$lod_rna[2], stan_data$lod_rna[2])
pred_dat$rna_hat_q1 <- replace(pred_dat$rna_hat_q1, pred_dat$rna_hat_q1 < stan_data$lod_rna[2], stan_data$lod_rna[2])
pred_dat$rna_hat_q99 <- replace(pred_dat$rna_hat_q99, pred_dat$rna_hat_q99 < stan_data$lod_rna[2], stan_data$lod_rna[2])

rna_corners <- tibble(
  id = rep(pred_dat$id, 3),
  time = c(
    mean(p$to_rna),
    mean(p$tp_rna),
    mean(p$tr_rna)
  ),
  time_q1 = c(
    quantile2(p$to_rna, 0.025),
    quantile2(p$tp_rna, 0.025),
    quantile2(p$tr_rna, 0.025)
  ),
  time_q99 = c(
    quantile2(p$to_rna, 0.975),
    quantile2(p$tp_rna, 0.975),
    quantile2(p$tr_rna, 0.975)
  ),
  rna_hat = c(
    mean(p$to_rna_hat),
    mean(p$tp_rna_hat),
    mean(p$tr_rna_hat)
  ),
  rna_hat_q1 = c(
    quantile2(p$to_rna_hat_mean, 0.025),
    quantile2(p$tp_rna_hat_mean, 0.025),
    quantile2(p$tr_rna_hat_mean, 0.025)
  ),
  rna_hat_q99 = c(
    quantile2(p$to_rna_hat_mean, 0.975),
    quantile2(p$tp_rna_hat_mean, 0.975),
    quantile2(p$tr_rna_hat_mean, 0.975)
  ),
  type = rep(c("o", "p", "r"), each = length(p$rna_hat))
)

rna_corners <- 
  rna_corners |>
  group_by(id, type) |>
  slice(1) |>
  ungroup()

# viral culture
pred_dat$pfu_hat <- mean(p$pfu_hat)
pred_dat$pfu_hat_q1 <- quantile2(p$pfu_hat, 0.025)
pred_dat$pfu_hat_q99 <- quantile2(p$pfu_hat, 0.975)
pred_dat$pfu_hat <- replace(pred_dat$pfu_hat, pred_dat$pfu_hat < stan_data$lod_pfu[2], stan_data$lod_pfu[2])
pred_dat$pfu_hat_q1 <- replace(pred_dat$pfu_hat_q1, pred_dat$pfu_hat_q1 < stan_data$lod_pfu[2], stan_data$lod_pfu[2])
pred_dat$pfu_hat_q99 <- replace(pred_dat$pfu_hat_q99, pred_dat$pfu_hat_q99 < stan_data$lod_pfu[2], stan_data$lod_pfu[2])

pfu_corners <- tibble(
  id = rep(pred_dat$id, 3),
  time = c(
    mean(p$to_pfu),
    mean(p$tp_pfu),
    mean(p$tr_pfu)
  ),
  time_q1 = c(
    quantile2(p$to_pfu, 0.025),
    quantile2(p$tp_pfu, 0.025),
    quantile2(p$tr_pfu, 0.025)
  ),
  time_q99 = c(
    quantile2(p$to_pfu, 0.975),
    quantile2(p$tp_pfu, 0.975),
    quantile2(p$tr_pfu, 0.975)
  ),
  pfu_hat = c(
    mean(p$to_pfu_hat),
    mean(p$tp_pfu_hat),
    mean(p$tr_pfu_hat)
  ),
  pfu_hat_q1 = c(
    quantile2(p$to_pfu_hat_mean, 0.025),
    quantile2(p$tp_pfu_hat_mean, 0.025),
    quantile2(p$tr_pfu_hat_mean, 0.025)
  ),
  pfu_hat_q99 = c(
    quantile2(p$to_pfu_hat_mean, 0.975),
    quantile2(p$tp_pfu_hat_mean, 0.975),
    quantile2(p$tr_pfu_hat_mean, 0.975)
  ),
  type = rep(c("o", "p", "r"), each = length(p$pfu_hat))
)

pfu_corners <- 
  pfu_corners |>
  group_by(id, type) |>
  slice(1) |>
  ungroup()

pred_dat$lfd_hat <- mean(p$lfd_hat)

plot_dat <- bind_rows(
  pred_dat,
  rna_corners,
  pfu_corners
) |>
  arrange(id, time) |>
  group_by(id) |>
  fill(pid, all_of(x_vars), source, .direction = "updown") |>
  ungroup()
  

# leave out 8, 39, and 45 for now from ATACCC
ata_plot_smp <- ata_plot_smp[!ata_plot_smp %in% c(8, 39, 45)]

g1 <-
  ggplot(plot_dat[plot_dat$source == 2, ]) +
  facet_wrap(~id) +#~factor(id, labels = unique(paste0("ID: ", id)))) +
  geom_rect(aes(
    xmin = time - 0.4,
    xmax = time + 0.4,
    ymin = 29,
    ymax = 32,
    fill = lfd_hat
  ),
  data = plot_dat[plot_dat$source == 2 &
                    !is.na(plot_dat$lfd_hat), ]) +
  geom_text(
    aes(
      x = time,
      y = 34,
      label = ifelse(lfd == 1, "X", "")
    ),
    data = plot_dat[plot_dat$source == 2 &
                      !is.na(plot_dat$lfd_hat), ],
    color = "black",
    size = 2.5
  ) +
  colorspace::scale_fill_continuous_sequential(
    name = "Probability LFD +",
    palette = "Emrld", 
    rev = FALSE
  ) +
  geom_point(
    aes(x = time, y = rna, group = id),
    color = "#4ca5ff",
    alpha = 0.75,
    shape = 16
  ) +
  geom_line(aes(x = time, y = rna_hat, group = id),
            data = plot_dat[plot_dat$source == 2 & !is.na(plot_dat$rna_hat),],
            color = "#4ca5ff") +
  geom_ribbon(
    aes(
      x = time,
      y = rna_hat,
      ymin = rna_hat_q99,
      ymax = rna_hat_q1,
      group = id
    ),
    data = plot_dat[plot_dat$source == 2 & !is.na(plot_dat$rna_hat),],
    alpha = 0.3,
    fill = "#4ca5ff"
  ) +
  geom_point(
    aes(x = time, y = pfu, group = id),
    color = "red",
    alpha = 0.5,
    shape = 16
  ) +
  geom_line(aes(x = time, y = pfu_hat, group = id),
            data = plot_dat[plot_dat$source == 2 &
                              !is.na(plot_dat$pfu_hat), ],
            color = "red") +
  geom_ribbon(
    aes(
      x = time,
      y = pfu_hat,
      ymin = pfu_hat_q1,
      ymax = pfu_hat_q99,
      group = id
    ),
    data = plot_dat[plot_dat$source == 2 &
                      !is.na(plot_dat$pfu_hat), ],
    alpha = 0.3,
    fill = "red"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal"
  ) +
  labs(
    x = "days from peak",
    y = "log count per ml"
  )

ggsave(
  filename = "3_figures/ataccc_fit.pdf",
  plot = g1,
  device = "pdf",
  width = 17,
  height = 10
)
g1


# NBA plot ----------------------------------------------------------------

g2 <- ggplot(pred_dat[pred_dat$source == 2, ], aes(x = time, y = rna, group = id)) +
  facet_wrap(~factor(id, labels = unique(paste0("ID: ", id)))) +
  geom_rect(aes(
    xmin = time - 0.4,
    xmax = time + 0.4,
    ymin = 27,
    ymax = 30,
    fill = lfd_hat
  )) +
  colorspace::scale_fill_continuous_sequential(
    name = "Probability LFD +",
    palette = "Emrld", 
    rev = FALSE
  ) +
  geom_point(color = "#4ca5ff",  alpha = 0.75, shape = 16) +
  geom_line(aes(x = time, y = rna_hat), color = "#4ca5ff") +
  geom_ribbon(
    aes(
      x = time,
      y = rna_hat,
      ymin = rna_hat_q99,
      ymax = rna_hat_q1
    ),
    alpha = 0.3,
    fill = "#4ca5ff"
  ) +
  geom_point(aes(x = time, y = pfu, group = id)) +
  geom_line(aes(x = time, y = pfu_hat), color = "red") +
  geom_ribbon(
    aes(
      x = time,
      y = pfu_hat,
      ymin = pfu_hat_q1,
      ymax = pfu_hat_q99
    ),
    alpha = 0.3,
    fill = "red"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal"
  ) +
  ylim(c(0, 30)) +
  labs(
    x = "days from peak",
    y = "log count per ml"
  )

ggsave(
  filename = "3_figures/nba_fit.pdf",
  plot = g2,
  device = "pdf",
  width = 17,
  height = 10
)

g2


# UIUC --------------------------------------------------------------------

g3 <-
  ggplot(plot_dat[plot_dat$source == 3, ]) +
  facet_wrap(~id) +#~factor(id, labels = unique(paste0("ID: ", id)))) +
  geom_rect(aes(
    xmin = time - 0.4,
    xmax = time + 0.4,
    ymin = 29,
    ymax = 32,
    fill = lfd_hat
  ),
  data = plot_dat[plot_dat$source == 3 &
                    !is.na(plot_dat$lfd_hat), ]) +
  geom_text(
    aes(
      x = time,
      y = 34,
      label = ifelse(lfd == 1, "X", "")
    ),
    data = plot_dat[plot_dat$source == 3 &
                      !is.na(plot_dat$lfd_hat), ],
    color = "black",
    size = 2.5
  ) +
  colorspace::scale_fill_continuous_sequential(
    name = "Probability LFD +",
    palette = "Emrld", 
    rev = FALSE
  ) +
  geom_point(
    aes(x = time, y = rna, group = id),
    color = "#4ca5ff",
    alpha = 0.75,
    shape = 16
  ) +
  geom_line(aes(x = time, y = rna_hat, group = id),
            data = plot_dat[plot_dat$source == 3 & !is.na(plot_dat$rna_hat),],
            color = "#4ca5ff") +
  geom_ribbon(
    aes(
      x = time,
      y = rna_hat,
      ymin = rna_hat_q99,
      ymax = rna_hat_q1,
      group = id
    ),
    data = plot_dat[plot_dat$source == 3 & !is.na(plot_dat$rna_hat),],
    alpha = 0.3,
    fill = "#4ca5ff"
  ) +
  geom_point(
    aes(x = time, y = pfu, group = id),
    color = "red",
    alpha = 0.5,
    shape = 16
  ) +
  geom_line(aes(x = time, y = pfu_hat, group = id),
            data = plot_dat[plot_dat$source == 3 &
                              !is.na(plot_dat$pfu_hat), ],
            color = "red") +
  geom_ribbon(
    aes(
      x = time,
      y = pfu_hat,
      ymin = pfu_hat_q1,
      ymax = pfu_hat_q99,
      group = id
    ),
    data = plot_dat[plot_dat$source == 3 &
                      !is.na(plot_dat$pfu_hat), ],
    alpha = 0.3,
    fill = "red"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal"
  ) +
  labs(
    x = "days from peak",
    y = "log count per ml"
  )

ggsave(
  filename = "3_figures/uiuc_fit.pdf",
  plot = g3,
  device = "pdf",
  width = 17,
  height = 10
)
g3


# HCT ---------------------------------------------------------------------


g4 <-
  ggplot(plot_dat[plot_dat$source == 4, ]) +
  facet_wrap(~id) +#~factor(id, labels = unique(paste0("ID: ", id)))) +
  geom_rect(aes(
    xmin = time - 0.4,
    xmax = time + 0.4,
    ymin = 29,
    ymax = 32,
    fill = lfd_hat
  ),
  data = plot_dat[plot_dat$source == 4 &
                    !is.na(plot_dat$lfd_hat), ]) +
  geom_text(
    aes(
      x = time,
      y = 34,
      label = ifelse(lfd == 1, "X", "")
    ),
    data = plot_dat[plot_dat$source == 4 &
                      !is.na(plot_dat$lfd_hat), ],
    color = "black",
    size = 2.5
  ) +
  colorspace::scale_fill_continuous_sequential(
    name = "Probability LFD +",
    palette = "Emrld", 
    rev = FALSE
  ) +
  geom_point(
    aes(x = time, y = rna, group = id),
    color = "#4ca5ff",
    alpha = 0.75,
    shape = 16
  ) +
  geom_line(aes(x = time, y = rna_hat, group = id),
            data = plot_dat[plot_dat$source == 4 & !is.na(plot_dat$rna_hat),],
            color = "#4ca5ff") +
  geom_ribbon(
    aes(
      x = time,
      y = rna_hat,
      ymin = rna_hat_q99,
      ymax = rna_hat_q1,
      group = id
    ),
    data = plot_dat[plot_dat$source == 4 & !is.na(plot_dat$rna_hat),],
    alpha = 0.3,
    fill = "#4ca5ff"
  ) +
  geom_point(
    aes(x = time, y = pfu, group = id),
    color = "red",
    alpha = 0.5,
    shape = 16
  ) +
  geom_line(aes(x = time, y = pfu_hat, group = id),
            data = plot_dat[plot_dat$source == 4 &
                              !is.na(plot_dat$pfu_hat), ],
            color = "red") +
  geom_ribbon(
    aes(
      x = time,
      y = pfu_hat,
      ymin = pfu_hat_q1,
      ymax = pfu_hat_q99,
      group = id
    ),
    data = plot_dat[plot_dat$source == 4 &
                      !is.na(plot_dat$pfu_hat), ],
    alpha = 0.3,
    fill = "red"
  ) +
  theme_minimal() +
  theme(
    legend.position = 'bottom',
    legend.direction = "horizontal"
  ) +
  labs(
    x = "days from peak",
    y = "log count per ml"
  )

ggsave(
  filename = "3_figures/hic_fit.pdf",
  plot = g4,
  device = "pdf",
  width = 17,
  height = 10
)
g4
