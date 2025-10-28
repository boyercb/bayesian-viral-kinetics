mod <- cmdstan_model("1_code/stan/double_triangle.stan")

# fit <- mod$sample(
#   data = list(
#     K = 4,
#     M = stan_data$N_ids,
#     N = stan_data$N,
#     lod_rna = stan_data$lod_rna,
#     lod_pfu = stan_data$lod_pfu,
#     id = stan_data$id,
#     time = stan_data$time,
#     rna = stan_data$rna,
#     pfu = stan_data$pfu,
#     source = stan_data$source,
#     rna_exist = stan_data$rna_exist,
#     pfu_exist = stan_data$pfu_exist
#   ),
#   seed = 134,
#   chains = 4
# )

smp <- filter(ataccc_dat) |>
  mutate(
    rna = replace(rna, rna_exist == 0, 0),
    pfu = replace(pfu, pfu_exist == 0, 0),
    lfd = replace(lfd, lfd_exist == 0, 0)
  )

smp2 <- filter(hct_dat) |>
  mutate(
    rna = replace(rna, rna_exist == 0, 0),
    pfu = replace(pfu, pfu_exist == 0, 0),
    lfd = replace(lfd, lfd_exist == 0, 0)
  )

smp3 <- filter(uiuc_dat) |>
  mutate(
    rna = replace(rna, rna_exist == 0, 0),
    pfu = replace(pfu, pfu_exist == 0, 0),
    lfd = replace(lfd, lfd_exist == 0, 0)
  )

smp_stack <- 
  bind_rows(
    smp,
    smp2,
    smp3,
    .id = "source"
  ) |>
  group_by(source, id) |>
  mutate(id = cur_group_id(),
         pfu_type = ifelse(source == 3, 2, 1)) |>
  ungroup() |>
  arrange(id)

fit <- mod$optimize(
  data = list(
    K = 3,
    M = c(
      length(unique(smp$id)),
      length(unique(smp2$id)),
      length(unique(smp3$id))
    ),
    N = c(
      nrow(smp),
      nrow(smp2),
      nrow(smp3)
    ),
    lod_rna = c(
      ct_to_rna(40, type = "ata") + 0.01,
      0,
      1 / (ct_to_rna(47 - 1/log(10), type = "uiuc-ct") - ct_to_rna(47, type = "uiuc-ct"))
      #
    ),
    lod_pfu = c(
      2.3,
      0,
      2.3
      #log(5)
    ),
    prior_fn = 0.05,
    prior_fp = 0.05,
    fp_mean = c(
      1 / (ct_to_rna(40 - 1/log(10), type = "ata") - ct_to_rna(40, type = "ata")),
      1 / 0.3,
      1 / (ct_to_rna(47 - 1/log(10), type = "uiuc-ct") - ct_to_rna(47, type = "uiuc-ct"))
    ),
    id = smp_stack$id,
    time = smp_stack$time,
    rna = smp_stack$rna,
    pfu = smp_stack$pfu,
    source = smp_stack$source,
    pfu_type = smp_stack$pfu_type,
    rna_exist = smp_stack$rna_exist,
    pfu_exist = smp_stack$pfu_exist,
    test_error = 1
  ),
  init = list(
    list(
      sigma_rna = 2,
      sigma_pfu = 3,
      fp = 0.9,
      fn = 0.9,
      alpha = c(0, 0),
      tau_tp = 0,
      tau_dp = 0,
      tau_wp = 0, 
      tau_wr = 0,
      tau0_tp = 0,
      tau0_dp = 0,
      tau0_wp = 0, 
      tau0_wr = 0,
      dp_mean_pfu = 10,
      wp_mean_pfu = 5,
      wr_mean_pfu = 5,
      tp_mean_pfu = 0,
      tp_i_pfu = rep(0, length(unique(smp_stack$id))),
      dp_i_pfu = rep(0, length(unique(smp_stack$id))),
      wp_i_pfu = rep(0, length(unique(smp_stack$id))),
      wr_i_pfu = rep(0, length(unique(smp_stack$id))),
      tp_i_rna = rep(0, length(unique(smp_stack$id))),
      dp_i_rna = rep(0, length(unique(smp_stack$id))),
      wp_i_rna = rep(0, length(unique(smp_stack$id))),
      wr_i_rna = rep(0, length(unique(smp_stack$id)))
    )
  ),
  iter = 10000
)

# fit2 <- mod$optimize(
#   data = list(
#     K = 2,
#     M = c(
#       length(unique(smp$id)),
#       length(unique(smp2$id))
#     ),
#     N = c(
#       nrow(smp),
#       nrow(smp2)
#     ),
#     lod_rna = c(
#       ct_to_rna(40, type = "ata") + 0.01,
#       0
#       #
#     ),
#     lod_pfu = c(
#       2.3,
#       0
#       #log(5)
#     ),
#     prior_fn = 0.05,
#     prior_fp = 0.05,
#     fp_mean = c(
#       1 / (ct_to_rna(40 - 1/log(10), type = "ata") - ct_to_rna(40, type = "ata")),
#       1 / 0.3
#     ),
#     id = smp_stack$id,
#     time = smp_stack$time,
#     rna = smp_stack$rna,
#     pfu = smp_stack$pfu,
#     source = smp_stack$source,
#     pfu_type = smp_stack$pfu_type,
#     rna_exist = smp_stack$rna_exist,
#     pfu_exist = smp_stack$pfu_exist,
#     test_error = 0
#   ),
#   init = list(
#     list(
#       sigma_rna = 2,
#       sigma_pfu = 3,
#       tau_tp = 0,
#       tau_dp = 0,
#       tau_wp = 0, 
#       tau_wr = 0,
#       dp_mean_pfu = 10,
#       wp_mean_pfu = 5,
#       wr_mean_pfu = 5,
#       tp_mean_pfu = 0,
#       tp_i_pfu = rep(0, length(unique(smp_stack$id))),
#       dp_i_pfu = rep(0, length(unique(smp_stack$id))),
#       wp_i_pfu = rep(0, length(unique(smp_stack$id))),
#       wr_i_pfu = rep(0, length(unique(smp_stack$id))),
#       tp_i_rna = rep(0, length(unique(smp_stack$id))),
#       dp_i_rna = rep(0, length(unique(smp_stack$id))),
#       wp_i_rna = rep(0, length(unique(smp_stack$id))),
#       wr_i_rna = rep(0, length(unique(smp_stack$id)))
#     )
#   ),
#   iter = 10000
# )

predict_individual_effects <- function(fit, id, type = "mle") {
  if (type == "mle") {
    dp_pfu <- fit$mle("dp_mean_pfu") * exp(fit$mle("dp_i_pfu"))
    wp_pfu <- fit$mle("wp_mean_pfu") * exp(fit$mle("wp_i_pfu"))
    wr_pfu <- fit$mle("wr_mean_pfu") * exp(fit$mle("wr_i_pfu"))
    tp_pfu <- fit$mle("tp_i_pfu")

    dp_rna <- (exp(fit$mle("tau0_dp")) + exp(fit$mle("tau_dp")) * dp_pfu) * exp(fit$mle("dp_i_rna"))
    wp_rna <- (exp(fit$mle("tau0_wp")) + exp(fit$mle("tau_wp")) * wp_pfu) * exp(fit$mle("wp_i_rna"))
    wr_rna <- (exp(fit$mle("tau0_wr")) + exp(fit$mle("tau_wr")) * wr_pfu) * exp(fit$mle("wr_i_rna"))
    tp_rna <- fit$mle("tau0_tp") + fit$mle("tau_tp") * tp_pfu + fit$mle("tp_i_rna")


  }

  data.frame(
    id,
    dp_pfu,
    wp_pfu,
    wr_pfu,
    tp_pfu,
    dp_rna,
    wp_rna,
    wr_rna,
    tp_rna,
    row.names = NULL
  )
}
# predict_individual_effects <- function(fit, id, type = "mle") {
#   if (type == "mle") {
#     dp_pfu <- fit$mle("dp_mean_pfu") + fit$mle("dp_i_pfu")
#     wp_pfu <- fit$mle("wp_mean_pfu") + fit$mle("wp_i_pfu")
#     wr_pfu <- fit$mle("wr_mean_pfu") + fit$mle("wr_i_pfu")
#     tp_pfu <- fit$mle("tp_i_pfu")
#     
#     dp_rna <- fit$mle("tau0_dp") + fit$mle("tau_dp") * dp_pfu + fit$mle("dp_i_rna")
#     wp_rna <- fit$mle("tau0_wp") + fit$mle("tau_wp") * wp_pfu + fit$mle("wp_i_rna")
#     wr_rna <- fit$mle("tau0_wr") + fit$mle("tau_wr") * wr_pfu + fit$mle("wr_i_rna")
#     tp_rna <- fit$mle("tau0_tp") + fit$mle("tau_tp") * tp_pfu + fit$mle("tp_i_rna")
#   
#     
#   }
#   
#   data.frame(
#     id,
#     dp_pfu,
#     wp_pfu,
#     wr_pfu,
#     tp_pfu,
#     dp_rna,
#     wp_rna,
#     wr_rna,
#     tp_rna,
#     row.names = NULL
#   )
# }
pred <- predict_individual_effects(fit, id = unique(smp_stack$id))

pred$lod_rna <-
  ifelse(pred$id <= 45, ct_to_rna(40, type = "ata") + 0.01, 0)

pred$lod_pfu <- 
  ifelse(pred$id <= 45, 2.3, 0)

# pred2$lod_rna <-
#   ifelse(pred$id <= 45, ct_to_rna(40, type = "ata") + 0.01, 0)
# 
# pred2$lod_pfu <-
#   ifelse(pred$id <= 45, 2.3, 0)


ggplot(pred, aes(group = id)) +
  facet_wrap(~id) +
  geom_point(
    data = smp_stack,
    aes(x = time, y = pfu),
    color = "#1F78B4",
    alpha = 0.4,
    shape = 16
  ) +
  geom_segment(aes(
    x = -Inf,
    y = lod_pfu,
    xend = tp_pfu - wp_pfu,
    yend = lod_pfu
  ),
  color = "#1F78B4") +
  geom_segment(aes(
    x = pmax(tp_pfu - wp_pfu, -10),
    y = lod_pfu,
    xend = tp_pfu,
    yend = dp_pfu
  ), color = "#1F78B4") +
  geom_segment(aes(
    x = tp_pfu,
    y = dp_pfu,
    xend = pmin(tp_pfu + wr_pfu, 20),
    yend = lod_pfu
  ), color = "#1F78B4") +
  geom_segment(aes(
    x = tp_pfu + wr_pfu,
    y = lod_pfu,
    xend = Inf,
    yend = lod_pfu
  ),
  color = "#1F78B4") +
  geom_point(
    data = smp_stack,
    aes(x = time, y = rna),
    color = "#FF7F00",
    alpha = 0.4,
    shape = 16
  ) +
  geom_segment(aes(
    x = -Inf,
    y = lod_rna,
    xend = tp_rna - wp_rna,
    yend = lod_rna
  ),
  color = "#FF7F00") +
  geom_segment(aes(
    x = pmax(tp_rna - wp_rna, -10),
    y = lod_rna,
    xend = tp_rna,
    yend = dp_rna
  ),
  color = "#FF7F00") +
  geom_segment(aes(
    x = tp_rna,
    y = dp_rna,
    xend = pmin(tp_rna + wr_rna, 20),
    yend = lod_rna
  ),
  color = "#FF7F00") +
  geom_segment(aes(
    x = tp_rna + wr_rna,
    y = lod_rna,
    xend = Inf,
    yend = lod_rna
  ),
  color = "#FF7F00") +
  # geom_segment(aes(
  #   x = -Inf,
  #   y = lod_pfu,
  #   xend = tp_pfu - wp_pfu,
  #   yend = lod_pfu
  # ), data = pred2,
  # color = "#1F78B4", linetype = "dashed") +
  # geom_segment(aes(
  #   x = pmax(tp_pfu - wp_pfu, -10),
  #   y = lod_pfu,
  #   xend = tp_pfu,
  #   yend = dp_pfu
  # ), data = pred2, color = "#1F78B4", linetype = "dashed") +
  # geom_segment(aes(
  #   x = tp_pfu,
  #   y = dp_pfu,
  #   xend = pmin(tp_pfu + wr_pfu, 20),
  #   yend = lod_pfu
  # ), data = pred2, color = "#1F78B4", linetype = "dashed") +
  # geom_segment(aes(
  #   x = tp_pfu + wr_pfu,
  #   y = lod_pfu,
  #   xend = Inf,
  #   yend = lod_pfu
  # ), data = pred2,
  # color = "#1F78B4", linetype = "dashed") +
  # geom_segment(aes(
  #   x = -Inf,
  #   y = lod_rna,
  #   xend = tp_rna - wp_rna,
  #   yend = lod_rna
  # ), data = pred2,
  # color = "#FF7F00", linetype = "dashed") +
  # geom_segment(aes(
  #   x = pmax(tp_rna - wp_rna, -10),
  #   y = lod_rna,
  #   xend = tp_rna,
  #   yend = dp_rna
  # ), data = pred2,
  # color = "#FF7F00", linetype = "dashed") +
  # geom_segment(aes(
  #   x = tp_rna,
  #   y = dp_rna,
  #   xend = pmin(tp_rna + wr_rna, 20),
  #   yend = lod_rna
  # ), data = pred2,
  # color = "#FF7F00", linetype = "dashed") +
  # geom_segment(aes(
  #   x = tp_rna + wr_rna,
  #   y = lod_rna,
  #   xend = Inf,
  #   yend = lod_rna
  # ), data = pred2,
  # color = "#FF7F00", linetype = "dashed") +
  xlim(c(-10, 20)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    axis.text.y = element_text(size = 8)
  ) 


ggplot(pred, aes(x = dp_pfu, y = dp_rna)) +
  geom_point() +
  geom_smooth(method = "lm") 

ggplot(pred |> mutate(dp_i_rna = dp_pfu * exp(fit$mle("dp_i_rna") + fit$mle("tau_dp"))), 
       aes(x = dp_pfu, y = dp_i_rna)) +
  geom_point() +
  geom_smo


ggplot(d, aes(x = counter, y = 40 - replace(saliva_ct_n, saliva_ct_n == 0, 40))) +
  facet_wrap(~user_id) +
  geom_point()






