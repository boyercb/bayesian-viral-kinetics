pefun <- function(t, tp, wp, wr, dp) {
  ifelse(
    t <= tp, 
    dp / wp * (t - (tp - wp)), 
    dp - dp / wr * (t - tp) 
  )
}

gendata <- function(n) {
  # define sample times and ids
  times <- seq(-6, 10, 0.5)
  id <- rep(1:n, each = length(times))
  x <- rep(times, n)
  
  # define test characteristics
  lod <- 1
  fnr <- 0.05
  fpr <- 0.05
  
  # draw individual parameters
  wp <- rnorm(n, 4, 1)
  wr <- rnorm(n, 7, 1)
  dp <- rnorm(n, 13, 1)
  
  # draw measurements 
  wp <- rep(wp, each = length(times))
  wr <- rep(wr, each = length(times))
  dp <- rep(dp, each = length(times))
  tp <- rep(0, length(times))
  
  # classic measurement error
  y <- rnorm(length(x), pefun(x, tp, wp, wr, dp), 1.5)
  
  # apply level of detection
  y <- ifelse(y < lod, lod, y)
  
  # add false positives and false negatives
  fp_raw <- rbinom(length(x), 1, fpr)
  fn_raw <- rbinom(length(x), 1, fnr)
  
  fp <- ifelse(fp_raw == 1 & y == lod, 1, 0)
  fn <- ifelse(fn_raw == 1 & y >= lod, 1, 0)
  
  y <- ifelse(fp == 1, lod + rexp(length(x), 1), y)
  y <- ifelse(fn == 1, lod, y)
  
  # return data frame of observed values
  data.frame(id, x, y, dp, wp, wr, tp, fp, fn)
}

df <- gendata(64)

mod <- cmdstan_model("1_code/stan/measurement_error.stan")


fit <- mod$sample(
  data = list(
    N = nrow(df),
    id = df$id,
    t = df$x,
    v_obs = df$y,
    lod = 1,
    mu = 1
  ),
  seed = 123,
  chains = 6,
  parallel_chains = 6,
  refresh = 100,
  iter_warmup = 1000,
  iter_sampling = 2000
)

# fit <- mod$optimize(
#   data = list(
#     N = nrow(df),
#     id = df$id,
#     t = df$x,
#     v_obs = df$y,
#     lod = 1,
#     fp = 0.1,
#     fn = 0.1
#   ),
#   iter = 10000
# )

ind_df <- df |> group_by(id) |> slice(1) |> ungroup()
ind_df$dp_hat <- E(as_draws_rvars(fit$draws(variables = "dp"))$dp)
ind_df$wp_hat <- E(as_draws_rvars(fit$draws(variables = "wp"))$wp)
ind_df$wr_hat <- E(as_draws_rvars(fit$draws(variables = "wr"))$wr)
ind_df$tp_hat <- E(as_draws_rvars(fit$draws(variables = "tp"))$tp)

ind_df$dp_lwr <- quantile2(as_draws_rvars(fit$draws(variables = "dp"))$dp, probs = 0.025)
ind_df$wp_lwr <- quantile2(as_draws_rvars(fit$draws(variables = "wp"))$wp, probs = 0.025)
ind_df$wr_lwr <- quantile2(as_draws_rvars(fit$draws(variables = "wr"))$wr, probs = 0.025)
ind_df$tp_lwr <- quantile2(as_draws_rvars(fit$draws(variables = "tp"))$tp, probs = 0.025)

ind_df$dp_upr <- quantile2(as_draws_rvars(fit$draws(variables = "dp"))$dp, probs = 0.975)
ind_df$wp_upr <- quantile2(as_draws_rvars(fit$draws(variables = "wp"))$wp, probs = 0.975)
ind_df$wr_upr <- quantile2(as_draws_rvars(fit$draws(variables = "wr"))$wr, probs = 0.975)
ind_df$tp_upr <- quantile2(as_draws_rvars(fit$draws(variables = "tp"))$tp, probs = 0.975)

ribbon_df <- 
  tibble(
    x = df$x,
    id = df$id,
    y = E(as_draws_rvars(fit$draws(variables = "v"))$v),
    y_lwr = quantile2(as_draws_rvars(fit$draws(variables = "v"))$v, probs = 0.025),
    y_upr = quantile2(as_draws_rvars(fit$draws(variables = "v"))$v, probs = 0.975)
  )

ribbon_df <- 
  mutate(
    ribbon_df,
    y_lwr = ifelse(y_lwr < 1, 1, y_lwr),
    y_upr = ifelse(y_upr < 1, 1, y_upr)
  )

ribbon_df_corners <- 
  tibble(
    x = ind_df$tp_hat,
    id = ind_df$id,
    y = E(as_draws_rvars(fit$draws(variables = "dp"))$dp),
    y_lwr = quantile2(as_draws_rvars(fit$draws(variables = "dp"))$dp, probs = 0.025),
    y_upr = quantile2(as_draws_rvars(fit$draws(variables = "dp"))$dp, probs = 0.975)
  )

ribbon_df <-
  bind_rows(
    ribbon_df,
    ribbon_df_corners
  )
# df$v_true <- pefun(df$x, df$tp, df$wp, df$wr, df$dp)
# df$v_true <- replace(df$v_true, df$v_true <= 1, 1)

ggplot(ind_df, aes(x = x, y = y, group = id)) +
  geom_point(data = df, aes(shape = factor(fn))) +
  geom_segment(aes(x=-Inf, y=1, xend=tp-wp, yend=1), linetype = "dashed") + 
  geom_segment(aes(x=tp-wp, y=1, xend=tp, yend=dp), linetype = "dashed") + 
  geom_segment(aes(x=tp, y=dp, xend=tp+wr, yend=1), linetype = "dashed") + 
  geom_segment(aes(x=tp+wr, y=1, xend=Inf, yend=1), linetype = "dashed") + 
  geom_ribbon(
    data = ribbon_df,
    aes(
      x = x,
      ymin = y_lwr,
      ymax = y_upr
    ),
    alpha = 0.2,
    fill = "black"
  ) +
  geom_segment(aes(x=-Inf, y=1, xend=tp_hat-wp_hat, yend=1)) + 
  geom_segment(aes(x=tp_hat-wp_hat, y=1, xend=tp_hat, yend=dp_hat)) + 
  geom_segment(aes(x=tp_hat, y=dp_hat, xend=tp_hat+wr_hat, yend=1)) + 
  geom_segment(aes(x=tp_hat+wr_hat, y=1, xend=Inf, yend=1)) + 
  facet_wrap(~id) +
  theme_minimal() +
  scale_x_continuous(limits=c(-6, 10)) +
  scale_shape_manual(values = c(1, 4)) +
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
  ) +
# theme(legend.position="bottom",
  #       panel.spacing=unit(1,"lines"),
  #       strip.background=element_blank(),
  #       strip.text=element_text(face="bold"),
  #       axis.text=element_text(size=8),axis.title=element_text(size=8),
  #       legend.title=element_text(size=6),legend.text=element_text(size=6),
  #       plot.tag=element_text(face="bold"),
  #       plot.background = element_rect(fill="white",color="white")) +
  labs(x = "t", y = bquote("log("*V[t]*")"))
