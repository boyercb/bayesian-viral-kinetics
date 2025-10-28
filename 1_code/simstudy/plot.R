
pred <- data.frame(
  id = 1:100,
  tp = fit$mle("tp"),
  dp = fit$mle("dp"),
  wp = fit$mle("wp"),
  wr = fit$mle("wr")
)

pred_tc <- map_dfr(1:100, function(x) {
  params <- c(
    b = fit$mle("b")[x],
    d = fit$mle("d")[x],
    c = fit$mle("c")[x],
    p = fit$mle("p")[x]
  )
  names(params) <- c("b", "d", "c", "p")
  states <- c(
    T = 8e8,
    I = 0,
    V = 1
  )
  
  if (fit$mle("t0")[x] < 0) {
    times <- c(fit$mle("t0")[x], seq(0, 13, 0.01))
  } else {
    times <- seq(fit$mle("t0")[x], 13, 0.01)
  }
  p <- ode(
    y = states,
    times = times,
    func = targetcell,
    parms = params
  )
  
  df <- as.data.frame(p)
  df$id <- x
  df$V <- ifelse(df$V < exp(2.3), exp(2.3), df$V)
  return(df)
})

pred_truth <- map_dfr(1:100, function(x) {
  params <- c(
    b = 0.000005788,
    d = 2.36,
    c = 974,
    p = 2.15
  )
  
  names(params) <- c("b", "d", "c", "p")
  states <- c(
    T = 8e8,
    I = 0,
    V = 1
  )
  
  p <- ode(
    y = states,
    times = seq(0, 13, 0.01),
    func = targetcell,
    parms = params
  )
  
  df <- as.data.frame(p)
  df$id <- x
  df$V <- ifelse(df$V < exp(2.3), exp(2.3), df$V)
  return(df)
})

ggplot(pred, aes(group = id)) +
  facet_wrap(~id) +
  geom_point(
    data = d,
    aes(x = time, y = V_obs),
    color = "#1F78B4",
    alpha = 0.4,
    shape = 16
  ) +
  geom_line(data = pred_truth, aes(x = time, y = log(V)), color = "#FF7F00") +
  geom_line(data = pred_tc, aes(x = time, y = log(V)), color = "#1F78B4", linetype = "dashed") +
  geom_segment(aes(
    x = -Inf,
    y = 2.3,
    xend = tp - wp,
    yend = 2.3
  ),
  color = "#1F78B4") +
  geom_segment(aes(
    x = pmax(tp - wp, -10),
    y = 2.3,
    xend = tp,
    yend = dp
  ), color = "#1F78B4") +
  geom_segment(aes(
    x = tp,
    y = dp,
    xend = pmin(tp + wr, 20),
    yend = 2.3
  ), color = "#1F78B4") +
  geom_segment(aes(
    x = tp + wr,
    y = 2.3,
    xend = Inf,
    yend = 2.3
  ),
  color = "#1F78B4") +
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

calc_corners <- function(tp, wp, wr, dp, lod) {
  to <- lod * wp / dp + (tp - wp)
  tr <- (dp - lod) * wr / dp + tp 
  
  return(data.frame(to = to, tp = tp, tr = tr))
}

pred_corners <- calc_corners(pred$tp, pred$wp, pred$wr, pred$dp, 2.3)

pred_result <- 
  pred |>
  mutate(
    wp = pred_corners$tp - pred_corners$to,
    wr = pred_corners$tr - pred_corners$tp
  )

pred_tc_result <- pred_tc |>
  group_by(id) |>
  summarise(
    dp = max(log(V)),
    tp = max(I(log(V) == max(log(V))) * time),
    wp = tp - min(time * I(log(V) > 2.3) + I(log(V) <= 2.3) * 99),
    wr = max(time * I(log(V) > 2.3)) - tp,
  )

pred_truth_result <- pred_truth |>
  group_by(id) |>
  summarise(
    dp = max(log(V)),
    tp = max(I(log(V) == max(log(V))) * time),
    wp = tp - min(time * I(log(V) > 2.3) + I(log(V) <= 2.3) * 99),
    wr = max(time * I(log(V) > 2.3)) - tp,
  )

left_join(
  pred_result,
  pred_tc_result,
  by = "id",
  suffix = c("_pe", "_tc")
) |>
  left_join(
    pred_truth_result,
    by = "id",
  ) |>
  summarise(
    tp_pe_bias = mean(tp_pe - tp),
    dp_pe_bias = mean(dp_pe - dp),
    wp_pe_bias = mean(wp_pe - wp),
    wr_pe_bias = mean(wr_pe - wr),
    tp_tc_bias = mean(tp_tc - tp),
    dp_tc_bias = mean(dp_tc - dp),
    wp_tc_bias = mean(wp_tc - wp),
    wr_tc_bias = mean(wr_tc - wr)
  )