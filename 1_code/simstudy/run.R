
d <- get_truth()

sims <- lapply(results, function(x, d) {
  rows <- lapply(1:length(x$fits), function(y) {
    if (x$fits[[y]]$return_codes() == 0) {
      params <- c(
        b = x$fits[[y]]$mle("b"),
        d = x$fits[[y]]$mle("d"),
        c = x$fits[[y]]$mle("c"),
        p = x$fits[[y]]$mle("p")
      )
      names(params) <- c("b", "d", "c", "p")
      
      states <- c(
        T = 8e8,
        I = 0,
        V = 1
      )
      
      if (x$fits[[y]]$mle("t0") < 0) {
        times <- seq(x$fits[[y]]$mle("t0"), 13, 0.01)
      } else {
        times <- seq(x$fits[[y]]$mle("t0"), 13, 0.01)
      }
      
      p <- ode(
        y = states,
        times = times,
        func = targetcell,
        parms = params
      )
      
      df <- as.data.frame(p)
      df$id <- y
      df$V <- ifelse(df$V < exp(2.3), exp(2.3), df$V)
      
      pred_corners_pe <- calc_corners(
        tp = x$fits[[y]]$mle("tp_pe"), 
        wp = x$fits[[y]]$mle("wp_pe"), 
        wr = x$fits[[y]]$mle("wr_pe"), 
        dp = x$fits[[y]]$mle("dp_pe"), 
        lod = 2.3
      )
      
      pred_corners_sm <- calc_corners(
        tp = x$fits[[y]]$mle("tp_sm"), 
        wp = x$fits[[y]]$mle("wp_sm"), 
        wr = x$fits[[y]]$mle("wr_sm"), 
        dp = x$fits[[y]]$mle("dp_sm"), 
        lod = 2.3
      )
      
      tibble(
        id = y,
        tp_pe = x$fits[[y]]$mle("tp_pe"),
        dp_pe = x$fits[[y]]$mle("dp_pe"),
        wp_pe = pred_corners_pe$tp - pred_corners_pe$to,
        wr_pe = pred_corners_pe$tr - pred_corners_pe$tp,
        tp_sm = x$fits[[y]]$mle("tp_sm"),
        dp_sm = x$fits[[y]]$mle("dp_sm"),
        wp_sm = pred_corners_sm$tp - pred_corners_sm$to,
        wr_sm = pred_corners_sm$tr - pred_corners_sm$tp,
        dp_tc = max(log(df$V)),
        tp_tc = max(I(log(df$V) == max(log(df$V))) * df$time),
        wp_tc = tp_tc - min(df$time * I(log(df$V) > 2.3) + I(log(df$V) <= 2.3) * 99),
        wr_tc = max(df$time * I(log(df$V) > 2.3)) - tp_tc,
        dp_truth = d$dp_truth,
        tp_truth = d$tp_truth,
        wp_truth = d$wp_truth,
        wr_truth = d$wr_truth,
        row.names = NULL
      )
    } else {
      tibble(
        id = y,
        tp_pe = NA,
        dp_pe = NA,
        wp_pe = NA,
        wr_pe = NA,
        tp_sm = NA,
        dp_sm = NA,
        wp_sm = NA,
        wr_sm = NA,
        dp_tc = NA,
        tp_tc = NA,
        wp_tc = NA,
        wr_tc = NA
      )
    }
    
    
  })
  bind_rows(rows)
}, d = d)


sims <- bind_rows(sims, .id = "scenario")

simsum <- sims |>
  group_by(scenario) |>
  summarise(
    bias_tp_pe = mean(tp_pe - tp_truth, na.rm = TRUE),
    bias_tp_sm = mean(tp_sm - tp_truth, na.rm = TRUE),
    bias_tp_tc = mean(tp_tc - tp_truth, na.rm = TRUE),
    bias_dp_pe = mean(dp_pe - dp_truth, na.rm = TRUE),
    bias_dp_sm = mean(dp_sm - dp_truth, na.rm = TRUE),
    bias_dp_tc = mean(dp_tc - dp_truth, na.rm = TRUE),
    bias_wp_pe = mean(wp_pe - wp_truth, na.rm = TRUE),
    bias_wp_sm = mean(wp_sm - wp_truth, na.rm = TRUE),
    bias_wp_tc = mean(wp_tc - wp_truth, na.rm = TRUE),
    bias_wr_pe = mean(wr_pe - wr_truth, na.rm = TRUE),
    bias_wr_sm = mean(wr_sm - wr_truth, na.rm = TRUE),
    bias_wr_tc = mean(wr_tc - wr_truth, na.rm = TRUE),
    bias_tt_pe = mean((wp_pe + wr_pe) - (wp_truth + wr_truth), na.rm = TRUE),
    bias_tt_tc = mean((wp_tc + wr_tc) - (wp_truth + wr_truth), na.rm = TRUE),
    bias_tt_sm = mean((wp_sm + wr_sm) - (wp_truth + wr_truth), na.rm = TRUE)
  ) |>
  pivot_longer(-scenario) |>
  separate(name, c("name", "type", "param"), sep = "_") |>
  pivot_wider(names_from = scenario, values_from = value)


ggplot(sims |> pivot_longer(
  -c(scenario, id),
  names_pattern = "(.*)_(.*)",
  names_to = c("var", "type")
),
aes(x =type, y= value)) +
  facet_wrap(scenario ~ var, scales = "free") +
  geom_violin(aes(fill = type), draw_quantiles = c(0.5)) +
  geom_boxplot(outlier.shape = NA, ) + 
  theme_bw()
