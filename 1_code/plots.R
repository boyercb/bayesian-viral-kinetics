library(tidyverse)
library(cmdstanr)
library(progressr)
library(deSolve)
library(extrafont) 
font_import(pattern = "lmroman10-regular") 
handlers("progress")
loadfonts()
source("1_code/simstudy/bin.R")

params <- c(
  b = 0.000005788,
  d = 2.36,
  c = 974,
  p = 2.15
)

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
df$V <- ifelse(log(df$V) < 2.3, 2.3, log(df$V))

p <- ggplot(df, aes(x = time, y = V)) +
  geom_line(color = "red", linewidth = 1.2) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  )  +
  labs(
    y = NULL
  )

ggsave(
  filename = "3_figures/tcl_example_1.pdf",
  plot = p,
  device = "pdf",
  width = 8/2.54,
  height = 6 / 2.54,
  bg = "transparent"
)


# target cell limited model ode system
targetcell2 <- function(t, states, params) {
  with(as.list(c(states, params)), {
    dT <- -b * V * T
    dI <- b * V * T - d * I
    dV <- p * I - c * V
    dR <- q * p * I - e * R
    
    list(c(dT, dI, dV, dR))
  })
}

params <- c(
  b = 0.000005788,
  d = 2.36,
  c = 974,
  p = 2.15,
  q = 5,
  e = 100
)

states <- c(
  T = 8e8,
  I = 0,
  V = 1,
  R = 0
)

p <- ode(
  y = states,
  times = seq(0, 13, 0.01),
  func = targetcell2,
  parms = params
)

df <- as.data.frame(p)
df$V <- ifelse(log(df$V) < 2.3, 2.3, log(df$V))
df$R <- ifelse(log(df$R) < 5, 5, log(df$R))

p <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = V), color = "red", linewidth = 1.2) +
  geom_line(aes(y = R), color = "darkblue", linewidth = 1.2) +
  coord_cartesian(expand = FALSE, clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "transparent",
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent")
  ) +
  labs(
    y = NULL
  )

ggsave(
  filename = "3_figures/tcl_example_2.pdf",
  plot = p,
  device = "pdf",
  width = 8 / 2.54,
  height = 6 / 2.54,
  bg = "transparent"
)

