table1a <- 
  stacked_dat |>
  group_by(id) |>
  slice(1) |> 
  group_by(source) |>
  summarise(
    across(all_of(x_vars_w_refs), sum, .names = "{.col}_n"),
    across(all_of(x_vars_w_refs), function(x) sum(x) / n(), .names = "{.col}_p")
  ) |>
  pivot_longer(-source) |>
  mutate(
    type = str_extract(name, "_([np])$", group = 1),
    name = str_remove(name, "_[np]$")
  ) |>
  pivot_wider(
    names_from = c(type, source),
    values_from = value
  ) |>
  mutate(
    name = x_labels[name],
    across(starts_with("p_"), ~ .x * 100)
  )

kable(
  x = table1a,
  format = "latex",
  col.names = c("Characteristic", rep(c("N", "\\%"), 5)),
  align = paste0("l", paste(rep("c", ncol(table1) - 1), collapse = "")),
  digits = 1,
  booktabs = TRUE,
  linesep = "",
  escape = FALSE
) |>
  kable_styling() |>
  add_header_above(c(
    " " = 1,
    "NBA" = 2,
    "ATACCC" = 2,
    "UIUC" = 2,
    "HCT" = 2,
    "Legacy" = 2
  ))

table1b <- 
  stacked_dat |>
  group_by(source) |>
  summarise(
    across(all_of(miss_vars), sum, .names = "{.col}_n"),
  ) |>
  pivot_longer(-source) |>
  mutate(
    type = str_extract(name, "_([n])$", group = 1),
    name = str_remove(name, "_[n]$")
  ) |>
  pivot_wider(
    names_from = c(type, source),
    values_from = value
  ) |>
  mutate(
    name = case_when(
      name == "rna_exist" ~ "Ct values",
      name == "pfu_exist" ~ "Viral cultures",
      name == "lfd_exist" ~ "Rapid LFD tests",
      name == "sym_exist" ~ "Symptom diaries"
    ),
    across(starts_with("p_"), ~ .x * 100)
  )

kable(
  x = table1b,
  format = "latex",
  col.names = c("Lab values", "NBA",
                "ATACCC",
                "UIUC",
                "HCT",
                "Legacy"),
  align = paste0("l", paste(rep("c", ncol(table1b) - 1), collapse = "")),
  digits = 1,
  booktabs = TRUE,
  linesep = "",
  escape = FALSE
) |>
  kable_styling()

bind_rows(
  table1a,
  table1b
)

library(bayesplot)
p <- mcmc_trace(kinetics, pars = paste0("beta_dp_rna[", 1:13, "]"))
ggsave(
  filename = "3_figures/trace_plots.pdf",
  plot = p,
  device = "pdf",
  width = 9,
  height = 5
)
