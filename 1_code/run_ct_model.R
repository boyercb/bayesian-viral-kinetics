# project 1 fit model to Ct values and determine

library(tidyverse)
library(cmdstanr)

ct_dat <- read_csv("0_data/ct_dat_refined.csv")
ct_mod <- cmdstan_model("1_code/ct_model_alt.stan")

# ids <- sample(ct_dat$InfectionEvent, )
# ct_dat <- filter(ct_dat, InfectionEvent %in% ids)

x_vars <- c(
  "age_[30,50)",
  "age_[50,100)",
  "recurrence",
  "delta",
  "omicron",
  "ba4ba5",
  "other",
  "vaccinated_boosted",
  "vaccinated_unboosted",
  "vaccinated_unreported",
  "unvaccinated_unreported",
  "unreported_primary"
)

refs <- c('age_[0,30)', 'alpha', 'unvaccinated_unboosted')

ct_dat <- ct_dat |>
  mutate(
    time = TestDateIndex,
    ct = CtT1,
    `age_[0,30)` = AgeGrp == "[0,30)",
    `age_[30,50)` = AgeGrp == "[30,50)",
    `age_[50,100)` = AgeGrp == "[50,100)",
    recurrence = InfNum >= 2,
    alpha = LineageBroad == "Alpha",
    delta = LineageBroad == "Delta",
    omicron = LineageBroad %in% c("BA.1", "BA.2"),
    ba4ba5 = LineageBroad %in% c("BA.4", "BA.5"),
    other = LineageBroad %in% c("None", "other"),
    unvaccinated_unboosted = VaccinationStatus == "Not Vaccinated" & BoosterStatus == "Not Boosted",
    vaccinated_boosted = VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Boosted",
    vaccinated_unboosted = VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Not Boosted",
    vaccinated_unreported = VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Not Reported",
    unvaccinated_unreported = VaccinationStatus == "Not Vaccinated" & BoosterStatus == "Not Reported",
    unreported_primary = VaccinationStatus == "Not Reported",
    drop = is.na(AgeGrp) | is.na(InfNum) | is.na(LineageBroad) | is.na(VaccinationStatus) | is.na(BoosterStatus),
    across(all_of(c(refs, x_vars, "drop")), as.numeric)
  ) |>
  filter(drop != 1) |>
  group_by(InfectionEvent) |>
  mutate(id = cur_group_id()) |>
  ungroup()

stan_data <- list(
  N = nrow(ct_dat),
  P = 12,
  lod = 40,
  kappa = 3.435,
  alpha = 1.418,
  sens = 0.99,
  fpmean = 1/log(10),
  id = ct_dat$id,
  time = ct_dat$time,
  ct = ct_dat$ct,
  x = ct_dat[, c(x_vars)],
  priorsd = 0.25,
  sigma_prior = c(0, 0.5),
  tp_prior = c(0, 2),
  dp_midpoint = 20,
  wp_midpoint = 5,
  wr_midpoint = 12
)

# ct_fit <- ct_mod$optimize(data = stan_data)

ct_fit <- ct_mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 100
)

dt <- ct_fit$summary()

