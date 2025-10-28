
x_vars <- c(
  "age_[30,50)",
  "age_[50,100)",
  "recurrence",
  "alpha",
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

miss_vars <- c(
  "rna_exist",
  "pfu_exist",
  "lfd_exist",
  "sym_exist"
)

refs <- c('age_[0,30)', 'first_infection', 'prealpha', 'unvaccinated_unboosted')

x_vars_w_refs <- c(
  refs[1],
  x_vars[1:2],
  refs[2],
  x_vars[3],
  refs[3],
  x_vars[4:8],
  refs[4],
  x_vars[9:13]
)

x_labels <- c(
  "age_[0,30)" = "Age: [0,30)",
  "age_[30,50)" = "Age: [30,50)",
  "age_[50,100)" = "Age: [50,100)",
  "first_infection" = "Recurrence: No",
  "recurrence" = "Recurrence: Yes",
  "prealpha" = "Variant: Pre-Alpha",
  "alpha" = "Variant: Alpha",
  "delta" = "Variant: Delta",
  "omicron" = "Variant: Omicron",
  "ba4ba5" = "Variant: BA.4/BA.5",
  "other" = "Variant: Other",
  "unvaccinated_unboosted" = "History: Unvaccinated",
  "vaccinated_boosted" = "History: Vaccinated boosted",
  "vaccinated_unboosted" = "History: Vaccinated unboosted",
  "vaccinated_unreported" = "History: Vaccinated unreported",
  "unvaccinated_unreported" = "History: Unreported",
  "unreported_primary" = "History: Boosted unreported primary" 
)

# clean NBA data ----------------------------------------------------------

nba_dat <- nba_dat |>
  mutate(
    person = PersonID,
    time = TestDateIndex,
    ct = CtT1,
    rna = ct_to_rna(ct, type = "nba"),
    pfu = NA,
    lfd = NA,
    sym = NA,
    sym_ever = NA,
    sym_onset = NA,
    pfu_type = 1,
    `age_[0,30)` = AgeGrp == "[0,30)",
    `age_[30,50)` = AgeGrp == "[30,50)",
    `age_[50,100)` = AgeGrp == "[50,100)",
    recurrence = InfNum >= 2,
    first_infection = 1 - recurrence,
    prealpha = LineageBroad == "None",
    alpha = LineageBroad == "Alpha",
    delta = LineageBroad == "Delta",
    omicron = LineageBroad %in% c("BA.1", "BA.2"),
    ba4ba5 = LineageBroad %in% c("BA.4", "BA.5"),
    other = LineageBroad %in% c("Other"),
    unvaccinated_unboosted = 
      VaccinationStatus == "Not Vaccinated" & BoosterStatus == "Not Boosted",
    vaccinated_boosted = 
      VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Boosted",
    vaccinated_unboosted = 
      VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Not Boosted",
    vaccinated_unreported = 
      VaccinationStatus == "Fully Vaccinated" & BoosterStatus == "Not Reported",
    unvaccinated_unreported = 
      VaccinationStatus == "Not Vaccinated" & BoosterStatus == "Not Reported",
    unreported_primary = VaccinationStatus == "Not Reported",
    rna_exist = !is.na(CtT1),
    pfu_exist = FALSE,
    lfd_exist = FALSE,
    sym_exist = FALSE,
    drop = is.na(AgeGrp) | is.na(InfNum) | is.na(LineageBroad) | 
      is.na(VaccinationStatus) | is.na(BoosterStatus),
    across(all_of(c(
      refs, x_vars, "drop", "rna_exist", "pfu_exist", "lfd_exist", "sym_exist"
    )), as.numeric)
  ) |>
  filter(drop != 1) |>
  filter(InfNum <= 2) |>
  group_by(InfectionEvent) |>
  mutate(id = cur_group_id(), pid = id) |>
  ungroup()


# clean ATACCC data --------------------------------------------------------

ataccc_dat <- 
  ataccc_dat |>
  mutate(
    id = participant,
    pid = participant,
    time = days_since_peak,
    rna = log(copy),
    ct = 40.0 - (rna - 3.435) * 1.418,
    ct = replace(ct, ct > 40, 40),
    rna = replace(rna, rna == 0, 3.435),
    pfu = replace(log(pfu), log(pfu) == 0, 2.3),
    pfu_type = 1,
    lfd = ifelse(LFD > 0, 1, 0),
    prealpha = WGS == "Pre-Alpha",
    alpha = WGS == "Alpha",
    delta = WGS == "Delta",
    vaccinated_unboosted = vaccinated,
    rna_exist = !is.na(copy),
    pfu_exist = !is.na(pfu),
    lfd_exist = !is.na(LFD),
    sym_exist = 0,
    recurrence = 0,
    first_infection = 1,
    omicron = 0,
    ba4ba5 = 0,
    other = 0,
    unvaccinated_unboosted = 1,
    vaccinated_boosted = 0,
    vaccinated_unreported = 0,
    unvaccinated_unreported = 0,
    unreported_primary = 0,
    `age_[0,30)` = 0,
    `age_[30,50)` = 1,
    `age_[50,100)` = 0,
    across(c(
      alpha,
      delta,
      vaccinated_unboosted,
      rna_exist,
      pfu_exist,
      lfd_exist,
      sym_exist,
    ),
    as.numeric)
  )

# add symptom data
ataccc_dat <- left_join(ataccc_dat, ataccc_sym, by = "participant")
ataccc_dat <-
  ataccc_dat |>
  group_by(id) |>
  mutate(
    sym_onset = ifelse(
      symptom_onset_relative_to_peak_vl == "Asymptomatic",
      max(days_since_peak),
      as.numeric(symptom_onset_relative_to_peak_vl)
    ),
    sym_ever = ifelse(symptom_onset_relative_to_peak_vl != "Asymptomatic", 1, 0),
    symptom_onset_relative_to_peak_vl = as.numeric(
      replace(symptom_onset_relative_to_peak_vl,
              symptom_onset_relative_to_peak_vl == "Asymptomatic",
              NA)
    ),
    sym = ifelse(
      is.na(symptom_onset_relative_to_peak_vl),
      NA,
      ifelse(
        days_since_peak < symptom_onset_relative_to_peak_vl,
        0,
        1
      )
    ),
    sym_exist = as.numeric(!is.na(sym_onset))
  ) |>
  ungroup()

# drop observations where everything is missing
ataccc_dat <-
  filter(ataccc_dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)

exclude_pfu_raw <- c(12, 18, 23, 25, 41, 56)
ataccc_dat <- filter(ataccc_dat, !id %in% exclude_pfu_raw)


# clean UIUC data ---------------------------------------------------------

uiuc_dat <- 
  uiuc_dat |>
  mutate(
    # id = Ind,
    pid = as.numeric(str_remove(Ind, " \\*+")),
    time = Time,
    ct = Saliva_Ct,
    cn = replace(Nasal_CN, Nasal_CN == 48, 45),
    rna_ct = ct_to_rna(ct, type = "uiuc-ct"),
    # rna = replace(rna, rna == 0, 3.435),
    rna_cn = ct_to_rna(cn, type = "uiuc-cn"),
    # rna_cn = replace(rna_cn, rna_cn == 0, 3.435),
    rna = ifelse(
      !is.na(rna_cn) & !is.na(rna_ct), 
      log(exp(rna_ct) / 2 + exp(rna_cn) / 2),
      ifelse(
        is.na(rna_cn),
        rna_ct,
        rna_cn
      )),
    pfu = replace(Virus_pos_days, Virus_pos_days == 5.9, 6),
    # pfu = -1 * (pfu - 7), #ifelse(Virus_pos_days != 5.9, 1, 0),
    pfu_type = 2,
    lfd = ifelse(Antigen == "Pos", 1, 0),
    recurrence = 0,
    first_infection = 1,
    prealpha = ifelse(!Lineage %in% c("B.1.1.7","P.1"), 1, 0),
    alpha = ifelse(Lineage == "B.1.1.7", 1, 0),
    delta = 0,
    omicron = 0,
    ba4ba5 = 0,
    other = ifelse(Lineage == "P.1", 1, 0),
    unvaccinated_unboosted = 1,
    vaccinated_unboosted = 0,
    vaccinated_boosted = 0,
    vaccinated_unreported = 0,
    unvaccinated_unreported = 0,
    unreported_primary = 0,
    `age_[0,30)` = ifelse(Age < 30, 1, 0),
    `age_[30,50)` = ifelse(Age >= 30 & Age < 50, 1, 0),
    `age_[50,100)` = ifelse(Age >= 50, 1, 0),
    ct_exist = as.numeric(!is.na(ct)),
    cn_exist = as.numeric(!is.na(cn)),
    rna_exist = as.numeric(!is.na(rna)),
    pfu_exist = as.numeric(!is.na(pfu)),
    lfd_exist = as.numeric(!is.na(lfd)),
    sym_exist = 0
  )

uiuc_dat <- uiuc_dat |>
  group_by(Ind) |>
  mutate(id = cur_group_id()) |>
  ungroup()

uiuc_dat <- 
  uiuc_dat |>
  group_by(pid) |>
  mutate(
    peak_day = sum(ifelse(rna == max(rna, na.rm = TRUE), time, NA), na.rm = TRUE),
    days_since_peak = time - peak_day,
    time = days_since_peak
  ) |>
  ungroup() 

uiuc_dat <- 
  uiuc_dat |>
  group_by(Ind) |>
  arrange(Ind, Time) |>
  mutate(
    obs = row_number()
  )

uiuc_sym <- 
  uiuc_sym |>
  group_by(Ind) |>
  arrange(Ind, Time) |>
  mutate(
    obs = row_number()
  )

uiuc_dat <- left_join(uiuc_dat, uiuc_sym, by = c("Ind", "Time", "obs"))

uiuc_dat <-
  uiuc_dat |>
  group_by(id) |>
  mutate(
    sym = ifelse(
      days_since_peak < sym_onset,
      0,
      1
    ),
    sym_exist = as.numeric(!is.na(sym_onset))
  ) |> 
  ungroup()

uiuc_dat <- filter(uiuc_dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)


# clean human challenge trial data ----------------------------------------

hct_dat <-
  hct_dat |>
  mutate(
    # id = pid,
    # time = Time,
    rna_ct = log(replace(qpcr_throat, qpcr_throat == 0, 1)),
    rna_cn = log(replace(qpcr_nose, qpcr_nose == 0, 1)),
    rna = ifelse(
      !is.na(rna_cn) & !is.na(rna_ct), 
      log(exp(rna_ct) / 2 + exp(rna_cn) / 2),
      ifelse(
        is.na(rna_cn),
        rna_ct,
        rna_cn
      )),
    pfu = log(replace(throat, throat == 0, 1)),
    pfu_type = 1,
    ffa_t = log(replace(ffa_throat, ffa_throat == 0, 1)),
    ffa_n = log(replace(ffa_nose, ffa_nose == 0, 1)),
    ffa = ifelse(
      !is.na(ffa_n) & !is.na(ffa_t), 
      log(exp(ffa_t) / 2 + exp(ffa_n) / 2),
      ifelse(
        is.na(ffa_n),
        ffa_t,
        ffa_n
      )),
    lfd = lfd,
    recurrence = 0,
    first_infection = 1,
    prealpha = 1,
    alpha = 0,
    delta = 0,
    omicron = 0,
    ba4ba5 = 0,
    other = 0,
    unvaccinated_unboosted = 1,
    vaccinated_unboosted = 0,
    vaccinated_boosted = 0,
    vaccinated_unreported = 0,
    unvaccinated_unreported = 0,
    unreported_primary = 0,
    `age_[0,30)` = 1,
    `age_[30,50)` = 0,
    `age_[50,100)` = 0,
    ct_exist = as.numeric(!is.na(rna_ct)),
    cn_exist = as.numeric(!is.na(rna_cn)),
    rna_exist = as.numeric(!is.na(rna)),
    pfu_exist = as.numeric(!is.na(pfu)),
    lfd_exist = as.numeric(!is.na(lfd))
  )

hct_dat <- 
  hct_dat |>
  group_by(pid) |>
  mutate(
    peak_day = sum(ifelse(rna == max(rna, na.rm = TRUE), day, NA), na.rm = TRUE),
    days_since_peak = day - peak_day,
    time = days_since_peak
  ) |>
  ungroup() 

hct_sym <- hct_sym |>
  select(-ends_with("_desc")) |>
  filter(
    !collection_method %in% c(
      "Paper diary card completed",
      "Not yet admitted",
      "Subject discharged"
    )
  )

hct_sym <- 
  hct_sym |>
  group_by(id, day) |>
  arrange(id, day, timepoint) |>
  summarise(
    across(runny_nose:score, \(x) max(x, na.rm = TRUE)),
    .groups = "drop"
  )
  
hct_sym <- 
  hct_sym |>
  group_by(id) |>
  arrange(id, day) |>
  mutate(
    sym = ifelse(cumsum(score) > 0, 1, 0),
    sym_ever = max(sym, na.rm = TRUE),
    x = sym * day + (1 - sym) * max(day),
    sym_onset = min(x, na.rm = TRUE)
  ) |>
  ungroup() 

hct_dat <- left_join(hct_dat, hct_sym, by = c("id", "day"))

hct_dat <- 
  hct_dat |>
  group_by(id) |>
  mutate(
    sym_exist = as.numeric(!is.na(sum(score, na.rm = TRUE))),
    sym_ever = first(sym_ever),
    sym_onset = first(sym_onset)
  ) |>
  ungroup()

hct_dat <- filter(hct_dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)


# clean Legacy data -------------------------------------------------------

legacy_dat <-
  legacy_dat |>
  mutate(
    id = id,
    pid = id,
    # time = Time,
    ct = ct_value,
    rna = ct_to_rna(ct_value, type = "nba") - 2,
    pfu = NA,
    lfd = NA,
    pfu_type = 1,
    sym = case_when(
      symptoms == "symptomatic" & swab_date < symptom_onset_date ~ 0,
      symptoms == "symptomatic" & swab_date >= symptom_onset_date ~ 1,
      symptoms == "asymptomatic" ~ 0
    ),
    sym_ever = as.numeric(symptoms == "symptomatic"),
    `age_[0,30)` = age_group == "20-34",
    `age_[30,50)` = age_group == "35-49",
    `age_[50,100)` = age_group == "50+",
    recurrence = infection_id >= 2,
    first_infection = 1 - recurrence,
    prealpha = VOC == "None",
    alpha = VOC == "Alpha",
    delta = VOC == "Delta",
    omicron = VOC %in% c("Omicron (BA.1)", "Omicron (BA.2)"),
    ba4ba5 = VOC %in% c("Omicron (BA.4)", "Omicron (BA.5)"),
    other = VOC %in% c("Other"),
    unvaccinated_unboosted = 0,
    vaccinated_unboosted = no_vaccines == 2,
    vaccinated_boosted = no_vaccines == 3,
    vaccinated_unreported = 0,
    unvaccinated_unreported = 0,
    unreported_primary = 0,
    rna_exist = as.numeric(!is.na(rna)),
    pfu_exist = 0,
    lfd_exist = 0,
    sym_exist = as.numeric(!is.na(sym)),
    across(c(
      sym,
      alpha,
      delta,
      omicron,
      ba4ba5,
      other, 
      vaccinated_unboosted,
      vaccinated_boosted,
      rna_exist,
      pfu_exist,
      lfd_exist,
      sym_exist,
    ),
    as.numeric)
  )

legacy_dat <- 
  legacy_dat |>
  group_by(pid) |>
  mutate(
    peak_day = sum(ifelse(rna == max(rna, na.rm = TRUE), t, NA), na.rm = TRUE),
    peak_date = max(ifelse(rna == max(rna, na.rm = TRUE), swab_date, NA), na.rm = TRUE),
    days_since_peak = t - peak_day,
    sym_onset = as.numeric(symptom_onset_date - peak_date),
    sym_onset = replace(sym_onset, sym_ever == 0, max(days_since_peak)),
    time = days_since_peak
  ) |>
  ungroup() 

legacy_dat <-
  filter(legacy_dat,
         (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0 & ct_type == "ct_value")

add_onset_obs <- 
  legacy_dat |>
  group_by(pid) |>
  slice(1) |>
  select(
    pid,
    id,
    time,
    rna,
    pfu_type,
    sym,
    sym_ever,
    sym_onset,
    sym_exist,
    symptom_onset_date,
    swab_date,
    all_of(x_vars_w_refs),
    all_of(miss_vars)
  ) |>
  filter(!is.na(sym_onset) & symptom_onset_date != swab_date) |>
  mutate(
    time = sym_onset,
    sym = 1,
    swab_date = NA,
    rna = NA,
    rna_exist = 0,
  ) |>
  ungroup()

legacy_dat <- 
  bind_rows(legacy_dat, add_onset_obs) |>
  arrange(pid, time) 

  
# create stacked dataset --------------------------------------------------

keep_vars <-
  c(
    "pid",
    "id",
    "time",
    "rna",
    "pfu",
    "pfu_type",
    "lfd",
    "sym",
    "sym_onset",
    "sym_ever",
    x_vars_w_refs,
    miss_vars
  )

stacked_dat <- bind_rows(
  select(nba_dat, all_of(keep_vars)),
  select(ataccc_dat, all_of(keep_vars)),
  select(uiuc_dat, all_of(keep_vars)),
  select(hct_dat, all_of(keep_vars)),
  select(legacy_dat, all_of(keep_vars)),
  .id = "source"
)

stacked_dat <- 
  stacked_dat |> 
  group_by(source, id) |>
  mutate(id = cur_group_id()) |>
  ungroup() 
  
stacked_dat <-
  stacked_dat |>
  mutate(
    rna = replace(rna, rna_exist == 0, 0),
    pfu = replace(pfu, pfu_exist == 0, 0),
    lfd = replace(lfd, lfd_exist == 0, 0),
    sym_onset = replace(sym_onset, sym_exist == 0, 99),
    sym_ever = replace(sym_ever, sym_exist == 0, 0)
  )


# create Stan list --------------------------------------------------------

stan_data <- list(
  K = 5,
  M = c(
    length(unique(nba_dat$id)),
    length(unique(ataccc_dat$id)),
    length(unique(uiuc_dat$id)),
    length(unique(hct_dat$id)),
    length(unique(legacy_dat$id))
  ),
  N = c(
    nrow(nba_dat),
    nrow(ataccc_dat),
    nrow(uiuc_dat),
    nrow(hct_dat),
    nrow(legacy_dat)
  ),
  lod_rna = c(
    ct_to_rna(40, type = "nba") + 0.01,
    ct_to_rna(40, type = "ata") + 0.01,
    ct_to_rna(47, type = "uiuc-ct") + 0.01,
    log(1000) + 0.01,
    ct_to_rna(40, type = "nba") + 0.01
  ),
  lod_pfu = c(
    0,
    2.3,
    2.3,
    log(5),
    0
  ),
  fp_mean = c(
    1 / (ct_to_rna(40 - 1/log(10), type = "nba") - ct_to_rna(40, type = "nba")),
    1 / (ct_to_rna(40 - 1/log(10), type = "ata") - ct_to_rna(40, type = "ata")),
    1 / (ct_to_rna(47 - 1/log(10), type = "uiuc-ct") - ct_to_rna(47, type = "uiuc-ct")),
    1 / 0.3,
    1 / (ct_to_rna(40 - 1/log(10), type = "nba") - ct_to_rna(40, type = "nba"))
  ),
  id = stacked_dat$id,
  time = stacked_dat$time,
  rna = stacked_dat$rna,
  pfu = stacked_dat$pfu,
  lfd = stacked_dat$lfd,
  sym_onset = stacked_dat$sym_onset,
  sym_ever = stacked_dat$sym_ever,
  source = as.numeric(stacked_dat$source),
  pfu_type = stacked_dat$pfu_type,
  P = 13,
  x = stacked_dat |>
    group_by(id) |>
    slice(1) |> 
    ungroup() |>
    select(all_of(x_vars)),
  rna_exist = stacked_dat$rna_exist,
  pfu_exist = stacked_dat$pfu_exist,
  lfd_exist = stacked_dat$lfd_exist,
  sym_exist = stacked_dat$sym_exist,
  ind_effects = 1,
  test_error = 1,
  adj_pfu = 0,
  adj_rna = 1,
  adj_lfd = 0,
  source_pfu = 0,
  source_rna = 0,
  source_lfd = 0,
  source_sym = 0,
  prior_fp = 0.05,
  prior_fn = 0.01,
  prior_dp_mean = 17,
  prior_dp_cv = 0.7,
  prior_wp_mean = 8,
  prior_wp_cv = 0.7,
  prior_wr_mean = 15,
  prior_wr_cv = 0.7,
  prior_sigma_sd = 5,
  prior_beta_sd = 1,
  prior_i_sd = 1,
  prior_k_sd = 1,
  prior_lfd_mean = 0.01,
  # dp_max = ct_to_rna(0, type = "nba"),
  # wp_max = 30,
  # wr_max = 30,
  # sigma_max = 10,
  additional = TRUE
)