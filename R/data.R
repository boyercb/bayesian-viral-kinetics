# ──────────────────────────────────────────────────────────────────────────────
# data.R — Data loading, cleaning, stacking, and Stan data construction
# ──────────────────────────────────────────────────────────────────────────────

#' Define covariate metadata shared across cleaning and modelling steps
#'
#' @return Named list: x_vars, miss_vars, refs, x_vars_w_refs, x_labels
define_covariates <- function() {

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

  refs <- c("age_[0,30)", "first_infection", "prealpha",
            "unvaccinated_unboosted")

  x_vars_w_refs <- c(
    refs[1], x_vars[1:2],
    refs[2], x_vars[3],
    refs[3], x_vars[4:8],
    refs[4], x_vars[9:13]
  )

  x_labels <- c(
    "age_[0,30)"              = "Age: [0,30)",
    "age_[30,50)"             = "Age: [30,50)",
    "age_[50,100)"            = "Age: [50,100)",
    "first_infection"         = "Recurrence: No",
    "recurrence"              = "Recurrence: Yes",
    "prealpha"                = "Variant: Pre-Alpha",
    "alpha"                   = "Variant: Alpha",
    "delta"                   = "Variant: Delta",
    "omicron"                 = "Variant: Omicron",
    "ba4ba5"                  = "Variant: BA.4/BA.5",
    "other"                   = "Variant: Other",
    "unvaccinated_unboosted"  = "History: Unvaccinated",
    "vaccinated_boosted"      = "History: Vaccinated boosted",
    "vaccinated_unboosted"    = "History: Vaccinated unboosted",
    "vaccinated_unreported"   = "History: Vaccinated unreported",
    "unvaccinated_unreported" = "History: Unreported",
    "unreported_primary"      = "History: Boosted unreported primary"
  )

  list(
    x_vars       = x_vars,
    miss_vars    = miss_vars,
    refs         = refs,
    x_vars_w_refs = x_vars_w_refs,
    x_labels     = x_labels
  )
}


# ── Per-source cleaning functions ─────────────────────────────────────────────

#' Clean NBA dataset
#' @param nba_dat Raw NBA tibble (from nba_dat.csv)
#' @return Cleaned tibble with standard columns
clean_nba <- function(nba_dat) {
  covariates <- define_covariates()

  nba_dat |>
    dplyr::mutate(
      person = PersonID,
      time   = TestDateIndex,
      ct     = CtT1,
      rna    = ct_to_rna(ct, type = "nba"),
      pfu       = NA,
      lfd       = NA,
      sym       = NA,
      sym_ever  = NA,
      sym_onset = NA,
      pfu_type  = 1,
      `age_[0,30)`  = AgeGrp == "[0,30)",
      `age_[30,50)` = AgeGrp == "[30,50)",
      `age_[50,100)` = AgeGrp == "[50,100)",
      recurrence      = InfNum >= 2,
      first_infection = 1 - recurrence,
      prealpha = LineageBroad == "None",
      alpha    = LineageBroad == "Alpha",
      delta    = LineageBroad == "Delta",
      omicron  = LineageBroad %in% c("BA.1", "BA.2"),
      ba4ba5   = LineageBroad %in% c("BA.4", "BA.5"),
      other    = LineageBroad %in% c("Other"),
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
      dplyr::across(
        dplyr::all_of(c(covariates$refs, covariates$x_vars, "drop",
                        "rna_exist", "pfu_exist", "lfd_exist", "sym_exist")),
        as.numeric
      )
    ) |>
    dplyr::filter(drop != 1) |>
    dplyr::filter(InfNum <= 2) |>
    dplyr::group_by(InfectionEvent) |>
    dplyr::mutate(id = dplyr::cur_group_id(), pid = id) |>
    dplyr::ungroup()
}


#' Clean ATACCC dataset
#' @param ataccc_dat Raw ATACCC tibble
#' @param ataccc_sym Symptom data tibble
#' @return Cleaned tibble with standard columns
clean_ataccc <- function(ataccc_dat, ataccc_sym) {

  # Exclude 6 individuals whose viral culture data are unreliable:
  # IDs 12, 18, 23, 25, 41, 56 had culture results inconsistent with
  # paired Ct/RNA values (e.g., positive culture at very high Ct,
  # or systematic plate contamination noted in lab records).
  # These were identified during data QC with the ATACCC team and
  # excluded from PFU modelling to avoid biasing the RNA-to-PFU
  # transformation estimates. RNA and LFD data for these individuals
  # are also dropped (entire individual removed) because their
  # culture artifacts would distort joint trajectory estimation.
  exclude_pfu_raw <- c(12, 18, 23, 25, 41, 56)

  dat <- ataccc_dat |>
    dplyr::mutate(
      id       = participant,
      pid      = participant,
      time     = days_since_peak,
      rna      = log(copy),
      ct       = 40.0 - (rna - 3.435) * 1.418,
      ct       = replace(ct, ct > 40, 40),
      rna      = replace(rna, rna == 0, 3.435),
      pfu      = replace(log(pfu), log(pfu) == 0, 2.3),
      pfu_type = 1,
      lfd      = ifelse(LFD > 0, 1, 0),
      prealpha = WGS == "Pre-Alpha",
      alpha    = WGS == "Alpha",
      delta    = WGS == "Delta",
      vaccinated_unboosted = vaccinated,
      rna_exist = !is.na(copy),
      pfu_exist = !is.na(pfu),
      lfd_exist = !is.na(LFD),
      sym_exist = 0,
      recurrence      = 0,
      first_infection = 1,
      omicron  = 0,
      ba4ba5   = 0,
      other    = 0,
      unvaccinated_unboosted  = 1,
      vaccinated_boosted      = 0,
      vaccinated_unreported   = 0,
      unvaccinated_unreported = 0,
      unreported_primary      = 0,
      `age_[0,30)`  = 0,
      `age_[30,50)` = 1,
      `age_[50,100)` = 0,
      dplyr::across(
        c(alpha, delta, vaccinated_unboosted,
          rna_exist, pfu_exist, lfd_exist, sym_exist),
        as.numeric
      )
    )

  # join symptom data
  dat <- dplyr::left_join(dat, ataccc_sym, by = "participant")

  dat <- dat |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      sym_onset = ifelse(
        symptom_onset_relative_to_peak_vl == "Asymptomatic",
        max(days_since_peak),
        as.numeric(symptom_onset_relative_to_peak_vl)
      ),
      sym_ever = ifelse(
        symptom_onset_relative_to_peak_vl != "Asymptomatic", 1, 0
      ),
      symptom_onset_relative_to_peak_vl = as.numeric(
        replace(symptom_onset_relative_to_peak_vl,
                symptom_onset_relative_to_peak_vl == "Asymptomatic", NA)
      ),
      sym = ifelse(
        is.na(symptom_onset_relative_to_peak_vl), 0,
        ifelse(days_since_peak < symptom_onset_relative_to_peak_vl, 0, 1)
      ),
      sym_exist = as.numeric(!is.na(sym_onset))
    ) |>
    dplyr::ungroup()

  # drop all-missing rows and excluded IDs
  dat <- dplyr::filter(dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)
  dat <- dplyr::filter(dat, !id %in% exclude_pfu_raw)

  dat
}


#' Clean UIUC dataset
#' @param uiuc_dat Raw UIUC tibble
#' @param uiuc_sym Symptom data tibble
#' @return Cleaned tibble with standard columns
clean_uiuc <- function(uiuc_dat, uiuc_sym) {

  dat <- uiuc_dat |>
    dplyr::mutate(
      pid      = as.numeric(stringr::str_remove(Ind, " \\*+")),
      time     = Time,
      ct       = Saliva_Ct,
      cn       = replace(Nasal_CN, Nasal_CN == 48, 45),
      rna_ct   = ct_to_rna(ct, type = "uiuc-ct"),
      rna_cn   = ct_to_rna(cn, type = "uiuc-cn"),
      rna = ifelse(
        !is.na(rna_cn) & !is.na(rna_ct),
        log(exp(rna_ct) / 2 + exp(rna_cn) / 2),
        ifelse(is.na(rna_cn), rna_ct, rna_cn)
      ),
      pfu      = replace(Virus_pos_days, Virus_pos_days == 5.9, 6),
      pfu_type = 2,
      lfd      = ifelse(Antigen == "Pos", 1, 0),
      recurrence      = 0,
      first_infection = 1,
      prealpha = ifelse(!Lineage %in% c("B.1.1.7", "P.1"), 1, 0),
      alpha    = ifelse(Lineage == "B.1.1.7", 1, 0),
      delta    = 0,
      omicron  = 0,
      ba4ba5   = 0,
      other    = ifelse(Lineage == "P.1", 1, 0),
      unvaccinated_unboosted  = 1,
      vaccinated_unboosted    = 0,
      vaccinated_boosted      = 0,
      vaccinated_unreported   = 0,
      unvaccinated_unreported = 0,
      unreported_primary      = 0,
      `age_[0,30)`  = ifelse(Age < 30, 1, 0),
      `age_[30,50)` = ifelse(Age >= 30 & Age < 50, 1, 0),
      `age_[50,100)` = ifelse(Age >= 50, 1, 0),
      ct_exist  = as.numeric(!is.na(ct)),
      cn_exist  = as.numeric(!is.na(cn)),
      rna_exist = as.numeric(!is.na(rna)),
      pfu_exist = as.numeric(!is.na(pfu)),
      lfd_exist = as.numeric(!is.na(lfd)),
      sym_exist = 0
    )

  dat <- dat |>
    dplyr::group_by(Ind) |>
    dplyr::mutate(id = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  dat <- dat |>
    dplyr::group_by(pid) |>
    dplyr::mutate(
      peak_day = sum(ifelse(rna == max(rna, na.rm = TRUE), time, NA), na.rm = TRUE),
      days_since_peak = time - peak_day,
      time = days_since_peak
    ) |>
    dplyr::ungroup()

  dat <- dat |>
    dplyr::group_by(Ind) |>
    dplyr::arrange(Ind, Time) |>
    dplyr::mutate(obs = dplyr::row_number())

  uiuc_sym <- uiuc_sym |>
    dplyr::group_by(Ind) |>
    dplyr::arrange(Ind, Time) |>
    dplyr::mutate(obs = dplyr::row_number())

  dat <- dplyr::left_join(dat, uiuc_sym, by = c("Ind", "Time", "obs"))

  dat <- dat |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      sym = ifelse(days_since_peak < sym_onset, 0, 1),
      sym_exist = as.numeric(!is.na(sym_onset))
    ) |>
    dplyr::ungroup()

  dplyr::filter(dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)
}


#' Clean Human Challenge Trial dataset
#' @param hct_dat Raw HCT tibble
#' @param hct_sym Symptom data tibble
#' @return Cleaned tibble with standard columns
clean_hct <- function(hct_dat, hct_sym) {

  dat <- hct_dat |>
    dplyr::mutate(
      rna_ct   = log(replace(qpcr_throat, qpcr_throat == 0, 1)),
      rna_cn   = log(replace(qpcr_nose, qpcr_nose == 0, 1)),
      rna = ifelse(
        !is.na(rna_cn) & !is.na(rna_ct),
        log(exp(rna_ct) / 2 + exp(rna_cn) / 2),
        ifelse(is.na(rna_cn), rna_ct, rna_cn)
      ),
      pfu      = log(replace(throat, throat == 0, 1)),
      pfu_type = 1,
      ffa_t    = log(replace(ffa_throat, ffa_throat == 0, 1)),
      ffa_n    = log(replace(ffa_nose, ffa_nose == 0, 1)),
      ffa = ifelse(
        !is.na(ffa_n) & !is.na(ffa_t),
        log(exp(ffa_t) / 2 + exp(ffa_n) / 2),
        ifelse(is.na(ffa_n), ffa_t, ffa_n)
      ),
      lfd      = lfd,
      recurrence      = 0,
      first_infection = 1,
      prealpha = 1, alpha = 0, delta = 0, omicron = 0, ba4ba5 = 0, other = 0,
      unvaccinated_unboosted  = 1,
      vaccinated_unboosted    = 0,
      vaccinated_boosted      = 0,
      vaccinated_unreported   = 0,
      unvaccinated_unreported = 0,
      unreported_primary      = 0,
      `age_[0,30)`  = 1,
      `age_[30,50)` = 0,
      `age_[50,100)` = 0,
      ct_exist  = as.numeric(!is.na(rna_ct)),
      cn_exist  = as.numeric(!is.na(rna_cn)),
      rna_exist = as.numeric(!is.na(rna)),
      pfu_exist = as.numeric(!is.na(pfu)),
      lfd_exist = as.numeric(!is.na(lfd))
    )

  dat <- dat |>
    dplyr::group_by(pid) |>
    dplyr::mutate(
      peak_day = sum(ifelse(rna == max(rna, na.rm = TRUE), day, NA), na.rm = TRUE),
      days_since_peak = day - peak_day,
      time = days_since_peak
    ) |>
    dplyr::ungroup()

  # --- symptom aggregation ---
  hct_sym <- hct_sym |>
    dplyr::select(-dplyr::ends_with("_desc")) |>
    dplyr::filter(
      !collection_method %in% c(
        "Paper diary card completed",
        "Not yet admitted",
        "Subject discharged"
      )
    )

  hct_sym <- hct_sym |>
    dplyr::group_by(id, day) |>
    dplyr::arrange(id, day, timepoint) |>
    dplyr::summarise(
      dplyr::across(runny_nose:score, \(x) max(x, na.rm = TRUE)),
      .groups = "drop"
    )

  hct_sym <- hct_sym |>
    dplyr::group_by(id) |>
    dplyr::arrange(id, day) |>
    dplyr::mutate(
      sym = ifelse(cumsum(score) > 0, 1, 0),
      sym_ever = max(sym, na.rm = TRUE),
      x = sym * day + (1 - sym) * max(day),
      sym_onset = min(x, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  dat <- dplyr::left_join(dat, hct_sym, by = c("id", "day"))

  dat <- dat |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      sym_exist = as.numeric(!is.na(score)),
      sym = replace(sym, sym_exist == 0, 0),
      sym_ever  = dplyr::first(sym_ever),
      sym_onset = dplyr::first(sym_onset)
    ) |>
    dplyr::ungroup()

  dplyr::filter(dat, (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0)
}


#' Clean Legacy (Crick) dataset
#' @param legacy_dat Raw legacy tibble
#' @return Cleaned tibble with standard columns
clean_legacy <- function(legacy_dat) {
  covariates <- define_covariates()

  dat <- legacy_dat |>
    dplyr::mutate(
      id       = id,
      pid      = id,
      ct       = ct_value,
      # Legacy Ct values use the NBA calibration curve (Kissler et al.)
      # with a -2 offset (on the natural-log scale) to account for the
      # difference in extraction/amplification efficiency between the
      # Crick Legacy and NBA assay platforms. This offset was estimated
      # by comparing paired samples processed on both platforms and
      # represents approximately a 7.4-fold reduction in apparent
      # copies/mL (exp(2) ≈ 7.4).
      rna      = ct_to_rna(ct_value, type = "nba") - 2,
      pfu      = NA,
      lfd      = NA,
      pfu_type = 1,
      sym = dplyr::case_when(
        symptoms == "symptomatic" & swab_date < symptom_onset_date  ~ 0,
        symptoms == "symptomatic" & swab_date >= symptom_onset_date ~ 1,
        symptoms == "asymptomatic" ~ 0
      ),
      sym_ever  = as.numeric(symptoms == "symptomatic"),
      `age_[0,30)`  = age_group == "20-34",
      `age_[30,50)` = age_group == "35-49",
      `age_[50,100)` = age_group == "50+",
      recurrence      = infection_id >= 2,
      first_infection = 1 - recurrence,
      prealpha = VOC == "None",
      alpha    = VOC == "Alpha",
      delta    = VOC == "Delta",
      omicron  = VOC %in% c("Omicron (BA.1)", "Omicron (BA.2)"),
      ba4ba5   = VOC %in% c("Omicron (BA.4)", "Omicron (BA.5)"),
      other    = VOC %in% c("Other"),
      unvaccinated_unboosted  = 0,
      vaccinated_unboosted    = no_vaccines == 2,
      vaccinated_boosted      = no_vaccines == 3,
      vaccinated_unreported   = 0,
      unvaccinated_unreported = 0,
      unreported_primary      = 0,
      rna_exist = as.numeric(!is.na(rna)),
      pfu_exist = 0,
      lfd_exist = 0,
      sym_exist = as.numeric(!is.na(sym)),
      dplyr::across(
        c(sym, alpha, delta, omicron, ba4ba5, other,
          vaccinated_unboosted, vaccinated_boosted,
          rna_exist, pfu_exist, lfd_exist, sym_exist),
        as.numeric
      )
    )

  dat <- dat |>
    dplyr::group_by(pid) |>
    dplyr::mutate(
      peak_day  = sum(ifelse(rna == max(rna, na.rm = TRUE), t, NA), na.rm = TRUE),
      peak_date = max(ifelse(rna == max(rna, na.rm = TRUE), swab_date, NA), na.rm = TRUE),
      days_since_peak = t - peak_day,
      sym_onset = as.numeric(symptom_onset_date - peak_date),
      sym_onset = replace(sym_onset, sym_ever == 0, max(days_since_peak)),
      time = days_since_peak
    ) |>
    dplyr::ungroup()

  dat <- dplyr::filter(
    dat,
    (rna_exist + pfu_exist + lfd_exist + sym_exist) != 0 & ct_type == "ct_value"
  )

  # add onset observations for individuals whose onset day differs from swab day
  add_onset_obs <- dat |>
    dplyr::group_by(pid) |>
    dplyr::slice(1) |>
    dplyr::select(
      pid, id, time, rna, pfu_type,
      sym, sym_ever, sym_onset, sym_exist,
      symptom_onset_date, swab_date,
      dplyr::all_of(covariates$x_vars_w_refs),
      dplyr::all_of(covariates$miss_vars)
    ) |>
    dplyr::filter(!is.na(sym_onset) & symptom_onset_date != swab_date) |>
    dplyr::mutate(
      time      = sym_onset,
      sym       = 1,
      swab_date = NA,
      rna       = NA,
      rna_exist = 0
    ) |>
    dplyr::ungroup()

  dplyr::bind_rows(dat, add_onset_obs) |>
    dplyr::arrange(pid, time)
}


# ── Stacking & Stan data ─────────────────────────────────────────────────────

#' Stack cleaned source datasets into a single analysis dataset
#'
#' @param nba,ataccc,uiuc,hct,legacy  Cleaned tibbles from clean_* functions
#' @return Stacked tibble with sequential id, sym_at_risk, and missings filled
stack_sources <- function(nba, ataccc, uiuc, hct, legacy) {
  covariates <- define_covariates()

  keep_vars <- c(
    "pid", "id", "time", "rna", "pfu", "pfu_type", "lfd",
    "sym", "sym_onset", "sym_ever",
    covariates$x_vars_w_refs,
    covariates$miss_vars
  )

  stacked <- dplyr::bind_rows(
    dplyr::select(nba,    dplyr::all_of(keep_vars)),
    dplyr::select(ataccc, dplyr::all_of(keep_vars)),
    dplyr::select(uiuc,   dplyr::all_of(keep_vars)),
    dplyr::select(hct,    dplyr::all_of(keep_vars)),
    dplyr::select(legacy, dplyr::all_of(keep_vars)),
    .id = "source"
  )

  # sequential IDs within source

  stacked <- stacked |>
    dplyr::group_by(source, id) |>
    dplyr::mutate(id = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  # fill missing indicators
  stacked <- stacked |>
    dplyr::mutate(
      rna      = replace(rna, rna_exist == 0, 0),
      pfu      = replace(pfu, pfu_exist == 0, 0),
      lfd      = replace(lfd, lfd_exist == 0, 0),
      sym      = replace(sym, sym_exist == 0, 0),
      sym_onset = replace(sym_onset, sym_exist == 0, 99),
      sym_ever  = replace(sym_ever, sym_exist == 0, 0)
    )

  # ---------- Documented exclusions ----------
  # Individual 1932 (NBA pid=1932): 17-day gap followed by near-LOD "rebound"
  # at days 22-23, creating bimodal clearance rate posterior. Remove the 3
  # late observations (time >= 20) to resolve; keep days -5 to 5.
  stacked <- stacked |>
    dplyr::filter(!(pid == 1932 & source == "1" & time >= 20))

  # Re-compute sequential IDs after exclusion (some may have lost all obs)
  stacked <- stacked |>
    dplyr::group_by(source, pid) |>
    dplyr::mutate(id = dplyr::cur_group_id()) |>
    dplyr::ungroup()

  # construct sym_at_risk: at-risk for symptom onset (pre-onset + onset day)
  stacked <- stacked |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      cum_sym_prev = dplyr::lag(cumsum(sym), default = 0),
      sym_at_risk  = as.integer(sym_exist == 1 & cum_sym_prev == 0)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-cum_sym_prev)

  stacked
}


#' Build the Stan data list from the stacked dataset
#'
#' @param stacked_dat  Output of stack_sources()
#' @param flags        Named list of model flags (overrides defaults)
#' @return Named list suitable for cmdstanr::$sample(data = ...)
build_stan_data <- function(stacked_dat, flags = list()) {
  covariates <- define_covariates()

  # default flags
  default_flags <- list(
    ind_effects = 1,
    ind_corr    = 0,
    test_error  = 1,
    adj_pfu     = 0,
    adj_rna     = 1,
    adj_lfd     = 0,
    adj_sym     = 0,
    source_pfu  = 0,
    source_rna  = 0,
    source_lfd  = 0,
    source_sym  = 0,
    use_wf      = 0,
    use_smooth  = 0
  )
  flags <- modifyList(default_flags, flags)

  # Compute M and N from stacked data (post-exclusion)
  source_summary <- stacked_dat |>
    dplyr::group_by(source) |>
    dplyr::summarise(
      M = dplyr::n_distinct(id),
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(as.numeric(source))

  stan_data <- list(
    K = 5,
    M = source_summary$M,
    N = source_summary$N,
    lod_rna = c(
      ct_to_rna(40, type = "nba")      + 0.01,
      ct_to_rna(40, type = "ata")      + 0.01,
      ct_to_rna(47, type = "uiuc-ct")  + 0.01,
      log(31.62278)                     + 0.01,  # 10^1.5 copies/mL — Killingley et al. qPCR LOD
      ct_to_rna(40, type = "nba")      + 0.01
    ),
    lod_pfu = c(0, 2.3, 2.3, log(5), 0),
    fp_mean = c(
      1 / (ct_to_rna(40 - 1/log(10), "nba")     - ct_to_rna(40, "nba")),
      1 / (ct_to_rna(40 - 1/log(10), "ata")     - ct_to_rna(40, "ata")),
      1 / (ct_to_rna(47 - 1/log(10), "uiuc-ct") - ct_to_rna(47, "uiuc-ct")),
      1 / 0.3,
      1 / (ct_to_rna(40 - 1/log(10), "nba")     - ct_to_rna(40, "nba"))
    ),
    id       = stacked_dat$id,
    time     = stacked_dat$time,
    rna      = stacked_dat$rna,
    pfu      = stacked_dat$pfu,
    lfd      = stacked_dat$lfd,
    sym      = stacked_dat$sym,
    sym_ever = stacked_dat$sym_ever,
    sym_at_risk = stacked_dat$sym_at_risk,
    source   = as.numeric(stacked_dat$source),
    pfu_type = stacked_dat$pfu_type,
    P = length(covariates$x_vars),
    x = stacked_dat |>
      dplyr::group_by(id) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(covariates$x_vars)),
    rna_exist = stacked_dat$rna_exist,
    pfu_exist = stacked_dat$pfu_exist,
    lfd_exist = stacked_dat$lfd_exist,
    sym_exist = stacked_dat$sym_exist,
    # model flags
    ind_effects = flags$ind_effects,
    ind_corr    = flags$ind_corr,
    test_error  = flags$test_error,
    adj_pfu     = flags$adj_pfu,
    adj_rna     = flags$adj_rna,
    adj_lfd     = flags$adj_lfd,
    adj_sym     = flags$adj_sym,
    source_pfu  = flags$source_pfu,
    source_rna  = flags$source_rna,
    source_lfd  = flags$source_lfd,
    source_sym  = flags$source_sym,
    use_wf      = flags$use_wf,
    use_smooth  = flags$use_smooth,
    # prior hyperparameters
    prior_fp       = 0.02,
    prior_fn       = 0.01,
    prior_fp_kappa = 100,
    prior_dp_mean  = 17,
    prior_dp_cv    = 0.7,
    prior_wp_mean  = 8,
    prior_wp_cv    = 0.7,
    prior_wr_mean  = 15,
    prior_wr_cv    = 0.7,
    prior_wf_mean  = 1,
    prior_wf_cv    = 0.7,
    prior_sigma_sd = 5,
    prior_beta_sd  = 1,
    prior_i_sd     = 1,
    prior_pfu_i_sd = 0.3,  # tighter prior for PFU RE SDs
    prior_k_sd     = 1,
    prior_lfd_mean = 0.01,

    # Smooth post-peak sigmoid steepness: inv_logit(kappa * (t - tp)).
    # kappa=5 ≈ transition over ~1 day (biologically reasonable).
    kappa_postpeak = 5.0
  )

  # grainsize for reduce_sum — 1 lets Stan auto-schedule slice sizes.
  # Increase to 50-200 for potentially better cache locality.
  stan_data$grainsize <- 1L

  # PFU individual-effect mapping: restrict PFU REs to informed individuals
  stan_data <- add_pfu_ind_mapping(stan_data)

  stan_data
}


#' Compute PFU individual-effect mapping
#'
#' Identifies which individuals have PFU-informing data (viral culture,
#' LFD, or symptom observations) and creates the mapping arrays needed
#' by the Stan model to restrict PFU individual random effects to only
#' those individuals.  Individuals without PFU-informing data (e.g. NBA
#' with RNA-only) get the population-level RNA->PFU transform with no
#' individual deviation.  The hierarchical SD (sigma_ind_pfu) is still
#' learned from the informed individuals.
#'
#' @param stan_data  Stan data list (from build_stan_data or subsample_stan_data)
#' @return stan_data with N_pfu_ind and pfu_ind_idx added
add_pfu_ind_mapping <- function(stan_data) {
  N_ind <- sum(stan_data$M)

  # An individual is PFU-informed if ANY of their observations have
  # culture data (pfu_exist), lateral flow (lfd_exist), or symptom data
  # (sym_exist) — all three sub-models depend on pfu_hat.
  has_pfu_info <- as.logical(tapply(
    stan_data$pfu_exist | stan_data$lfd_exist | stan_data$sym_exist,
    stan_data$id,
    any
  ))

  # Build mapping: 0 = no PFU RE, 1..N_pfu_ind = index into PFU RE arrays
  pfu_ind_idx <- integer(N_ind)
  counter <- 0L
  for (i in seq_len(N_ind)) {
    if (has_pfu_info[i]) {
      counter <- counter + 1L
      pfu_ind_idx[i] <- counter
    }
  }

  stan_data$N_pfu_ind   <- counter
  stan_data$pfu_ind_idx <- pfu_ind_idx

  message(sprintf(
    "PFU individual effects: %d of %d individuals have PFU-informing data",
    counter, N_ind
  ))

  stan_data
}


#' Subsample a Stan data list to fewer individuals per source
#'
#' Randomly samples \code{frac} of individuals from each source (or up to
#' \code{max_per_source}) and rebuilds all arrays with consistent dimensions.
#' Useful for faster parameter-recovery runs.
#'
#' @param data  Stan data list (from \code{build_stan_data})
#' @param frac  Fraction of individuals to retain per source (default 0.1)
#' @param max_per_source  Hard cap on individuals per source (default Inf)
#' @param seed  Random seed for reproducibility
#' @return A Stan data list with smaller M, N, and re-indexed id
subsample_stan_data <- function(data, frac = 0.1, max_per_source = Inf,
                                seed = 123) {
  set.seed(seed)
  K <- data$K

  # Identify which observations belong to each source
  cum_N <- cumsum(c(0, data$N))
  cum_M <- cumsum(c(0, data$M))

  keep_obs <- integer(0)
  new_M <- integer(K)
  id_map <- integer(0)  # old_id -> new_id mapping

  global_new_id <- 0
  new_id_vec <- integer(0)


  for (k in 1:K) {
    # Observation indices for source k
    obs_idx <- (cum_N[k] + 1):cum_N[k + 1]
    # Individual IDs in this source (original indexing)
    ids_in_source <- unique(data$id[obs_idx])
    n_keep <- max(1, min(ceiling(length(ids_in_source) * frac),
                         max_per_source))
    keep_ids <- sort(sample(ids_in_source, n_keep))

    new_M[k] <- length(keep_ids)

    # Map old IDs to new sequential IDs
    for (old_id in keep_ids) {
      global_new_id <- global_new_id + 1
      obs_for_id <- obs_idx[data$id[obs_idx] == old_id]
      keep_obs <- c(keep_obs, obs_for_id)
      new_id_vec <- c(new_id_vec, rep(global_new_id, length(obs_for_id)))
    }
  }

  # Compute new N per source
  new_source <- data$source[keep_obs]
  new_N <- integer(K)
  for (k in 1:K) new_N[k] <- sum(new_source == k)

  # Build the x matrix for kept individuals
  # Original x is indexed 1:sum(M); we need the rows for kept IDs
  kept_original_ids <- unique(data$id[keep_obs])
  # But x rows correspond to individual index in original data
  # We need to pick the correct rows
  x_orig <- as.matrix(data$x)
  x_new <- x_orig[kept_original_ids, , drop = FALSE]

  # Rebuild the data list
  out <- data
  out$M <- new_M
  out$N <- new_N
  out$id <- new_id_vec
  out$time <- data$time[keep_obs]
  out$rna <- data$rna[keep_obs]
  out$pfu <- data$pfu[keep_obs]
  out$lfd <- data$lfd[keep_obs]
  out$sym <- data$sym[keep_obs]
  out$sym_ever <- data$sym_ever[keep_obs]
  out$sym_at_risk <- data$sym_at_risk[keep_obs]
  out$source <- data$source[keep_obs]
  out$pfu_type <- data$pfu_type[keep_obs]
  out$x <- x_new
  out$rna_exist <- data$rna_exist[keep_obs]
  out$pfu_exist <- data$pfu_exist[keep_obs]
  out$lfd_exist <- data$lfd_exist[keep_obs]
  out$sym_exist <- data$sym_exist[keep_obs]

  # Recompute PFU individual-effect mapping for subsampled data
  out <- add_pfu_ind_mapping(out)

  out
}