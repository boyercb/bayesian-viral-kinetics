# ──────────────────────────────────────────────────────────────────────────────
# utils.R — General-purpose utility functions
# ──────────────────────────────────────────────────────────────────────────────

#' Format a number to k decimal places (trimmed)
specd <- function(x, k) trimws(format(round(x, k), nsmall = k))

#' Piecewise-linear viral kinetics model with optional flat-top
#'
#' @param t   Time (scalar or vector)
#' @param tp  Time of peak
#' @param wp  Proliferation width
#' @param wr  Clearance width
#' @param dp  Peak viral load (log-scale)
#' @param wf  Flat-top duration at peak (days); 0 = sharp peak (default)
#' @return Viral load at time t (log-scale)
pefun <- function(t, tp, wp, wr, dp, wf = 0) {
  ifelse(
    t <= tp,
    dp / wp * (t - (tp - wp)),
    ifelse(
      t <= tp + wf,
      dp,
      dp - dp / wr * (t - tp - wf)
    )
  )
}

#' Smooth trajectory approximation with optional flat-top
#'
#' Log-sum-exp approximation to the piecewise-linear trajectory.
#' The falling arm is shifted right by wf, creating an approximate
#' plateau of width ~wf. When wf = 0, reduces to the standard
#' smooth two-arm envelope.
#'
#' @param t   Time (scalar or vector)
#' @param tp  Time of peak
#' @param wp  Proliferation width
#' @param wr  Clearance width
#' @param dp  Peak viral load (log-scale)
#' @param wf  Flat-top duration at peak (days); 0 = no flat-top (default)
#' @return Viral load at time t (log-scale)
smfun <- function(t, tp, wp, wr, dp, wf = 0) {
  a <- dp / wp      # proliferation rate
  b <- dp / wr      # clearance rate
  arg1 <- pmin(-a * (t - tp), 50)
  arg2 <- pmin( b * (t - (tp + wf)), 50)
  raw <- dp + log((a + b) / (b * exp(arg1) + a * exp(arg2)))
  # Soft-cap at dp: scaled soft_min(raw, dp) with k=10.
  # Prevents overshoot when wf > 0; negligible for wf=0 (~0.07 log-units).
  k_cap <- 10
  delta <- k_cap * (raw - dp)
  # Numerically stable softplus: log1p(exp(x)) ≈ x for x >> 0
  sp <- ifelse(delta > 30, delta, log1p(exp(delta)))
  raw - sp / k_cap
}

#' Trajectory dispatcher: calls pefun or smfun based on use_smooth flag
#'
#' @param t   Time (scalar or vector)
#' @param tp  Time of peak
#' @param wp  Proliferation width
#' @param wr  Clearance width
#' @param dp  Peak viral load (log-scale)
#' @param wf  Flat-top duration (days); 0 = no flat-top
#' @param use_smooth  Logical; if TRUE, use smfun; otherwise pefun
#' @return Viral load at time t (log-scale)
traj_fun <- function(t, tp, wp, wr, dp, wf = 0, use_smooth = FALSE) {
  if (use_smooth) {
    smfun(t, tp, wp, wr, dp, wf)
  } else {
    pefun(t, tp, wp, wr, dp, wf)
  }
}

#' Repeat vector as rows of a matrix
rep_row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

#' Repeat vector as columns of a matrix
rep_col <- function(x, n) {

  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

# Backward-compatible aliases (used in prior_predictive)
rep.row <- rep_row
rep.col <- rep_col

#' Convert Ct value to log-RNA copies/mL using study-specific calibration
#'
#' @param x    Ct value(s)
#' @param type Study identifier: "nba", "ata", "uiuc-cn", "uiuc-ct", "hct-cn", "hct-ct"
#' @return log-RNA on the natural log scale (except "ata" which returns log10)
ct_to_rna <- function(x, type) {

  if (type == "nba") {
    # Kissler et al. (2021) NEJMc2102507
    rna <- (x - 40.93733) / (-3.60971) + log10(250)
    rna <- rna * log(10)

  } else if (type == "ata") {
    # Singanayagam et al. (2022) Lancet ID
    rna <- (40 - x) / 1.418 + 3.435

  } else if (type == "uiuc-cn") {
    # Ke et al. (2022) Nat Microbiol — nasal swabs
    rna <- (11.35 - 0.25 * x)
    rna <- rna * log(10)

  } else if (type == "uiuc-ct") {
    # Ke et al. (2022) Nat Microbiol — saliva swabs
    rna <- (14.24 - 0.28 * x)
    rna <- rna * log(10)

  } else if (type == "hct-cn") {
    # Killingley et al. (2022) Nature Medicine — nasal swabs
    # HCT data are provided as log copies/mL directly (not Ct values),
    # so no calibration conversion is needed. The raw values from the
    # challenge study dataset are already on the natural-log scale.
    rna <- x

  } else if (type == "hct-ct") {
    # Killingley et al. (2022) Nature Medicine — throat swabs
    # Same as nasal: data are already log copies/mL, no conversion needed.
    rna <- x
  }

  rna
}

#' Calculate onset / peak / resolution times from triangle parameters
#'
#' @param tp   Peak time
#' @param wp   Proliferation width
#' @param wr   Clearance width
#' @param dp   Peak viral load
#' @param lod  Limit of detection
#' @param wf   Flat-top duration (days, default 0)
#' @return data.frame with columns to (onset), tp (peak), tr (resolution)
calc_corners <- function(tp, wp, wr, dp, lod, wf = 0) {
  to <- lod * wp / dp + (tp - wp)
  tr <- (dp - lod) * wr / dp + tp + wf
  data.frame(to = to, tp = tp, tr = tr)
}

#' Build a 3-row data frame of onset / peak / resolution corner points
add_corners <- function(id, tp, wp, wr, dp, lod, wf = 0, name = "rna") {
  corner_times <- calc_corners(tp, wp, wr, dp, lod, wf)
  df <- data.frame(
    id = rep(id, 3),
    time = unlist(corner_times)
  )
  df[[name]] <- c(lod, dp, lod)
  rownames(df) <- NULL
  df
}

#' Inverse-logit
expit <- function(x) exp(x) / (1 + exp(x))

#' Logit
logit <- function(x) log(x / (1 - x))

# rvars-compatible wrappers (for posterior package)
rvar_pefun <- posterior::rfun(pefun)
rvar_smfun <- posterior::rfun(smfun)
rvar_traj_fun <- posterior::rfun(traj_fun)
rvar_calc_corners <- posterior::rfun(calc_corners)
