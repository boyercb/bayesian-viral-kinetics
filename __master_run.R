
# packages ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(posterior)
library(kableExtra)
library(truncnorm)
library(mvtnorm)


# globals -----------------------------------------------------------------

RUN_MCMC <- FALSE
RUN_TEST <- FALSE

# load data ---------------------------------------------------------------

# NBA dataset
nba_dat <- read_csv("0_data/nba_dat.csv")

# ATACCC dataset
ataccc_dat <- read_csv("0_data/ataccc_dat.csv")
ataccc_sym <- read_csv("0_data/ataccc_sx_dat.csv")

# college students
uiuc_dat <- read_csv("0_data/uiuc_dat.csv")
uiuc_sym <- read_csv("0_data/uiuc_sx_dat.csv")

# human challenge trial
hct_dat <- read_csv("0_data/hct_dat.csv")
hct_sym <- read_csv("0_data/hct_sx_dat.csv")

# legacy cohort 
legacy_dat <- read_csv("0_data/legacy_dat.csv")
  

# if running in test mode sample just a few observations
if (RUN_TEST) {
  n_samples <- 1
  
  # NBA datset
  nba_smp <- sample(nba_dat$InfectionEvent, n_samples)
  nba_dat <- filter(nba_dat, InfectionEvent %in% nba_smp)
  
} 


# run code ----------------------------------------------------------------

source("1_code/functions.R")

source("1_code/clean_data.R")

source("1_code/prior_predictive.R")

if (RUN_MCMC) {
  source("1_code/fit_model.R")
} else {
  kinetics_model <- read_rds("1_code/stan/kinetics_model.rds")
}

source("1_code/check_model.R")

source("1_code/model_summaries.R")

source("1_code/prediction.R")