## 01_preprocess_behavioural.R --------------------------------------------
## Load raw behavioural data, construct derived variables, and save
## cleaned dataset for downstream scripts.

source("R/00_setup.R")

## 1. Load raw data --------------------------------------------------------

d <- read.csv(
  here::here("data", "behaviour", "raw_lf.csv"),
  header = TRUE
)

## 2. Generic recoding --------------------------------------------

d <- d %>%
  mutate(
    env = case_when(
      env_bin ==  1 ~ "rich",
      env_bin == -1 ~ "poor",
      TRUE          ~ NA_character_
    ),
    env = factor(env, levels = c("poor", "rich")),
    
    mag_factor = paste0(mag, "-points")
  )


## 3. Policy-related variables --------------------------------------------

d <- d %>%
  mutate(
    policy_switch = case_when(
      is.na(prev_policy)            ~ NA_real_,
      response == prev_policy       ~ 0,
      response != prev_policy       ~ 1
    ),
    
    pursue_switch = ifelse(response == 1 & prev_policy == 0, 1, 0),
    reject_switch = ifelse(response == 0 & prev_policy == 1, 1, 0),
    
    switch_type = case_when(
      pursue_switch == 1 ~ "pursue-switch",
      reject_switch == 1 ~ "reject-switch",
      TRUE               ~ "none"
    ),
    
    switch_congruence = case_when(
      pursue_switch == 1 & env_bin == -1 & mag == 10 ~ 1,
      pursue_switch == 1 & env_bin ==  1 & mag == 10 ~ 0,
      reject_switch == 1 & env_bin ==  1 & mag == 10 ~ 1,
      reject_switch == 1 & env_bin == -1 & mag == 10 ~ 0,
      TRUE                                            ~ NA_real_
    ),
    
    switch_type_bin = case_when(
      pursue_switch == 1 ~  1,
      reject_switch == 1 ~ -1,
      TRUE               ~ NA_real_
    )
  )

## 4. Ave. Value PE variable ----------------------------------------------

d <- d %>%
  mutate(
    mu_val_pe_z = mag_z - mu_val_z
  )

## 5. Save cleaned dataset -------------------------------------------------

saveRDS(
  d,
  here::here("data", "behaviour", "clean_lf.rds")
)
