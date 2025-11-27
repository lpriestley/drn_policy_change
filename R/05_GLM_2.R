## 05_GLM_2.R --------------------------------------------
## - GLM2.* 

## Setup -------------------------------------------------------------------

source("R/00_setup.R")

d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

## 1. Descriptive statistics: policy switches ------------------------------
n_switches <- d %>%
  filter(mag == 10, policy_switch == 1) %>%
  count(subject, mag, policy_switch, name = "n") %>%
  summarise(
    mu    = mean(n, na.rm = TRUE),
    sem   = sd(n, na.rm = TRUE) / sqrt(sum(!is.na(n))),
    total = sum(n, na.rm = TRUE),
    .groups = "drop"
  )

print(n_switches)

## 2. GLM2.1: Option-specific policy  ---------------------

GLM2.1 <- fit_glmer_binom(
  formula = response ~ mag_z + prev_policy + trial_z +
    (mag_z + prev_policy | subject),
  data    = d,
  name    = "GLM_2.1"
)

## 3. GLM2.2: Action-history  ---------------------

GLM2.2 <- fit_glmer_binom(
  formula = response ~ mag_z + prev_action + trial_z +
    (prev_action | subject),
  data    = d,
  name    = "GLM_2.2"
)


## 4. GLM2.3: Frequency of policy switches by magnitude --------------------

GLM2.3a <- fit_glmer_binom(
  formula = policy_switch ~ factor(mag) + trial_n_block_z +
    (factor(mag) + trial_n_block_z | subject),
  data    = d %>% filter(mag != 5),
  name    = "GLM_2.3a"
)

GLM2.3b <- fit_glmer_binom(
  formula = policy_switch ~ factor(mag) + trial_n_block_z +
    (factor(mag) + trial_n_block_z | subject),
  data    = d %>% filter(mag != 50),
  name    = "GLM_2.3b"
)


## 5. GLM2.4: Policy-switch types in rich-vs-poor env -----------

GLM2.4a <- fit_glmer_binom(
  formula = policy_switch ~ prev_policy + (prev_policy | subject),
  data    = d %>% filter(mag == 10, env == "poor"),
  name    = "GLM_2.4a"
)

GLM2.4b <- fit_glmer_binom(
  formula = policy_switch ~ prev_policy + (prev_policy | subject),
  data    = d %>% filter(mag == 10, env == "rich"),
  name    = "GLM_2.4b"
)


## 6. GLM2.5: Pursue-vs-reject switches in rich-vs-poor env ---------------------

GLM2.5a <- fit_glmer_binom(
  formula = pursue_switch ~ env_bin + (env_bin | subject),
  data    = d %>% filter(mag == 10, prev_policy == 0),
  name    = "GLM_2.5a"
)

GLM2.5b <- fit_glmer_binom(
  formula = reject_switch ~ env_bin + (env_bin | subject),
  data    = d %>% filter(mag == 10, prev_policy == 1),
  name    = "GLM_2.5b"
)


## 7. GLM2.6: Policy-switch vs average value -------------------

GLM2.6a <- fit_glmer_binom(
  formula = pursue_switch ~ mu_val_z + (mu_val_z | subject),
  data    = d %>% filter(mag == 10, prev_policy == 0),
  name    = "GLM_2.6a"
)

GLM2.6b <- fit_glmer_binom(
  formula = reject_switch ~ mu_val_z + (mu_val_z | subject),
  data    = d %>% filter(mag == 10, prev_policy == 1),
  name    = "GLM_2.6b"
)

GLM2.6c <- fit_glmer_binom(
  formula = pursue_switch ~ mu_val_z + trial_z + (1 | subject),
  data    = d %>% filter(mag == 5, prev_policy == 0),
  name    = "GLM_2.6c"
)

GLM2.6d <- fit_glmer_binom(
  formula = reject_switch ~ mu_val_z + trial_z + (1 | subject),
  data    = d %>% filter(mag == 5, prev_policy == 1),
  name    = "GLM_2.6d"
)

GLM2.6e <- fit_glm_binom(
  formula = pursue_switch ~ mu_val_z + trial_z,
  data    = d %>% filter(mag == 50, prev_policy == 0),
  name    = "GLM_2.6e"
)

GLM2.6f <- fit_glm_binom(
  formula = reject_switch ~ mu_val_z + trial_z,
  data    = d %>% filter(mag == 50, prev_policy == 1),
  name    = "GLM_2.6f"
)

## 2. Switch-timing  -----------------------------------------
congr_switches <- d %>%
  filter(switch_congruence == 1)

## Per-subject mean timing within each environment
change_loc <- congr_switches %>%
  group_by(subject, env_bin) %>%
  summarise(m = mean(trial_n_block, na.rm = TRUE), .groups = "drop") %>%
  group_by(env_bin) %>%
  summarise(
    mean_time = mean(m, na.rm = TRUE),
    sd_time   = sd(m, na.rm = TRUE),
    .groups   = "drop"
  )

print(change_loc)

## Overall SD across all congruent switches
congr_switches_sd <- sd(congr_switches$trial_n_block, na.rm = TRUE)
print(congr_switches_sd)

GLM2.7 <- lm(
  trial_n_block_z ~ env_bin,
  data = congr_switches
)

saveRDS(
  GLM2.7,
  here::here("results", "behaviour", "GLM_2.7.rds")
)

sink(here::here("results", "behaviour", "GLM_2.7.txt"))
print(summary(GLM2.7))
sink()
