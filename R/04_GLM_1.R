## 04_GLM_1.R ----------------------------------------------
## - GLM1.*

## Setup ------------------------------------------------------------------

source("R/00_setup.R")

# Load data
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))


## 1. Descriptive statistics ----------------------------------------------

subj_means <- d %>%
  group_by(env_bin, subject) %>%
  summarise(subj_means = mean(response, na.rm = TRUE), .groups = "drop")

subj_means_summary <- subj_means %>%
  group_by(env_bin) %>%
  summarise(
    mu  = 100 * mean(subj_means, na.rm = TRUE),
    sem = 100 * (sd(subj_means, na.rm = TRUE) /
                   sqrt(sum(!is.na(subj_means)))),
    .groups = "drop"
  ) %>% 
  print()


## 2. GLM1.1: Environment × magnitude effects -----------------------------

GLM1.1a <- fit_glmer_binom(
  formula = response ~ mag_z * env_bin + trial_z + (mag_z:env_bin | subject),
  data    = d,
  name    = "GLM_1.1a"
)

## GLM1.1b mid / low / high

GLM1.1b_mid <- fit_glmer_binom(
  formula = response ~ env_bin + trial_z + (env_bin + trial_z | subject),
  data    = d %>% filter(mag == 10),
  name    = "GLM_1.1b_mid"
)

GLM1.1b_low <- fit_glmer_binom(
  formula = response ~ env_bin + trial_z + (env_bin + trial_z | subject),
  data    = d %>% filter(mag == 5),
  name    = "GLM_1.1b_low"
)

GLM1.1b_high <- fit_glmer_binom(
  formula = response ~ env_bin + trial_z + (1 | subject),
  data    = d %>% filter(mag == 50),
  name    = "GLM_1.1b_high"
)


## 3. GLM1.2: Environment × time-in-environment ---------------------------

GLM1.2_mid <- fit_glmer_binom(
  formula = response ~ env_bin * trial_n_block_z + (1 | subject),
  data    = d %>% filter(mag == 10),
  name    = "GLM_1.2_mid"
)

GLM1.2_low <- fit_glmer_binom(
  formula = response ~ env_bin * trial_n_block_z + (1 | subject),
  data    = d %>% filter(mag == 5),
  name    = "GLM_1.2_low"
)

GLM1.2_high <- fit_glmer_binom(
  formula = response ~ env_bin * trial_n_block_z + (1 | subject),
  data    = d %>% filter(mag == 50),
  name    = "GLM_1.2_high"
)

## Within-env time effects for mid-option only

GLM1.2_mid_poor <- fit_glmer_binom(
  formula = response ~ trial_n_block_z + (1 | subject),
  data    = d %>% filter(mag == 10, env_bin == -1),
  name    = "GLM_1.2_mid_poor"
)

GLM1.2_mid_rich <- fit_glmer_binom(
  formula = response ~ trial_n_block_z + (1 | subject),
  data    = d %>% filter(mag == 10, env_bin == 1),
  name    = "GLM_1.2_mid_rich"
)


## 4. GLM1.3: Average value effects ---------------------------------------

GLM1.3_mid <- fit_glmer_binom(
  formula = response ~ mu_val_z + trial_z + (mu_val_z | subject),
  data    = d %>% filter(mag == 10),
  name    = "GLM_1.3_mid"
)

GLM1.3_low <- fit_glmer_binom(
  formula = response ~ mu_val_z + trial_z + (mu_val_z | subject),
  data    = d %>% filter(mag == 5),
  name    = "GLM_1.3_low"
)

GLM1.3_high <- fit_glmer_binom(
  formula = response ~ mu_val_z + trial_z + (mu_val_z | subject),
  data    = d %>% filter(mag == 50),
  name    = "GLM_1.3_high"
)


## 5. Acceptance-rate vs total reward -------------------------------------

total_rewards <- d %>%
  group_by(subject) %>%
  summarise(total_reward = sum(reward, na.rm = TRUE), .groups = "drop") %>%
  mutate(total_reward = as.numeric(scale(total_reward)))

acc_rates <- d %>%
  group_by(subject, mag, env) %>%
  summarise(acc = mean(response, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = c(env, mag),
    names_sep   = "_",
    values_from = acc
  ) %>%
  left_join(total_rewards, by = "subject") %>%
  ungroup()

GLM1.4 <- lm(total_reward ~ poor_10 * rich_10, data = acc_rates)

# Save regression
saveRDS(
  GLM1.4,
  here::here("results", "behaviour", "GLM_1.4.rds")
)

sink(here::here("results", "behaviour", "GLM_1.4_summary.txt"))
sink()
