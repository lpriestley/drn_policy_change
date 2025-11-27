## 10_fig2.R ---------------------------------------------------------------
##  Figure 2

source("R/00_setup.R")
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

# ==========================================================================
# FIGURE 2A
# ==========================================================================

fig_2_A_data <- d %>%
  group_by(subject, mag, env) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  mutate(m_percent = m * 100)

fig_2_A <- ggplot(fig_2_A_data,
                aes(x = factor(mag), fill = factor(env),
                    color = factor(env), y = m_percent)
) +
  stat_summary(fun = "mean", geom = "bar",
               color = "black", width = 0.85, alpha = 0.75,
               position = position_dodge(width = 1)) +
  geom_point(
    size = 1.5, colour = "black", shape = 21, alpha = 0.40,
    position = position_dodge(width = 1)
  ) +
  scale_y_continuous(name = "Pursue-rate [%]") +
  scale_x_discrete(name = "Rw. value [points]", labels = c("5", "10", "50")) +
  scale_fill_manual(values = c("lightskyblue1", "darkblue"), guide = "none") +
  scale_color_manual(values = c("lightskyblue1", "darkblue"), guide = "none") +
  theme(
    aspect.ratio = 1
  )

save_source_data(fig_2_A_data, "fig_2_A.csv")
save_figure(fig_2_A, "fig_2_A.pdf", width = 3, height = 3)



# ==========================================================================
# FIGURE 2B
# ==========================================================================

## mag == 10
fig_2_B_i_data <- d %>%
  filter(mag == 10) %>%
  group_by(env, trial_n_block, mag_factor) %>%
  summarise(m = mean(response, na.rm = TRUE) * 100, .groups = "drop")

fig_2_B <- ggplot(fig_2_B_i_data,
                 aes(x = trial_n_block, y = m, group = env, fill = env)
) +
  stat_summary(fun = "mean", geom = "point",
               size = 2, shape = 21) +
  geom_smooth(method = "lm", alpha = 0.35, color = "black", linewidth = 0.5) +
  scale_y_continuous(name = "Pursue-rate [%]",
                     breaks = seq(0, 100, 20), limits = c(40, 80)) +
  scale_x_continuous(name = "Time-in-env. [Trials]", limits = c(0, 15)) +
  scale_fill_manual(
    values = c("skyblue", "darkblue"),
    labels = c("Poor", "Rich")
  ) +
  facet_grid(
    ~ mag_factor,
    labeller = labeller(
      mag_factor = c("10-points" = "(i) 10 points")
    )
  ) +
  theme(
    legend.position = "none",
    aspect.ratio = 1,
    strip.background = element_blank()
  )

save_source_data(fig_2_B_i_data, "fig_2_B_i.csv")
save_figure(fig_2_B, "fig_2_B_i.pdf", width = 3, height = 3)


## mag == 5
fig_2_B_ii_data <- d %>%
  filter(mag == 5) %>%
  group_by(env, trial_n_block, mag_factor) %>%
  summarise(m = mean(response, na.rm = TRUE) * 100, .groups = "drop")

fig_2_B_ii <- ggplot(fig_2_B_ii_data,
                  aes(x = trial_n_block, y = m, group = env, fill = env)
) +
  stat_summary(fun = "mean", geom = "point",
               size = 2, shape = 21) +
  geom_smooth(method = "lm", alpha = 0.35, color = "black", linewidth = 0.5) +
  scale_y_continuous(
    name = "Pursue-rate [%]",
    breaks = c(0, 10, 20), limits = c(-5, 21)
  ) +
  scale_x_continuous(
    name = "Time-in-env. [Trials]",
    limits = c(0, 15), 
    breaks = seq(0, 15, 5)
  ) +
  scale_fill_manual(
    values = c("skyblue", "darkblue"),
    labels = c("Poor", "Rich")
  ) +
  facet_grid(
    ~ mag_factor,
    labeller = labeller(
      mag_factor = c("5-points" = "(ii) 5 points")
    )
  ) + 
  theme(
    legend.position = "none",
    aspect.ratio = 0.25,
    strip.background = element_blank()
  )

save_source_data(fig_2_B_ii_data, "fig_2_B_ii.csv")
save_figure(fig_2_B_ii, "fig_2_B_ii.pdf", width = 3, height = 3)


## mag == 50
fig_2_B_iii_data <- d %>%
  filter(mag == 50) %>%
  group_by(env, trial_n_block, mag_factor) %>%
  summarise(m = mean(response, na.rm = TRUE) * 100, .groups = "drop")

fig_2_B_iii <- ggplot(fig_2_B_iii_data,
                   aes(x = trial_n_block, y = m, group = env, fill = env)
) +
  stat_summary(fun = "mean", geom = "point",
               size = 2, shape = 21) +
  geom_smooth(method = "lm", alpha = 0.35, color = "black", linewidth = 0.5) +
  scale_y_continuous(
    name = "Pursue-rate [%]",
    breaks = c(90), limits = c(90, 102)
  ) +
  scale_x_continuous(
    name = "Time-in-env. [Trials]",
    limits = c(0, 15), 
    breaks = seq(0, 15, 5)
  ) +
  scale_fill_manual(
    values = c("skyblue", "darkblue"),
    labels = c("Poor", "Rich")
  ) +
  facet_grid(
    ~ mag_factor,
    labeller = labeller(
      mag_factor = c("50-points" = "(iii) 50 points")
    )
  ) + 
  theme(
    legend.position = "none",
    aspect.ratio = 0.25,
    strip.background = element_blank()
  )

save_source_data(fig_2_B_iii_data, "fig_2_B_iii.csv")
save_figure(fig_2_B_iii, "fig_2_B_iii.pdf", width = 3, height = 3)



# ==========================================================================
# FIGURE 2C
# ==========================================================================

acc_rates <- d %>%
  group_by(subject, mag, env) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = c(env, mag),
    names_sep = "_",
    values_from = m
  )

total_rewards <- d %>%
  group_by(subject) %>%
  summarise(total_reward = sum(reward, na.rm = TRUE), .groups = "drop") %>%
  mutate(total_reward = as.numeric(scale(total_reward)))

fig_2_C_data <- acc_rates %>%
  left_join(total_rewards, by = "subject") %>%
  mutate(env_diff = poor_10 - rich_10)

fig_2_C <- ggplot(fig_2_C_data,
                aes(x = env_diff * 100, y = total_reward)
) +
  geom_vline(xintercept = 0, linetype = 3, size = 0.5, color = "red") +
  geom_point(
    shape = 21, color = "black", fill = "lightblue1",
    size = 2
  ) +
  geom_smooth(method = "lm", se = FALSE,
              color = "darkblue", linetype = 2, linewidth = 1) +
  scale_x_continuous(name = TeX("$\\Delta$%Purs. [Poor - Rich]")) +
  scale_y_continuous(name = "Total reward [Z]") +
  theme(
    aspect.ratio = 1.3
  )

save_source_data(fig_2_C_data, "fig_2_C.csv")
save_figure(fig_2_C, "fig_2_C.pdf", width = 3.2, height = 3)