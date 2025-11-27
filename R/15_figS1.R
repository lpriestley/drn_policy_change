## 15_figS1.R ---------------------------------------------------------------
## Supplementary Figure S1

source("R/00_setup.R")
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

# ==========================================================================
# FIGURE S1B
# ==========================================================================

fig_S1_B_data <- d %>%
  group_by(subject, env_bin, block_count) %>%
  summarise(max_trials = max(trial_n_block, na.rm = TRUE), .groups = "drop") %>%
  filter(max_trials > 1) %>%   # preserves original filtering
  group_by(subject, env_bin) %>%
  summarise(m_trials = mean(max_trials, na.rm = TRUE), .groups = "drop") %>%
  mutate(env_factor = factor(env_bin,
                             levels = c(-1, 1),
                             labels = c("Poor", "Rich"))
  )

fig_S1_B <- ggplot(fig_S1_B_data,
                 aes(x = env_factor, y = m_trials)) +
  geom_point(
    aes(fill = env_factor),
    size = 2,
    shape = 21,
    alpha = 0.40
  ) +
  stat_summary(
    aes(fill = env_factor),
    fun.data = 'mean_se', 
    geom = 'pointrange', 
    shape = 22, 
    size = 1, 
    alpha = 1) + 
  scale_fill_manual(values = c("skyblue2", "darkblue")) +
  scale_y_continuous(name = "Ave. block length [trials]") +
  scale_x_discrete(name = "Environment-type") +
  theme_fig() +
  theme(
    aspect.ratio = 2,
    legend.position = "none"
  )

save_source_data(fig_S1_B_data, "fig_S1_B.csv")
save_figure(fig_S1_B, "fig_S1_B.pdf", width = 3, height = 3.5)



# ==========================================================================
# FIGURE S1C
# ==========================================================================

fig_S1_C_data <- d %>%
  group_by(subject, mag, block_count, env) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  group_by(subject, mag, env) %>%
  summarise(mu = mean(m, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = env, values_from = mu) %>%
  mutate(env_diff = (poor - rich) * 100)   # express as % for plotting

fig_S1_C <- ggplot(fig_S1_C_data,
                 aes(x = factor(mag), y = env_diff, fill = factor(mag))) +
  geom_point(shape = 21, alpha = 0.40, size = 2) +
  stat_summary(fun.data = "mean_se", geom = "pointrange",
               color = "black", shape = 21, alpha = 0.75, size = 1) +
  stat_summary(fun = "mean", geom = "line",
               color = "blue", linetype = 2, group = 1) +
  scale_fill_manual(values = c("azure1", "skyblue2", "darkblue")) +
  scale_y_continuous(name = TeX("$\\Delta$%Purs. [Poor - Rich]"),
                     limits = c(-35, 60)) +
  scale_x_discrete(name = "Rw.-mag.[points]") +
  theme_fig() +
  theme(aspect.ratio = 1.75,
        legend.position = "none")

save_source_data(fig_S1_C_data, "fig_S1_C.csv")
save_figure(fig_S1_C, "fig_S1_C.pdf", width = 4.0, height = 3.5)



# ==========================================================================
# FIGURE S1D
# ==========================================================================

fig_S1_D_data <- d %>%
  drop_na(mu_val_z) %>%
  mutate(mu_bin = cut(mu_val_z,
                      breaks = unique(quantile(mu_val_z,
                                               probs = seq(0, 1, length.out = 11))),
                      labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_bin) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject, mag, mu_val_graph) %>%
  summarise(pursue_rate = 100 * mean(response, na.rm = TRUE),
            .groups = "drop")

fig_S1_D <- ggplot(fig_S1_D_data,
                 aes(x = mu_val_graph, y = pursue_rate,
                     fill = factor(mag))) +
  stat_summary(fun = "mean", geom = "point",
               size = 2.5, shape = 21,
               alpha = 0.75, color = "black") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.5, alpha = 0.40) +
  scale_fill_manual(values = c("azure1", "skyblue2", "darkblue")) +
  scale_x_continuous(name = "Ave. value [Z]") +
  scale_y_continuous(name = "Pursue-rate [%]",
                     breaks = seq(0, 100, 20)) +
  facet_wrap(~ factor(mag, labels = c("(i) 5-point", "(ii) 10-point", "(iii) 50-point")),
             scales = "free") +
  theme_fig() +
  theme(aspect.ratio = 1.2,
        strip.background = element_blank(), 
        legend.position = "none")

save_source_data(fig_S1_D_data, "fig_S1_D.csv")
save_figure(fig_S1_D, "fig_S1_D.pdf", width = 9, height = 3.5)
