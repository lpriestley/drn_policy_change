## 08_fig1.R ---------------------------------------------------------------
##  Figure 1

source("R/00_setup.R")
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

# ==========================================================================
# FIGURE 1B: Environment-specific frequencies
# ==========================================================================

fig_1_B_data <- d %>%
  filter(subject == 103) %>%
  mutate(
    p_low  = ifelse(env_bin ==  1, 1/6 + 0.001,
                    ifelse(env_bin == -1, 3/6 + 0.001, NA)),
    p_mod  = 2/6,
    p_high = ifelse(env_bin == -1, 1/6 - 0.001,
                    ifelse(env_bin ==  1, 3/6 - 0.001, NA))
  ) %>%
  select(env_bin, p_low, p_mod, p_high) %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "option",
    values_to = "p_opt"
  ) %>%
  mutate(
    option_factor = factor(option,
                           levels = c("p_low", "p_mod", "p_high"),
                           labels = c("5", "10", "50")
    ),
    env_factor = factor(env_bin, levels = c(-1, 1), labels = c("Poor", "Rich"))
  )

fig_1_B <- ggplot(fig_1_B_data,
                aes(x = option_factor, y = p_opt, fill = option_factor)
) +
  stat_summary(fun = "mean", geom = "bar",
               color = "black", alpha = 0.65, width = 0.65) +
  scale_x_discrete(name = "Rw-magnitude") +
  scale_y_continuous(name = "p(Option)", limits = c(0, 0.55)) +
  scale_fill_manual(
    values = c("#CD7F32", "#C0C0C0", "#DAA520"),
    guide = "none"
  ) +
  facet_grid(~ env_factor) +
  coord_cartesian(ylim = c(0.15, 0.55)) +
  theme_fig() +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank()
  )

save_source_data(fig_1_B_data, "fig_1_B.csv")
save_figure(fig_1_B, "fig_1_B.pdf", width = 4.5, height = 3.5)



# ==========================================================================
# FIGURE 1C: Example environments (subject 106)
# Programmatic transitions + original cropping
# ==========================================================================

example_session <- d %>%
  filter(subject == 106) %>%
  mutate(env = factor(env, levels = c("poor", "rich")))

env_transitions <- which(diff(example_session$env_bin) != 0) + 1

fig_1_C_data <- example_session %>%
  filter(block_count %in% c(9, 10)) %>%
  mutate(trial = 1:n())

fig_1_C <- ggplot() +
  geom_vline(
    xintercept = 14,
    color      = "black",
    linetype   = 3,
    size       = 0.50
  ) +
  geom_point(
    data  = fig_1_C_data,
    aes(x = trial, y = mag, fill = factor(mag)),
    shape = 21,
    alpha = 0.65,
    size  = 4
  ) +
  geom_smooth(
    data   = fig_1_C_data,
    aes(x = trial, y = mu_val, linetype = "Ave.-value"),  # mapped for legend
    color    = "red",
    alpha    = 0.65,
    linewidth = 0.5,
    method   = "loess",
    se       = FALSE,
    span     = 0.4
  ) +
  scale_x_continuous(
    name   = "Time [trials]",
    breaks = seq(0, 40, 10)
  ) +
  scale_y_continuous(
    name   = "Rw.-value",
    breaks = seq(0, 40, 20)
  ) +
  # Hide fill legend (points)
  scale_fill_manual(
    values = c("#CD7F32", "#C0C0C0", "#DAA520"),
    guide  = "none"
  ) +
  scale_linetype_manual(
    name   = NULL,                     
    values = c("Ave.-value" = 2)       
  ) +
  coord_cartesian(ylim = c(-10, 60)) +
  theme(
    legend.position      = c(0.5, 1),
    legend.justification = c(0.5, 1),
    legend.direction     = "horizontal",
    aspect.ratio         = 0.45
  )

save_source_data(fig_1_C_data, "fig_1_C.csv")
save_figure(fig_1_C, "fig_1_C.pdf", width = 5, height = 3)


