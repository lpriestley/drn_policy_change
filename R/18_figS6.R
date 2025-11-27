## 18_figS6.R ---------------------------------------------------------------
## Figure S6 

source("R/00_setup.R")

## -------------------------------------------------------------------------
## Load data
## -------------------------------------------------------------------------

peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")
d_peak    <- readRDS(peak_path)

tc_path <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
d_tc    <- readRDS(tc_path)

beh_path <- here::here("data", "behaviour", "clean_lf.rds")
d_beh    <- readRDS(beh_path)

## ========================================================================
## Figure S6A 
## ========================================================================

fig_S6_A_i_data <- d_tc %>%
  filter(
    glm      == "GLM_0S6.1",
    contrast == "contr_1",
    region   == "MBD"
  )

fig_S6_A_i <- ggplot(
  fig_S6_A_i_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.10) +
  geom_ribbon(alpha = 0.75, fill = "grey90") +
  geom_line(colour = "black", linewidth = 0.60) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff [a.u.]")
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ factor(region, levels = "MBD", labels = "Pursue-vs-reject")) +
  theme(
    legend.position  = "none",
    aspect.ratio     = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines")
  )

save_source_data(fig_S6_A_i_data, "fig_S6_A_i.csv")
save_figure(fig_S6_A_i, "fig_S6_A_i.pdf",
            width = 3, height = 3)

fig_S6_A_ii_data <- d_peak %>%
  filter(
    region   == "MBD",
    contrast == "contr1",
    GLM      == "GLM_0S6.1"
  )

fig_S6_A_ii <- ggplot(
  fig_S6_A_ii_data,
  aes(x = GLM, y = peak, fill = GLM)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    colour   = "black",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0),
    fill     = "white"
  ) +
  geom_point(
    alpha = 0.50,
    shape = 21,
    colour = "black",
    size = 2.5,
    fill = "white"
  ) +
  scale_y_continuous(
    name     = TeX("$\\beta$-coeff [a.u.]"),
    position = "right",
    labels   = scales::number_format(accuracy = 0.01), 
    breaks = seq(-0.10, 0.20, 0.10)
  ) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(
    axis.title.x    = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    aspect.ratio    = 2.5,
    strip.background = element_blank()
  ) + 
  border(color = 'black')

save_source_data(fig_S6_A_ii_data, "fig_S6_A_ii.csv")
save_figure(fig_S6_A_ii, "fig_S6_A_ii.pdf",
            width = 2.0, height = 2)

## ========================================================================
## Figure S6B
## ========================================================================

# Compute value-difference regressor
d_beh <- d_beh %>%
  mutate(mu_val_pe_z = mag_z - mu_val_z)

fig_S6_B_data <- d_beh %>%
  drop_na(mu_val_pe_z) %>%
  mutate(
    mu_val_pe_ind = cut(
      mu_val_pe_z,
      breaks = unique(quantile(mu_val_pe_z,
                               probs = seq(0, 1, length.out = 11))),
      labels = FALSE,
      include.lowest = TRUE
    )
  ) %>%
  group_by(mu_val_pe_ind) %>%
  mutate(mu_val_pe_graph = mean(mu_val_pe_z, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject, mag, prev_policy, mu_val_pe_graph) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop")

fig_S6_B <- ggplot(fig_S6_B_data, aes(x = mu_val_pe_graph, y = m * 100)) +
  geom_smooth(method = "lm",
              fill = "grey90",
              colour = "black",
              linewidth = 0.5) +
  stat_summary(
    fun = "mean",
    geom     = "point",
    fill     = "grey90",
    shape    = 21,
    colour   = "black",
    size     = 2.5,
    linewidth = 0.5
  ) +
  theme_fig() +
  scale_y_continuous(
    name   = "Pursue-rate [%]",
    breaks = seq(0, 100, 20)
  ) +
  scale_x_continuous(
    name = "Value-difference [Z]"
  ) +
  theme(
    aspect.ratio    = 1.2,
    legend.position = "none"
  )

save_source_data(fig_S6_B_data, "fig_S6_B.csv")
save_figure(fig_S6_B, "fig_S6_B.pdf",
            width = 3, height = 3)

## ========================================================================
## Figure S6C 
## ========================================================================

fig_S6_C_i_data <- d_tc %>%
  filter(
    glm      == "GLM_0S6.2",
    contrast == "contr_1",
    region   == "MBD"
  )

fig_S6_C_i <- ggplot(
  fig_S6_C_i_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.10) +
  geom_ribbon(alpha = 0.60, fill = "grey90") +
  geom_line(colour = "black", linewidth = 0.60) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff [a.u.]")
  ) +
  scale_y_continuous(
    breaks = seq(-0.03, 0.06, 0.03)
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ factor(region, levels = "MBD", labels = "Value-difference")) +
  theme(
    legend.position  = "none",
    aspect.ratio     = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines")
  )

save_source_data(fig_S6_C_i_data, "fig_S6_C_i.csv")
save_figure(fig_S6_C_i, "fig_S6_C_i.pdf",
            width = 3, height = 3)


fig_S6_C_ii_data <- d_peak %>%
  filter(
    region   == "MBD",
    contrast == "contr1",
    GLM      == "GLM_0S6.2"
  )

fig_S6_C_ii <- ggplot(
  fig_S6_C_ii_data,
  aes(x = GLM, y = peak, fill = GLM)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    colour   = "black",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0),
    fill     = "white"
  ) +
  geom_point(
    alpha = 0.50,
    shape = 21,
    colour = "black",
    size = 2.5,
    fill = "white"
  ) +
  scale_y_continuous(
    name     = TeX("$\\beta$-coeff [a.u.]"),
    labels   = scales::number_format(accuracy = 0.01),
    breaks   = seq(-0.10, 0.20, 0.10), 
    position = "right"
  ) +
  scale_x_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(
    axis.title.x    = element_blank(),
    axis.ticks.x    = element_blank(),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    aspect.ratio    = 2.5,
    strip.background = element_blank()
  ) + 
  border(color = 'black')

save_source_data(fig_S6_C_ii_data, "fig_S6_C_ii.csv")
save_figure(fig_S6_C_ii, "fig_S6_C_ii.pdf",
            width = 2.0, height = 2.0)

