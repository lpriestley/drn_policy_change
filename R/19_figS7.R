## 19_figS7.R ---------------------------------------------------------------
## Figure S7 

source("R/00_setup.R")

## -------------------------------------------------------------------------
## Load data
## -------------------------------------------------------------------------

peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")
d_peak    <- readRDS(peak_path)

tc_path <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
d_tc    <- readRDS(tc_path)

## ========================================================================
## Figure S7A
## ========================================================================

fig_S7A_data <- d_peak %>%
  filter(
    region   %in% c("DRN", "MBD", "MS", "HB", "LC"),
    contrast == "contr1",
    GLM      == "GLM_0S7.1"
  )

fig_S7A <- ggplot(
  fig_S7A_data,
  aes(x = region, y = peak)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    colour   = "black",
    fill     = "grey90",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0)
  ) +
  geom_point(
    alpha  = 0.50,
    shape  = 21,
    fill   = "white",
    size   = 2.5
  ) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_discrete(name = "ROI") +
  theme_fig() +
  theme(
    aspect.ratio    = 0.75,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

save_source_data(fig_S7A_data, "fig_S7_A.csv")
save_figure(fig_S7A, "fig_S7_A.pdf",
            width = 4.5, height = 3.0)

## ========================================================================
## Figure S7B 
## ========================================================================

fig_S7B_data <- d_tc %>%
  filter(
    glm      == "GLM_0S7.1",
    contrast == "contr_1",
    region   == "MS"
  )

fig_S7B <- ggplot(
  fig_S7B_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.10) +
  geom_ribbon(fill = "grey90", alpha = 0.60) +
  geom_line(aes(group = contrast), colour = "black", linewidth = 0.60) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff [a.u.]")
  ) +
  scale_y_continuous(
    breaks = c(-0.05, 0.00, 0.05, 0.10)
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ factor(region, levels = "MS", labels = "MS")) +
  theme(
    legend.position  = "none",
    aspect.ratio     = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines")
  )

save_source_data(fig_S7B_data, "fig_S7_B.csv")
save_figure(fig_S7B, "fig_S7_B.pdf",
            width = 3, height = 3)

## ========================================================================
## Figure S7C
## ========================================================================

exploratory_peaks <- d_peak %>%
  filter(
    region   %in% c("DRN", "MS"),
    contrast == "contr1",
    GLM      == "GLM_0S7.1"
  )

environment_peaks <- d_peak %>%
  filter(
    region   %in% c("DRN", "MS"),
    contrast == "contr2",
    GLM      == "GLM_0S7.2"
  )

fig_S7C_data <- bind_rows(exploratory_peaks, environment_peaks)

fig_S7C <- ggplot(
  fig_S7C_data,
  aes(x = region, y = peak)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    colour   = "black",
    fill = 'grey90',
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0)
  ) +
  geom_point(
    alpha = 0.50,
    shape = 21,
    fill  = 'white',
    size  = 2.5
  ) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    labels = scales::number_format(accuracy = 0.01),
    breaks = seq(-0.15, 0.30, 0.15)
  ) +
  scale_x_discrete(name = "ROI") +
  theme_fig() +
  facet_grid(
    ~ factor(
      GLM,
      levels = c("GLM_0S7.1", "GLM_0S7.2"),
      labels = c("Exploratory", "Env.-specific")
    )
  ) +
  theme(
    aspect.ratio    = 1.2,
    strip.background = element_blank(),
    panel.spacing.x  = unit(1, "lines"), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

save_source_data(fig_S7C_data, "fig_S7_C.csv")
save_figure(fig_S7C, "fig_S7_C.pdf",
            width = 6.0, height = 3.0)

