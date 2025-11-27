## 17_figS5.R ---------------------------------------------------------------
## Figure S5

source("R/00_setup.R")

## -------------------------------------------------------------------------
## Load data
## -------------------------------------------------------------------------

peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")
d_peak    <- readRDS(peak_path)

tc_path <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
d_tc    <- readRDS(tc_path)

## ========================================================================
## Figure  S5A
## ========================================================================

fig_S5_A_data <- d_peak %>%
  filter(
    region   %in% c("DRN", "fourth_ventricle", "MRN"),
    contrast == "contr1",
    GLM      == "GLM_04.1_supp"
  )

fig_S5_A <- ggplot(
  fig_S5_A_data,
  aes(x = region, y = peak, fill = region)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    color    = "black",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0)
  ) +
  geom_point(alpha = 0.50, shape = 21, size = 2.5) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_discrete(
    name   = "Region",
    labels = c("DRN", "Vent.", "MRN")
  ) +
  scale_fill_manual(
    name   = "Region",
    values = c("DRN" = "deeppink4",
               "fourth_ventricle" = "white",
               "MRN" = "white")
  ) +
  scale_colour_manual(
    name   = "Region",
    values = c("DRN" = "deeppink4",
               "fourth_ventricle" = "black",
               "MRN" = "black")
  ) +
  theme_fig() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank()
  )

save_source_data(fig_S5_A_data, "fig_S5_A.csv")
save_figure(fig_S5_A, "fig_S5_A.pdf",
            width = 4.5, height = 3.0)

## ========================================================================
## Figure S5B
## ========================================================================


fig_S5_B_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    # Only GLM_0S5.1_win* or GLM_0S5.2_win*
    str_detect(GLM, "^GLM_0S5\\.(1|2)_win")
  ) %>%
  mutate(
    switch_type = case_when(
      str_detect(GLM, "^GLM_0S5\\.1_win") ~ "congruent_poor",
      str_detect(GLM, "^GLM_0S5\\.2_win") ~ "congruent_rich",
      TRUE                                ~ NA_character_
    ),
    window_start = str_extract(GLM, "(?<=win)\\d+") %>% as.integer()
  ) %>%
  filter(!is.na(switch_type)) %>%
  group_by(switch_type, window_start) %>%
  summarise(
    mean_peak = mean(peak, na.rm = TRUE),
    se        = sd(peak,  na.rm = TRUE) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  ) %>%
  filter(window_start <= 10) %>%
  group_by(switch_type) %>%
  mutate(
    n_comparisons = n(),
    alpha_bonf    = 0.05 / n_comparisons,
    z_crit        = qnorm(1 - alpha_bonf / 2),
    ci_bonf       = z_crit * se
  ) %>%
  ungroup() %>%
  mutate(
    significant = (mean_peak - ci_bonf > 0) | (mean_peak + ci_bonf < 0),
    switch_type = factor(
      switch_type,
      levels = c("congruent_poor", "congruent_rich"),
      labels = c("(i) Congruent poor", "(ii) Congruent rich")
    )
  )

fig_S5_B <- ggplot(
  fig_S5_B_data,
  aes(x = window_start, y = mean_peak,
      fill = switch_type, color = switch_type)
) +
  geom_hline(yintercept = 0, color = "red", linetype = 2) +
  geom_line(color = "black", alpha = 0.30) +
  geom_errorbar(
    aes(ymin = mean_peak - ci_bonf,
        ymax = mean_peak + ci_bonf,
        alpha = significant),
    width = 0,
    color = "black"
  ) +
  geom_point(
    aes(alpha = significant),
    size  = 4,
    color = "black",
    shape = 21
  ) +
  scale_x_continuous(
    name   = "Time-window",
    breaks = c(1, 5, 10),
    labels = c("(t=1)\n:(t=6)", "(t=5)\n:(t=10)", "(t=10)\n:(t=15)")
  ) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    breaks = seq(-0.10, 0.10, 0.10)
  ) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3)) +
  scale_fill_manual(
    name   = "Switch-type",
    values = c("(i) Congruent poor" = "deeppink4",
               "(ii) Congruent rich" = "pink")
  ) +
  scale_color_manual(
    name   = "Switch-type",
    values = c("(i) Congruent poor" = "deeppink4",
               "(ii) Congruent rich" = "pink")
  ) +
  facet_wrap(~ switch_type, nrow = 1) +
  theme_fig() +
  theme(
    legend.position = 'none',
    axis.title.x    = element_text(margin = margin(t = 5)),
    aspect.ratio    = 0.50,
    strip.background = element_blank(),
    panel.spacing.x  = unit(2, "lines")
  )

save_source_data(fig_S5_B_data, "fig_S5_B.csv")
save_figure(fig_S5_B, "fig_S5_B.pdf",
            width = 9, height = 3.0)

## ========================================================================
## Figure S5C 
## ========================================================================

fig_S5_C_i_data <- d_tc %>%
  filter(
    glm      == "GLM_0S5.3",
    contrast == "contr_1",
    region   == "DRN"
  )

fig_S5_C_i <- ggplot(
  fig_S5_C_i_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_ribbon(alpha = 0.65, fill = "grey90") +
  geom_line(linewidth = 0.60, colour = "black") +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.25) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(
    breaks = c(-0.05, 0.00, 0.05, 0.10),
    limits = c(-0.05, 0.06)
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ factor(region, labels = c("Trial-after"))) +
  theme(
    legend.position = "none",
    aspect.ratio    = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines")
  )

save_source_data(fig_S5_C_i_data, "fig_S5_C_i.csv")
save_figure(fig_S5_C_i, "fig_S5_C_i.pdf",
            width = 3.5, height = 3.5)

fig_S5_C_ii_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM      == "GLM_0S5.3"
  )

fig_S5_C_ii <- ggplot(
  fig_S5_C_ii_data,
  aes(x = region, y = peak, fill = region)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    color    = "black",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0)
  ) +
  geom_point(alpha = 0.50, shape = 21, fill = NA, size = 2.5) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    labels = scales::number_format(accuracy = 0.01), 
    position = "right" 
  ) +
  scale_x_discrete(
    name   = "ROI",
    labels = c("DRN")
  ) +
  scale_fill_manual(
    name   = "ROI",
    values = c("DRN" = "white")
  ) +
  theme_fig() +
  theme(
    legend.position  = "none",
    aspect.ratio     = 3,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) + 
  border(color = 'black')

save_source_data(fig_S5_C_ii_data, "fig_S5_C_ii.csv")
save_figure(fig_S5_C_ii, "fig_S5_C_ii.pdf",
            width = 2.5, height = 2.5)

## ========================================================================
## Figure S5D
## ========================================================================

fig_S5_D_i_data <- d_tc %>%
  filter(
    glm      == "GLM_0S5.4",
    contrast == "contr_1",
    region   == "DRN"
  )

fig_S5_D_i <- ggplot(
  fig_S5_D_i_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_ribbon(alpha = 0.65, fill = "grey90") +
  geom_line(linewidth = 0.60, colour = "black") +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.25) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.25) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(
    breaks = c(-0.05, 0.00, 0.05, 0.10),
    limits = c(-0.05, 0.06)
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ factor(region, labels = c("Encounter-after"))) +
  theme(
    legend.position = "none",
    aspect.ratio    = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines")
  )

save_source_data(fig_S5_D_i_data, "fig_S5_D_i.csv")
save_figure(fig_S5_D_i, "fig_S5_D_i.pdf",
            width = 3.5, height = 3.5)

fig_S5_D_ii_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM      == "GLM_0S5.4"
  )

fig_S5_D_ii <- ggplot(
  fig_S5_D_ii_data,
  aes(x = region, y = peak, fill = region)
) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = "red") +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    color    = "black",
    width    = 0.5,
    alpha    = 0.75,
    position = position_dodge2(width = 0)
  ) +
  geom_point(alpha = 0.50, shape = 21, fill = NA, size = 2.5) +
  scale_y_continuous(
    name   = TeX("$\\beta$-coeff [a.u.]"),
    labels = scales::number_format(accuracy = 0.01), 
    position = 'right'
  ) +
  scale_x_discrete(
    name   = "ROI",
    labels = c("DRN")
  ) +
  scale_fill_manual(
    name   = "ROI",
    values = c("DRN" = "white")
  ) +
  theme_fig() +
  theme(
    legend.position  = "none",
    aspect.ratio     = 3,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"), 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  ) + 
  border(color = 'black')

save_source_data(fig_S5_D_ii_data, "fig_S5_D_ii.csv")
save_figure(fig_S5_D_ii, "fig_S5_D_ii.pdf",
            width = 2.5, height = 2.5)
