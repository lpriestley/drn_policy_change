## 14_fig6.R ---------------------------------------------------------------
## Figure 6

source("R/00_setup.R")

## -------------------------------------------------------------------------
## Load analysed fMRI timecourse data
## -------------------------------------------------------------------------

tc_path <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
d_tc    <- readRDS(tc_path)

## -------------------------------------------------------------------------
## Figure 6 – DRN–dACC PPI 
## -------------------------------------------------------------------------

fig_6_data <- d_tc %>%
  filter(
    glm      == "GLM_04.3a",
    contrast == "contr_1",
    region   == "DRN"
  ) %>%
  mutate(region_label = "DRN-dACC")

fig_6 <- ggplot(
  fig_6_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se)
) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.10) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.10) +
  geom_ribbon(aes(fill = contrast, group = contrast),
              alpha = 0.60, fill = "skyblue1") +
  geom_line(aes(group = contrast), colour = "black", linewidth = 0.60) +
  labs(
    x = "Time [s]",
    y = TeX("$\\beta$-coeff. [a.u.]")
  ) +
  scale_y_continuous(
    breaks   = seq(-0.06, 0.03, 0.03),
    position = "right"
  ) +
  scale_x_continuous(
    limits = c(-2, 8),
    breaks = c(0, 4, 8)
  ) +
  theme_fig() +
  facet_grid(~ region_label) +
  theme(
    aspect.ratio    = 1,
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"),
    legend.position  = "none"
  )

save_source_data(fig_6_data, "fig_6.csv")
save_figure(fig_6, "fig_6.pdf",
            width = 3, height = 3)
