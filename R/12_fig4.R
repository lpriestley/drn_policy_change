## 12_fig4.R ---------------------------------------------------------------
## Figure 4 
source("R/00_setup.R")

library(dplyr)
library(ggplot2)
library(latex2exp)
library(stringr)

## -------------------------------------------------------------------------
## Load data
## -------------------------------------------------------------------------

tc_path   <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")

d_tc   <- readRDS(tc_path)
d_peak <- readRDS(peak_path)

out_dir <- here::here("results", "fmriTC")  # optional, if you want to use it

## -------------------------------------------------------------------------
## FIGURE 4B
## -------------------------------------------------------------------------

fig_4_B_data <- d_peak %>%
  filter(
    region %in% c("DRN", "MBD", "LC", "MS", "HB"),
    contrast == "contr1",
    GLM %in% c("GLM_04.1", "GLM_04.1b")
  ) %>%
  mutate(region = factor(region, levels = c("DRN", "MBD", "MS", "LC", "HB")))

fig_4_B <- ggplot(fig_4_B_data, aes(x = region, y = peak, fill = GLM, group = GLM)) +
  stat_summary(
    fun   = "mean",
    geom  = "bar",
    colour = "black",
    width  = 0.5,
    position = position_dodge(width = 0.8, preserve = "total")
  ) +
  geom_point(
    shape = 21,
    alpha  = 0.25,
    position = position_dodge(width = 0.8, preserve = "total"),
    size = 2
  ) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  labs(
    x = "Region",
    y = TeX("$\\beta$-coeff [a.u.]")
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("GLM_04.1" = "deeppink4", "GLM_04.1b" = "azure1"),
    labels = NULL
  ) +
  theme_fig() +
  theme(
    legend.position = "none",
    aspect.ratio    = 0.2,
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
save_source_data(fig_4_B_data, "fig_4_B.csv")
save_figure(fig_4_B, "fig_4_B.pdf",
            width = 8, height = 3)

# ==========================================================================
# FIGURE 4C 
# ==========================================================================

fig_4_C_data <- d_tc %>%
  filter(
    glm %in% c("GLM_04.1", "GLM_04.1b"),
    contrast == "contr_1",
    region   == "DRN"
  ) %>%
  mutate(
    facet_label = dplyr::case_when(
      glm == "GLM_04.1"  ~ "(i) Policy switch [DRN]",
      glm == "GLM_04.1b" ~ "(ii) Action switch [DRN]"
    )
  )

fig_4_C <- ggplot(
  fig_4_C_data,
  aes(x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = facet_label)
) +
  geom_ribbon(alpha = 0.65) +
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
  scale_fill_manual(
    guide  = "none",
    values = c(
      "(i) Policy switch [DRN]"   = "deeppink4",
      "(ii) Action switch [DRN]"   = "azure2"
    )
  ) +
  theme_fig() +
  facet_wrap(~ facet_label, nrow = 2) +
  theme(
    legend.position  = "none",
    strip.clip       = "off",                  
    aspect.ratio     = 1.5,
    strip.text       = element_text(
      size = 14,
      hjust = 0.5,                              
      vjust = 0.5
    ),
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.2, "lines"),
    panel.spacing.y  = unit(1.2, "lines") 
  )

save_source_data(fig_4_C_data, "fig_4_C.csv")
save_figure(fig_4_C, "fig_4_C.pdf", width = 3.15, height = 6)

## -------------------------------------------------------------------------
## FIGURE 4D: Further analysis of policy switches in DRN
## -------------------------------------------------------------------------

## ---- helper  function ----------------------------------------------------

plot_DRN_pair <- function(data, fill_vals, show_y_label = TRUE) {
  ggplot(data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 1, colour = "red") +
    stat_summary(
      fun      = "mean",
      geom     = "bar",
      colour   = "black",
      width    = 0.5,
      alpha    = 0.75,
      position = position_dodge2(width = 0)
    ) +
    geom_point(alpha = 0.50, shape = 21, size = 2.5) +
    scale_y_continuous(
      name   = if (show_y_label) TeX("$\\beta$-coeff [a.u.]") else NULL,
      labels = scales::number_format(accuracy = 0.01),
      limits = c(-0.30, 0.50),
      breaks = seq(-0.25, 0.50, 0.25)
    ) +
    scale_x_discrete(labels = NULL) +
    scale_fill_manual(values = fill_vals) +
    theme_fig() +
    theme(
      axis.title.y      = if (show_y_label) element_text(size = 18) else element_blank(),
      axis.title.x      = element_blank(),
      axis.ticks.x      = element_blank(),      # no x ticks
      axis.line.x       = element_blank(),      # no x axis line
      axis.text         = element_text(size = 14),
      axis.text.y       = if (show_y_label) element_text(size = 14) else element_blank(),
      legend.position   = "none",
      legend.background = element_blank(),
      aspect.ratio      = 2.5,
      strip.text.x      = element_text(size = 18),
      strip.background  = element_blank(), 
      panel.spacing.x   = unit(0.2, "lines")
    )
}

## -------------------------- Fig 4D-i -------------------------------------

fig_4_D_i_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM %in% c("GLM_04.2a", "GLM_04.2c")
  ) %>%
  mutate(facet_label = "(i)")

fig_4_D_i <- plot_DRN_pair(
  fig_4_D_i_data,
  fill_vals     = c("GLM_04.2a" = "deeppink4", "GLM_04.2c" = "pink"),
  show_y_label  = TRUE
) +
  facet_wrap(~ facet_label) +
  theme(strip.text.x = element_text(size = 18, hjust = 0))  # left-align label

save_source_data(fig_4_D_i_data, "fig_4_D_i.csv")
save_figure(fig_4_D_i, "fig_4_D_i.pdf",
            width = 2.5, height = 4)

## -------------------------- Fig 4D-ii ------------------------------------

fig_4_D_ii_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM %in% c("GLM_04.2a", "GLM_04.2b")
  ) %>%
  mutate(facet_label = "(ii)")

fig_4_D_ii <- plot_DRN_pair(
  fig_4_D_ii_data,
  fill_vals     = c("GLM_04.2a" = "deeppink4", "GLM_04.2b" = "azure2"),
  show_y_label  = TRUE    # keep TRUE so panel size matches 4D-i
) +
  facet_wrap(~ facet_label) +
  theme(strip.text.x = element_text(size = 18, hjust = 0))

save_source_data(fig_4_D_ii_data, "fig_4_D_ii.csv")
save_figure(fig_4_D_ii, "fig_4_D_ii.pdf",
            width = 2.5, height = 4)

## -------------------------- Fig 4D-iii -----------------------------------

fig_4_D_iii_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM %in% c("GLM_04.2c", "GLM_04.2d")
  ) %>%
  mutate(facet_label = "(iii)")

fig_4_D_iii <- plot_DRN_pair(
  fig_4_D_iii_data,
  fill_vals     = c("GLM_04.2c" = "pink", "GLM_04.2d" = "azure2"),
  show_y_label  = TRUE
) +
  facet_wrap(~ facet_label) +
  theme(strip.text.x = element_text(size = 18, hjust = 0))

save_source_data(fig_4_D_iii_data, "fig_4_D_iii.csv")
save_figure(fig_4_D_iii, "fig_4_D_iii.pdf",
            width = 2.5, height = 4)

## -------------------------- Fig 4D-iv ------------------------------------

fig_4_D_iv_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM %in% c("GLM_04.2a", "GLM_04.2a_low")
  ) %>%
  mutate(facet_label = "(iv)")

fig_4_D_iv <- plot_DRN_pair(
  fig_4_D_iv_data,
  fill_vals     = c("GLM_04.2a" = "deeppink4", "GLM_04.2a_low" = "azure2"),
  show_y_label  = TRUE
) +
  facet_wrap(~ facet_label) +
  theme(strip.text.x = element_text(size = 18, hjust = 0))

save_source_data(fig_4_D_iv_data, "fig_4_D_iv.csv")
save_figure(fig_4_D_iv, "fig_4_D_iv.pdf",
            width = 2.5, height = 4)

## -------------------------- Fig 4D-v -------------------------------------

fig_4_D_v_data <- d_peak %>%
  filter(
    region   == "DRN",
    contrast == "contr1",
    GLM %in% c("GLM_04.2c", "GLM_04.2c_low")
  ) %>%
  mutate(facet_label = "(v)")

fig_4_D_v <- plot_DRN_pair(
  fig_4_D_v_data,
  fill_vals     = c("GLM_04.2c" = "pink1", "GLM_04.2c_low" = "azure2"),
  show_y_label  = TRUE
) +
  facet_wrap(~ facet_label) +
  theme(strip.text.x = element_text(size = 18, hjust = 0))

save_source_data(fig_4_D_v_data, "fig_4_D_v.csv")
save_figure(fig_4_D_v, "fig_4_D_v.pdf",
            width = 2.5, height = 4)
