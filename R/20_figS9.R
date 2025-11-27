## 20_figS9_RSA_supp.R ----------------------------------------------------
## Supplementary Figure S9 

source("R/00_setup.R")

## ------------------------------------------------------------------------
## Load RSA data (cosine + Pearson)
## ------------------------------------------------------------------------

rsa_cosine_path  <- here::here("data", "fMRI", "rsa", "rsa_cosine_clean_lf.rds")
rsa_pearson_path <- here::here("data", "fMRI", "rsa", "rsa_pearson_clean_lf.rds")

d_rsa_cosine  <- readRDS(rsa_cosine_path)
d_rsa_pearson <- readRDS(rsa_pearson_path)

# Ensure 'constant' column is present for faceting
d_rsa_cosine$constant  <- 1
d_rsa_pearson$constant <- 1

## ========================================================================
## Figure S9A
## ========================================================================

## S9A-i – dACC H1 vs behaviour (COSINE)

fig_S9_A_i_data <- d_rsa_cosine %>%
  dplyr::select(subject, p_accept_10_diff, dACC_H1, constant)

fig_S9_A_i <- ggplot(
  fig_S9_A_i_data,
  aes(x = p_accept_10_diff * 100, y = dACC_H1)
) +
  geom_point(shape = 21, fill = "red2", size = 2.5) +
  geom_smooth(
    method    = "lm",
    colour    = "red4",
    linewidth = 0.65,
    alpha     = 0.50,
    se        = FALSE,
    linetype  = 2
  ) +
  scale_y_continuous(
    name   = "H1(d3 - d1)",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_continuous(
    name = TeX("$\\Delta$%Purs. [Poor - Rich]")
  ) +
  theme_fig() +
  facet_grid(~ factor(constant, labels = "dACC")) +
  theme(
    aspect.ratio     = 1,
    strip.background = element_blank()
  )

save_source_data(fig_S9_A_i_data, "fig_S9_A_i.csv")
save_figure(fig_S9_A_i, "fig_S9_A_i.pdf",
            width = 3.5, height = 3.0)

## S9A-ii – AI H1 vs behaviour (COSINE)

fig_S9_A_ii_data <- d_rsa_cosine %>%
  dplyr::select(subject, p_accept_10_diff, AI_H1, constant)

fig_S9_A_ii <- ggplot(
  fig_S9_A_ii_data,
  aes(x = p_accept_10_diff * 100, y = AI_H1)
) +
  geom_point(shape = 21, fill = "lightblue1", size = 2.5) +
  geom_smooth(
    method    = "lm",
    colour    = "blue4",
    linewidth = 0.65,
    alpha     = 0.50,
    se        = FALSE,
    linetype  = 2
  ) +
  scale_y_continuous(
    name   = "H1(d3 - d1)",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_continuous(
    name = TeX("$\\Delta$%Purs. [Poor - Rich]")
  ) +
  theme_fig() +
  facet_grid(~ factor(constant, labels = "AI")) +
  theme(
    aspect.ratio     = 1,
    strip.background = element_blank()
  )

save_source_data(fig_S9_A_ii_data, "fig_S9_A_ii.csv")
save_figure(fig_S9_A_ii, "fig_S9_A_ii.pdf",
            width = 3.5, height = 3.0)


## S9B – H1 distance-change across ROIs (COSINE)

fig_S9_B_data <- d_rsa_cosine %>%
  tidyr::pivot_longer(
    cols      = tidyselect::ends_with("_H1"),
    names_to  = "region",
    values_to = "H1_distance_change"
  ) %>%
  dplyr::select(subject, region, H1_distance_change)

fig_S9_B <- ggplot(
  fig_S9_B_data,
  aes(x = region, y = H1_distance_change)
) +
  stat_summary(
    fun   = "mean",
    geom  = "bar",
    color = "black",
    fill  = "white",
    width = 0.85
  ) +
  geom_point(
    fill  = "white",
    alpha = 0.50,
    color = "black",
    shape = 21
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = 2,
    color      = "red"
  ) +
  scale_x_discrete(
    labels = c(
      "AI_H1"    = "AI",
      "dACC_H1"  = "dACC",
      "dmPFC_H1" = "dmPFC",
      "latFP_H1"  = "lFPC",
      "pgACC_H1" = "pgACC",
      "vpgACC_H1" = "sgACC"
    )
  ) +
  scale_y_continuous(name = "H1 [d3 - d1]") +
  theme_fig() +
  theme(
    aspect.ratio = 0.75,
    axis.text.x  = element_text(size = 14, vjust = 1, angle = 30),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

save_source_data(fig_S9_B_data, "fig_S9_B.csv")
save_figure(fig_S9_B, "fig_S9_B.pdf",
            width = 6.0, height = 3.0)

## S9C – H2 distance-change across ROIs (COSINE)

fig_S9_C_data <- d_rsa_cosine %>%
  tidyr::pivot_longer(
    cols      = tidyselect::ends_with("_H2"),
    names_to  = "region",
    values_to = "H2_distance_change"
  ) %>%
  dplyr::select(subject, region, H2_distance_change)

fig_S9_C <- ggplot(
  fig_S9_C_data,
  aes(x = region, y = H2_distance_change)
) +
  stat_summary(
    fun   = "mean",
    geom  = "bar",
    color = "black",
    fill  = "white",
    width = 0.85
  ) +
  geom_point(
    fill  = "white",
    alpha = 0.50,
    color = "black",
    shape = 21
  ) +
  geom_hline(
    yintercept = 0,
    linetype   = 2,
    color      = "red"
  ) +
  scale_x_discrete(
    labels = c(
      "AI_H2"    = "AI",
      "dACC_H2"  = "dACC",
      "dmPFC_H2" = "dmPFC",
      "latFP_H2"  = "lFPC",
      "pgACC_H2" = "pgACC",
      "vpgACC_H2" = "sgACC"
    )
  ) +
  scale_y_continuous(name = "H2 [d4 - d2]") +
  theme_fig() +
  theme(
    aspect.ratio = 0.75,
    axis.text.x  = element_text(size = 14, vjust = 1, angle = 30),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

save_source_data(fig_S9_C_data, "fig_S9_C.csv")
save_figure(fig_S9_C, "fig_S9_C.pdf",
            width = 6.0, height = 3.0)


## S9D – Pearson-based analogues of Fig 4C-i and 4D-i

## S9D-i – dACC H2 vs behaviour (PEARSON)

fig_S9_D_i_data <- d_rsa_pearson %>%
  dplyr::select(subject, p_accept_10_diff, dACC_H2, constant)

fig_S9_D_i <- ggplot(
  fig_S9_D_i_data,
  aes(x = p_accept_10_diff * 100, y = dACC_H2)
) +
  geom_point(
    shape = 21,
    fill  = "red2",
    size  = 2.5
  ) +
  geom_smooth(
    method    = "lm",
    colour    = "red4",
    linewidth = 0.65,
    alpha     = 0.50,
    se        = FALSE,
    linetype  = 2
  ) +
  scale_y_continuous(
    name   = "H2(d4 - d2)",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_continuous(
    name = TeX("$\\Delta$%Purs. [Poor - Rich]")
  ) +
  theme_fig() +
  theme(
    aspect.ratio     = 1,
    strip.background = element_blank()
  ) +
  facet_grid(~ factor(constant, labels = "dACC"))

save_source_data(
  fig_S9_D_i_data,
  "fig_S9_D_i.csv"
)
save_figure(
  fig_S9_D_i,
  "fig_S9_D_i.pdf",
  width = 3.5, height = 3.0
)

## S9D-ii – AI H2 vs behaviour (PEARSON)

fig_S9_D_ii_data <- d_rsa_pearson %>%
  dplyr::select(subject, p_accept_10_diff, AI_H2, constant)

fig_S9_D_ii <- ggplot(
  fig_S9_D_ii_data,
  aes(x = p_accept_10_diff * 100, y = AI_H2)
) +
  geom_point(
    shape = 21,
    fill  = "lightblue1",
    size  = 2.5
  ) +
  geom_smooth(
    method    = "lm",
    colour    = "blue4",
    linewidth = 0.65,
    alpha     = 0.50,
    se        = FALSE,
    linetype  = 2
  ) +
  scale_y_continuous(
    name   = "H2(d4 - d2)",
    labels = scales::number_format(accuracy = 0.01)
  ) +
  scale_x_continuous(
    name = TeX("$\\Delta$%Purs. [Poor - Rich]")
  ) +
  theme_fig() +
  theme(
    aspect.ratio     = 1,
    strip.background = element_blank()
  ) +
  facet_grid(~ factor(constant, labels = "AI"))

save_source_data(
  fig_S9_D_ii_data,
  "fig_S9_D_ii.csv"
)
save_figure(
  fig_S9_D_ii,
  "fig_S9_D_ii.pdf",
  width = 3.5, height = 3.0
)
