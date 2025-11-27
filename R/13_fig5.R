## 13_fig5.R ---------------------------------------------------------------
## Figure 5 

source("R/00_setup.R")

## -------------------------------------------------------------------------
## Load data
## -------------------------------------------------------------------------

peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")
d_peak    <- readRDS(peak_path)

rsa_path <- here::here("data", "fMRI", "rsa", "rsa_cosine_clean_lf.rds")
d_rsa    <- readRDS(rsa_path)

## ------------------------------------------------------------------------
## Generic helper function: construct similarity matrix for heatmaps
## ------------------------------------------------------------------------

create_similarity_matrix <- function(data, d1_col, d2_col, d3_col, d4_col, env_label) {
  expand.grid(
    subject = unique(data$subject),
    opt_1   = c("low", "mid", "high"),
    opt_2   = c("low", "mid", "high")
  ) %>%
    dplyr::left_join(
      data %>%
        dplyr::select(
          subject,
          p_accept_10_diff,
          mean_split,
          dplyr::all_of(c(d1_col, d2_col, d3_col, d4_col))
        ),
      by = "subject"
    ) %>%
    dplyr::mutate(
      d_val = dplyr::case_when(
        opt_1 == opt_2 ~ NA_real_,
        (opt_1 == "low" & opt_2 == "mid") |
          (opt_1 == "mid" & opt_2 == "low") ~ .data[[d3_col]],
        (opt_1 == "mid" & opt_2 == "high") |
          (opt_1 == "high" & opt_2 == "mid") ~ .data[[d4_col]],
        (opt_1 == "low" & opt_2 == "high") |
          (opt_1 == "high" & opt_2 == "low") ~ NA_real_,
        TRUE ~ NA_real_
      ),
      env = env_label
    ) %>%
    dplyr::select(subject, opt_1, opt_2, d_val, p_accept_10_diff, env, mean_split)
}


## ------------------------------------------------------------------------
## Figure 5C -  RSA results for dACC
## ------------------------------------------------------------------------

dACC_rich <- create_similarity_matrix(
  d_rsa, "dACC_d1", "dACC_d2", "dACC_d3", "dACC_d4", "rich"
)

dACC_poor <- create_similarity_matrix(
  d_rsa, "dACC_d1", "dACC_d2", "dACC_d1", "dACC_d2", "poor"
)

d_rsa_dACC <- dplyr::bind_rows(dACC_rich, dACC_poor)

d_rsa_dACC_summary <- d_rsa_dACC %>%
  dplyr::group_by(mean_split, opt_1, opt_2, env) %>%
  dplyr::summarise(
    mean_d_val = mean(d_val, na.rm = TRUE),
    .groups    = "drop"
  )

dACC_global_stats <- d_rsa_dACC_summary %>%
  dplyr::summarise(
    global_min    = min(mean_d_val, na.rm = TRUE),
    global_max    = max(mean_d_val, na.rm = TRUE),
    global_median = median(mean_d_val, na.rm = TRUE)
  )

# Fig. 5C-i
fig_5_C_i_data <- d_rsa %>%
  dplyr::select(subject, p_accept_10_diff, dACC_H2, constant)

fig_5_C_i <- ggplot(
  fig_5_C_i_data,
  aes(x = p_accept_10_diff * 100, y = dACC_H2)
) +
  geom_point(shape = 21, fill = "red2", size = 2.5) +
  geom_smooth(
    method    = "lm",
    color     = "red4",
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
    name = TeX("$\\Delta$%acc. [Poor - Rich]")
  ) +
  theme_fig() +
  facet_grid(~ factor(constant, labels = "dACC")) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank()
  )

save_source_data(fig_5_C_i_data, "fig_5_C_i.csv")
save_figure(fig_5_C_i, "fig_5_C_i.pdf",
            width = 3.5, height = 3)

# Fig 5C-ii
## Adaptive subjects
fig_5_C_ii_adaptive_data <- d_rsa_dACC_summary %>%
  dplyr::filter(mean_split == "adaptive")

na_tiles <- fig_5_C_ii_adaptive_data %>%
  dplyr::filter(is.na(mean_d_val) | is.nan(mean_d_val))

fig_5_C_ii_adaptive <- ggplot(
  fig_5_C_ii_adaptive_data,
  aes(
    x    = factor(opt_1, levels = c("high", "mid", "low")),
    y    = factor(opt_2, levels = c("low", "mid", "high")),
    fill = mean_d_val
  )
) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(
    name      = "Cosine distance",
    low       = "#B2182B",
    mid       = "white",
    high      = "#2166AC",
    na.value  = "grey90",
    midpoint  = dACC_global_stats$global_median,
    limits    = c(dACC_global_stats$global_min - 0.01,
                  dACC_global_stats$global_max + 0.01),
    breaks    = round(
      c(
        dACC_global_stats$global_min - 0.01,
        (dACC_global_stats$global_min + dACC_global_stats$global_max) / 2,
        dACC_global_stats$global_max + 0.01
      ),
      2
    ),
    space = "Lab"
  ) +
  geom_point(
    data = na_tiles,
    aes(
      x = factor(opt_1, levels = c("high", "mid", "low")),
      y = factor(opt_2, levels = c("low", "mid", "high"))
    ),
    inherit.aes = FALSE,
    shape  = 4,      # 'x' shape; use 3 for '+' if you prefer
    size   = 11,
    stroke = 0.5,
    colour = "black"
  ) +
  # Axis titles + relabelled ticks
  scale_x_discrete(
    name   = "Rw. magnitude",
    labels = c("50", "10", "5")    # levels: high, mid, low
  ) +
  scale_y_discrete(
    name   = "Rw. magnitude",
    labels = c("5", "10", "50")    # levels: low, mid, high
  ) +
  theme_fig() +
  theme(
    legend.position  = "right",
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    aspect.ratio     = 1
  ) +
  facet_grid(~ factor(env, levels = c("poor", "rich"),
                      labels = c("Poor", "Rich"))) +
  coord_fixed()

save_source_data(fig_5_C_ii_adaptive_data, "fig_5_C_ii_adaptive.csv")
save_figure(fig_5_C_ii_adaptive, "fig_5_C_ii_adaptive.pdf",
            width = 6.0, height = 3.0)

## Non-adaptive subjects
fig_5_C_ii_nonadaptive_data <- d_rsa_dACC_summary %>%
  dplyr::filter(mean_split == "non-adaptive")

na_tiles <- fig_5_C_ii_nonadaptive_data %>%
  dplyr::filter(is.na(mean_d_val) | is.nan(mean_d_val))

fig_5_C_ii_nonadaptive <- ggplot(
  fig_5_C_ii_nonadaptive_data,
  aes(
    x    = factor(opt_1, levels = c("high", "mid", "low")),
    y    = factor(opt_2, levels = c("low", "mid", "high")),
    fill = mean_d_val
  )
) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(
    name      = "Cosine distance",
    low       = "#B2182B",
    mid       = "white",
    high      = "#2166AC",
    na.value  = "grey90",
    midpoint  = dACC_global_stats$global_median,
    limits    = c(dACC_global_stats$global_min - 0.01,
                  dACC_global_stats$global_max + 0.01),
    breaks    = round(
      c(
        dACC_global_stats$global_min - 0.01,
        (dACC_global_stats$global_min + dACC_global_stats$global_max) / 2,
        dACC_global_stats$global_max + 0.01
      ),
      2
    ),
    space = "Lab"
  ) +
  geom_point(
    data = na_tiles,
    aes(
      x = factor(opt_1, levels = c("high", "mid", "low")),
      y = factor(opt_2, levels = c("low", "mid", "high"))
    ),
    inherit.aes = FALSE,
    shape  = 4,
    size   = 11,
    stroke = 0.5,
    colour = "black"
  ) +
  # Axis titles + relabelled ticks
  scale_x_discrete(
    name   = "Rw. magnitude",
    labels = c("50", "10", "5")    # levels: high, mid, low
  ) +
  scale_y_discrete(
    name     = "Rw. magnitude",
    labels   = c("5", "10", "50"), # levels: low, mid, high
    position = "right"
  ) +
  theme_fig() +
  theme(
    legend.position  = "right",
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    aspect.ratio     = 1
  ) +
  facet_grid(~ factor(env, levels = c("poor", "rich"),
                      labels = c("Poor", "Rich"))) +
  coord_fixed()

save_source_data(fig_5_C_ii_nonadaptive_data, "fig_5_C_ii_nonadaptive.csv")
save_figure(fig_5_C_ii_nonadaptive, "fig_5_C_ii_nonadaptive.pdf",
            width = 6.0, height = 3.0)


## ------------------------------------------------------------------------
## Figure 5D -  RSA results for AI
## ------------------------------------------------------------------------

AI_rich <- create_similarity_matrix(
  d_rsa, "AI_d1", "AI_d2", "AI_d3", "AI_d4", "rich"
)

AI_poor <- create_similarity_matrix(
  d_rsa, "AI_d1", "AI_d2", "AI_d1", "AI_d2", "poor"
)

d_rsa_AI <- dplyr::bind_rows(AI_rich, AI_poor)

d_rsa_AI_summary <- d_rsa_AI %>%
  dplyr::group_by(mean_split, opt_1, opt_2, env) %>%
  dplyr::summarise(
    mean_d_val = mean(d_val, na.rm = TRUE),
    .groups    = "drop"
  )

AI_global_stats <- d_rsa_AI_summary %>%
  dplyr::summarise(
    global_min    = min(mean_d_val, na.rm = TRUE),
    global_max    = max(mean_d_val, na.rm = TRUE),
    global_median = median(mean_d_val, na.rm = TRUE)
  )

# Fig 5D-i 

fig_5_D_i_data <- d_rsa %>%
  dplyr::select(subject, p_accept_10_diff, AI_H2, constant)

fig_5_D_i <- ggplot(
  fig_5_D_i_data,
  aes(x = p_accept_10_diff * 100, y = AI_H2)
) +
  geom_point(shape = 21, fill = "lightblue1", size = 2.5) +
  geom_smooth(
    method    = "lm",
    color     = "blue4",
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
    name = TeX("$\\Delta$%acc. [Poor - Rich]")
  ) +
  theme_fig() +
  facet_grid(~ factor(constant, labels = "AI")) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank()
  )

save_source_data(fig_5_D_i_data, "fig_5_D_i.csv")
save_figure(fig_5_D_i, "fig_5_D_i.pdf",
            width = 3.5, height = 3)

## Fig 5D-ii – AI similarity heatmaps
## Adaptive subjects
fig_5_D_ii_adaptive_data <- d_rsa_AI_summary %>%
  dplyr::filter(mean_split == "adaptive")

na_tiles_AI_adapt <- fig_5_D_ii_adaptive_data %>%
  dplyr::filter(is.na(mean_d_val) | is.nan(mean_d_val))

fig_5_D_ii_adaptive <- ggplot(
  fig_5_D_ii_adaptive_data,
  aes(
    x    = factor(opt_1, levels = c("high", "mid", "low")),
    y    = factor(opt_2, levels = c("low", "mid", "high")),
    fill = mean_d_val
  )
) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(
    name      = "Cosine distance",
    low       = "#B2182B",
    mid       = "white",
    high      = "#2166AC",
    na.value  = "grey90",
    midpoint  = AI_global_stats$global_median,
    limits    = c(AI_global_stats$global_min - 0.01,
                  AI_global_stats$global_max + 0.01),
    breaks    = round(
      c(
        AI_global_stats$global_min - 0.01,
        (AI_global_stats$global_min + AI_global_stats$global_max) / 2,
        AI_global_stats$global_max + 0.01
      ),
      2
    ),
    space = "Lab"
  ) +
  geom_point(
    data = na_tiles_AI_adapt,
    aes(
      x = factor(opt_1, levels = c("high", "mid", "low")),
      y = factor(opt_2, levels = c("low", "mid", "high"))
    ),
    inherit.aes = FALSE,
    shape  = 4,   # 'x' shape
    size   = 11,
    stroke = 0.5,
    colour = "black"
  ) +
  # Axis titles + relabelled ticks
  scale_x_discrete(
    name   = "Rw. magnitude",
    labels = c("50", "10", "5")      # levels: high, mid, low
  ) +
  scale_y_discrete(
    name   = "Rw. magnitude",
    labels = c("5", "10", "50")      # levels: low, mid, high
  ) +
  theme_fig() +
  theme(
    legend.position  = "right",
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    aspect.ratio     = 1
  ) +
  facet_grid(~ factor(env, levels = c("poor", "rich"),
                      labels = c("Poor", "Rich"))) +
  coord_fixed()

save_source_data(fig_5_D_ii_adaptive_data, "fig_5_D_ii_adaptive.csv")
save_figure(fig_5_D_ii_adaptive, "fig_5_D_ii_adaptive.pdf",
            width = 6.0, height = 3.0)

## Non-adaptive subjects
fig_5_D_ii_nonadaptive_data <- d_rsa_AI_summary %>%
  dplyr::filter(mean_split == "non-adaptive")

na_tiles_AI_nonadapt <- fig_5_D_ii_nonadaptive_data %>%
  dplyr::filter(is.na(mean_d_val) | is.nan(mean_d_val))

fig_5_D_ii_nonadaptive <- ggplot(
  fig_5_D_ii_nonadaptive_data,
  aes(
    x    = factor(opt_1, levels = c("high", "mid", "low")),
    y    = factor(opt_2, levels = c("low", "mid", "high")),
    fill = mean_d_val
  )
) +
  geom_tile(color = "black", linewidth = 0.5) +
  scale_fill_gradient2(
    name      = "Cosine distance",
    low       = "#B2182B",
    mid       = "white",
    high      = "#2166AC",
    na.value  = "grey90",
    midpoint  = AI_global_stats$global_median,
    limits    = c(AI_global_stats$global_min - 0.03,
                  AI_global_stats$global_max + 0.03),
    breaks    = round(
      c(
        AI_global_stats$global_min - 0.01,
        (AI_global_stats$global_min + AI_global_stats$global_max) / 2,
        AI_global_stats$global_max + 0.01
      ),
      2
    ),
    space = "Lab"
  ) +
  geom_point(
    data = na_tiles_AI_nonadapt,
    aes(
      x = factor(opt_1, levels = c("high", "mid", "low")),
      y = factor(opt_2, levels = c("low", "mid", "high"))
    ),
    inherit.aes = FALSE,
    shape  = 4,
    size   = 11,
    stroke = 0.5,
    colour = "black"
  ) +
  # Axis titles + relabelled ticks
  scale_x_discrete(
    name   = "Rw. magnitude",
    labels = c("50", "10", "5")      # levels: high, mid, low
  ) +
  scale_y_discrete(
    name     = "Rw. magnitude",
    labels   = c("5", "10", "50"),   # levels: low, mid, high
    position = "right"
  ) +
  theme_fig() +
  theme(
    legend.position  = "right",
    strip.text.y     = element_blank(),
    strip.background = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    aspect.ratio     = 1
  ) +
  facet_grid(~ factor(env, levels = c("poor", "rich"),
                      labels = c("Poor", "Rich"))) +
  coord_fixed()

save_source_data(fig_5_D_ii_nonadaptive_data, "fig_5_D_ii_nonadaptive.csv")
save_figure(fig_5_D_ii_nonadaptive, "fig_5_D_ii_nonadaptive.pdf",
            width = 6.0, height = 3.0)


## -------------------------------------------------------------------------
## Figure 5E – DRN / AI / dACC peaks from GLM4.2a and 4.2c
## -------------------------------------------------------------------------

fig_5_E_data <- d_peak %>%
  filter(
    region   %in% c("DRN", "AI", "dACC"),
    contrast == "contr1",
    GLM      %in% c("GLM_04.2a_cortex", "GLM_04.2c_cortex")
  )

fig_5_E <- ggplot(fig_5_E_data, aes(x = GLM, y = peak, fill = GLM, group = GLM)) +
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
    labels = scales::number_format(accuracy = 0.01),
    limits = c(-0.30, 0.60),
    breaks = seq(-0.25, 0.50, 0.25)
  ) +
  scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = c("GLM_04.2a_cortex" = "deeppink4",
                               "GLM_04.2c_cortex" = "pink")) +
  theme_fig() +
  facet_wrap(
    ~ region,
    ncol = 1,
    labeller = labeller(
      region = c(
        "DRN"  = "(iii) DRN",
        "AI"   = "(i) AI",
        "dACC" = "(ii) dACC"
      )
    )
  ) + 
  theme(
    legend.position = "none", 
    strip.background = element_blank(),
    panel.spacing.x  = unit(0.1, "lines"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    aspect.ratio = 1.75
  )

save_source_data(fig_5_E_data, "fig_5_E.csv")
save_figure(fig_5_E, "fig_5_E.pdf",
            width = 3, height = 7)




