## 16_figS2.R ---------------------------------------------------------------
## Supplementary Figure S2

source("R/00_setup.R")
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

# A clean decile-binning helper (inline for transparency)
bin_mu <- function(v) {
  bins <- quantile(v, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  cut(v, breaks = unique(bins), include.lowest = TRUE, labels = FALSE)
}

# ===========================================================================
# Figure S2A 
# ==========================================================================

## Shared data wrangling helper
make_S2_data <- function(d, switch_var, facet_label_text) {
  d %>%
    filter(mag == 5) %>%
    drop_na(mu_val_z) %>%
    mutate(mu_bin = bin_mu(mu_val_z)) %>%
    group_by(mu_bin) %>%
    mutate(mu_val_graph = mean(mu_val_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(subject, mu_val_graph) %>%
    summarise(
      switch_rate = 100 * mean(.data[[switch_var]], na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    mutate(facet_label = facet_label_text)
}

## Panel (i): Pursue-switch [5-pt]
fig_S2_A_i_data <- make_S2_data(
  d,
  switch_var       = "pursue_switch",
  facet_label_text = "(i) Pursue-sw. [5-pt]"
)

## Panel (ii): Reject-switch [5-pt]
fig_S2_A_ii_data <- make_S2_data(
  d,
  switch_var       = "reject_switch",
  facet_label_text = "(ii) Reject-sw. [5-pt]"
)

## Combine for faceting
fig_S2_A_data <- bind_rows(fig_S2_A_i_data, fig_S2_A_ii_data)

fig_S2_A <- ggplot(
  fig_S2_A_data,
  aes(x = mu_val_graph, y = switch_rate)
) +
  geom_smooth(
    method    = "lm",
    color     = "black",
    fill      = "deeppink4",
    linewidth = 0.5
  ) +
  stat_summary(
    fun   = "mean",
    geom  = "point",
    shape = 21,
    fill  = "deeppink4",
    color = "black",
    size  = 2
  ) +
  scale_y_continuous(
    name   = "Switch-rate [%]",
    breaks = seq(0, 20, 10)
  ) +
  scale_x_continuous(name = "Ave. value [Z]") +
  theme_fig() +
  facet_wrap(~ facet_label, ncol = 1) +
  theme(
    aspect.ratio = 0.25, 
    strip.background = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 30))

save_source_data(fig_S2_A_data, "fig_S2_A.csv")
save_figure(fig_S2_A, "fig_S2_A.pdf", width = 6, height = 3)

# ===========================================================================
# Figure S2B 
# ==========================================================================

## Shared data wrangling helper for S2B
make_S2_B_data <- function(d, switch_var, facet_label_text) {
  d %>%
    filter(mag == 50) %>%
    drop_na(mu_val_z) %>%
    mutate(mu_bin = bin_mu(mu_val_z)) %>%
    group_by(mu_bin) %>%
    mutate(mu_val_graph = mean(mu_val_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(subject, mu_val_graph) %>%
    summarise(
      switch_rate = 100 * mean(.data[[switch_var]], na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    mutate(facet_label = facet_label_text)
}

## Panel (i): Pursue-switch [50-pt]
fig_S2_B_i_data <- make_S2_B_data(
  d,
  switch_var       = "pursue_switch",
  facet_label_text = "(i) Pursue-sw. [50-pt]"
)

## Panel (ii): Reject-switch [50-pt]
fig_S2_B_ii_data <- make_S2_B_data(
  d,
  switch_var       = "reject_switch",
  facet_label_text = "(ii) Reject-sw. [50-pt]"
)

## Combine for faceting
fig_S2_B_data <- bind_rows(fig_S2_B_i_data, fig_S2_B_ii_data)

fig_S2_B <- ggplot(
  fig_S2_B_data,
  aes(x = mu_val_graph, y = switch_rate)
) +
  geom_smooth(
    method    = "lm",
    color     = "black",
    fill      = "deeppink4",
    linewidth = 0.5
  ) +
  stat_summary(
    fun   = "mean",
    geom  = "point",
    shape = 21,
    fill  = "deeppink4",
    color = "black",
    size  = 2
  ) +
  scale_y_continuous(
    name   = "Switch-rate [%]",
    breaks = seq(0, 20, 10)
  ) +
  scale_x_continuous(name = "Ave. value [Z]") +
  theme_fig() +
  facet_wrap(~ facet_label, ncol = 1) +
  theme(
    aspect.ratio     = 0.25,
    strip.background = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 20))

save_source_data(fig_S2_B_data, "fig_S2_B.csv")
save_figure(fig_S2_B, "fig_S2_B.pdf", width = 6, height = 3)

# ===========================================================================
# Figure S2C 
# ===========================================================================

fig_S2_C_data <- d %>%
  filter(switch_congruence == 1) %>%
  mutate(
    env_bin = factor(
      env_bin,
      levels = c(-1, 1),
      labels = c("Poor", "Rich")
    )
  ) %>%
  group_by(subject, env_bin) %>%
  summarise(
    mean_trial_switch = mean(trial_n_block, na.rm = TRUE),
    .groups = "drop"
  )

fig_S2_C <- ggplot(
  fig_S2_C_data,
  aes(x = env_bin, y = mean_trial_switch, fill = env_bin)
) +
  # group-level mean line across subjects
  stat_summary(
    fun      = "mean",
    geom     = "line",
    color    = "black",
    linetype = 2,
    linewidth = 0.5,
    group    = 1
  ) +
  # subject-level points
  geom_point(
    alpha = 0.30,
    shape = 21
  ) +
  # group-level mean ± CI across subjects
  stat_summary(
    fun.data = "mean_ci",
    geom     = "point",
    color    = "black",
    size     = 7,
    shape    = 22
  ) +
  scale_fill_manual(values = c("deeppink4", "pink")) +
  scale_y_continuous(name = "Trial-in-block") +
  scale_x_discrete(name = "Environment") +
  theme_fig() +
  theme(
    aspect.ratio   = 1.5,
    legend.position = "none",
    axis.line.x    = element_blank(),
    axis.ticks.x   = element_blank()
  )

# save subject-level source data and figure
save_source_data(fig_S2_C_data, "fig_S2_C.csv")
save_figure(fig_S2_C, "fig_S2_C.pdf", width = 3, height = 3)