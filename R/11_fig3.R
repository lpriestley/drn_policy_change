## 11_fig3.R ---------------------------------------------------------------
## Figure 3

source("R/00_setup.R")
d <- readRDS(here::here("data", "behaviour", "clean_lf.rds"))

# ======================================================================
# FIGURE 3B 
# ======================================================================

## Policy data
fig_3_B_policy_data <- d %>%
  drop_na(prev_policy) %>%
  group_by(subject, prev_policy) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    prev_choice = factor(
      prev_policy,
      levels = c(0, 1),
      labels = c("Rej.", "Purs.")
    ),
    m_percent   = m * 100,
    facet_label = "(ii) Policy"
  )

## Action data
fig_3_B_action_data <- d %>%
  drop_na(prev_action) %>%
  group_by(subject, prev_action) %>%
  summarise(m = mean(response, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    prev_choice = factor(
      prev_action,
      levels = c(0, 1),
      labels = c("Rej.", "Purs.")
    ),
    m_percent   = m * 100,
    facet_label = "(i) Action"
  )

## Combine
fig_3_B_data <- bind_rows(fig_3_B_policy_data, fig_3_B_action_data)

fig_3_B <- ggplot(
  fig_3_B_data,
  aes(
    x    = prev_choice,
    y    = m_percent,
    fill = interaction(facet_label, prev_choice, sep = ":")
  )
) +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    color    = "black",
    alpha    = 0.75,
    width    = 0.8,
    position = "dodge"
  ) +
  geom_point(
    size  = 2.5,
    alpha = 0.5,
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    guide  = "none",
    values = c(
      "(i) Action:Rej."   = "azure1",
      "(i) Action:Purs."  = "skyblue3",
      "(ii) Policy:Rej."  = "azure1",
      "(ii) Policy:Purs." = "deeppink4"
    )
  ) +
  scale_x_discrete(name = "Previous behaviour") +
  scale_y_continuous(
    name   = "Pursue-rate [%]",
    limits = c(0, 100)
  ) +
  facet_grid(~ facet_label) +
  theme_fig() +
  theme(
    aspect.ratio     = 1.5,
    legend.position  = "none",
    strip.background = element_blank(),
    panel.spacing    = unit(1.2, "lines")
  )

save_source_data(fig_3_B_data, "fig_3_B.csv")
save_figure(fig_3_B, "fig_3_B.pdf", width = 4.5, height = 3)


save_source_data(fig_3_B_data, "fig_3_B.csv")
save_figure(fig_3_B, "fig_3_B.pdf", width = 4.5, height = 3)


# ==========================================================================
# FIGURE 3C
# ==========================================================================

fig_3_C_data <- d %>%
  drop_na(policy_switch) %>%
  mutate(policy_switch_factor = factor(policy_switch)) %>%
  count(subject, mag, policy_switch_factor, .drop = FALSE) %>%
  filter(policy_switch_factor == 1)

fig_3_C <- ggplot(
  fig_3_C_data,
  aes(x = factor(mag), y = n, fill = factor(mag))
) +
  stat_summary(
    fun   = "mean",
    geom  = "bar",
    color = "black",
    alpha = 0.75
  ) +
  geom_point(
    size  = 2.5,
    alpha = 0.40,
    shape = 21,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("azure1", "deeppink4", "white")
  ) +
  scale_x_discrete(name = "Reward mag.") +
  scale_y_continuous(name = "Policy switches [N]") +
  theme_fig() +
  theme(
    aspect.ratio    = 1,
    legend.position = "none"   # remove legend
  )


save_source_data(fig_3_C_data, "fig_3_C.csv")
save_figure(fig_3_C, "fig_3_C.pdf", width = 3, height = 3)


# ==========================================================================
# FIGURE 3D
# ==========================================================================

fig_3_D_data <- d %>%
  filter(env %in% c("poor", "rich")) %>%
  group_by(subject, mag, env) %>%
  summarise(
    m_pursue = mean(pursue_switch, na.rm = TRUE),
    m_reject = mean(reject_switch, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  pivot_longer(
    cols      = starts_with("m_"),
    names_to  = "switch_type",
    values_to = "switch_rate"
  ) %>%
  filter(mag == 10) %>%
  mutate(
    switch_type        = factor(switch_type, labels = c("Purs.", "Rej.")),
    switch_rate_percent = switch_rate * 100,
    facet_label        = ifelse(env == "poor", "(i) Poor env.", "(ii) Rich env.")
  )

fig_3_D <- ggplot(
  fig_3_D_data,
  aes(x = switch_type, y = switch_rate_percent, fill = switch_type)
) +
  stat_summary(
    fun   = "mean",
    geom  = "bar",
    color = "black",
    alpha = 0.75,
    width = 0.75
  ) +
  geom_point(
    size  = 2.5,
    shape = 21,
    alpha = 0.40
  ) +
  scale_fill_manual(values = c("deeppink4", "azure1")) +
  scale_y_continuous(
    name   = "Switch-rate [%]",
    limits = c(0, 40),
    breaks = seq(0, 30, 10)
  ) +
  scale_x_discrete(name = "Switch-type") +
  facet_grid(~ facet_label) +
  theme_fig() +
  theme(
    aspect.ratio    = 1.25,
    legend.position = "none",
    strip.background = element_blank()
  )

# Save combined source data and figure
save_source_data(fig_3_D_data, "fig_3_D.csv")
save_figure(fig_3_D, "fig_3_D.pdf", width = 4.5, height = 3)


# ==========================================================================
# FIGURE 3E
# ==========================================================================

fig_3_E_data <- d %>%
  filter(mag == 10) %>%
  drop_na(mu_val_z) %>%
  mutate(mu_bin = cut(mu_val_z,
                      breaks = unique(quantile(mu_val_z, probs = seq(0, 1, length.out = 11))),
                      labels = FALSE, include.lowest = TRUE)) %>%
  group_by(mu_bin) %>%
  mutate(mu_val_graph = mean(mu_val_z, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(subject, prev_policy, mu_val_graph) %>%
  summarise(
    pursue = 100 * mean(pursue_switch, na.rm = TRUE),
    reject = 100 * mean(reject_switch, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c("pursue", "reject"),
               names_to = "switch_type", values_to = "switch_rate") %>%
  mutate(switch_type = recode(switch_type,
                              pursue = "(i) Pursue switch",
                              reject = "(ii) Reject switch"))

fig_3_E <- ggplot(fig_3_E_data,
                aes(x = mu_val_graph, y = switch_rate,
                    fill = switch_type)) +
  geom_smooth(method = "lm", color = "black",
              linewidth = 0.5) +
  stat_summary(fun = 'mean', geom = 'point', shape = 21, color = 'black', size = 3, linewidth = 0.5) + 
  scale_fill_manual(values = c("deeppink4", "azure1")) +
  scale_y_continuous(name = "Switch-rate [%]") +
  scale_x_continuous(name = "Ave. value [Z]") +
  facet_grid(~ switch_type) +
  theme_fig() +
  theme(aspect.ratio = 1,
        legend.position = "none", 
        strip.background = element_blank())

save_source_data(fig_3_E_data, "fig_3_E.csv")
save_figure(fig_3_E, "fig_3_E.pdf", width = 5, height = 3)

