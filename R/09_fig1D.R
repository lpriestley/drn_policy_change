## 09_fig1D.R ----------------------------------------------------
## Figure 1D 

source("R/00_setup.R")

set.seed(123)  # for reproducibility


## ------------------------------------------------------------------------
## Task parameters
## ------------------------------------------------------------------------

low_mag  <- 5
mod_mag  <- 10
high_mag <- 50

p_rich <- c(1/6, 2/6, 3/6)
p_poor <- c(3/6, 2/6, 1/6)

block_length <- 60      # trials per block
n_blocks     <- 10       # number of blocks to simulate (matches legacy)

block_dur <- 4.5 * 60    # seconds per block (4.5 min)

iti_dur      <- 4.5
gocue_dur    <- 4.5
response_dur <- 1.5
ao_dur       <- 4.5
rw_dur       <- 4.5

## Convenience for fixed frame times
t_reject <- iti_dur + gocue_dur + response_dur
# For accepted trials we add a variable decision time plus ao + rw


## ------------------------------------------------------------------------
## Agent parameters
## ------------------------------------------------------------------------

p_acc_low_rich  <- 0.05
p_acc_high_rich <- 1.00

p_acc_low_poor  <- 0.05
p_acc_high_poor <- 1.00

# Grid over mid-option acceptance probabilities in rich vs poor
p_acc_mod_rich <- seq(0.20, 1.00, 0.10)
p_acc_mod_poor <- seq(0.20, 1.00, 0.10)

n_reps <- 1e3  # number of Monte-Carlo repetitions

## ------------------------------------------------------------------------
## Vectorised block simulation
## ------------------------------------------------------------------------

simulate_block <- function(env_label, p_mod_poor, p_mod_rich) {
  # Choose environment-specific probabilities
  if (env_label == "rich") {
    p_vec <- p_rich
  } else {
    p_vec <- p_poor
  }
  
  # Draw the sequence of offers for this block
  options <- rep(
    x     = c(low_mag, mod_mag, high_mag),
    times = p_vec * block_length
  ) |>
    sample()
  
  n_trials <- length(options)
  
  if (env_label == "rich") {
    prob_acc <- dplyr::case_when(
      options == low_mag  ~ p_acc_low_rich,
      options == mod_mag  ~ p_mod_rich,
      options == high_mag ~ p_acc_high_rich,
      TRUE                ~ 0
    )
  } else {
    prob_acc <- dplyr::case_when(
      options == low_mag  ~ p_acc_low_poor,
      options == mod_mag  ~ p_mod_poor,
      options == high_mag ~ p_acc_high_poor,
      TRUE                ~ 0
    )
  }
  
  # Sample responses for all trials at once
  response <- rbinom(n = n_trials, size = 1, prob = prob_acc)
  
  # Decision times only for accepted trials
  decision_time <- numeric(n_trials)
  n_accept      <- sum(response == 1)
  if (n_accept > 0) {
    decision_time[response == 1] <- runif(
      n    = n_accept,
      min  = 0.15,
      max  = response_dur
    )
  }
  
  durations <- numeric(n_trials)
  durations[response == 0] <- t_reject
  durations[response == 1] <- iti_dur + gocue_dur +
    decision_time[response == 1] + ao_dur + rw_dur
  
  cum_times         <- cumsum(durations)
  n_completed_block <- sum(cum_times <= block_dur)
  
  if (n_completed_block == 0) {
    return(list(total_rw_block = 0, n_trials_block = 0))
  }
  
  # Only count reward from trials that occur within the block duration
  valid_idx <- seq_len(n_completed_block)
  total_rw_block <- sum(options[valid_idx][response[valid_idx] == 1])
  
  list(
    total_rw_block = total_rw_block,
    n_trials_block = n_completed_block
  )
}

## ------------------------------------------------------------------------
## Simulation function: single “session” for given mid-option policies
## ------------------------------------------------------------------------

simulate_session <- function(p_mod_poor, p_mod_rich) {
  total_rw     <- 0
  total_trials <- 0
  
  for (block in seq_len(n_blocks)) {
    
    env_label <- if (block %% 2 == 1) "rich" else "poor"
    
    block_out <- simulate_block(
      env_label  = env_label,
      p_mod_poor = p_mod_poor,
      p_mod_rich = p_mod_rich
    )
    
    total_rw     <- total_rw + block_out$total_rw_block
    total_trials <- total_trials + block_out$n_trials_block
  }
  
  list(
    total_rw     = total_rw,
    total_trials = total_trials
  )
}


## ------------------------------------------------------------------------
## Run simulation grid over p(accept | mag=10, poor) × p(accept | mag=10, rich)
## ------------------------------------------------------------------------

n_poor  <- length(p_acc_mod_poor)
n_rich  <- length(p_acc_mod_rich)
n_grid  <- n_poor * n_rich * n_reps

# Preallocate results (much faster than growing a list of tibbles)
sim_results <- tibble::tibble(
  rep             = integer(n_grid),
  p_acc_mod_poor  = numeric(n_grid),
  p_acc_mod_rich  = numeric(n_grid),
  total_rw        = numeric(n_grid),
  total_trials    = integer(n_grid)
)

idx <- 1L

for (poor_val in p_acc_mod_poor) {
  for (rich_val in p_acc_mod_rich) {
    for (rep in seq_len(n_reps)) {
      
      sim <- simulate_session(
        p_mod_poor = poor_val,
        p_mod_rich = rich_val
      )
      
      sim_results$rep[idx]            <- rep
      sim_results$p_acc_mod_poor[idx] <- poor_val
      sim_results$p_acc_mod_rich[idx] <- rich_val
      sim_results$total_rw[idx]       <- sim$total_rw
      sim_results$total_trials[idx]   <- sim$total_trials
      
      idx <- idx + 1L
    }
  }
}


## ------------------------------------------------------------------------
## Summarise across repetitions and plot heatmap (Figure 1D)
## ------------------------------------------------------------------------

heatmap_summary <- sim_results %>%
  dplyr::group_by(p_acc_mod_poor, p_acc_mod_rich) %>%
  dplyr::summarise(
    mean_total_rw     = mean(total_rw,     na.rm = TRUE),
    mean_total_trials = mean(total_trials, na.rm = TRUE),
    .groups           = "drop"
  )

## ------------------------------------------------------------------------
## Figure 1D – styled to match RSA similarity matrices
## ------------------------------------------------------------------------

library(grid)  # for unit()

step <- diff(sort(unique(heatmap_summary$p_acc_mod_poor)))[1]

x_min <- min(heatmap_summary$p_acc_mod_poor) - step / 2
x_max <- max(heatmap_summary$p_acc_mod_poor) + step / 2
y_min <- min(heatmap_summary$p_acc_mod_rich) - step / 2
y_max <- max(heatmap_summary$p_acc_mod_rich) + step / 2

legend_min <- 2900
legend_max <- 3300

fig_1D <- ggplot(
  heatmap_summary,
  aes(
    x    = p_acc_mod_poor,
    y    = p_acc_mod_rich,
    fill = mean_total_rw
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  annotate(
    "rect",
    xmin     = x_min,
    xmax     = x_max,
    ymin     = y_min,
    ymax     = y_max,
    colour   = "black",
    fill     = NA,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    name   = "p(Pursue | Poor)",
    breaks = seq(0.2, 1.0, 0.2),
    limits = c(x_min, x_max),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "p(Pursue | Rich)",
    breaks = seq(0.2, 1.0, 0.2),
    limits = c(y_min, y_max),
    expand = c(0, 0)
  ) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(9, "YlOrRd"),  # invert palette here
    name    = "Total points",
    limits  = c(legend_min, legend_max),
    breaks  = c(legend_min, legend_max),
    labels  = c(as.integer(legend_min), as.integer(legend_max))
  ) + 
  guides(
    fill = guide_colorbar(
      title.position  = "top",
      title.hjust     = 0.5,
      direction       = "horizontal",
      label.position  = "bottom",
      barwidth        = unit(3.5, "cm"),
      barheight       = unit(0.4, "cm")
    )
  ) +
  coord_fixed(expand = FALSE) +
  theme_fig() +
  theme(
    legend.position    = "top",
    legend.direction   = "horizontal",
    legend.box         = "vertical",
    legend.box.spacing = unit(0, "pt"), 
    strip.background   = element_blank(),
    panel.border       = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.x        = element_text(margin = margin(t = 2)),
    axis.text.y        = element_text(margin = margin(r = 2)),
    plot.margin        = margin(2, 5, 5, 5)
  )

save_source_data(heatmap_summary, "fig_1_D.csv")

save_figure(
  fig_1D,
  "fig_1_D.pdf",
  width  = 3.5,
  height = 3.5
)
