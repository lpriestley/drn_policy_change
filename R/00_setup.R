## 00_setup.R -------------------------------------------------------------
## Shared setup: packages, options, generic functions

## Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(ggpubr)
  library(latex2exp)
  library(afex)
  library(patchwork)
})

## Global figure theme ---------------------------------------------------

theme_fig <- function(
    base_size = 18,
    axis_text_size = 14,
    strip_text_size = 16,
    legend_text_size = 14
) {
  theme_pubr(base_size = base_size) +
    theme(
      text            = element_text(size = base_size),
      axis.text       = element_text(size = axis_text_size),
      strip.text      = element_text(size = strip_text_size),
      legend.text     = element_text(size = legend_text_size),
      legend.background = element_rect(fill = "grey99")
    )
}
theme_set(theme_fig())

## Generic functions --------------------------------------------------------

bin_decile <- function(x, probs = seq(0, 1, length.out = 11)) {
  breaks <- unique(quantile(x, probs = probs, na.rm = TRUE))
  
  # Assign each observation to bin
  idx <- cut(
    x,
    breaks = breaks,
    labels = FALSE,
    include.lowest = TRUE
  )
  
  # Mean of x within each bin
  means <- tapply(x, idx, function(z) mean(z, na.rm = TRUE))
  
  # Build output vector of same length
  out <- rep(NA_real_, length(x))
  ok  <- !is.na(idx)
  
  # tapply names the result by bin index; coerce to character for indexing
  out[ok] <- means[as.character(idx[ok])]
  
  out
}

save_source_data <- function(data, filename) {
  readr::write_csv(
    data,
    here::here("source_data", filename)
  )
}

save_figure <- function(plot, filename, width, height) {
  ggplot2::ggsave(
    filename = here::here("figs", filename),
    plot     = plot,
    width    = width,
    height   = height
  )
}

fit_glm_binom <- function(formula, data, name = NULL, save = TRUE) {
  m <- stats::glm(
    formula = formula,
    data    = data,
    family  = "binomial"
  )
  
  if (!is.null(name) && isTRUE(save)) {
    saveRDS(
      m,
      here::here("results", "behaviour", paste0(name, ".rds"))
    )
    
    summary_path <- here::here("results", "behaviour", paste0(name, "_summary.txt"))
    sink(summary_path)
    on.exit(sink(), add = TRUE)
    print(summary(m))
  }
  
  m
}

fit_glmer_binom <- function(formula, data, name = NULL, save = TRUE) {
  m <- lme4::glmer(
    formula = formula,
    data    = data,
    family  = "binomial",
    control = lme4::glmerControl(
      optimizer = "optimx",
      optCtrl   = list(method = "nlminb")
    )
  )
  
  if (!is.null(name) && isTRUE(save)) {
    # Save model object
    saveRDS(
      m,
      here::here("results", "behaviour", paste0(name, ".rds"))
    )
    
    # Save human-readable summary
    summary_path <- here::here("results", "behaviour", paste0(name, "_summary.txt"))
    sink(summary_path)
    on.exit(sink(), add = TRUE)
    print(summary(m))
  }
  
  m
}
