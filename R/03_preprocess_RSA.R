## 03_preprocess_rsa.R --------------------------------------------------------
## Preprocess RSA outputs (Pearson r and cosine)

source("R/00_setup.R")

## ---------------------------------------------------------------------------
## Paths
## ---------------------------------------------------------------------------

rsa_dir   <- here::here("data", "fMRI", "rsa")
out_dir   <- here::here("data", "fMRI", "rsa")

pearson_path <- file.path(rsa_dir, "rsa_pearson_raw.csv")
cosine_path  <- file.path(rsa_dir, "rsa_cosine_raw.csv")

## ---------------------------------------------------------------------------
## Helper: load + preprocess a single RSA file
## ---------------------------------------------------------------------------

preprocess_rsa <- function(path, metric_label) {
  df <- readr::read_csv(path, show_col_types = FALSE) %>%
    # constant column for facet labels (dACC / AI etc.)
    mutate(
      constant = 1
    )
  
  # Behavioural split: adaptive vs non-adaptive
  mean_behav_change <- mean(df$p_accept_10_diff, na.rm = TRUE)
  
  df <- df %>%
    mutate(
      mean_split = if_else(
        p_accept_10_diff > mean_behav_change,
        "adaptive",
        "non-adaptive"
      ),
      mean_split = factor(mean_split,
                          levels = c("non-adaptive", "adaptive")),
      metric = metric_label
    )
  
  attr(df, "mean_behav_change") <- mean_behav_change
  
  df
}

## ---------------------------------------------------------------------------
## Preprocess Pearson r and cosine RSA outputs
## ---------------------------------------------------------------------------

d_rsa_pearson <- preprocess_rsa(pearson_path, metric_label = "pearson_r")
d_rsa_cosine  <- preprocess_rsa(cosine_path,  metric_label = "cosine")

## Combined long-format object (both metrics)
d_rsa_all <- bind_rows(d_rsa_pearson, d_rsa_cosine)

## ---------------------------------------------------------------------------
## Save preprocessed data
## ---------------------------------------------------------------------------

saveRDS(d_rsa_pearson,
        file = file.path(out_dir, "rsa_pearson_clean_lf.rds"))

saveRDS(d_rsa_cosine,
        file = file.path(out_dir, "rsa_cosine_clean_lf.rds"))

saveRDS(d_rsa_all,
        file = file.path(out_dir, "rsa_all_lf.rds"))
