## 07_RSAl.R ------------------------------------------------------
## Inferential statistics for RSA data

source("R/00_setup.R")

## Setup ---------------------------------------------------------------------

rsa_preproc_dir <- here::here("data", "fMRI", "rsa")

d_rsa_cosine  <- readRDS(file.path(rsa_preproc_dir, "rsa_cosine_clean_lf.rds"))
d_rsa_pearson <- readRDS(file.path(rsa_preproc_dir, "rsa_pearson_clean_lf.rds"))

out_dir <- here::here("results", "rsa")

## Generic function to perform stats on RSA dataset ----------------------

run_rsa_stats <- function(d_rsa) {
  
  ## dACC / AI H2: adaptive vs non-adaptive -------------------------------
  
  dACC_H2_t <- t.test(
    d_rsa$dACC_H2[d_rsa$mean_split == "adaptive"],
    d_rsa$dACC_H2[d_rsa$mean_split == "non-adaptive"]
  )
  
  AI_H2_t <- t.test(
    d_rsa$AI_H2[d_rsa$mean_split == "adaptive"],
    d_rsa$AI_H2[d_rsa$mean_split == "non-adaptive"]
  )
  
  ## Correlations with behavioural change (Δ% accept) ---------------------
  
  dACC_H2_cor <- cor.test(d_rsa$dACC_H2, d_rsa$p_accept_10_diff)
  AI_H2_cor   <- cor.test(d_rsa$AI_H2,  d_rsa$p_accept_10_diff)
  
  ## Supplementary H1 correlations ----------------------------------------
  
  dACC_H1_cor <- cor.test(d_rsa$dACC_H1, d_rsa$p_accept_10_diff)
  AI_H1_cor   <- cor.test(d_rsa$AI_H1,  d_rsa$p_accept_10_diff)
  
  list(
    dACC_H2_t   = dACC_H2_t,
    AI_H2_t     = AI_H2_t,
    dACC_H2_cor = dACC_H2_cor,
    AI_H2_cor   = AI_H2_cor,
    dACC_H1_cor = dACC_H1_cor,
    AI_H1_cor   = AI_H1_cor
  )
}

## 1. Cosine RSA stats -------------------------------------------------------

RSA_cosine_stats <- run_rsa_stats(d_rsa_cosine)

# Save R object
saveRDS(
  RSA_cosine_stats,
  file = file.path(out_dir, "RSA_cosine_stats.rds")
)

# Human-readable summary
sink(file.path(out_dir, "RSA_cosine_stats_summary.txt"))

cat("RSA inferential statistics – cosine metric\n")
cat("==========================================\n\n")

cat("dACC_H2: adaptive vs non-adaptive (t-test)\n")
print(RSA_cosine_stats$dACC_H2_t)
cat("\n\n")

cat("AI_H2: adaptive vs non-adaptive (t-test)\n")
print(RSA_cosine_stats$AI_H2_t)
cat("\n\n")

cat("Correlation: dACC_H2 ~ p_accept_10_diff\n")
print(RSA_cosine_stats$dACC_H2_cor)
cat("\n\n")

cat("Correlation: AI_H2 ~ p_accept_10_diff\n")
print(RSA_cosine_stats$AI_H2_cor)
cat("\n\n")

cat("Supplementary: dACC_H1 ~ p_accept_10_diff\n")
print(RSA_cosine_stats$dACC_H1_cor)
cat("\n\n")

cat("Supplementary: AI_H1 ~ p_accept_10_diff\n")
print(RSA_cosine_stats$AI_H1_cor)
cat("\n\n")

sink()

## 2. Pearson r RSA stats ----------------------------------------------------

RSA_pearson_stats <- run_rsa_stats(d_rsa_pearson)

# Save R object
saveRDS(
  RSA_pearson_stats,
  file = file.path(out_dir, "RSA_pearson_stats.rds")
)

# Human-readable summary
sink(file.path(out_dir, "RSA_pearson_stats_summary.txt"))

cat("RSA inferential statistics – Pearson r metric\n")
cat("============================================\n\n")

cat("dACC_H2: adaptive vs non-adaptive (t-test)\n")
print(RSA_pearson_stats$dACC_H2_t)
cat("\n\n")

cat("AI_H2: adaptive vs non-adaptive (t-test)\n")
print(RSA_pearson_stats$AI_H2_t)
cat("\n\n")

cat("Correlation: dACC_H2 ~ p_accept_10_diff\n")
print(RSA_pearson_stats$dACC_H2_cor)
cat("\n\n")

cat("Correlation: AI_H2 ~ p_accept_10_diff\n")
print(RSA_pearson_stats$AI_H2_cor)
cat("\n\n")

cat("Supplementary: dACC_H1 ~ p_accept_10_diff\n")
print(RSA_pearson_stats$dACC_H1_cor)
cat("\n\n")

cat("Supplementary: AI_H1 ~ p_accept_10_diff\n")
print(RSA_pearson_stats$AI_H1_cor)
cat("\n\n")

sink()
