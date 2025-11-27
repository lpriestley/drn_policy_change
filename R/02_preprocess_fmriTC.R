## 02_preprocess_fmriTC.R ----------------------------------------------------
## Preprocess fMRI timecourse and peak outputs from MATLAB

source("R/00_setup.R")

raw_dir   <- here::here("data", "fMRI", "timecourse_outputs")
out_dir   <- here::here("data", "fMRI", "timecourse_preprocessed")

# ============================================================================
# Helper: Recode ROI names
# ============================================================================

recode_roi <- function(x) {
  dplyr::recode(
    x,
    "dACCsphere7" = "dACC",
    "AIsphere7" = "AI",
    "DRN_custom" = "DRN",
    "HB" = "Hb",
    "LC_Pauli" = "LC",
    "MDB" = "MDB",
    "PP" = "PPN",
    "BF" = "MS",
    .default = x
  )
}

# ============================================================================
# Load + prepare timecourse CSVs
# ============================================================================

tc_files <- list.files(raw_dir, pattern = "*_timecourse.csv", full.names = TRUE)

timecourse_list <- lapply(tc_files, function(f) {
  
  df <- read.csv(f, header = TRUE)
  
  # Standardise column names
  df <- df %>%
    rename(
      m       = tc_data1,
      se      = tc_data2,
      region  = tc_data3,
      glm     = tc_data4,
      contrast = tc_data5
    ) %>%
    mutate(
      region = recode_roi(region),
      timepoint = 1:n(),
      timepoint = scales::rescale(timepoint, to = c(-2, 8))
    )
  
  df
})

timecourses <- bind_rows(timecourse_list)

saveRDS(timecourses, file = file.path(out_dir, "timecourses.rds"))

# ============================================================================
# Load + prepare peak CSVs
# ============================================================================

peak_files <- list.files(raw_dir, pattern = "*peaks.csv", full.names = TRUE)

peak_list <- lapply(peak_files, function(f) read.csv(f))

peaks_raw <- bind_rows(peak_list)

# wide → long
peaks <- peaks_raw %>%
  tidyr::pivot_longer(
    cols = c(
      "MBD", "DRN_custom", "LC_Pauli", "BF", "HB",
      "AIsphere7", "dACCsphere7",
      "fourth_ventricle", "MRN"
    ),
    names_to = "region",
    values_to = "peak"
  ) %>%
  mutate(region = recode_roi(region))

saveRDS(peaks, file = file.path(out_dir, "peaks.rds"))

