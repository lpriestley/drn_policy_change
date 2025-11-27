## 06_GLM4.R --------------------------------------------------------
## GLM4.*

source("R/00_setup.R")

tc_path   <- here::here("data", "fMRI", "timecourse_preprocessed", "timecourses.rds")
peak_path <- here::here("data", "fMRI", "timecourse_preprocessed", "peaks.rds")

d_tc   <- readRDS(tc_path)
d_peak <- readRDS(peak_path)

out_dir <- here::here("results", "fmriTC")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# run t-test
do_ttest <- function(x, y = NULL, paired = TRUE, name) {
  tt <- t.test(x, y, paired = paired)
  tibble(
    test = name,
    estimate  = tt$estimate[[1]],
    statistic = tt$statistic,
    df        = tt$parameter,
    p_value   = tt$p.value,
    CI_low    = tt$conf.int[1],
    CI_high   = tt$conf.int[2]
  )
}

# save stats
write_stats <- function(df, fname) {
  readr::write_csv(df, file.path(out_dir, fname))
}

# Bonferroni-Holm correction 

apply_holm <- function(df) {
  df %>% mutate(p_holm = p.adjust(p_value, method = "holm"))
}

# Extractor -------------------------------------------------------

get_peak <- function(glm, region_name, contr="contr1") {
  d_peak %>% 
    filter(GLM == glm, region == region_name, contrast == contr) %>% 
    pull(peak)
}

# -----------------------------------------------------------------------------
# GLM4.1
# -----------------------------------------------------------------------------

# 4.1a: effect of policy switch in all subcortical ROIs
GLM_4.1a_regions <- c("DRN", "MBD", "LC", "Hb", "MS")
GLM_4.1a_peaks <- d_peak %>%
  filter(GLM == "GLM_04.1a",
         contrast == "contr1",
         region %in% GLM_4.1a_regions) %>%
  group_by(region) %>%
  summarise(peak_values = list(peak), .groups = "drop")

GLM_4.1a_raw <- GLM_4.1a_peaks %>%
  mutate(
    test = map2(
      peak_values, region,
      ~ do_ttest(
        x = .x,
        y = NULL,
        paired = FALSE,
        name = paste0(.y, '_vs_0', sep='')
      )
    )
  ) %>%
  select(test) %>%
  unnest(test)

GLM_4.1a_bh <- apply_holm(GLM_4.1a_raw) %>% print()

write_stats(GLM_4.1a_bh, "GLM_4.1a.csv")

# GLM4.1a follow-up: effect of policy-switch in DRN vs {MBD, LC, Hb, MS}
regions_other <- c("MBD", "LC", "Hb", "MS")

GLM_4.1a_drn_vs_other_raw <- map_dfr(regions_other, function(rg) {
  do_ttest(
    x = get_peak(glm="GLM_04.1a", region="DRN"),
    y = get_peak(glm="GLM_04.1a", region=rg),
    paired = TRUE,
    name = paste0("policy_switch_DRN_vs_", rg)
  )
})

GLM_4.1a_drn_vs_other_bh <- apply_holm(GLM_4.1a_drn_vs_other_raw) %>% print()

write_stats(GLM_4.1a_drn_vs_other_bh, "GLM_4.1a_drn_vs_other.csv")

# 4.1b: effect of action switch in all subcortical ROIs
GLM_4.1b_regions <- c("DRN", "MBD", "LC", "Hb", "MS")
GLM_4.1b_peaks <- d_peak %>%
  filter(GLM == "GLM_04.1b",
         contrast == "contr1",
         region %in% GLM_4.1b_regions) %>%
  group_by(region) %>%
  summarise(peak_values = list(peak), .groups = "drop")

GLM_4.1b_raw <- GLM_4.1b_peaks %>%
  mutate(
    test = map2(
      peak_values, region,
      ~ do_ttest(
        x = .x,
        y = NULL,
        paired = FALSE,
        name = paste0(.y, '_vs_0', sep='')
      )
    )
  ) %>%
  select(test) %>%
  unnest(test)

GLM_4.1b_bh <- apply_holm(GLM_4.1b_raw) %>% print()

write_stats(GLM_4.1a_bh, "GLM_4.1b.csv")

# GLM4.1a vs GLM4.1b: effect of policy-switch vs action-switch in DRN

GLM_4.1.a_vs_4.1b_drn <- do_ttest(
  x      = get_peak("GLM_04.1a", "DRN"),
  y      = get_peak("GLM_04.1b", "DRN"),
  paired = TRUE,
  name   = "DRN_GLM_4.1a_vs_4.1b"
) %>% print()

write_stats(GLM_4.1.a_vs_4.1b_drn, "GLM_4.1a_vs_4.1b_drn.csv")

# GLM4.1a_supp: effect of policy-switch in DRN vs {MRN, 4th ventricle}
GLM_4.1a_supp_regions <- c("DRN", "MRN", "fourth_ventricle")
GLM_4.1a_supp_peaks <- d_peak %>%
  filter(GLM == "GLM_04.1a_supp",
         contrast == "contr1",
         region %in% GLM_4.1a_supp_regions) %>%
  group_by(region) %>%
  summarise(peak_values = list(peak), .groups = "drop")

GLM_4.1a_supp_raw <- GLM_4.1a_supp_peaks %>%
  mutate(
    test = map2(
      peak_values, region,
      ~ do_ttest(
        x = .x,
        y = NULL,
        paired = FALSE,
        name = paste0(.y, '_vs_0', sep='')
      )
    )
  ) %>%
  select(test) %>%
  unnest(test)

GLM_4.1a_supp_bh <- apply_holm(GLM_4.1a_supp_raw) %>% print()

write_stats(GLM_4.1a_supp_bh, "GLM_4.1a_supp.csv")

GLM_4.1a_supp_drn_vs_other_raw <- map_dfr(c('MRN', 'fourth_ventricle'), function(rg) {
  do_ttest(
    x = get_peak(glm="GLM_04.1a_supp", region="DRN"),
    y = get_peak(glm="GLM_04.1a_supp", region=rg),
    paired = TRUE,
    name = paste0("DRN_vs_", rg)
  )
})

GLM_4.1a_supp_drn_vs_other_bh <- apply_holm(GLM_4.1a_supp_drn_vs_other_raw) %>% print()

write_stats(GLM_4.1a_supp_drn_vs_other_bh, "GLM_4.1a_supp_drn_vs_other.csv")

# -----------------------------------------------------------------------------
# GLM4.2
# -----------------------------------------------------------------------------

# GLM_04.2a: congruent (i.e. pursue) policy-change in poor environments (DRN)
GLM_4.2a_drn<- do_ttest(
  x      = get_peak("GLM_04.2a", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.2a_drn)
write_stats(GLM_4.2a_drn, "GLM_4.2a.csv")

# GLM_04.2b: incongruent (i.e. pursue) policy-change in rich environments (DRN)
GLM_4.2b_drn <- do_ttest(
  x      = get_peak("GLM_04.2b", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.2b_drn)
write_stats(GLM_4.2b_drn, "GLM_4.2b.csv")

# GLM_04.2c: congruent (i.e. reject) policy-change in rich environments (DRN)
GLM_4.2c_drn <- do_ttest(
  x      = get_peak("GLM_04.2c", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.2c_drn)
write_stats(GLM_4.2c_drn, "GLM_4.2c.csv")

# GLM_04.2d: incongruent (i.e. reject) policy-change in poor environments (DRN)
GLM_4.2d_drn <- do_ttest(
  x      = get_peak("GLM_04.2d", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.2d_drn)
write_stats(GLM_4.2d_drn, "GLM_4.2d.csv")

# GLM_04.2a vs GLM_04.2b: congruent-vs-incongruent pursue policy-changes (DRN)
GLM_4.2a_vs_4.2b_drn <- do_ttest(
  x      = get_peak("GLM_04.2a", "DRN"),
  y      = get_peak("GLM_04.2b", "DRN"),
  paired = TRUE,
  name   = "DRN_4.2a_vs_4.2b"
)

print(GLM_4.2a_vs_4.2b_drn)
write_stats(GLM_4.2a_vs_4.2b_drn, "GLM_4.2a_vs_4.2b.csv")

# GLM_04.2c vs GLM_04.2d: congruent-vs-incongruent reject policy-changes (DRN)
GLM_4.2c_vs_4.2d_drn <- do_ttest(
  x      = get_peak("GLM_04.2c", "DRN"),
  y      = get_peak("GLM_04.2d", "DRN"),
  paired = TRUE,
  name   = "DRN_4.2c_vs_4.2d"
)

print(GLM_4.2c_vs_4.2d_drn)
write_stats(GLM_4.2c_vs_4.2d_drn, "GLM_4.2c_vs_4.2d.csv")

# GLM_04.2a_low: congruent pursue change for 5-pt option (DRN)
GLM_4.2a_low_drn <- do_ttest(
  x      = get_peak("GLM_04.2a_low", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_4.2a_low_vs_0"
)

print(GLM_4.2a_low_drn)
write_stats(GLM_4.2a_low_drn, "GLM_4.2a_low.csv")

# GLM_04.2c_low: congruent reject change for 5-pt option (DRN)
GLM_4.2c_low_drn <- do_ttest(
  x      = get_peak("GLM_04.2c_low", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_4.2c_low_vs_0"
)

print(GLM_4.2c_low_drn)
write_stats(GLM_4.2c_low_drn, "GLM_4.2c_low.csv")

# GLM_04.2a vs GLM_04.2a_low: congruent pursue-changes for 10-pt vs 5-pt option (DRN)
GLM_4.2a_vs_4.2a_low_drn <- do_ttest(
  x      = get_peak("GLM_04.2a", "DRN"),
  y      = get_peak("GLM_04.2a_low", "DRN"),
  paired = TRUE,
  name   = "DRN_4.2a_vs_4.2a_low"
)

print(GLM_4.2a_vs_4.2a_low_drn)
write_stats(GLM_4.2a_vs_4.2a_low_drn, "GLM_4.2a_vs_4.2a_low.csv")

# GLM_04.2c vs GLM_04.2c_low: congruent reject-changes for 10-pt vs 5-pt option (DRN)
GLM_4.2c_vs_4.2c_low_drn <- do_ttest(
  x      = get_peak("GLM_04.2c", "DRN"),
  y      = get_peak("GLM_04.2c_low", "DRN"),
  paired = TRUE,
  name   = "DRN_4.2c_vs_4.2c_low"
)

print(GLM_4.2c_vs_4.2c_low_drn)
write_stats(GLM_4.2c_vs_4.2c_low_drn, "GLM_4.2c_vs_4.2c_low.csv")

# GLM_04.2a_supp: congruent (i.e. pursue) policy-change in poor environments (DRN vs LC)
GLM_4.2a_supp<- do_ttest(
  x      = get_peak("GLM_04.2a_supp", "DRN"),
  y      = get_peak("GLM_04.2a_supp", "LC"),
  paired = T,
  name   = "DRN_vs_LC"
)
print(GLM_4.2a_supp)
write_stats(GLM_4.2a_supp, "GLM_4.2a_supp.csv")

# GLM_04.2c_supp: congruent (i.e. pursue) policy-change in rich environments (DRN vs LC)
GLM_4.2c_supp<- do_ttest(
  x      = get_peak("GLM_04.2c_supp", "DRN"),
  y      = get_peak("GLM_04.2c_supp", "LC"),
  paired = T,
  name   = "DRN_vs_LC"
)
print(GLM_4.2c_supp)
write_stats(GLM_4.2c_supp, "GLM_4.2c_supp.csv")

# GLM_04.2a_cortex: congruent policy-change in poor environemnts (dACC, AI)

GLM_4.2a_cortex_regions <- c("dACC", "AI")

GLM_4.2a_cortex <- map_dfr(GLM_4.2a_cortex_regions, function(rg) {
  do_ttest(
    x      = get_peak("GLM_04.2a_cortex", rg),
    y      = NULL,
    paired = FALSE,
    name   = paste0(rg, "_vs_0")
  )
})

GLM_4.2a_cortex_bh <- apply_holm(GLM_4.2a_cortex)

print(GLM_4.2a_cortex_bh)

write_stats(GLM_4.2a_cortex_bh, "GLM_4.2a_cortex.csv")

# GLM_04.2c_cortex: congruent policy-change in rich environemnts (dACC, AI)

GLM_4.2c_cortex_regions <- c("dACC", "AI")

GLM_4.2c_cortex <- map_dfr(GLM_4.2c_cortex_regions, function(rg) {
  do_ttest(
    x      = get_peak("GLM_04.2c_cortex", rg),
    y      = NULL,
    paired = FALSE,
    name   = paste0(rg, "_vs_0")
  )
})

GLM_4.2c_cortex_bh <- apply_holm(GLM_4.2c_cortex)

print(GLM_4.2c_cortex_bh)

write_stats(GLM_4.2c_cortex_bh, "GLM_4.2c_cortex.csv")

# -----------------------------------------------------------------------------
# GLM4.3
# -----------------------------------------------------------------------------

# GLM_04.3a: DRN-to-dACC PPI 
GLM_4.3a <- do_ttest(
  x      = get_peak("GLM_04.3a", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.3a)
write_stats(GLM_4.3a, "GLM_4.3a.csv")

# GLM_04.3b: DRN-to-AI PPI 
GLM_4.3b <- do_ttest(
  x      = get_peak("GLM_04.3b", "DRN"),
  y      = NULL,
  paired = FALSE,
  name   = "DRN_vs_0"
)
print(GLM_4.3b)
write_stats(GLM_4.3b, "GLM_4.3b.csv")

# -----------------------------------------------------------------------------
# GLM4.4
# -----------------------------------------------------------------------------

GLM_4.4_regions <- c("dACC", "AI")

GLM_4.4 <- map_dfr(GLM_4.4_regions, function(rg) {
  do_ttest(
    x      = get_peak("GLM_04.5", rg),
    y      = NULL,
    paired = FALSE,
    name   = paste0(rg, "_vs_0")
  )
})

print(GLM_4.4)

GLM_4.4_bh <- apply_holm(GLM_4.4)
print(GLM_4.4_bh)

write_stats(GLM_4.4_bh, "GLM_4.4.csv")

# -----------------------------------------------------------------------------
# GLM4.5
# -----------------------------------------------------------------------------
GLM_4.5_regions_cortex <- c("AI", "dACC")

GLM_4.5_drn_vs_cortex <- map_dfr(GLM_4.5_regions_cortex, function(rg) {
  do_ttest(
    x      = get_peak("GLM_04.5", "DRN"),
    y      = get_peak("GLM_04.5", rg),
    paired = TRUE,
    name   = paste0("DRN_vs_", rg)
  )
})

GLM_4.5_drn_vs_cortex_bh <- apply_holm(GLM_4.5_drn_vs_cortex)

print(GLM_4.5_drn_vs_cortex_bh)

write_stats(GLM_4.5_drn_vs_cortex_bh, "GLM_4.5_DRN_vs_cortex.csv")

# -----------------------------------------------------------------------------
# GLMS5
# -----------------------------------------------------------------------------

# GLMS5.1: congruent pursue-change in moving time windows (DRN)

win_indices <- 1:10

GLM_S5.1 <- map_dfr(win_indices, function(i) {
  
  glm_name <- sprintf("GLM_0S5.1_win%d", i)
  
  do_ttest(
    x      = get_peak(glm_name, "DRN"),
    y      = NULL,
    paired = FALSE,
    name   = paste0("window", i, "_DRN_vs_0")
  )
  
})

GLM_S5.1_bh <- apply_holm(GLM_S5.1) %>% print()
write_stats(GLM_S5.1_bh, "GLM_S5.1.csv")

# GLMS5.2: congruent pursue-change in moving time windows (DRN)

win_indices <- 1:10

GLM_S5.2 <- map_dfr(win_indices, function(i) {
  
  glm_name <- sprintf("GLM_0S5.2_win%d", i)
  
  do_ttest(
    x      = get_peak(glm_name, "DRN"),
    y      = NULL,
    paired = FALSE,
    name   = paste0("window", i, "_DRN_vs_0")
  )
  
})

GLM_S5.2_bh <- apply_holm(GLM_S5.2) %>% print()
write_stats(GLM_S5.2_bh, "GLM_S5.2.csv")

# -----------------------------------------------------------------------------
# GLMS6
# -----------------------------------------------------------------------------

# GLMS6.1: action initiation, i.e., pursue-vs-reject (MBD)
GLM_S6.1 <- do_ttest(
  x      = get_peak("GLM_0S6.1", "MBD"),
  y      = NULL,
  paired = FALSE,
  name   = "MBD_vs_0"
)
print(GLM_S6.1)
write_stats(GLM_S6.1, "GLM_S6.1.csv")

# GLMS6.2: value-difference i.e., rw_val[t] - ave_val[t] (MBD)
GLM_S6.2 <- do_ttest(
  x      = get_peak("GLM_0S6.2", "MBD"),
  y      = NULL,
  paired = FALSE,
  name   = "MBD_vs_0"
)
print(GLM_S6.2)
write_stats(GLM_S6.2, "GLM_S6.2.csv")

# -----------------------------------------------------------------------------
# GLMS7
# -----------------------------------------------------------------------------

# GLMS7.1: effect of incongruent changes 
GLM_S7.1_regions <- c("DRN", "MBD", "LC", "Hb", "MS")
GLM_S7.1_peaks <- d_peak %>%
  filter(GLM == "GLM_0S7.1",
         contrast == "contr1",
         region %in% GLM_S7.1_regions) %>%
  group_by(region) %>%
  summarise(peak_values = list(peak), .groups = "drop")

GLM_S7.1 <- GLM_S7.1_peaks %>%
  mutate(
    test = map2(
      peak_values, region,
      ~ do_ttest(
        x = .x,
        y = NULL,
        paired = FALSE,
        name = paste0(.y, '_vs_0', sep='')
      )
    )
  ) %>%
  select(test) %>%
  unnest(test)

GLM_S7.1_bh <- apply_holm(GLM_S7.1) %>% print()

write_stats(GLM_S7.1_bh, "GLM_S7.1.csv")

# GLMSS7.2: compare effects of incongruent and congruent changes (DRN & MS)
anova_data <- d_peak %>%
  filter(
    GLM %in% c("GLM_0S7.2"),
    region %in% c("DRN", "MS"),
    contrast %in% c("contr1", "contr2")
  ) %>%
  select(GLM, region, contrast, peak) %>%
  drop_na(peak) %>%
  group_by(GLM, contrast, region) %>%
  mutate(ID = row_number()) %>%
  ungroup() %>%
  mutate(
    ID       = factor(ID),
    region   = factor(region,   levels = c("DRN", "MS")),
    contrast = factor(contrast, levels = c("contr1", "contr2"))
  )

anova_res <- afex::aov_ez(
  id     = "ID",
  dv     = "peak",
  data   = anova_data,
  within = c("region", "contrast"),
  type   = 3
)

anova_table <- as.data.frame(anova_res$anova_table) %>%
  rownames_to_column("effect")

write_stats(anova_table, "GLM_S7.2_ANOVA.csv")

GLM_S7.2_drn_vs_ms <- do_ttest(
  x      = get_peak("GLM_0S7.2", "DRN", "contr2"),
  y      = get_peak("GLM_0S7.2", "MS", "contr2"),
  paired = TRUE,
  name   = "DRN_vs_MS_contr2"
)
print(GLM_S7.2_drn_vs_ms)
write_stats(GLM_S7.2_drn_vs_ms, "GLM_S7.2_drn_vs_ms.csv")



