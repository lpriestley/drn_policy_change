# README

Analysis code for *ACTIVITY IN HUMAN DORSAL RAPHE NUCLEUS SIGNALS CHANGES IN BEHAVIOURAL POLICY *.

This repository contains all code required to:

- preprocess behavioural and fMRI timecourse data
- compute GLM and RSA-based statistics
- reproduce all main and supplementary figures

R-based analysis scripts live under `rev_3/R/` and can be run individually or via a single orchestration script (`run_all.R`).

MATLAB-based analysis scripts live under `rev_3/MATLAB/`. 

---

## 1. Getting started

### 1.1. Requirements

- R (≥ 4.0 recommended)
- R packages:
  - Core: `tidyverse`, `here`, `lme4`, `lmerTest`
  - Plotting: `ggplot2`, `ggpubr`, `latex2exp`, `patchwork`, `gghalves`, `ggdist`
  - Misc: `scales`, `reshape2`, `rstatix`, `grid`
- All packages are loaded via `rev_3/R/00_setup.R`.  
  If you source that file and get “package not found” errors, just `install.packages("<name>")` as needed.

### 1.2. Directory structure

Layout (simplified):

```text
project_root/
  rev_3/
    R/
      00_setup.R
      01_preprocess_behaviour.R
      02_preprocess_fmriTC.R
      03_preprocess_RSA.R
      04_GLM_1.R
      05_GLM_2.R
      06_GLM_3.R
      07_RSA.R
      08_fig1.R
      09_fig1D.R
      10_fig2.R
      11_fig3.R
      12_fig4.R
      13_fig5.R
      14_fig6.R
      15_figS1.R
      16_figS2.R
      17_figS5.R
      18_figS6.R
      19_figS7.R
      20_figS9.R
      run_all.R
  data/
    behaviour/
      clean_lf.rds              # produced by 01_preprocess_behaviour.R
      ...
    fMRI/
      timecourse_preprocessed/
        peaks.rds
        timecourses.rds
      rsa/
        rsa_cosine_clean_lf.rds
        rsa_pearson_clean_lf.rds
        ...
  rev_3/results/
    ...
  rev_3/results/rsa/
    ...                         # RSA inferential stats
  (figure / source-data dirs – configured in 00_setup.R via save_figure(), save_source_data())
