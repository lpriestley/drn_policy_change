source("R/00_setup.R")

scripts <- c(
  "R/01_preprocess_behaviour.R",
  "R/02_preprocess_fmriTC.R",
  "R/03_preprocess_RSA.R",    
  "R/04_GLM_1.R", 
  "R/05_GLM_2.R", 
  "R/06_GLM_4.R", 
  "R/07_RSA.R", 
  "R/08_fig1.R",
  #"R/09_fig1D.R", # Simulation will take several minutes to run
  "R/10_fig2.R",
  "R/11_fig3.R",
  "R/12_fig4.R",
  "R/13_fig5.R",
  "R/14_fig6.R",
  "R/15_figS1.R",
  "R/16_figS2.R",
  "R/17_figS5.R",
  "R/18_figS6.R",
  "R/19_figS7.R",
  "R/20_figS9.R"
)

for (s in scripts) {
  message("Running: ", s)
  source(s, echo = TRUE, max.deparse.length = Inf)
}
