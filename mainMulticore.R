# ------------------------------------------------------------------------------
# Setup
# ==============================================================================
rm(list=ls())
# Please change below working directory to where the unzipped folder locates
setwd("~")
source("preparePackages.R")
registerDoFuture()
plan(multisession)

# ------------------------------------------------------------------------------
# Estimation
# ==============================================================================
backtest = FALSE
startAge = 50
endAge   = 100
for(pp in c("A","B","C")) for(ss in 1:2) foreach(dd = 1:19) %dopar% source("allCountries.R", local = TRUE)
for(pp in c("A")) for(ss in 1:2) foreach(dd = 1:19) %dopar% source("allCountriesSlopeConstant.R", local = TRUE)
# Sweden backtest
for(pp in c("A")) for(backtest in c(1900,1950)) for(dd in 16) for(ss in 1:2) source("allCountries.R", local = TRUE)

# ------------------------------------------------------------------------------
# Metropolis Hastings Test
# ==============================================================================
backtest = FALSE
source("metropolisHastingsTest.R")
foreach(dd = 1:38) %dopar% {
  print(dd)
  metropolisHastingsTest(dd, backtest = backtest)
}
# Sweden backtest
for(backtest in c(1900,1950)) for(dd in c(16,16+19)) metropolisHastingsTest(dd = dd, backtest = backtest)

# ------------------------------------------------------------------------------
# Estimate CIs
# ==============================================================================
backtest = FALSE
robust = FALSE
mumax = 2/3
source("calculateCI.R")
# Sweden backtest
for(backtest in c(1900,1950)) source("calculateCI.R")

# ------------------------------------------------------------------------------
# Calculate Life Expectencies
# ==============================================================================
source("cohortLifeExpectencies.R")


# ------------------------------------------------------------------------------
# Robustness
# ==============================================================================
pp = "A"
startAge = 50
endAge   = 110
for(ss in 1:2) foreach (dd = 1:19) %dopar% source("allCountries.R", local = TRUE)

dd = 19
metropolisHastingsTest(dd = dd, backtest = backtest)

robust = TRUE
mumax = 2/3
source("calculateCI.R")
mumax = 1
source("calculateCI.R")
# ------------------------------------------------------------------------------
# Output
# ==============================================================================
# Figures
rmarkdown::render("Output/figure2-4.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure5.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure6.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure7.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure8.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure9.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS1-S2.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS3-S4.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS5-S7.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS8-S14.Rmd", knit_root_dir = "../")
# Tables
rmarkdown::render("Output/table1.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/table2.Rmd", knit_root_dir = "../")

