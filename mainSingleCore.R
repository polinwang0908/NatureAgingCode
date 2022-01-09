# ------------------------------------------------------------------------------
# Setup
# ==============================================================================
rm(list=ls())
# Please change below working directory to where the unzipped folder locates
setwd("~")
source("preparePackages.R")

# ------------------------------------------------------------------------------
# Estimation
# ==============================================================================
backtest = FALSE
for(pp in c("A","B","C")) for(ss in 1:2) for(dd in 1:19) source("allCountries.R", local = TRUE)
for(pp in c("A")) for(ss in 1:2) for(dd in 1:19) source("allCountriesSlopeConstant.R", local = TRUE)
# Sweden backtest
for(pp in c("A")) for(backtest in c(1900,1950)) for(dd in 16) for(ss in 1:2) source("allCountries.R", local = TRUE)

# ------------------------------------------------------------------------------
# Metropolis Hastings Test
# ==============================================================================
source("metropolisHastingsTest.R")
backtest = FALSE
for(dd in 1:38) metropolisHastingsTest(dd, backtest = backtest)
# Sweden backtest
for(backtest in c(1900,1950)) for(dd in c(16,16+19)) metropolisHastingsTest(dd = dd, backtest = backtest)

# ------------------------------------------------------------------------------
# Estimate CIs
# ==============================================================================
backtest = FALSE
source("calculateCI.R")
# Sweden backtest
for(backtest in c(1900,1950)) source("calculateCI.R")

# ------------------------------------------------------------------------------
# Calculate Life Expectencies
# ==============================================================================
source("cohortLifeExpectencies.R")

# ------------------------------------------------------------------------------
# Output
# ==============================================================================
# Figures
rmarkdown::render("Output/figure2.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure3.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure4.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS1S2S3.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figure5S4S5.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS6S7.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/figureS8S9S10.Rmd", knit_root_dir = "../")
# Tables
rmarkdown::render("Output/tableS1.Rmd", knit_root_dir = "../")
rmarkdown::render("Output/tableM1.Rmd", knit_root_dir = "../")

