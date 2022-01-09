This code estimates the model in the paper, "Longevity records have not increased in 25
years. So has the human lifespan reached its limit?", and produces figures and tables 
based off the estimation.

==========================================================================================
* There are four folders where the data and results are stored:

	1. Data
	2. Estimates
	3. MH
	4. Output

* It contains the following scripts:

	1. allCountries.R: 
	This script estimates the Bayesian cohort-based model with variable DRA (the rate of
	increase in mortality rates by age) and stores the estimates and model results in the
	"Estimates" folder. There are folders for three different priors, A, B, and, C.

	2. allCountriesSlopeConstant.R: 
	This script estimates the Bayesian cohort-based model with fixed DRA (the rate of
	increase in mortality rates by age) and stores the estimates and model results in the
	"Estimates" folder. The fixed DRA case is calculated for prior A.

	3. calculateCI.R: 
	This script takes the MH results and estimates to calculate the 95% confidence
	intervals for the age at which cohort mortality hazard first reaches 2/3. The results 
	are stored in the "Estimates/CI" folder.

	4. cohortLifeExpectencies.R: 
	This script takes the mortality data and calculates the actual changes in remaining
	cohort base life expectancy at age 50 over ten-year birth cohorts.

	5. funCohortModel.R: 
	This script contains all functions used in "allCountries.R".

	6. funCohortModelSlopeConstant.R: 
	This script contains all functions used in "allCountriesSlopeConstant.R".

	7. main.R: 
	This is the main model execution.
	
	8. mainSingleCore.R:
	This is the non-parallel version of main.R, which does exactly the same thing but
	using only single core.

	9. metropolisHastingsTest.R: 
	This script takes the model estimates and runs Metropolis-Hastings algorithm to sample
	the posterior distribution.

	10. plotCI.R: 
	This script plots the CIs figures used in output.

	11. pareparePackages.R: 
	This script parepares all packages needed in this set of scripts.

	12. readMortalityData.R: 
	This script reads off the mortality data and organizes it in the format that's used.

	13. sampleHM.R: 
	This script set up the parameters of MH samples. The paper uses 100,000 MH draws in
	which we have 1,000 burn-in and take every 99th draw. This gives a sample of 1,000
	draws.

	14. setCountries.R: 
	The paper uses 19 rich industrialized countries that are specified here.

	15. setup.R: 
	This script sets up the prior distribution for the Bayesian cohort-based model.

* Use of the codes. All execution is in "main.R", including model estimation, MH sampling,
  calculating CIs, and producing tables and figures in the paper. Please see the following
  blocks of execution:

	1. Setup: 
	This sets up the working directory, prepares all packages needed and loads them.
	Please specify the directory to the local machine.

	2. Estimation: 
	This uses the "allCountries.R" and "allCountriesSlopeConstant.R" to estimate the
	mortality model. The results are stored in the "Estimates" folder. The program is set
	up to run all 19 countries and both sexes, and the back test for Sweden. To run a
	specific prior, sex, and country, one can specify "pp" (prior that takes string values
	of "A", "B", or "C"), "ss" (integer 1 is males and 2 is females), and "dd" (integer 
	1-19, see countries in "setCountries.R" that have them in alphabetical order). Once
	the parameters are set, execute "allCountries.R" or "allCountriesSlopeConstant.R". For 
	example, to run USA females with variable DRA using prior A, one can run the following
	line of code.
	# pp = "A"; ss = 2; dd = 19; source("allCountries.R")
	
	"allCountries.R" takes approximately 15 minutes to estimate one prior, one country, 
	and one sex. "allCountriesSlopeConstant.R" takes approximately 10 minutes to estimate 
	one prior, one country, and one sex.

	3. Metropolis Hastings Test: 
	This takes the model estimates and samples the prior distribution.
	
	Approximate run time: 19 hours for one sex, one country with 100,000 draws.

	4. Estimate CIs: 
	This takes the MH results and estimates the confidence intervals for the age at which
	cohort mortality hazard first reaches 2/3.
	
	Approximate run time: 20 minutes for one sex, one country.

	5. Calculate Life Expectancies: 
	This uses the mortality data to calculate the life expectancies.

	6. Output: 
	All figures are tables in the paper are produced here using above calculations and
	stored in the "Output" folder. Output contains word documents including Figure 2, 3,
	4, 5, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, and Table M1 and S1.

* Note: Some of the executions are written in "main.R" to run in parallel for the sake of 
        speed. An alternative is provided as "mainSingleCore.R" that works exactly the
        same with a single core. Approximate run time provided above can vary 
        significantly due to the length of data available, hence the number of parameters 
        to estimate (e.g. Sweden takes the longest). The executions have no run time 
        approximated above usually take seconds to finish. 
        
        The run time approximation is based on the following computer specs:

			OS Name:                   Microsoft Windows 10 Enterprise
			OS Version:                10.0.18363 N/A Build 18363
			Processor(s):              Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz 3.00 GHz
			System Type:			   64-bit Operating System, x64-based processor
			Installed Memory (RAM):    64.0 GB
