#' Execute Simulated Exact Sequential Analysis in Multi-Site Observational Database Studies
#'
#' SequentialDesign is designed for planning observational database studies that 
#' use exact sequential analysis. It is designed to be used in conjunction with 
#' the R package Sequential. This package is appropriate to use when one is
#' performing a multi-site observational database study (i.e., an epidemiologic study) and 
#' planning to use sequential statistical analysis. This package supports two 
#' types of observational study designs: 
#' \itemize{
#' 	\item a self-controlled risk interval design which creates binomial data, and 
#' 	\item a current v. historical design which creates Poisson data. 
#' }
#' The goal of this package is to allow the investigator to plan for the optimal study.
#' 
#' The simulation has the following steps: 
#' \enumerate{
#' 	\item Sample Size Calculations for the study using the R Sequential package 
#' 	\item Given these sample size calculations and an exposure uptake function, 
#' calculate new exposure accrual in calendar time for the exposures of interest.
#' 	\item Given the simulated exposure information, generate adverse events of 
#' interest according to a pre-specified effect size.
#' 	\item Perform sequential analysis on these simulated data. 
#' 	\item Generate calendar time descriptive statistics with respect to stopping points. 
#' }
#' These steps will be discussed in more detail. \cr
#' First, the investigator should work with the R package Sequential in order to 
#' calculate design parameters for their study. These are the statistical 
#' parameters that govern stopping points in statistical analysis. The relevant 
#' ones required for this analysis are: maxSampleSize, totalAlpha, minEvents, 
#' AlphaSpendType, AlphaParameter. If binomial data is being used for sequential 
#' analysis of a self-controlled risk interval design, then MatchRatio is also 
#' needed. \cr
#' Second, this function will generate incident exposure to a simulated 
#' study population based on the parameters of an exposure accrual function. \cr
#' Third, with incremental exposure accrual information, new adverse events 
#' will be assigned based on user-specified characteristics. 
#' This function also allows for outcome misclassification so true positive 
#' adverse events, false positive adverse events, and false negative adverse 
#' events are all simualted. \cr
#' Fourth, Sequential analysis is implemented on these simulated data using 
#' function in R Sequential. \cr 
#' Fifth, the investigator is able to generate descriptive statistics in 
#' calendar time to enable the investigator to plan for their analysis. \cr
#' Simulating sequential analysis in observational data requires many parameter 
#' inputs about 
#' \itemize{
#' 	\item the parameters that control the epidemiologic study design, 
#' 	\item the parameters that describe the characteristics of the databases, and 
#' 	\item the parameters of the simulation. 
#' }
#' In addition to the parameter inputs, there are many sub-functions that are 
#' needed to perform different steps in the simulation. These sub-functions 
#' are not intended to be run as stand-alone functions but rather always 
#' in the sequence specified in this function.
#' @param seed Seed used for randomization
#' @param N Number of simulations to be created. Because adverse event assignment is stochastic, this number is usually at least 10,000.
#' @param t0 Initial time point, a number in units of either days, weeks, months, or years. It is important to be consistent.
#' @param tf Final time point, a number in units of either days, weeks, months, or years.
#' @param NStrata Number of strata in the observational study design, where a "stratum" can be defined by age categories, sex, and any other defining characteristics. Event rate of the adverse event of interest is also segregated by strata and database population size is also segregated by strata. For example, a single strata might 0-17 year old females.
#' @param strataRatio Ratio of individuals within a single strata for exposed and unexposed individuals.  The number of elements in this list should be 2*NStrata.
#' @param EventRate Rate of event accrual given in events /person-time where the time constant is the same constant being used throughout the study. Additionally, the number of elements in the EventRate matrix should be equal to the NStrata.
#' @param sensitivity True sensitivity of the outcome of interest.  sensitivity = (true positive case) / (true positive case + false negative case). 
#' @param PPVest True positive predictive value of outcome in the unexposed group.  PPV = (true positive case) / (true positive case + false positive case).
#' @param RR Intended relative risk to detect (and therefore to simulate) in the dataset.
#' @param MatchRatio Single numeric value. In a self-controlled risk interval design, it is the ratio of the length of the control window to the length of the risk window.
#' @param maxSampleSize Maximum number of events before sequential analysis is ended or the upper limit on sample size expressed in terms of total number of events. This is the same variable as N from R Sequential.
#' @param maxTest Number of tests to perform on simulation data.
#' @param totalAlpha Total amount of alpha available to spend.
#' @param minEvents Minimum number of events needed before the null hypothesis can be rejected. Represented as M in R Sequential.
#' @param AlphaSpendType Method of alpha spenditure.  Available values are "Wald" or "power-type". This is the same as AlphaSpending R Sequential.
#' @param AlphaParameter Rho parameter for power-type alpha spending function. This is the same as rho in R Sequential.
#' @param rate Rate of exposure/cohort accrual.
#' @param offset Offset for exposure/cohort accrual.
#' @param address Output folder where Sequential TXT files are to be stored. These should be preserved between runs, as detailed within the Sequential package.
#' @param ... additional arguments to be passed to or from methods.
#' @keywords sequential 
#' @export
#' @import Sequential
#' @examples
#' #paramtest <- initialize.data(seed=8768, N=1, t0=0, tf=2, NStrata=2, 
#' #strataRatio=c(0.2, 0.3, 0.3, 0.2), EventRate=c(0.4, 0.5), sensitivity=0.9, PPVest=0.9, RR=3.0, 
#' #MatchRatio=1, maxSampleSize=200, maxTest=1, totalAlpha=0.05, minEvents=5, AlphaSpendType="Wald",
#' #AlphaParameter=0.5, address=getwd(), rate=20, offset=30)
#' #exposed1 <- create.exposure(paramtest)
#' #exposed2 <- sim.exposure(exposed1, paramtest)
#' #SCRI.seq(exposed2, paramtest)


seq_wrapper <- function(seed, N, t0, tf, NStrata, strataRatio, 
		EventRate, sensitivity, PPVest, RR, 
		MatchRatio, maxSampleSize, maxTest, totalAlpha, minEvents,
		AlphaSpendType, AlphaParameter, rate, offset, address, ...) {
		
		setwd(address)
		
		params <- initialize.data(seed=seed, N=N, t0=t0, tf=tf, 
						NStrata=NStrata, strataRatio=strataRatio, 
						EventRate=EventRate, sensitivity=sensitivity, PPVest=PPVest, RR=RR, 
						MatchRatio=MatchRatio, maxSampleSize=maxSampleSize, 
						maxTest=maxTest, totalAlpha=totalAlpha, minEvents=minEvents,
						AlphaSpendType=AlphaSpendType, AlphaParameter=AlphaParameter, 
						address=address, rate=rate, offset=offset)
		
		set.seed(params$seed)
		
		exposure <- create.exposure(params)
		
		exposure2 <- sim.exposure(exposure, params)
		
		SCRI.results <- SCRI.seq(exposure2, params)
		
		return(SCRI.results)
	
}