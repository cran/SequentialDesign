#' Create Simulated Sequential Data Parameter Data Frame
#'
#' The function creates a data frame with all the needed parameters for simulation and 
#' initializes the simulation problem. Do not run initialize.data as a stand-alone function.
#' @param seed Seed used for randomization. 
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
#' @param AlphaSpendType Method of alpha expenditure.  Available values are "Wald" or "power-type". This is the same as AlphaSpending R Sequential.
#' @param AlphaParameter Rho parameter for power-type alpha spending function. This is the same as rho in R Sequential.
#' @param address File directory where data for sequential analysis is stored for future tests.
#' @param rate Rate of exposure/cohort accrual.
#' @param offset Offset for exposure/cohort accrual.
#' @keywords sequential
#' @export
#' @examples
#' initialize.data(seed=8768, N=1, t0=0, tf=2, NStrata=2, strataRatio=c(0.2, 0.3, 0.3, 0.2),
#' EventRate=c(0.4, 0.5), sensitivity=0.9, PPVest=0.9, RR=3.0, MatchRatio=1, maxSampleSize=200, 
#' maxTest=1, totalAlpha=0.05, minEvents=5, AlphaSpendType="Wald", AlphaParameter=0.5, address=getwd(),
#' rate=20, offset=30)

initialize.data <- function(seed, N, t0, tf, NStrata, strataRatio,
		EventRate, sensitivity, PPVest, RR,
		MatchRatio, maxSampleSize, maxTest, totalAlpha, minEvents,
		AlphaSpendType, AlphaParameter, address, rate, offset) {

	params <- list(seed=seed, 
					N=N,
					t0=t0,
					tf=tf,
					NStrata = NStrata,
					strataRatio=strataRatio,
					EventRate=EventRate,
					sensitivity=sensitivity,
					PPVest=PPVest,
					RR=RR,
					MatchRatio=MatchRatio,
					maxSampleSize=maxSampleSize,
					maxTest=maxTest,
					totalAlpha=totalAlpha,
					minEvents=minEvents,
					AlphaSpendType=AlphaSpendType,
					AlphaParameter=AlphaParameter,
					address=address,
					rate=rate,
					offset=offset
					)
	return(params)
}
