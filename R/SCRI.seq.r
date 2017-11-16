#' Perform sequential analysis on true and misclassified binomial data 
#'
#' This function performs a prespecified number of binomial sequential analyses on real and misclassified binomial data that is designed to simulate a self-controlled risk interval design. Do not run SCRI.seq
#' as a stand-alone function.
#' @param data Output from sim.exposure that contains real and observed data.
#' @param params Output from initialize.data function.
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

SCRI.seq <- function(data, params) {
  
	x <- c(params$t0:params$tf)
	timeN <- length(x)
	
	z <- params$MatchRatio
	p <- 1/(1+z)
	
	result <- list()
	real.result <- list()
	real.result.str <- list()
	real.result.sim <- list()
	test.result <- list()
	test.result.str <- list()
	test.result.sim <- list()
	
	### Check if text files already exist
	
	## create txt files in new folder, not where package is defined

	
	for (k in 1:params$N) {
		for (j in 1:(params$NStrata*2)) {
		  
		  realname <- paste("datareal", j, "_", k, sep="")
		  testname <- paste("datatest", j, "_", k, sep="")
		  
		  AnalyzeSetUp.Binomial(name=realname,N=params$maxSampleSize, alpha=params$totalAlpha, pp=p, M=params$minEvents, AlphaSpendType=params$AlphaSpendType, rho=params$AlphaParameter, address=params$address)
		  AnalyzeSetUp.Binomial(name=testname,N=params$maxSampleSize, alpha=params$totalAlpha, pp=p, M=params$minEvents, AlphaSpendType=params$AlphaSpendType, rho=params$AlphaParameter, address=params$address)
		  
			for (i in 1:params$maxTest) {
		 		  
				result <- Analyze.Binomial(name=realname, test=i, p=p, cases=data[[1]][i,j,k], controls=data[[2]][i,j,k])
				real.result[[i]] <- result
				
				
				result <- Analyze.Binomial(name=testname, test=i, p=p, cases=data[[3]][i,j,k], controls=data[[4]][i,j,k])
				test.result[[i]] <- result
			}
			real.result.str[[j]] <- real.result
			test.result.str[[j]] <- test.result
		}
		real.result.sim[[k]] <- real.result.str
		test.result.sim[[k]] <- test.result.str
	}
	
	results <- list(real.results=real.result.sim, test.results=test.result.sim)
	return(results)
}