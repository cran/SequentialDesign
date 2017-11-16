#' Create Exposure Data Matrix 
#'
#' This function create.exposure is a sub-function used in conjunction with the 
#' initialize.data function and creates an exposure matrix. The columns 
#' represent strata for the observational data and the rows represent new 
#' exposures in unit time. It take cumulative data and segregates it by time 
#' period. Do not run create.exposure as a stand-alone function.
#' @param params This is a set of parameters from the initialize.data function that allows for simulation of a sequence of sequential exposures in unit time.
#' @keywords sequential 
#' @export
#' @examples
#' paramtest <- initialize.data(seed=8768, N=1, t0=0, tf=2, NStrata=2, 
#' strataRatio=c(0.2, 0.3, 0.3, 0.2), EventRate=c(0.4, 0.5), sensitivity=0.9, PPVest=0.9, RR=3.0, 
#' MatchRatio=1, maxSampleSize=200, maxTest=1, totalAlpha=0.05, minEvents=5, AlphaSpendType="Wald",
#' AlphaParameter=0.5, address=getwd(), rate=20, offset=30)
#' create.exposure(paramtest)

create.exposure <- function(params) {
	x <- c(params$t0:params$tf)
	timeN <- length(x)

	out.exposure <- array(NA, c(timeN, params$NStrata*2, params$N))
	rownames(out.exposure) <- x
	colnames(out.exposure) <- paste(c(rep("T", params$NStrata), rep("C", params$NStrata)), 1:params$NStrata, sep="")
	dimnames(out.exposure)[[3]] <- c(1:params$N)
	
	exposed <- array(NA, c(1,length(x),params$N))
	for (i in 1:params$N) {
		exposed[,,i] <- fun.exposure(params$rate, params$offset, x)+1
	}
	
	for (i in 1:params$N) {
		out.exposure[,,i]<-apply(matrix(params$strataRatio), 1, "*", exposed[,,i])+1
	}

	decum.exposed <- fun.decum(out.exposure)
	return(decum.exposed)
}