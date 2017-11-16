#' Create Simulated Exposure Matrix for Real and Observed Data 
#'
#' This function creates an exposure matrix with real and observed data 
#' after taking into account true positive, false negative, and false 
#' positive rates. The columns represent strata for the observational data #' and the rows represent new events in unit time. Do not run sim.exposure
#' as a stand-alone function.
#' @param exposed.matrix Output exposure matrix from create.exposure function.
#' @param params Output from initialize.data function.
#' @keywords sequential
#' @export
#' @import stats
#' @examples
#' paramtest <- initialize.data(seed=8768, N=1, t0=0, tf=2, NStrata=2, 
#' strataRatio=c(0.2, 0.3, 0.3, 0.2), EventRate=c(0.4, 0.5), sensitivity=0.9, PPVest=0.9, RR=3.0, 
#' MatchRatio=1, maxSampleSize=200, maxTest=1, totalAlpha=0.05, minEvents=5, AlphaSpendType="Wald",
#' AlphaParameter=0.5, address=getwd(), rate=20, offset=30)
#' exposed1 <- create.exposure(paramtest)
#' sim.exposure(exposed1, paramtest)

sim.exposure <- function(exposed.matrix, params) {
	EventRate <- rep(params$EventRate, each=2)
	Priskcase <- params$RR / (params$RR+params$MatchRatio)
	Priskperson <- 1/ (params$MatchRatio + 1)
	CaseRate <- (params$RR + params$MatchRatio) * EventRate / (params$MatchRatio + 1)
	TP.rate <- params$sensitivity*CaseRate
	FN.rate <- (1-params$sensitivity) * CaseRate

  
	PPVrisk <- params$RR*params$PPVest*(1-EventRate)/(params$RR*(params$PPVest-EventRate)+1-params$PPVest)
	
	################!!##
	# specificity <- (1-CaseRate)/(1-params$sensitivity * CaseRate)
	# PPVrisk2 <- params$PPVest*(params$sensitivity*EventRate + (1-specificity)*(1-EventRate))*params$RR/(params$sensitivity*EventRate*params$RR + (1-specificity)*(1-EventRate*params$RR))

	
	#######################
	
	
	FP.risk.rate <- (params$RR*EventRate*params$sensitivity/PPVrisk)-(params$RR*EventRate*params$sensitivity)
	FP.comp.rate <- (EventRate*params$sensitivity/params$PPVest)-(EventRate*params$sensitivity)
	
	TP.decum <- array(NA, dim(exposed.matrix))
	FN.decum <- array(NA, dim(exposed.matrix))
	FP.decum.risk <- array(NA, dim(exposed.matrix))
	FP.decum.comp <- array(NA, dim(exposed.matrix))

	for (i in 1:dim(exposed.matrix)[2]) {
	  TP.decum[,i,] <- round(exposed.matrix[,i,]*TP.rate[i])
	  FN.decum[,i,] <- round(exposed.matrix[,i,]*FN.rate[i])
	  FP.decum.risk[,i,] <- round(exposed.matrix[,i,]*FP.risk.rate[i])
	  FP.decum.comp[,i,] <- round(exposed.matrix[,i,]*FP.comp.rate[i])
	}
	
	TP.draw <- array(NA, dim(TP.decum))
	for (j in 1:dim(TP.decum)[3]) {
		for (i in 1:dim(TP.decum)[2]){
			TP.draw[,i,j] <- apply(matrix(TP.decum[,i,j]), 1, function(x) rpois(1, x))
		}
	}

	TP.risk <- array(NA, dim(TP.decum))
	TP.comp <- TP.risk
	for (j in 1:dim(TP.draw)[3]) {
		for (i in 1:dim(TP.draw)[2]) {
			TP.risk[,i,j] <- apply(matrix(TP.decum[,i,j]), 1, function(x) rbinom(1, x, Priskcase))
			TP.comp[,i,j] <- TP.decum[,i,j] - TP.risk[,i,j]
		}
	}

	FN.draw <- array(NA, dim(FN.decum))
	for (j in 1:dim(FN.decum)[3]) {
		for (i in 1:dim(FN.decum)[2]){
			FN.draw[,i,j] <- apply(matrix(FN.decum[,i,j]), 1, function(x) rpois(1, x))
		}
	}

	FN.risk <- array(NA, dim(FN.decum))
	FN.comp <- FN.risk
	for (j in 1:dim(FN.draw)[3]) {
		for (i in 1:dim(FN.draw)[2]) {
			FN.risk[,i,j] <- apply(matrix(FN.decum[,i,j]), 1, function(x) rbinom(1, x, Priskcase))
			FN.comp[,i,j] <- FN.decum[,i,j] - FN.risk[,i,j]
		}
	}

	FP.risk <- array(NA, dim(FP.decum.risk))
	for (j in 1:dim(FP.decum.risk)[3]) {
		for (i in 1:dim(FP.risk)[2]){
			FP.risk[,i,j] <- apply(matrix(FP.decum.risk[,1,j]), 1, function(x) rpois(1, x))
		}
	}

	FP.comp <- array(NA, dim(FP.decum.comp))
	for (j in 1:dim(FP.decum.comp)[3]) {
		for (i in 1:dim(FP.comp)[2]){
			FP.comp[,i,j] <- apply(matrix(FP.decum.comp[,1,j]), 1, function(x) rpois(1, x))
		}
	}
	
	real.risk <- TP.risk + FN.risk
	real.comp <- TP.comp + FN.comp
	
	test.risk <- TP.risk + FP.risk
	test.comp <- TP.comp + FP.comp
	
	real.risk.acc <- array(NA, dim(real.risk))
	real.comp.acc <- array(NA, dim(real.comp))
	test.risk.acc <- array(NA, dim(test.risk))
	test.comp.acc <- array(NA, dim(test.comp))
	for (i in 1:dim(real.risk)[2]) {
		real.risk.acc[,i,] <- cumsum(real.risk[,i,])
		real.comp.acc[,i,] <- cumsum(real.comp[,i,])
		test.risk.acc[,i,] <- cumsum(test.risk[,i,])
		test.comp.acc[,i,] <- cumsum(test.comp[,i,])
	}

	sim.data <- list(real.risk=real.risk.acc, 
				real.comp=real.comp.acc, 
				test.risk=test.risk.acc, 
				test.comp=test.comp.acc)
	
	return(sim.data)
}