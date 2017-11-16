#' Define Exposure Accumulation Function
#'
#' This function creates a dataset simulating the accumulation of individuals exposed to 
#' treatment under the self-controlled risk interval design.
#' @param rate Rate of accumulation.
#' @param offset Initial exposed population.
#' @param t Time at which individuals are exposed.
#' @keywords sequential exposure simulation
#' @export
#' @examples
#' fun.exposure(rate=100, offset=20, t=20)

fun.exposure <- function(rate, offset, t) {
	exp.data <- offset + rate*t
	return(exp.data)
}