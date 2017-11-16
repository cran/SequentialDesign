#' Decumulate Values in a Matrix
#'
#' The function takes a matrix of values that have accumulated over rows and
#' returns a matrix of the incremental increase between each row. Do not run fun.decum
#' as a stand-alone function
#' @param matrix This is a matrix of values to be decumulated where the cumulation occurs over rows within the same column.
#' @keywords sequential
#' @export
#' @examples
#' testarray <- array(NA, dim=c(5,2,2))
#' testarray[,,1] <- cbind(c(1:5), c(1:5)*2)
#' testarray[,,2] <- cbind(c(1:5)*1.5, c(1:5)*3)
#' fun.decum(testarray)

fun.decum <- function(matrix) {
  decum.matrix <- array(NA, dim(matrix))
  decum.matrix[1,,] <- matrix[1,,]
  for (i in dim(matrix)[1]:2) {
    decum.matrix[i,,] <- matrix[i,,]-matrix[i-1,,]
  }
  return(decum.matrix)
}