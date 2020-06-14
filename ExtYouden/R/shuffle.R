#' shuffle
#'
#' This script generates a randome shuffle within each row of a matrix
#'
#' @param x matrix of a design, b rows and k columns
#' @return Shuffled matrixx
#' @export

shuffle <- function(x) {
  a = matrix(NA, nrow = NROW(x), ncol = dim(x)[2])
  for (i in 1:NROW(x)) {
    a[i, ] = sample(x[i, ])
  }
  return(a)

}
