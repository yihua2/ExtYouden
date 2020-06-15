#' shuffle
#'
#' This script generates a randome shuffle within each row of a matrix
#'
#' @param x matrix of a design, b rows and k columns
#' @return Shuffled matrixx
#' @export
#' @examples
#' X = matrix(c(0,1,2,2,0,1), byrow = TRUE,nrow = 2)
#' shuffle(X)

shuffle <- function(x) {
  a = matrix(NA, nrow = NROW(x), ncol = dim(x)[2])
  for (i in 1:NROW(x)) {
    a[i, ] = sample(x[i, ])
  }
  return(a)

}
