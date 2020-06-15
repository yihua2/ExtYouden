#' sdr_exist
#'
#' This script checks if the remaining design has SDR for each block
#'
#' @param x matrix
#' @param block_by_row indication of orientation of the finding
#' @return TRUE/FALSE if the remaining design has SDR
#' @export
#' @examples
#' X = matrix(c(0,1,2,2,0,1), byrow = TRUE,nrow = 2)
#' sdr_exist(X)


sdr_exist <- function(x, block_by_row = TRUE){
  if (block_by_row!=TRUE){
    x = t(x)
  }
  if(length(unique(as.vector(x)))<dim(x)[1]){
    return(FALSE)
  }
  else return(TRUE)
}
