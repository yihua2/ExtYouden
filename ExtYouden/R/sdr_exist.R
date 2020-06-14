#' sdr_exist
#'
#' This script checks if the remaining design has SDR for each block
#'
#' @param x matrix
#' @param block_by_row indication of orientation of the finding
#' @return TRUE/FALSE if the remaining design has SDR
#' @export


###
sdr_exist <- function(x, block_by_row = T){
  if (block_by_row!=T){
    x = t(x)
  }
  if(length(unique(as.vector(x)))<dim(x)[1]){
    return(F)
  }
  else return(T)
}
