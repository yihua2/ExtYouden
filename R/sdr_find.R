#' sdr_find
#'
#' This script finds a Set of Distinct Representative of a matrix
#'
#' @param x matrix
#' @param block_by_row indication of orientation of the finding
#' @return A Set of Distinct Representative of the matrix
#' @export
#' @examples
#' X = matrix(c(0,1,2,2,0,1), byrow = TRUE,nrow = 2)
#' sdr_find(X)

sdr_find <- function(x, block_by_row = TRUE){
  if(!block_by_row){
    x = t(x)
  }
  nr = dim(x)[1]
  nc = dim(x)[2]

  # Build a record for number of remaining available counts in the play
  sym = sort(unique(c(x)))
  counts = list()
  for (i in sym){
    counts[[as.character(i)]] = nc
  }


  fill = NULL
  ct_update = counts
  cand = x[1,]
  i=1
  while (i <= (nr)){
    prob = NULL
    for (j in 1:length(cand)){prob[j] = ct_update[[as.character(cand[j])]]}
    tmp = cand[sample(prob = prob/(sum(prob)), x = length(cand), size = 1)]
    if (i<nr){
      sdr_remain = length(setdiff(c(x[(i+1):nr,]), c(fill,tmp)))>=nr-i
    }else{
      sdr_remain = TRUE
    }


    if (sdr_remain && (! (tmp%in% fill))){
      fill[i] = tmp

      # update counts
      ct_update[[as.character(tmp)]] = ct_update[[as.character(tmp)]] - 1
      i = i+1
      if (i<=nr){cand = x[i,]} else cand = NULL

      next
    }else if (!sdr_remain){
      fill = NULL
      i = 1
      ct_update = counts
      cand = x[1,]
    }else{
      cand = setdiff(cand,tmp)

    }
    if (length(cand)==0){
      cand = x[1,]
      fill = NULL

      ct_update = counts
      i = 1
    }

  }
  return(fill)

}
