library(data.table)
library(reshape2)
library(MASS)
#' efficiency
#'
#' This script calculates the relative efficiency of design with given optimality criterion.
#'
#' @param design matrix of a design, b rows and k columns
#' @param v number of treatments and controls
#' @param b number of blocks
#' @param k block size
#' @param ctrl Number of controls required in each block, 1 for 1 control, 0 for no control
#' @param type Type of efficiency
#' @importFrom data.table fsetdiff data.table
#' @importFrom reshape2 melt
#' @importFrom MASS ginv
#' @return The relative efficiency of the given design
#' @export
#' @examples
#' design = generate_eyd(6,5,5,1)
#' efficiency(design,6,5,5,1, type = "D")



efficiency<-function(design,v,b,k,ctrl, type = "D"){

  if(type=="D"){
    info = informat(design,v,b,k,ctrl)
    return(det(info$Md)/det(info$upper_M))
  }else if(type=="A"){
    info =informat(design,v,b,k,ctrl)
    return(sum(diag(info$Md))/sum(diag(info$upper_M)))

  }else if(type=="E"){
    info = informat(design,v,b,k,ctrl)
    return(max(diag(info$Md))/max(diag(info$upper_M)))
  }


}

