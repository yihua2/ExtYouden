#' informat
#'
#' This script calculates the information matrix of the design
#'
#' @param design matrix of a design, b rows and k columns
#' @param v number of treatments and controls
#' @param b number of blocks
#' @param k block size
#' @param ctrl Number of controls required in each block, 1 for 1 control, 0 for no control
#' @return The information matrix of the design
#' @export
#' @examples
#' \donttest{
#'    design = generate_eyd(6,5,5,1)
#'    informat(design,6,5,5,1)
#' }




informat <- function(design,v,b,k,ctrl){
  if (ctrl ==0){
    r = b*k/v
    k0 = k
    v0 = v
  }else{
    k0 = k - 1
    v0 = v - 1
    r = b*k0/v0
  }
  long = melt(t(design))

  long_b = long[,2:3]
  long_r = long[,c(1,3)]

  N_b = t(table(long_b)) # incidence matrix for block(plate)
  N_r = t(table(long_r)) # incidence matrix for row

  # Calculate  information matrix for beta
  if (ctrl ==0){
    Cd = diag(c(rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(rep(r,v0)),ncol = 1)%*%matrix(c(rep(r,v0)),ncol = v)/b/k
  }else{
    Cd = diag(c(b,rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(b,rep(r,v0)),ncol = 1)%*%matrix(c(b,rep(r,v0)),ncol = v)/b/k

  }
  # Calculate upper bound
  Td = matrix(0,nrow = b*k, ncol = v)
  U = matrix(0,nrow = b*k, ncol = b)
  P = matrix(0,nrow = b*k, ncol = k)
  for(i in 1:NROW(long_b)){
    Td[i,long_b[i,]$value+1] = 1
    U[i,ceiling(i/k)]=1
    P[i,ifelse(i%%k==0,k,i%%k)] =1
  }

  upper = t(Td)%*%(diag(rep(1,b*k))-U%*%ginv(t(U)%*%U)%*%t(U))%*%Td
  Md = Cd[2:v,2:v]
  upper_M = upper[2:v,2:v]

  return(list(Cd = Cd,upper = upper,Md = Md,upper_M = upper_M))
}
