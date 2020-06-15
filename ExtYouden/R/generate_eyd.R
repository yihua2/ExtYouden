library(data.table)
library(ibd)

#' generate_eyd
#'
#' This script generates an Extended Youden design with ctrl controls required in each block, default 1
#'
#' @param v number of treatments and controls
#' @param b number of blocks
#' @param k block size
#' @param ctrl number of controls required in each block, 1 for 1 control, 0 for no control
#' @return An Extended Youden Design
#' @importFrom data.table fsetdiff
#' @export
#' @examples
#' \donttest{
#'    generate_eyd(6,20,3,0)
#' }

generate_eyd <- function(v,b,k, ctrl){
    #adding 1 control to each block

    if (ctrl == 0) {
      r = b * k / v
      try(if (((b * k) %% v > 0) |
              ((r * (k - 1)) %% (v - 1) > 0) |
              (b < v))
        stop("No valid design available"))
      k0 = k
      v0 = v
      t = r %% k
      lambda = r * (k - 1) / (v - 1)

      D = ibd::bibd(v, b, r, k, lambda, pbar = FALSE)$design
      D1 = D
      s = (k * v - t * v) / k # assume starting from a bib design


      ### add dummy columns
      dummy = NULL
      for (i in 1:v0) {
        dummy = c(dummy, rep(i, k - t))
      }
      D2 = rbind(D1, matrix(dummy, nrow = s))


      ### re-label trtments/control
      D3 = D2

      for (j in 1:v0) {
        # 0 no control
        # shuffle the indices of original blocks in which treatment j exist
        ind  = sample(which((
          apply(D3 == j, FUN = sum, MARGIN = 1) * c(1:(b + s)) > 0
        ) * (
          apply(D3 == j, FUN = sum, MARGIN = 1) * c(1:(b + s)) < b + 1
        ) == 1))

        if (j > 0) {
          m = floor(r / k)
          r1 = r + k - t
          tt = t
        } else {
          m = floor(b / k)
          r1 = b + k - t0
          tt = t0
        }

        if (m == 1) {
          if (tt > 0) {
            D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(m * k + 1):min((m + 1) * k,length(ind))], ] == j,
              j + m * v ,
              D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ]
            )

          }
          D3[(b + 1):(b + s), ] = ifelse(
            D3[(b + 1):(b + s), ] == j,
            j + v ,
            D3[(b + 1):(b + s), ]
          )
        } else if (m > 1) {
          for (i in 1:(m - 1)) {
            D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ] == j,
              j + i * v ,
              D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ]
            )
          }
          if (tt > 0) {
            D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(m *k + 1):min((m + 1) * k, length(ind))], ] == j,
              j + m * v ,
              D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ]
            )

          }
          D3[(b + 1):(b + s), ] = ifelse(
            D3[(b + 1):(b + s), ] == j,
            j + m *v ,
            D3[(b + 1):(b + s), ]
          )
        }

      }
      #fill row
      new_fill = NULL
      j = 1
      remain = D3

      while (j <= k) {
        new_fill = cbind(new_fill, sdr_find(remain)) # in sdr_find_new
        remain = matrix(NA, nrow = b + s, ncol = k - dim(new_fill)[2])
        for (i in 1:(b + s)) {
          remain[i, ] = c(
            fsetdiff(data.table(D3[i, ]), data.table(new_fill[i, ]), all = TRUE),
            recursive = TRUE,
            use.names = FALSE
          )
        }

        if (j < k && sdr_exist(remain) == TRUE) {
          #sdr_exist in code_bioassay
          j = j + 1
        }
        else if (j == k)
          break
        else  if (j < k &&  sdr_exist(remain) == FALSE) {
          j = 1
          new_fill = NULL
        }
      }


      #transform back to original treatment label
      a = new_fill
      if (r / k > 0) {
        #|| b/k>0){
        a = (a %% v)[1:b, ]
      }

    } else {
      # 1 control
      k0 = k - 1
      v0 = v - 1
      r = b * k0 / v0
      try(if (((b * k0) %% v0 > 0) |
              ((r * (k0 - 1)) %% (v0 - 1) > 0) |
              (b < v0))
        stop("No valid design available"))

      t = r %% k
      t0 = b %% k

      lambda = r * (k0 - 1) / (v0 - 1)
      D = ibd::bibd(v0, b, r, k0, lambda, pbar = FALSE)$design

      #add control
      D1 = cbind(rep(0, b), D)
      s = (k * v - t * v0 - b %% k) / k # assume starting from a bib design

      #add dummy columns
      dummy = rep(0, k - b %% k)
      for (i in 1:v0) {
        dummy = c(dummy, rep(i, k - t))
      }
      D2 = rbind(D1, matrix(dummy, nrow = s))


      #re-label trtments/control
      D3 = D2
      for (j in 0:v0) {
        # shuffle the indices of original blocks in which treatment j exist
        ind  = sample(which((
          apply(D3 == j, FUN = sum, MARGIN = 1) * c(1:(b + s)) > 0
        ) * (
          apply(D3 == j, FUN = sum, MARGIN = 1) * c(1:(b + s)) < b + 1
        ) == 1))
        if (j > 0) {
          m = floor(r / k)
          r1 = r + k - t
          tt = t
        } else {
          m = floor(b / k)
          r1 = b + k - t0
          tt = t0
        }

        if (m == 1) {
          if (tt > 0) {
            D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] == j,
              j + m * v ,
              D3[ind[(m * k +1):min((m + 1) * k, length(ind))], ]
            )
          }
          D3[(b + 1):(b + s), ] = ifelse(
            D3[(b + 1):(b + s), ] == j,
            j + v ,
            D3[(b +1):(b + s), ]
          )
        } else if (m > 1) {
          for (i in 1:(m - 1)) {
            D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ] == j,
              j + i * v ,
              D3[ind[(i * k + 1):min((i + 1) * k, length(ind))], ]
            )
          }
          if (tt > 0) {
            D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] = ifelse(
              D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ] == j,
              j + m * v ,
              D3[ind[(m * k + 1):min((m + 1) * k, length(ind))], ]
            )
          }
          D3[(b + 1):(b + s), ] = ifelse(
            D3[(b + 1):(b + s), ] == j,
            j + m * v ,
            D3[(b + 1):(b + s), ]
          )
        }

      }

      #fill row
      new_fill = NULL
      j = 1
      remain = D3

      while (j <= k) {
        new_fill = cbind(new_fill, sdr_find(remain)) # using new sdr_find function
        remain = matrix(NA, nrow = b + s, ncol = k - dim(new_fill)[2])
        for (i in 1:(b + s)) {
          remain[i, ] = c(
            fsetdiff(data.table(D3[i, ]), data.table(new_fill[i, ]), all = TRUE),
            recursive = TRUE,
            use.names = FALSE
          )
        }

        if (j < k && sdr_exist(remain) == TRUE) {
          j = j + 1
        }
        else if (j == k)
          break
        else  if (j < k &&  sdr_exist(remain) == FALSE) {
          j = 1
          new_fill = NULL
        }
      }


      #transform back to original treatment label
      a = new_fill
      if (r / k > 0 || b / k > 0) {
        a = (a %% v)[1:b, ]
      }
    }


    return(a)
  }

