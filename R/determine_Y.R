#' Determine experiment result vector
#'
#' Determine the experiment result vector (either in binary or Ct-Values) from
#' testing matrix M and input vector (samples) X. Only works for constant
#' weight codes at the moment.
#'
#' @param M: Testing Matrix of 1s and 0s
#' @param X: Vector with either positives as 1s or with individual Viral Loads
#' @param mode: 'binary' returns binary outcome vector (test positive or negative),
#' 'Ct' returns Ct-Valued outcome vector
#' @param Se: Test Sensitivity (not implemented for 'Ct)
#' @param Sp: Test Specificity (not implemented for 'Ct)
#' @param spread: The maximum deviation of a Ct-Value from its expected value.
#'
#' @return Returns test result vector in form specified via input
#'
#' @export
#'

determine_Y <- function(M, X, mode = 'Ct',Se = 1, Sp = 1, spread = 4) {
  n = ncol(M)
  if (mode == 'binary') {
    Y <- vector(mode = "integer",length = nrow(M))
    for (j in 1:nrow(M)) {
      for (i in 1:n) {
        if (M[j,i] == 1 & X[i] == 1) { #item positive and in the test
          Y[j] <-  1
        }
      }
    }
    # add some measurement noise
    #Y[which(Y == 1)] <- (runif(length(Y[which(Y == 1)])) <= Se)
    #Y[which(Y == 0)] <- (runif(length(Y[which(Y == 0)])) >  Sp)
    return(Y)
  } else if (mode == 'Ct') {
    Y <- vector(mode = "double",length = nrow(M))
    Y <- M %*% cbind(X)
    for (i in 1:nrow(M)) {
      w <- sum(M[i,])#length(which(M[i,] != 0))
      Y[i] <- Y[i]/w #devide each test VL by weight of test to get accurate dilution
    }

    #Y <- GT.Testbed::VL_to_Ct(Y) # added GT.Testbed::... here, does it make any difference in R packages whether I do or not?
    # dropped vectorized version for now because of 0s.
    for (i in 1:nrow(Y)) {
      Y[i] <- VL_to_Ct(Y[i],spread = spread, mode = 'dist')
    }

    #### Add false positives here #####
    # take some negative results with probability q and turn them positive by assigning a very high Ct-Value (38-40)
    ###################################

    #### Add false negatives here #####
    # simply all Ct's > 40
    ###################################

    return(Y)
  }
}
