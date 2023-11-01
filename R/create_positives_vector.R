#' Create Vector of positive Items
#'
#' Takes in number of items and either a probability of an item being positive
#' or a total number of positives. Returns a length n 0-1 or 0-VL vector of negatives
#' and positives.
#'
#' Functionalities for changing VL generation to be implemented later.
#'
#' @import stats
#'
#' @param n: Number of items
#' @param mode: 'iid' for i.i.d. positives or 'fixed_d' for a fixed number of
#' positives in random positions
#' @param output: Either a 'binary' vector or a vector with Viral loads ('VL')
#' @param p: Probability of item being positive. For i.i.d. mode p is taken as
#' probability and for fixed d mode d is taken to be ceil(p * n)
#'
#' @return Returns a vector of 0s and either 1s ('binary') or viral loads ('vL')
#'
#' @export
#'

create_positives_vector <- function(n, mode = 'iid', output = 'VL', p = 0.01) {
  X <- vector(mode = "integer", length = n)
  if (mode == 'iid') {
    true_positives <- (stats::runif(n) <= p) # i.i.d.
  } else if (mode == 'fixed_d') {
    d <- (p*n)
    true_positives <- sample(n,d)
  }
  if (output == 'binary') {
    X[true_positives] <- 1
  } else if (output == 'VL') {
    #idx <- sample(length(VL_Pool),length(X[true_positives]),replace = TRUE)
    X[true_positives] <- VL_generator(length(X[true_positives]), 'Cleary', lower_bound = 100, upper_bound = 1e8)#VL_generator(length(X[true_positives]), 'Cleary')
  }
  return(X)
}
