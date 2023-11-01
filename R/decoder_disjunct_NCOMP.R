#' Decodes disjunct matrices using NCOMP
#'
#' Returns the positive items using NCOMP decoder \[1\], where positive items are defined by being in at
#' least d-NCOMP.tolerance tests (where d is the disjunctness of the
#' test matrix).
#'
#' @param n: Number of subjects
#' @param k: disjunctness
#' @param M: Test Matrix
#' @param Y: Output vector
#' @param NCOMP.tolerance: Tolerance of decoder determining number of tests
#' an item is in that are allowed to be negative while still declaring the item
#' positive. Increases false positive but decreases false negative rate.
#'
#' @return Returns a vector of inferred positive items (1 := item is positive, 0 := item is negative)
#'
#' @references Chan, CL et al. Non-adaptive probabilistic group testing with
#' noisy measurements: Near-optimal bounds with efficient algorithms. In:
#' Forty-Ninth Annual Allerton Conference Allerton House, UIUC, Illinois, USA
#' September 28 - 30, 2011.
#'
#' @export
#'


decoder_disjunct_NCOMP <- function(n, d, M, Y, NCOMP.tolerance = 0) {
  k <- d + 1
  postests <- vector(mode = "integer", length = n)
  inferred.positives <- vector(mode = "integer", length = n)
  for (i in 1:n) {
    postests[i] <- 0
    for (j in 1:nrow(M)) {
      if ((M[j,i] == 1) & (Y[j] == 1)) {
        postests[i] <- postests[i] + 1
      }
    }
    if (postests[i] > (sum(M[,i]) - 1 - NCOMP.tolerance)) {#(k-1-NCOMP.tolerance)) {
      inferred.positives[i] <- 1
    } else {
      inferred.positives[i] <- 0
    }
  }
  return(inferred.positives)
}
