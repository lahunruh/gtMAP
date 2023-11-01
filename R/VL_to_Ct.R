#' Viral Load to Ct-Value
#'
#' Converts a Viral Load into a Ct-value, infusing noise either from
#' a uniform spread or a distribution. Default values are loosely taken from
#' data in \[1\].
#'
#' @import stats
#'
#' @param VL: viral Load number.
#' @param mode: 'range' (for choosing noise from uniform range) or 'dist' for
#' choosing noise via a distribution.
#' @param spread: spread of Ct-values, either a range for mode == 'range'
#' or a standard deviation for mode == 'dist'.
#' @param baseline_Ct: baseline Ct at baseline Viral Load as starting point
#' for dilution based approximation
#' @param baseline_VL: Baseline Viral load generating baseline Ct as starting point
#' for dilution based approximation
#'
#' @return Returns a Ct-Value
#'
#' @references \[1\] Chantal BF Vogels et al. Analytical sensitivity and efficiency
#' comparisons of SARS CoV 2 RT qPCR primer probe sets. In:
#' Naturemicrobiology 5.10 (2020), pp. 1299 to 1305.
#'
#' @note Doesn't allow for false positives.
#'
#' @export
#'

VL_to_Ct <- function(VL, mode = 'range', spread = 2, baseline_Ct = 37, baseline_VL = 1,det=FALSE) {
  if (VL == 0) {
    return(0)
  }
  if (det) {
    return(baseline_Ct - log(VL/baseline_VL, base = 2))
  } else {
    if (mode == 'range') {
      Ct <- baseline_Ct - log(VL/baseline_VL, base = 2)
      # Just as note, if vector with 0s is given these of course go to Inf
      noise <- stats::runif(length(VL)) * spread - (spread/2)
      return(Ct + noise)
    } else if (mode == 'dist') {
      Ct <- baseline_Ct - log(VL/baseline_VL, base = 2)
      noise <- stats::rnorm(1,mean = 0, sd = spread)
      return(Ct + noise)
    }
  }
}
