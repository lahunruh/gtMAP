#' Generate fake Viral Loads
#'
#' Generates fake viral loads randomly either from a uniform distribution or
#' by the process described by Cleary \[1\].
#'
#' @import stats
#'
#' @param num_samples: Number of positive swap samples to be taken
#' @param distribution: 'Uniform' or 'Cleary'
#' @param max_time: last time a sample can be taken
#' @param alpha...: See Cleary et al., 2020
#' @param lower_bound: Lower bound for uniform distribution
#' @param upper_bound: Upper bound for uniform distribution
#'
#' @return Returns vector of viral loads of length n.
#'
#' @references \[1\] Brian Cleary et al. \(2020\), Using viral load and epidemic dynamics to optimize
#' pooled testing in resource constrained settings. medRxiv, doi: 10.1101/2020.05.01.20086801,
#' version from October 6, 2020
#'
#' @export
#'

VL_generator <- function(num_samples, distribution,
                         alpha = 7.98, t_w = 15.1, t_inc = 5.52, t_g = 0, t_p = 3,
                         t_inf = 0, lower_bound = 0, upper_bound = 1e8){
 # choose distributions of VLs
  if (distribution == 'Uniform') {
    return(stats::runif(num_samples) * (upper_bound-lower_bound) + lower_bound)
  } else if (distribution == "Cleary") {
    samples <- vector(mode='double',length = num_samples)
    for (i in 1:num_samples) {
      time_of_sample <- stats::runif(1) * (t_inc-t_p+t_w)#t_w
      samples[i] <- viral_load_by_cleary(time_of_sample,alpha,t_inf,t_g,t_p,t_w,t_inc)
    }
    return(samples)
  }

}

viral_load_by_cleary <- function(sample_time,alpha,t_inf,t_g,t_p,t_w,t_inc) {
  if ((sample_time >= t_inf) & (sample_time <= t_g)) { #redundant for us as t_g will be 0
    return(0)
  } else if ((sample_time > t_g) & (sample_time <= t_p)) {
    return(10^(alpha/t_p * sample_time))
  } else if ((sample_time > t_p) & (sample_time <= (t_inc-t_p+t_w))) {
    return( 10^(alpha - (alpha/(t_inc-t_p+t_w)*(sample_time-t_p))) )
  }
}
