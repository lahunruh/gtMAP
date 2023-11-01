#' Optimizes Ct-Values of given status vector
#'
#' Optimizes Ct-Values of given status vector using a stan optimizer with a
#' given stan model using LL and a prior over Ct-Values.
#'
#'
#' @param stan_model: Compiled Stan model to be used
#' @param status_vec: Binary status vector of positive and negative items
#' @param M: Test Matrix
#' @param Y: Output vector with Ct-Values
#' @param fittedpars: The parameters of the Ct-Value prior distribution
#' @param sigma: The standard variation we assume Ct-Values can have. Preset to .5
#'
#' @return Returns the squared average error of the optimization process
#'
#' @export
#'

LL_plus_prior <- function(gt_model, status_vec, M, Y, fittedpars,sigma=0.5, out_val = 'SE',ctbase = 50, ctmult = 1) {
    stan_data <- list(N = ncol(M),
                      M = M,
                      Ct = as.vector(t(Y)),
                      T = nrow(M),
                      n = ncol(M),
                      sigma = sigma,
                      alpha = 1000,
                      beta = 19000,
                      a = fittedpars$estimate[1],
                      b = fittedpars$estimate[2],
                      m = fittedpars$estimate[3],
                      s = fittedpars$estimate[4],
                      w = status_vec,
		      ctbase = ctbase,
		      ctmult = ctmult)
    gt_out_optim <- rstan::optimizing(stanmodels$GT_bayesian_logsumexp_optim, stan_data,hessian = F, verbose = FALSE)
    if (out_val == 'SE') {
      return(gt_out_optim$par['total'])
    } else if (out_val == 'Post') {
      return(abs(gt_out_optim$value))
    }
}
