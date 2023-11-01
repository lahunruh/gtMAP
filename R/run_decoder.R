#' Optimizes Ct-Values of given status vector
#'
#' Optimizes Ct-Values of given status vector using a stan optimizer with a
#' given stan model using LL and a prior over Ct-Values.
#'
#' @import rstan
#'
#' @param Y: Group Test Ct-values
#' @param M: Design Matrix
#' @param n: Number of samples
#' @param l: Disjunctness of design matrix
#' @param fittedpars: The parameters of the Ct-Value prior distribution. Provide as list of (mean, sd, mean, sd). Preset to list(28.50433862,6.664942792,28.50433862,6.664942792)
#' @param sigma_est: The standard variation we assume Ct-Values can have. Preset to .5
#' @param out_val: Error used in simulated annealing. Default 'SE' (squared error)
#' @param sum_iterations: Number of times decoder is run. Default: 20
#' @param sum_threshold: Decision Threshold as intege. Default: 10
#' @param ctbase: Ct-Value generation Baseline. Default: 50
#' @param ctmult: Multiplier for VLs in Ct-Value generation. Default: 1
#' @param ctmindetect: Ct-value for potentially hidden samples. Default: 25
#'
#' @return Returns the indeces of the positive sampless
#'
#'
#' @export
#'

run_decoder <- function(Y, M, n, l, fittedpars = list(28.50433862,6.664942792,28.50433862,6.664942792), sigma_est = 0.5, out_val = 'SE', sum_iterations = 20, sum_threshold = 10, ctbase = 50, ctmult = 1, ctmindetect = 25){
  gtstan <- ''
  fitpars <- c()
  fitpars$estimate <- fittedpars
  fittedpars <- fitpars
  Y_binary <- c(Y>0)
  M_bin <- M
  inferred_pos_bin <- decoder_disjunct_NCOMP(n, l, M_bin, Y_binary, NCOMP.tolerance = 0)

  if (sum(inferred_pos_bin) > l) {

    M_truncated  <- data.frame(M[Y > 0,])
    M_truncated  <- (M_truncated[,inferred_pos_bin >0])
    Y <- Y[Y >0]
    if (sum(inferred_pos_bin) == 1) {
      out_vec <- decoder_SA(gt_stan, length(M_truncated), M_truncated, Y, fittedpars, prev_est = prev, iter_in = 10,iter_out = 50, temp_red = 'quad_mult', sigma = sigma_est, out_val = out_val)
    } else {
      sum_out_vec <- vector('integer',sum(inferred_pos_bin))
      for (cc in 1:20) {
        gc()
        out_vec <- decoder_SA(gt_stan, ncol(M_truncated), M_truncated, Y, fittedpars, prev_est = prev, iter_in = 10,iter_out = 100, temp_red = 'quad_mult', sigma = sigma_est, out_val = 'SE', ctbase = ctbase, ctmult=ctmult, ctmindetect = ctmindetect)
        sum_out_vec <- sum_out_vec + out_vec
      }
    }
    out_vec <- as.integer(sum_out_vec > sum_threshold)

    M_inf <- M_truncated[,out_vec]
    nam <- names(M_truncated)[which(out_vec>0)]
    for (cc in 1:length(nam)) {
      nam[cc] <- sub('.', '', nam[cc])
    }
    inf_pos <- strtoi(nam)

    X_rest <- X[(inferred_pos_bin >0)]
    X_rest <- X_rest>0
    for (idx in 1:sum(inferred_pos_bin)) {
      if (X_rest[idx] == T) {
        tp_cc[count]  <- sum_out_vec[idx]
        count <- count + 1
      } else {
        fp_cc[count2] <- sum_out_vec[idx]
        count2 <- count2 + 1
      }
    }
    positives <- inf_pos
  } else {
    positives <- list(which(inferred_pos_bin>0))
  }

  return(positives)
}
