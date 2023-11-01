#' Decodes with Simulated Annealing Approach
#'
#' Returns the positive items and negative items using the simulated annealing.
#' Parameter is the status vector. At each step the algorithm chooses a random
#' new status vector with Hamming distance of 1 to the current status vector.
#' The Ct-Values of the positive items in either status vector are optimized via
#' a supplied function. The error returned from this determines the evaluation
#' of either status vector, which is then accepted with a probability determined
#' by a sigmoid
#'
#'
#' @param n: Number of subjects
#' @param M: Test Matrix
#' @param Y: Output vector with Ct-Values
#' @param fittedpars: The parameters of the Ct-Value prior distribution
#' @param eval_fun: Handle of the optimization function used
#' @param iter_out: Outer loop iterations
#' @param iter_in: Inner loop iterations
#' @param t: Starting temperature
#' @param alpha: Parameter determining speed of temperature reduction
#' @param temp_red: Determining temperature reduction function. Either 'linear', 'geometric' or 'slow-decrease'. If another value is provided the sequence of temperatures provided via t_seq is used instead.
#' @param t_seq: Sequence of temperature values
#' @param sigma: The standard variation we assume Ct-Values can have. Preset to .5
#'
#' @return Returns a vector of inferred positive items (1 := item is positive, 0 := item is negative)
#'
#' @export
#'

decoder_SA <- function(gt_model, n, M, Y, fittedpars, eval_fun = LL_plus_prior,
                       iter_out = 50, iter_in = 10, t_init = 1000, alpha = 1e9,
                       temp_red = 'linear',t_seq = c(), sigma = 0.5, out_val = 'SE', init = 'random', ctbase = 37, ctmult = 100, ctmindetect = 30) {

    #initialize status vector as all 0s
    status_vec <- vector("logical",n)
    init_num <- sample(n,1)
    init_1s <- sample(n,init_num)
    status_vec[init_1s] <- 1

    #keep track of lowest eval status vector
    lowest_stat_vec <- status_vec
    min_eval <- 1e10

    # outer loop over temperature
    for (t_cc in 1:iter_out) {

      if (t_cc == 1) {
        accept_flag = FALSE
        flag <- FALSE
        t_tmp <- t_init
        cc <- 1
        deltas <- c()
        for (i in 1:100) {
          # generate starting vector
          status_vec <- vector("logical",n)
          init_num <- sample(n,1)
          init_1s <- sample(n,init_num)
          status_vec[init_1s] <- 1
          # choose random entry of status vector and invert
          rand_neighbour <- sample(n,1)
          status_vec_tmp <- status_vec
          status_vec_tmp[rand_neighbour] <- !status_vec[rand_neighbour]
          # get eval of current and proposed status vectors via evaluation function
          cur_eval <- eval_fun(gt_model, status_vec, M, Y,fittedpars,sigma, out_val,ctbase = ctbase, ctmult = ctmult)
          next_eval <- eval_fun(gt_model, status_vec_tmp, M, Y,fittedpars,sigma, out_val,ctbase = ctbase, ctmult = ctmult)
          deltas[i] <-  next_eval - cur_eval # difference in evaluation function
          #
        }
        while (!flag) {
          #################
          accept <-  0
          gc()
          for (i in 1:100) {
          #  # choose acceptance
            accept_prob <- exp(-deltas[i]/t_tmp) # 1/(1+exp(-deltas[i]/t_tmp))
            if (runif(1) <= accept_prob) {
              accept <-  accept + 1
            }
          }
          ##############
          accept_p <- accept/100
          if ((accept_p >= 0.8) & (accept_p <= 0.9)) {
            t0 <- t_tmp
            flag <- TRUE
          } else if (accept_p < 0.8) {
            t_tmp <- t_tmp * 2
          } else {
            t_tmp <- t_tmp / 3
          }
          cc <- cc + 1
        }
        t <- t0

        if (temp_red == 'quad_mult') {
          alpha <- (t0-1)/((iter_out)^2)
        } else if (temp_red == 'linear') {
          alpha <- t0/iter_out
        } else if (temp_red == 'exponential') {
          alpha <- nthroot(1/t0, iter_out)
        } else if (temp_red == 'slow-decrease') {
          alpha <- 1
        } else {
          t <- t_seq[t_cc]
        }
        #return(c(t,alpha))
        if (init == 'all_0'){
          status_vec <- vector("logical",n)
        }

      } else {
        # choose temperature decline
        if (temp_red == 'linear') {
          t <- t0 - (alpha * t_cc)
          if (t < 1) {t <- 1}
        } else if (temp_red == 'exponential') {
          t <- t0 * alpha^t_cc
        } else if (temp_red == 'slow-decrease') {
          t <- t/(1 + alpha * t)
        } else if (temp_red == 'quad_mult') {
          t <- t0 / (1 + alpha * t_cc^2)
        } else {
          t <- t_seq[t_cc]
        }
      }

      # inner loop over possibilities
      for (i in 1:iter_in) {

        #only in first iteration evaluate current status vector
        if ((i == 1) & (t_cc == 1)) {
          cur_eval <- eval_fun(gt_model, status_vec, M, Y,fittedpars,sigma, out_val,ctbase = ctbase, ctmult = ctmult)
        } else if (accept_flag) {
          cur_eval <- next_eval
          accept_flag <- FALSE
        }

        # choose random entry of status vector and invert
        rand_neighbour <- sample(n,1)
        status_vec_tmp <- status_vec
        status_vec_tmp[rand_neighbour] <- !status_vec[rand_neighbour]

        # get eval of current and proposed status vectors via evaluation function
        if (cur_eval < min_eval) {
          lowest_stat_vec <- status_vec
          min_eval <- cur_eval
        }        #track status vector optimum
        next_eval <- eval_fun(gt_model, status_vec_tmp, M, Y,fittedpars,sigma, out_val,ctbase = ctbase, ctmult = ctmult)
        if (i == iter_in) {
          if (next_eval < min_eval) {
            lowest_stat_vec <- status_vec
            min_eval <- next_eval
          }
        }
        delta_E <-  next_eval - cur_eval # cur_eval -next_eval # difference in evaluation function

        # choose acceptance
        accept_prob <- exp(-delta_E/t) # 1/(1+exp(-delta_E/t))
        if (!is.infinite(accept_prob)) {
          if (runif(1) <= accept_prob) {
            status_vec <- status_vec_tmp
            accept_flag <- TRUE
          }
        }
      }
    }


    stan_data <- list(N = ncol(M),
                      M = M,
                      Ct = as.vector(t(Y)),
                      T = nrow(M),
                      n = ncol(M), #n
                      sigma = sigma,#1,
                      alpha = 1000,
                      beta = 19000,#
                      a = fittedpars$estimate[1],
                      b = fittedpars$estimate[2],
                      m = fittedpars$estimate[3],
                      s = fittedpars$estimate[4],
                      w = status_vec,
		      ctbase = ctbase,
		      ctmult = ctmult)
    gt_out_optim <- rstan::optimizing(stanmodels$GT_bayesian_logsumexp_optim, stan_data,hessian = F, verbose = FALSE)#list(list(X = rep(1000,stan_data$N))), chains = 1)

    ctind <- c()
    for (entry in 1:length(status_vec)) {
      if (status_vec[entry] > 0) {
        ctind[entry] <- gt_out_optim$par[paste('Ct_ind[',entry,']',sep='')]
      } else {
        ctind[entry] <- 0
      }
    }

    t <- length(as.vector(t(Y)))
    N <- ncol(M)

    ctpred = vector('double',length = t)


    for (i in 1:t) {
      for (j in 1:N) {
        if (M[i,j] > 0) {
          if (status_vec[j] > 0) {
            ctpred[i] = ctpred[i] + M[i,j] * 2^(ctbase - ctind[j]) * ctmult # only now added ctbase and ctmult
          }
        }
      }
    }
    for (i in 1:t) {
      if (ctpred[i] > 0) {
        ctpred[i] = ctbase - log(ctpred[i]/(ctmult * sum(M[i,])), base = 2)
      }
    }

    new_stat_vec <- status_vec
    for (entry in 1:length(ctind)) {
      if (status_vec[entry] == 0) {
        status_vec_tmp <- status_vec
        status_vec_tmp[entry] <- 1
        ctind_bu <- ctind
        ctind[entry] <- ctmindetect

        ctpred_new = vector('double',length = t)

        for (i in 1:t) {
          for (j in 1:N) {
            if (M[i,j] > 0) {
              if (status_vec_tmp[j] > 0) {
                ctpred_new[i] = ctpred_new[i] + M[i,j] * 2^(ctbase - ctind[j]) * ctmult
              }
            }
          }
        }
        for (i in 1:t) {
          if (ctpred_new[i] > 0) {
              ctpred_new[i] = ctbase - log(ctpred_new[i]/(ctmult * sum(M[i,])), base = 2)
          }
        }

        outflag <- TRUE
        for (output in 1:length(ctpred)) {
          if (abs(ctpred_new[output] - ctpred[output]) > sigma) {
            outflag <- FALSE
          }
        }
        if (outflag) {
          new_stat_vec[entry] <- 1
        }
        ctind <- ctind_bu
      }
    }
    return(new_stat_vec)
}
