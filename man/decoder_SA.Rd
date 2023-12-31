% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decoder_SA.R
\name{decoder_SA}
\alias{decoder_SA}
\title{Decodes with Simulated Annealing Approach}
\usage{
decoder_SA(
  gt_model,
  n,
  M,
  Y,
  fittedpars,
  eval_fun = LL_plus_prior,
  iter_out = 50,
  iter_in = 10,
  t_init = 1000,
  alpha = 1e+09,
  temp_red = "linear",
  t_seq = c(),
  sigma = 0.5,
  out_val = "SE",
  init = "random",
  ctbase = 37,
  ctmult = 100,
  ctmindetect = 30
)
}
\arguments{
\item{n:}{Number of subjects}

\item{M:}{Test Matrix}

\item{Y:}{Output vector with Ct-Values}

\item{fittedpars:}{The parameters of the Ct-Value prior distribution}

\item{eval_fun:}{Handle of the optimization function used}

\item{iter_out:}{Outer loop iterations}

\item{iter_in:}{Inner loop iterations}

\item{t:}{Starting temperature}

\item{alpha:}{Parameter determining speed of temperature reduction}

\item{temp_red:}{Determining temperature reduction function. Either 'linear', 'geometric' or 'slow-decrease'. If another value is provided the sequence of temperatures provided via t_seq is used instead.}

\item{t_seq:}{Sequence of temperature values}

\item{sigma:}{The standard variation we assume Ct-Values can have. Preset to .5}
}
\value{
Returns a vector of inferred positive items (1 := item is positive, 0 := item is negative)
}
\description{
Returns the positive items and negative items using the simulated annealing.
Parameter is the status vector. At each step the algorithm chooses a random
new status vector with Hamming distance of 1 to the current status vector.
The Ct-Values of the positive items in either status vector are optimized via
a supplied function. The error returned from this determines the evaluation
of either status vector, which is then accepted with a probability determined
by a sigmoid
}
