#' Runs a test for the decoder after installation
#'
#' Creates a test design, vectors of positives and group test results and runs the decoder.
#'
#' @export
#'

run_test <- function() {

  n=121
  p = 0.1
  testing.mat = 'STD'
  mat.params = list(d = 3)
  X.output = 'VL'
  Y.mode = 'Ct'
  sd = 8581285
  M=c()
  set.seed(sd) #set seed

  # Create Testing Matrix
  message('Creating M')
  if (testing.mat == 'STD') {
    M <- generate_STD_testing_protocol(n, mat.params$d)
  } else if (testing.mat == 'provided'){
    M <- M
  }  else if (testing.mat != 'STD') {
    stop('Not yet implemented.')
  }

  message('Creating X')
  # Generate Viral Loads
  X <- create_positives_vector(n,p=p,output=X.output)

  message('Creating Y')
  # Determine test results
  Y <- determine_Y(M,X,mode = Y.mode)

  message('Running Decoder')
  run_decoder(Y, M, n, mat.params$d)

  message('Finished Running Successfully!')

  return()
}
