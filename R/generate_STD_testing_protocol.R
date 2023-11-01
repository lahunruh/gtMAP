#' Shifted Transversal Design - Generate Test Matrix - ADAPTED FOR LAB TESTING
#'
#' Generates the test matrix for the shifted transversal design \[1\].
#' Columns are subjects and rows are tests. Requires the number of subjects
#' to be a square number. If not, the next lowest square number will be taken
#' for the design and all the left over subjects will be grouped in one "control"
#' group at the end.
#'
#' @param n:        number of subjects/items (preferably square number)
#' @param d:        desired disjunctiveness
#' @param control:  T or F, turns automatically to TRUE if n isn't a square.
#'
#' @return Returns a d-disjunct Test Matrix from STD design with n columns and (d+1) * sqrt(n) rows (if control == F).
#'
#' @references \[1\] A new pooling strategy for high-throughput screening:
#' the Shifted Transversal Design. Thierry-Mieg, N. (2006), BMC Bio Inf, 7:28,
#' doi: 10.1186/1471-2105-7-28
#'
#' @export
#'

generate_STD_testing_protocol <- function(n, d, control = FALSE, sq =  c(7,0,1,4)) {
    if (n<1) {
      stop("Please provide positive n")
    }
    if (!is.wholenumber(sqrt(n))) { # check if n is square
      control = TRUE
    }
    k <- d + 1
    if (!control){ # if n is square
      q <- sqrt(n)
      M <- matrix(0,nrow = q*k, ncol = n)
      gamma <- 1 # just a reminder for now
      sequence <- c(q,0:(k-2))
      if (d == 3) { # customise for our test protocol
        sequence <-sq
      }
      for (j in 0:(k-1)) {#
        seq <- sequence[j+1]
        for (i in 0:(n-1)) {
          # assign rotated submatrix to M
          M[(q*j + 1):(q*j + q),i+1] <- rbind(rotate(q,s_i_j(q,i,seq)))
        }
      }
    } else { # if n isn't square
      q <- floor(sqrt(n))
      control.subs <- n - q^2
      M <- matrix(0,nrow = q*k + 1, ncol = n)
      gamma <- 1 # just a reminder for now
      for (j in 0:(k-1)) {
        for (i in 0:(n-1-control.subs)) {
          # assign rotated submatrix to M
          M[(q*j + 1):(q*j + q),i+1] <- rbind(rotate(q,s_i_j(q,i,j)))
        }
      }
      # add one control group with leftover subjects
      M[q*k + 1, ] <-  rbind(integer(n))
      M[q*k + 1, (q^2+1):n] <- 1
    }
    return(M)
}

s_i_j <- function(q,i,j,gamma = 1) {
  if (j<q) {
    s = i
    for (cc in 1:gamma) {
      s = s + j^cc * floor(i/q^cc)
    }
  } else {
    s = floor(i/q^gamma)
  }
  return(s)
}

rotate <- function(q,s) {
  C <- integer(q)
  C[s%%q+1] <- 1
  return(C)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}
