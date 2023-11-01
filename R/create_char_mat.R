#' Create Character Matrix
#'
#' Creates a matrix whose characters correspond to the labeling of the laboratory system.
#'
#' @param nrows: Number of rows
#' @param ncols: Number of columns
#'
#' @return Returns character matrix
#'
#' @export
#'

#################################################################
#                         FUNCTIONS                             #
#################################################################

create.char.mat <- function(nrows,ncols) {
  char.mat <- matrix( nrow = nrows, ncol = ncols)
  for (i in 1:nrow(char.mat)) {
    for (j in 1:ncol(char.mat)) {
      char.mat[i,j] <- paste(intToUtf8((64+i)),as.character(j),sep="")
    }
  }
  return(char.mat)
}
