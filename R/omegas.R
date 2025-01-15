#' Generate a dataset
#'
#'[omegas()] creates a matrix with specified diagonal and off-diagonal values.
#'
#'
#'
#' @param diag A value for the diagonal
#' @param offdiag A value for the off-diagonal
#' @param n_om The size of the matrix
#'
#' @returns A matrix with specified diagonal and off-diagonal values.
#' @export
#' @examples
#' #test


omegas <- function(diag,offdiag,n_om) { # Creates a matrix with specified diagonal and off-diagonal values.
  ## a simple utility function
  ## returns a matrix of size n_om, with constant diag and offdiag values
  m <- matrix(offdiag,n_om,n_om)
  diag(m) <- diag
  m
}
