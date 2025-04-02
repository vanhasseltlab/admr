#' Create a covariance matrix with specified diagonal and off-diagonal values
#'
#' @description
#' `omegas` creates a covariance matrix with specified diagonal and off-diagonal values.
#' This function is useful for creating initial or fixed covariance matrices for
#' random effects in pharmacometric models.
#'
#' @param diag Value for the diagonal elements of the matrix. This represents the
#'             variance of each random effect.
#' @param offdiag Value for the off-diagonal elements of the matrix. This represents
#'                the covariance between random effects.
#' @param n_om Size of the matrix (number of random effects).
#'
#' @returns A symmetric matrix of size `n_om` x `n_om` with:
#' \itemize{
#'   \item Diagonal elements equal to `diag`
#'   \item Off-diagonal elements equal to `offdiag`
#' }
#'
#' @details
#' The function creates a symmetric matrix where:
#' \itemize{
#'   \item All diagonal elements are set to `diag`
#'   \item All off-diagonal elements are set to `offdiag`
#' }
#'
#' This is particularly useful for creating:
#' \itemize{
#'   \item Initial covariance matrices for optimization
#'   \item Fixed covariance matrices for simulation
#'   \item Simple covariance structures for testing
#' }
#'
#' @examples
#' # Create a 3x3 matrix with diagonal = 0.09 and off-diagonal = 0
#' omega1 <- omegas(0.09, 0, 3)
#'
#' # Create a 5x5 matrix with diagonal = 0.09 and off-diagonal = 0.01
#' omega2 <- omegas(0.09, 0.01, 5)
#'
#' @export
omegas <- function(diag, offdiag, n_om) { # Creates a matrix with specified diagonal and off-diagonal values.
  ## a simple utility function
  ## returns a matrix of size n_om, with constant diag and offdiag values
  m <- matrix(offdiag,n_om,n_om)
  diag(m) <- diag
  m
}
