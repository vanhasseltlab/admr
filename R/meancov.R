#' Compute mean and covariance of a matrix
#'
#' @description
#' `meancov` computes the mean and covariance of a matrix, optionally with weights.
#' This function is used to convert raw data into aggregate form (mean and covariance)
#' for use in aggregate data modeling.
#'
#' @param m A numeric matrix or data frame containing the observations.
#'          Each row represents an individual, and each column represents a time point.
#' @param wt Optional vector of weights for each observation. If not provided,
#'           all observations are weighted equally.
#'
#' @returns A list containing:
#' \itemize{
#'   \item `E`: Vector of means for each time point
#'   \item `V`: Covariance matrix of the observations
#' }
#'
#' @details
#' The function computes:
#' \itemize{
#'   \item The mean of each column (time point) using `colMeans`
#'   \item The covariance matrix using `cov.wt` with maximum likelihood estimation
#' }
#'
#' If weights are provided, they are used in both the mean and covariance calculations.
#' The function uses maximum likelihood estimation for the covariance matrix, which
#' is appropriate for aggregate data modeling.
#'
#' @examples
#' # Create a matrix of observations
#' m <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'
#' # Compute mean and covariance without weights
#' result <- meancov(m)
#'
#' # Compute mean and covariance with weights
#' weights <- runif(10)
#' result_weighted <- meancov(m, weights)
#'
#' @export

meancov <- function(m,wt) {
  if (missing(wt)) {
    with(cov.wt(m,method="ML"),list(E=center,V=cov))
  } else {
    with(cov.wt(m,wt,method="ML"),list(E=center,V=cov))
  }
}
