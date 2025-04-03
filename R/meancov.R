#' Compute mean and covariance of a matrix
#'
#' @description
#' `meancov` computes the mean and covariance of a matrix, optionally with weights.
#' This function is used to convert raw data into aggregate form (mean and covariance)
#' for use in aggregate data modeling.
#'
#' @param m A numeric matrix or data frame containing the observations.
#'          Each row represents an individual, and each column represents a time point.
#'          For pharmacometric data, columns typically represent concentration measurements
#'          at different time points.
#' @param wt Optional vector of weights for each observation. If not provided,
#'           all observations are weighted equally. Weights can be used to account for
#'           different sample sizes or reliability of different data sources.
#'
#' @returns A list containing:
#' \itemize{
#'   \item `E`: Vector of means for each time point (population typical values)
#'   \item `V`: Covariance matrix representing the variability between individuals
#' }
#'
#' @details
#' The function computes:
#' \itemize{
#'   \item The mean of each column (time point) using `colMeans` for unweighted data
#'         or weighted means for weighted data
#'   \item The covariance matrix using `cov.wt` with maximum likelihood estimation,
#'         which provides unbiased estimates of the population covariance
#' }
#'
#' The maximum likelihood estimation method is used because:
#' \itemize{
#'   \item It provides unbiased estimates of the covariance matrix
#'   \item It is appropriate for aggregate data modeling where we want to estimate
#'         population parameters
#'   \item It handles both balanced and unbalanced designs through optional weights
#' }
#'
#' Key features:
#' \itemize{
#'   \item Handles missing data automatically through the underlying `cov.wt` function
#'   \item Provides numerically stable computations
#'   \item Can be used with both raw PK data and simulated data
#'   \item Supports weighted calculations for meta-analysis or combined analysis
#' }
#'
#' @examples
#' # Create a matrix of concentration measurements
#' # 10 subjects measured at 10 time points
#' m <- matrix(rnorm(100), nrow = 10, ncol = 10)
#'
#' # Compute unweighted mean and covariance
#' # Useful for single-study analysis
#' result <- meancov(m)
#'
#' # Compute weighted mean and covariance
#' # Useful for meta-analysis or when combining studies
#' weights <- runif(10)  # weights could represent study sizes
#' result_weighted <- meancov(m, weights)
#'
#' @export

meancov <- function(m, wt) {
  # Check if weights are provided
  if (missing(wt)) {
    # Unweighted case: compute ML estimates of mean and covariance
    # - center: contains column means (typical values at each time point)
    # - cov: contains ML estimate of covariance matrix
    with(cov.wt(m, method="ML"), list(E=center, V=cov))
  } else {
    # Weighted case: incorporate observation weights
    # - Weights could represent:
    #   * Study sizes in meta-analysis
    #   * Precision of measurements
    #   * Relative importance of observations
    with(cov.wt(m, wt, method="ML"), list(E=center, V=cov))
  }
}
