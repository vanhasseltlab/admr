#' Generate options for optimization
#'
#'[meancov()] computes the mean and covariance of a matrix, optionally with weights.
#'
#' @param m The matrix
#' @param wt The weights
#'
#' @returns A list with the mean and covariance
#' @export
#' @examples
#' #test

meancov <- function(m,wt) {
  if (missing(wt)) {
    with(cov.wt(m,method="ML"),list(E=center,V=cov))
  } else {
    with(cov.wt(m,wt,method="ML"),list(E=center,V=cov))
  }
}
