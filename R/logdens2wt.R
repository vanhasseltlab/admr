#' @noRd

logdens2wt <- function(x) { # Converts log-densities into normalized weights for importance sampling.
  x2 <- exp(x-max(x)) ## standardize to avoid floating point errors
  x2[!is.finite(x2)] <- max(x2[is.finite(x2)])
  x2/sum(x2)
}
