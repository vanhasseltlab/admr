#' @noRd

uptolow <- function(m) { # Copies the upper triangular part of a matrix to the lower triangular and vice versa.
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}
