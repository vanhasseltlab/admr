#' @noRd
#'
lowtoup <- function(m) { # Copies the upper triangular part of a matrix to the lower triangular
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  m
}
