#' @noRd

covlogit <- function(x) { # Performs logit and inverse logit transformations for values constrained to the [-1,1] domain.
  x2 <- (x+1)/2
  log(x2/(1-x2))
}
