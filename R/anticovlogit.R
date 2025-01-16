#' @noRd

anticovlogit <- function(x) { # Performs logit and inverse logit transformations for values constrained to the [-1,1] domain.
  2*exp(x)/(exp(x)+1)-1
}
