#' @noRd

gen_bi2 <- function(Omega,biseq) { # Generates random effects bi based on the covariance matrix in the options, with optional expansion.
  ## error handling for the  special case of zero omegas
  if (any(diag(Omega)==0))
    diag(Omega)[diag(Omega)==0] <- 1e-10
  biseq %*% chol(Omega)
}
