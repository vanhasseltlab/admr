#' @noRd

gen_bi <- function(opts,expanded=TRUE) { # Generates random effects bi based on the covariance matrix in the options, with optional expansion.
  Omega <- opts$p$Omega
  if (expanded) Omega <- Omega*opts$omega_expansion
  gen_bi2(Omega,opts$biseq)
}
