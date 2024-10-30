#' @noRd
#'
samplogdensfun <- function(bi,p,omega_expansion) { # Computes the log-density of random effects bi under an expanded and unexpanded covariance matrix.
  map_dbl(1:nrow(bi),function(i) {
    proposal <- dmnorm(bi[i,],mean=rep(0,nrow(p$Omega)),sigma=p$Omega*omega_expansion,log=TRUE)$den
    true <- dmnorm(bi[i,],mean=rep(0,nrow(p$Omega)),sigma=p$Omega,log=TRUE)$den
    true-proposal
  })
}
