#' @noRd

nllfun <- function(obsEV,predEV,invpredV,n=1) {   # Computes the negative log-likelihood given observed and predicted expected values and their covariance.
  if (missing(invpredV)) invpredV <- solve(predEV$V)
  resids <- c(obsEV$E-predEV$E) ## force to simple vector
  c(1/2*n*(log(det(predEV$V))+sum(diag(obsEV$V %*% invpredV))+t(resids) %*% invpredV %*% resids))
}
