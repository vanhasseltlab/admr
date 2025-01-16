#' @noRd

MCapprEV <- function(opts) { # Performs a Monte Carlo approximation to compute the expected values (mean and covariance) of the model outputs based on random effects.
  p <- opts$p
  ## p is normal-scale parameters
  bi <- gen_bi(opts)
  theta_i <- g_iter(opts,bi)
  m <- opts$f(opts$time,theta_i)
  if (opts$omega_expansion==1) {
    r <- meancov(m)
  } else {
    wt <- samplogdensfun(bi,p,opts$omega_expansion)
    r <- meancov(m,logdens2wt(wt))
  }
  opts$h(r,opts$p)
}
