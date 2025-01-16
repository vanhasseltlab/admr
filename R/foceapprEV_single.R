#' @noRd
#'
#'
foceapprEV_single <- function(opts,bi,biseq # Uses the First-Order Conditional Estimation (FOCE) method to approximate expected values for a single individual.
) {
  if (missing(bi)) bi <- matrix(0,1,ncol(opts$p$Omega))
  theta_i <- g_iter(opts,bi)
  Eind <- opts$f(opts$xt,theta_i)
  Etyp <- opts$f(opts$xt,opts$p$beta)
  d_f_d_bi <- jacobiann(function(et) opts$f(opts$xt,g_iter(opts,et)),bi)
  E <- t(c(Eind)-d_f_d_bi %*% t(bi))
  V <- d_f_d_bi %*% opts$p$Omega %*% t(d_f_d_bi)
  if (opts$interact) V <- opts$h(list(V=V,E=Eind),opts$p)$V else V <- opts$h(list(V=V,E=Etyp),opts$p)$V
  res <- list(E=E,V=V,Etyp=Etyp,Eind=Eind,d_f_d_bi=d_f_d_bi,weights=weights)
  res
}
