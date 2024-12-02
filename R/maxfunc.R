#' @noRd
#'
#'

maxfunc <- function(opts) {
  bi <- gen_bi(opts)
  origbeta <- opts$p$beta
  rawpreds <- opts$f(opts$xt,g_iter(opts,bi))
  NAfilters <- apply(rawpreds,1,function(i) !any(is.na(i)))
  adjust <- !any(is.na(opts$single_betas))
  propdens <- dmnorm(bi,mean=rep(0,nrow(opts$p$Omega)),sigma=opts$p$Omega*opts$omega_expansion,log=TRUE)$den
  function(pnew,returnEV=FALSE) {
    pneww <- opts$ptrans(pnew)
    newdens <- opts$p_thetai(pneww,origbeta,bi)
    wttot <- logdens2wt(newdens[NAfilters]-propdens[NAfilters])
    EVnow <- with(cov.wt(rawpreds[NAfilters,],wttot,method="ML"),list(E=center,V=cov))
    EVnow <- opts$h(EVnow,pneww)#' @noRd
#'
#'

maxfunc <- function(opts) {
  bi <- gen_bi(opts)
  origbeta <- opts$p$beta
  rawpreds <- opts$f(opts$xt,g_iter(opts,bi))
  NAfilters <- apply(rawpreds,1,function(i) !any(is.na(i)))
  adjust <- !any(is.na(opts$single_betas))
  propdens <- dmnorm(bi,mean=rep(0,nrow(opts$p$Omega)),sigma=opts$p$Omega*opts$omega_expansion,log=TRUE)$den
  function(pnew,returnEV=FALSE) {
    pneww <- opts$ptrans(pnew)
    newdens <- opts$p_thetai(pneww,origbeta,bi)
    wttot <- logdens2wt(newdens[NAfilters]-propdens[NAfilters])
    EVnow <- with(cov.wt(rawpreds[NAfilters,],wttot,method="ML"),list(E=center,V=cov))
    EVnow <- opts$h(EVnow,pneww)
    if (adjust) {
      kappa <- opts$f(opts$xt,opts$g(ifelse(opts$single_betas,pneww$beta,origbeta))) - rawpreds[1,]
      EVnow$E <- EVnow$E + kappa
    }
    if (returnEV) EVnow else nllfun(opts$obs,EVnow,n=opts$n)
  }
}

    if (adjust) {
      kappa <- opts$f(opts$xt,opts$g(ifelse(opts$single_betas,pneww$beta,origbeta))) - rawpreds[1,]
      EVnow$E <- EVnow$E + kappa
    }
    if (returnEV) EVnow else nllfun(opts$obs,EVnow,n=opts$n)
  }
}
