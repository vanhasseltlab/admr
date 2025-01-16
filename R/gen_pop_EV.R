#' @noRd

gen_pop_EV <- function(opts) { # Computes the expected population-level mean and covariance of the model predictions.
  if (opts$fo_appr) {
    bi <- gen_bi(opts)
    logweights <- samplogdensfun(bi,opts$p,opts$omega_expansion)
    weights <- logweights %>% logdens2wt()
    rawres <- map(1:nrow(bi),~foceapprEV_single(opts,bi[.,,drop=FALSE],opts$biseq[.,,drop=FALSE]
    ))
    E_pred <- Reduce('+',map2(rawres,weights,~.x$E*.y))
    V_pred <- Reduce('+',map2(rawres,weights,~.x$V*.y))
  } else {
    res <- MCapprEV(opts)
    E_pred <- res$E
    V_pred <- res$V
  }
  return(list(E=E_pred,V=V_pred))
}
