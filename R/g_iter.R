#' @noRd

g_iter <- function(opts,bi) { # Applies the model g to generate individual-specific parameter values theta_i based on random effects bi.
  if (!is.matrix(bi)) {
    opts$g(opts$p$beta,bi,opts$ai)
  } else {
    t(apply(bi,1,function(i) opts$g(opts$p$beta,i,opts$ai)))
  }
}
