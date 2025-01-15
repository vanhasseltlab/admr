#' Generate a fitting function for optimization
#'
#'[genfitfunc()] generates a fitting function for optimization, taking observed data and model options as inputs.
#'
#'
#' @param opts The model options
#' @param obs The observed data
#' @returns A fitting function
#' @export
#' @examples
#' #test
#'

genfitfunc <- function(opts,obs) { # Generates a fitting function for optimization, taking observed data and model options as inputs.
  ### Generate sequence of covariates ai, and initial values of random effects biseq
  if (!(is.matrix(obs) && nrow(obs)==opts$nsim && ncol(obs)==opts$time) && !all(names(obs)== c("E","V")))
    stop("obs argument should give the data in aggregate form or as an nsim x Xi matrix!")
  if (missing(obs)) {## if observed data is missing, use expected data.
    message("genfitfunc message: Obs not supplied, using expected data as obs...")
    obs <- gen_pop_EV(opts)[1:2]
  }
  if (is.matrix(obs)) {
    message("genfitfunc message: Converting obs from raw data to aggregate E and V...")
    obs  <- meancov(obs)
  }
  ## make sure the obs$E is a normal vector
  obs$E <- c(obs$E)
  ## append obs to opts
  opts$obs <- obs
  function(pp,givedetails=FALSE,opts_overrides) {
    if (!missing(opts_overrides)) opts <- opts_overrides
    if (missing(pp)) pp <- opts$pt
    if (!is.list(pp)) {
      opts$pt <- pp
      opts$p <- opts$ptrans(pp)
    } else {
      opts$p <- pp
      opts$pt <- p_to_optim(pp)$values
    }
    EV <- gen_pop_EV(opts)
    invV <- tryCatch(solve(EV[[2]]),error=function(e) return(NA))
    if (any(is.na(invV))) stop("Problems inverting V, V=",paste(signif(EV[[2]],5),collapse="\n"),
                               "\n solve(V)=",paste(signif(invV,5),collapse="\n"))
    var_nll <- nllfun(opts$obs,EV,invV,opts$n)
    if (givedetails) {
      attr(var_nll,"EV") <- EV
      attr(var_nll,"obs") <- obs
      attr(var_nll,"nllfun") <- function(EV) nllfun(obs,EV,n=opts$n)
      attr(var_nll,"opts") <- opts
    }
    return(var_nll)
  }
}
