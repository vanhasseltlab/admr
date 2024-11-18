#' Fitting aggregate data 2
#'
#'[fitEM2()] implements the Expectation-Maximization(EM) algorithm for parameter
#'estimation of the given aggregate data model, iterating over maximum
#'likelihood updates with weighted MC updates. This version of the function uses nloptr instead of optimx. TOL = 1e-10
#'
#' @param p0 initial parameter values
#' @param opts options
#' @param obs observed data
#' @param maxiter maximum number of iterations
#' @param convcrit_nll convergence criterion for the negative log-likelihood
#' @param nomap single model or multiple models
#' @returns A fitted model
#' @export
#' @examples
#' #test

fitEM2 <- function(p0,opts,obs,maxiter=100,convcrit_nll=0.001,nomap=TRUE) {
  # Initialization of model options and observed data (setup for E-step)
  if (nomap) {
    opts <- opts %>% p2opts(p0) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts),function(i) opts[[i]] %>% p2opts(p0) %>% obs2opts(obs[[i]]))
  }
  # Storage for tracking parameter values, NLL values, and time across iterations
  res <- tibble(p=vector("list",maxiter),
                nll=NA,
                appr_nll=NA,
                time=Sys.time(),
                iter=1)
  res$p[[1]] <- p0
  res$time[1] <- Sys.time()
  # Initial evaluation of negative log-likelihood (E-step)
  if (nomap) {
    res$nll[1] <- maxfunc(opts)(p0)
  } else {
    res$nll[1] <- Reduce('+',map(opts,~maxfunc(.)(p0)))
  }
  res$appr_nll[1] <- res$nll[1]
  message(paste0("iteration ",1,", nll=",res$nll[1]))
  pvals <- rep(0,length(p0)+1) # Initialize p-values for convergence checks

  # Begin the iterative process
  for (i in 2:maxiter) {
    # E-step: Construct the expectation based on the current parameter estimates
    if (nomap) {
      ff <- maxfunc(p2opts(opts,p0))
    } else {
      ffs <- map(opts,~maxfunc(p2opts(.,p0)))
      ff <- function(p) Reduce('+',map(ffs,~.(p)))
    }

    ff_nloptr <- function(params) {
      ff(params)  # assuming ff is your objective function
    }

    # M-step: Maximize the log-likelihood with respect to the parameters
    m0 <- nloptr::nloptr(
      x0 = p0,
      eval_f = ff_nloptr,
      lb = p0-1,
      ub = p0+1,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",  # equivalent to BOBYQA in nloptr
        xtol_rel = 1e-10,                # tolerance, adjust as needed
        maxeval = 100000,                    # max number of evaluations
        check_derivatives = F       # similar to setting KKT to FALSE
      )
    )

    # Update parameter estimates (M-step result)
    p0 <- m0$solution
    res$p[[i]] <- p0

    # E-step: Recalculate the negative log-likelihood with updated parameters
    res$time[i] <- Sys.time()
    if (nomap) {
      res$nll[i] <- maxfunc(p2opts(opts,p0))(p0)
    } else {
      res$nll[i] <- Reduce('+',map(opts,~maxfunc(p2opts(.,p0))(p0)))
    }
    res$appr_nll[i] <- m0$objective
    res$iter[i] <- i
    message(paste0("iteration ",i,", nll=",res$nll[i]))

    # Convergence check: stop if parameters are stationary (change in OFV or parameters)
    if (i>10) {
      pset <- do.call(rbind,res$p) %>% cbind(res$nll[!is.na(res$nll)]) %>% tail(10)
      pvals <- map(1:ncol(pset),function(colN) {
        xx <- 1:10
        yy <- pset[,colN]
        summary(lm(yy~xx))$coef[2,4]
      })
    }
    if (all(pvals>0.05)) {message("should break now due to stationary ofv+parameters");break()}
    if (abs(res$nll[i]-res$appr_nll[i]) <convcrit_nll) {message("should break now due to no difference between OFV and appr OFV");break()}
  }
  res[!is.na(res$nll),]
}
