#' Fitting aggregate data
#'
#'[timedEM()] implements the Expectation-Maximization(EM) algorithm for parameter
#'estimation of the given aggregate data model, iterating over maximum
#'likelihood updates with weighted MC updates. It is used to compare the performance
#'of the old and new implementation of aggregate data modelling.
#'
#'
#' @param init initial parameter values
#' @param opts options
#' @param obs observed data
#' @param maxiter maximum number of iterations
#' @param convcrit_nll convergence criterion for the negative log-likelihood
#' @param nomap single model or multiple models
#' @returns A fitted model
#' @export
#' @examples
#' #test

timedEM <- function(init,opts,obs,maxiter=100,convcrit_nll=0.0005,nomap=TRUE) { # Implements the Expectation-Maximization (EM) algorithm for parameter estimation, iterating over maximum likelihood updates.
  if (nomap) {
    opts <- opts %>% p2opts(init) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts),function(i) opts[[i]] %>% p2opts(init) %>% obs2opts(obs[[i]]))
  }
  res <- tibble(p=vector("list",maxiter),
                nll=NA,
                appr_nll=NA,
                time=Sys.time(),
                iter=1)
  res$p[[1]] <- init
  res$time[1] <- Sys.time()
  if (nomap) {
    res$nll[1] <- maxfunc(opts)(init)
  } else {
    res$nll[1] <- Reduce('+',map(opts,~maxfunc(.)(init)))
  }
  res$appr_nll[1] <- res$nll[1]
  message(paste0("iteration ",1,", nll=",res$nll[1]))
  pvals <- rep(0,length(init)+1)
  for (i in 2:maxiter) {
    if (nomap) {
      ff <- maxfunc(p2opts(opts,init))
    } else {
      ffs <- map(opts,~maxfunc(p2opts(.,init)))
      ff <- function(p) Reduce('+',map(ffs,~.(p)))
    }

    ff_nloptr <- function(params) {
      ff(params)
    }

    m0 <- nloptr::nloptr(
      x0 = init,
      eval_f = ff_nloptr,
      lb = init - 2,
      ub = init + 2,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",
        ftol_rel=sqrt(.Machine$double.eps),
        maxeval = 2000
      )
    )

    init <- m0$solution
    res$p[[i]] <- init
    res$time[i] <- Sys.time()
    if (nomap) {
      res$nll[i] <- maxfunc(p2opts(opts,init))(init)
    } else {
      res$nll[i] <- Reduce('+',map(opts,~maxfunc(p2opts(.,init))(init)))
    }
    res$appr_nll[i] <- m0$objective
    res$iter[i] <- i
    message(paste0("iteration ",i,", nll=",res$nll[i]))
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

