#' Fitting aggregate data
#'
#'[timedbobyqa()] implements the bobyqa algorithm for parameter estimation of the
#'given aggregate data model, iterating over maximum likelihood updates with MC
#'updates. Each iteration creates new MC samples and updates the parameter values.
#'It is used to compare the performance of the old and new implementation
#'of aggregate data modelling.
#'
#' @param init initial parameter values
#' @param opts options
#' @param obs observed data
#' @returns A fitted model
#' @export
#' @examples
#' #test

timedbobyqa <- function(init, opts, obs, nomap = TRUE) {

  if (nomap) {
    opts <- opts %>% p2opts(init) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts),function(i) opts[[i]] %>% p2opts(init) %>% obs2opts(obs[[i]]))
  }

  res <- tibble(p=vector("list",1e4),nll=NA,time=Sys.time(),iter=NA)

  i <- 1

  if (nomap) {
    fitfun <- genfitfunc(opts, obs)
  } else {
    fitfuns <- map(seq_along(opts), ~ genfitfunc(opts[[.x]], obs[[.x]]))
    fitfun <- function(p) Reduce('+', map(fitfuns, ~ .(p)))
  }

  fitfun2 <- function(p) {
    nllNow <- fitfun(p)
    res$p[[i]] <<- p
    res$time[i] <<- Sys.time()
    res$nll[i] <<- nllNow
    res$iter[i] <<- i

    if (i %% 50 == 0) {
      cat("Iteration:", i, "- NLL:", nllNow, "\n")
    }
    i <<- i + 1
    nllNow
  }

  est <- nloptr::nloptr(init,
                        fitfun2,
                        lb=init-2,
                        ub=init+2,
                        opts=list(
                          algorithm="NLOPT_LN_BOBYQA",
                          ftol_rel=.Machine$double.eps^2,
                          maxeval = 5000,                    # max number of evaluations
                          check_derivatives = F))
  est$solution
  res[!is.na(res$nll),]
}
