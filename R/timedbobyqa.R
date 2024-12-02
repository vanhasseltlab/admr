#' Fitting aggregate data
#'
#'[timedbobyqa()] implements the bobqa algorithm for parameter estimation of the
#'given aggregate data model, iterating over maximum likelihood updates with MC
#'updates. Each iteration creates new MC samples and updates the parameter values.
#'It is used to compare the performance of the old and new implementation
#'of aggregate data modelling.
#'
#' @param p0 initial parameter values
#' @param opts options
#' @param obs observed data
#' @returns A fitted model
#' @export
#' @examples
#' #test

timedbobyqa <- function(p0, opts, obs) {
  res <- tibble(p=vector("list",5e3),nll=NA,time=Sys.time(),iter=NA)
  i <- 1
  fitfun <-   fitfun <- genfitfunc(opts,obs)
  fitfun2 <- function(p) {
    nllNow <- fitfun(p)
    res$p[[i]] <<- p
    res$time[i] <<- Sys.time()
    res$nll[i] <<- nllNow
    res$iter[i] <<- i
    i <<- i+1
    nllNow
  }
  est <- nloptr::nloptr(p0,fitfun2,lb=p0-3,ub=p0+3,opts=list("algorithm"="NLOPT_LN_BOBYQA", "ftol_rel"=sqrt(.Machine$double.eps),
                                                             maxeval = 2000,                    # max number of evaluations
                                                             check_derivatives = F))
  est$solution
  res[!is.na(res$nll),]
}
