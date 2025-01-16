#' Generate a dataset
#'
#'[gendataset()] generates a dataset based on the model structure and random effects, with optional residual error.
#'
#'
#' @param opts A list of options
#' @param seed A seed for the random number generator
#' @param reserr A logical indicating whether residual error should be added
#' @param nlmixrform A logical indicating whether the dataset should be in nlmixr format
#' @returns A dataset
#' @export
#' @examples
#' #test


gendataset <- function(opts,seed=1,reserr=TRUE,nlmixrform=FALSE) {
  set.seed(seed)
  bi <- gen_bi(opts,FALSE)
  nsimobs <- length(opts$xt)*opts$nsim
  theta_i <- g_iter(opts,bi)
  res <- opts$f(opts$xt,theta_i)
  if (!reserr) return(res)
  if (!is.null(opts$p$Sigma_prop))
    res <- res*(1+rnorm(nsimobs,sd=sqrt(opts$p$Sigma_prop)))
  if (!is.null(opts$p$Sigma_add))
    res <- res+rnorm(nsimobs,sd=sqrt(opts$p$Sigma_add))
  if (!nlmixrform) {
    return(res)
  } else {
    tibble(dv=c(t(cbind(NA,res))),time=rep(c(0,opts$xt),times=opts$nsim),
           id=rep(seq_len(opts$nsim),each=1+length(opts$xt)),
           amt=rep(c(100,rep(NA,length(opts$xt))),opts$nsim),
           evid=101*as.integer(!is.na(.data$amt)),
           cmt=ifelse(is.na(.data$amt),2,1))
  }
}
