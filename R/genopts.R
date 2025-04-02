#' Generate options for optimization
#'
#'[genopts()] initializes and generates core options and settings for modeling
#'and optimization, including random effects, simulation settings, and
#'likelihood approximations.
#'
#'
#' @param f The likelihood function
#' @param time The data
#' @param p The parameter settings
#' @param h The error function
#' @param nsim The number of simulations
#' @param n The number of simulations
#' @param adist The distribution of the random effects
#' @param interact FOCEI interaction
#' @param fo_appr Whether to use first-order approximation
#' @param biseq The sequence of random effects
#' @param omega_expansion The expansion factor for the omega matrix
#' @param single_betas The beta matrix
#' @returns A fitted model
#' @export
#' @examples
#' #test


genopts <- function(f,time,p,h,nsim=1,n=30,adist=NULL,
                    interact=TRUE,fo_appr=(nsim<10),biseq=NA,
                    omega_expansion=1,single_betas=NA,
                    p_thetai = function(p,origbeta,bi) {
                      dmnorm(bi,mean=log(p$beta/origbeta),
                             sigma=p$Omega,log=TRUE)$den},
                    g = function(beta,bi=rep(0,length(beta)),ai) {
                      beta*exp(bi)
                    }) {
  if (is.null(n)) stop("n is null, breaking early")
  if (missing(h)) {
    h <- function(EV,p) {
      if (!is.null(p$Sigma_add))
        diag(EV$V) <- diag(EV$V)+p$Sigma_add
      ## assuming the proper E is given as input depending on whether opts$interact==TRUE
      if (!is.null(p$Sigma_prop))
        diag(EV$V) <- diag(EV$V)+p$Sigma_prop*EV$E^2
      EV
    }
  }
  if (nsim==1 & !fo_appr) stop("Error! nsim=1 but FO_appr is false!")
  totseq <- sobol(nsim,nrow(p$Omega)+length(adist))
  if (!is.matrix(totseq)) totseq <- matrix(totseq,nsim,nrow(p$Omega)+length(adist))
  p_res <- p_to_optim(p)
  ptrans <- p_res$backtransformfunc
  pderiv <- p_res$d_psi_d_psitrans_long
  pt <- p_res$values
  p2 <- ptrans(p_res$values) ## to remove fixed parameters as strings
  ## attempt to generate derivative functions of g
  d_g_d_beta_expr <- tryCatch(D(body(g)[[2]],"beta"),error=function(e) NA)
  d_g_d_beta <- tryCatch(function(beta,bi,...) eval(d_g_d_beta_expr),error=function(e) NA)
  d_g_d_bi_expr <- tryCatch(D(body(g)[[2]],"bi"),error=function(e) NA)
  d_g_d_bi <- tryCatch(function(beta,bi,...) eval(d_g_d_bi_expr),error=function(e) NA)
  if (any(is.na(biseq))) {
    biseq <- apply(totseq[,1:nrow(p$Omega),drop=F],2,qnorm)
    biseq <- matrix(biseq,nsim,nrow(p$Omega))
  }
  ai <- NULL
  if (!is.null(adist)) {
    aiseq <- totseq[,nrow(p$Omega)+seq_along(adist),drop=F]
    if (!is.list(adist)) adist <- list(adist)
    ai <- do.call(cbind,map(seq_along(adist),~adist[[.]](aiseq[,.])))
  }
  list(f=f,g=g,time=time,p=p2,h=h,ptrans=ptrans,pderiv=pderiv,
       nsim=nsim,n=n,interact=interact,fo_appr=fo_appr,
       ai=ai,biseq=biseq,pt=pt,
       d_g_d_beta=d_g_d_beta,d_g_d_bi=d_g_d_bi,
       omega_expansion=omega_expansion,
       p_thetai=p_thetai,single_betas=single_betas,
       d_bi_d_omega=p_res$d_bi_d_omega,d_omega_d_Omega=p_res$d_omega_d_Omega)
}
