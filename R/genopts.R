#' Generate options for aggregate data modeling
#'
#' @description
#' `genopts` initializes and generates core options and settings for aggregate data modeling
#' and optimization. It creates a comprehensive options object that contains all necessary
#' information for model fitting, including random effects, simulation settings, and
#' likelihood approximations.
#'
#' @param f The prediction function that simulates the model output given parameters and time points.
#' @param time Vector of time points at which to evaluate the model predictions.
#' @param p List containing initial parameter values and structure:
#'          \itemize{
#'            \item `beta`: Vector of population parameters
#'            \item `Omega`: Covariance matrix for random effects
#'            \item `Sigma_prop`: Proportional error variance (optional)
#'            \item `Sigma_add`: Additive error variance (optional)
#'          }
#' @param h The error function that computes the variance of the predictions.
#' @param nsim Number of Monte Carlo samples per iteration. Default is 1.
#' @param n Number of individuals in the dataset. Used for OFV, AIC, and BIC calculation.
#'          Default is 30.
#' @param adist Distribution of random effects. Default is NULL (normal distribution).
#' @param interact Logical indicating whether to use FOCEI interaction. Default is TRUE.
#' @param fo_appr Logical indicating whether to use first-order approximation.
#'                Default is TRUE if nsim < 10, FALSE otherwise.
#' @param biseq Sequence of random effects. Default is NA (generated internally).
#' @param omega_expansion Factor by which to expand the covariance matrix during estimation.
#'                       Default is 1.
#' @param single_betas Matrix of beta parameters for multiple models. Default is NA.
#' @param p_thetai Function to compute the log-density of random effects. Default is a
#'                 multivariate normal density.
#' @param g Function to transform population parameters to individual parameters. Default is
#'          exponential transformation.
#'
#' @returns A list containing:
#' \itemize{
#'   \item `f`: The prediction function
#'   \item `time`: Time points for evaluation
#'   \item `p`: Parameter structure and initial values
#'   \item `h`: The error function
#'   \item `nsim`: Number of Monte Carlo samples
#'   \item `n`: Number of individuals
#'   \item `adist`: Distribution of random effects
#'   \item `interact`: FOCEI interaction flag
#'   \item `fo_appr`: First-order approximation flag
#'   \item `biseq`: Random effects sequence
#'   \item `omega_expansion`: Covariance expansion factor
#'   \item `single_betas`: Beta parameters for multiple models
#'   \item `p_thetai`: Random effects density function
#'   \item `g`: Parameter transformation function
#'   \item `pt`: Transformed initial parameters
#' }
#'
#' @details
#' This function is the main entry point for setting up aggregate data modeling. It creates
#' an options object that contains all necessary information for model fitting, including
#' the prediction function, parameter structure, and various settings for the optimization
#' algorithm.
#'
#' @examples
#' # Define prediction function
#' predder <- function(time, theta_i, dose = 100) {
#'   # ... prediction function implementation ...
#' }
#'
#' # Create options
#' opts <- genopts(
#'   f = predder,
#'   time = c(.1, .25, .5, 1, 2, 3, 5, 8, 12),
#'   p = list(
#'     beta = c(cl = 5, v1 = 10, v2 = 30, q = 10, ka = 1),
#'     Omega = matrix(c(0.09, 0, 0, 0, 0,
#'                     0, 0.09, 0, 0, 0,
#'                     0, 0, 0.09, 0, 0,
#'                     0, 0, 0, 0.09, 0,
#'                     0, 0, 0, 0, 0.09), nrow = 5, ncol = 5),
#'     Sigma_prop = 0.04
#'   ),
#'   nsim = 2500,
#'   n = 500
#' )
#'
#' @export
genopts <- function(f, time, p, h, nsim = 1, n = 30, adist = NULL,
                    interact = TRUE, fo_appr = (nsim < 10), biseq = NA,
                    omega_expansion = 1, single_betas = NA,
                    p_thetai = function(p, origbeta, bi) {
                      dmnorm(bi, mean = log(p$beta/origbeta),
                             sigma = p$Omega, log = TRUE)$den
                    },
                    g = function(beta, bi = rep(0, length(beta)), ai) {
                      beta * exp(bi)
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
