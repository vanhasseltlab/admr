#' Convert parameters to optimizable form
#'
#' @description
#' `p_to_optim` converts a parameter list into a form suitable for optimization by transforming
#' the parameters and computing their derivatives. This function handles the transformation of
#' both fixed and random effect parameters.
#'
#' @param p A list containing the parameter structure:
#'          \itemize{
#'            \item `beta`: Vector of population parameters
#'            \item `Omega`: Covariance matrix for random effects
#'            \item `Sigma_prop`: Proportional error variance (optional)
#'            \item `Sigma_add`: Additive error variance (optional)
#'          }
#'
#' @returns A list containing:
#' \itemize{
#'   \item `values`: Transformed parameter values
#'   \item `backtransformfunc`: Function to back-transform parameters
#'   \item `d_psi_d_psitrans_long`: Derivatives of the transformation
#'   \item `d_bi_d_omega`: Derivatives of random effects with respect to Omega
#'   \item `d_omega_d_Omega`: Derivatives of omega with respect to Omega
#' }
#'
#' @details
#' The function performs several transformations:
#' \itemize{
#'   \item Log-transforms positive parameters (beta and diagonal elements of Omega)
#'   \item Transforms correlation parameters using the inverse hyperbolic tangent
#'   \item Computes derivatives for the transformation
#'   \item Handles fixed parameters (specified as strings)
#' }
#'
#' The back-transformation function can be used to convert the optimized parameters back to
#' their original scale. The derivatives are used in the optimization process to compute
#' gradients and Hessians.
#'
#' @examples
#' # Define parameters
#' p <- list(
#'   beta = c(cl = 5, v1 = 10, v2 = 30, q = 10, ka = 1),
#'   Omega = matrix(c(0.09, 0, 0, 0, 0,
#'                   0, 0.09, 0, 0, 0,
#'                   0, 0, 0.09, 0, 0,
#'                   0, 0, 0, 0.09, 0,
#'                   0, 0, 0, 0, 0.09), nrow = 5, ncol = 5),
#'   Sigma_prop = 0.04
#' )
#'
#' # Convert to optimizable form
#' p_optim <- p_to_optim(p)
#'
#' # Back-transform parameters
#' p_orig <- p_optim$backtransformfunc(p_optim$values)
#'
#' @export
p_to_optim <- function(p) {
  ### a master function to feed parameters to optim
  if (any(names(p) == "Sigma_exp")) {
    message("Message from p_to_optim:\n Sigma_exp will be transformed to Sigma_prop!")
    if (any(names(p) == "Sigma_prop")) {
      p$Sigma_prop <- p$Sigma_prop + p$Sigma_exp
    } else {
      p$Sigma_prop <- p$Sigma_exp
    }
    p$Sigma_exp <- NULL
  }

  # Ensure names are preserved during transformations
  p2 <- map(p, ~if (is.list(.) || is.matrix(.)) {
    p_transform_mat(.)
  } else {
    transformed <- p_transform_vec(.)
    if (is.vector(.)) names(transformed[[1]]) <- names(.)  # Preserve names
    transformed
  })

  # Extract transformed values while keeping names
  p3 <- map(p2, ~{
    res <- .[[1]]
    if (!is.null(names(res))) names(res) <- names(.[[1]])  # Ensure names persist
    res
  })

  # Clean up names to remove 'beta.' prefix
  valuess <- do.call(c, p3)
  names(valuess) <- gsub("^beta\\.", "", names(valuess))  # Remove 'beta.' prefix

  backtransformfunc <- function(pp) {
    if (any(is.na(pp))) stop("Parameters contain NA values!!")
    if (missing(pp)) pp <- valuess
    p5 <- relist(pp, p3)
    p5 <- map(seq_along(p5), ~p2[[.]][[2]](p5[[.]]))
    names(p5) <- names(p)
    p5
  }

  d_psi_d_psitrans_long <- function(pp) {
    if (missing(pp)) pp <- valuess
    res <- do.call(c, map(names(p3), function(n) {
      r <- p2[[n]]$d_psi_d_psitrans_long(pp[grepl(n, names(valuess))])
      r
    }))
    names(res) <- rep(names(p3), map_dbl(p3, length))
    res
  }

  tottransfunclist <- do.call(c, map(p2, ~.$transfunclist))
  d_psi_d_psitrans_short <- function(pp) {  ## new
    if (missing(pp)) pp <- valuess
    pp <- split(unlist(pp), factor(str_extract(names(tottransfunclist), ".*[^[:digit:]]")))
    res <- map2(p2, pp, ~.x$d_psi_d_psitrans_short(.y)) %>% unlist()
    names(res) <- names(tottransfunclist)
    res
  }

  list(
    values = valuess,
    backtransformfunc = backtransformfunc,
    d_psi_d_psitrans_long = d_psi_d_psitrans_long,
    d_psi_d_psitrans_short = d_psi_d_psitrans_short,
    d_bi_d_omega = p2$Omega$d_bi_d_omega,
    d_omega_d_Omega = p2$Omega$d_omega_d_Omega
  )
}
