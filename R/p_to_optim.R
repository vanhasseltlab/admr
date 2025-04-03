#' Convert parameters to optimizable form
#'
#' @description
#' `p_to_optim` converts a parameter list into a form suitable for optimization by transforming
#' the parameters and computing their derivatives. This function handles the transformation of
#' both fixed and random effect parameters.
#'
#' @param p A list containing the parameter structure:
#'          \itemize{
#'            \item `beta`: Vector of population parameters (fixed effects)
#'                         e.g., clearance (CL), volume (V), absorption rate (ka)
#'            \item `Omega`: Covariance matrix for random effects (between-subject variability)
#'                          diagonal elements are variances, off-diagonal are covariances
#'            \item `Sigma_prop`: Proportional error variance (optional)
#'                               represents CV² of residual error
#'            \item `Sigma_add`: Additive error variance (optional)
#'                              represents constant error magnitude
#'          }
#'
#' @returns A list containing:
#' \itemize{
#'   \item `values`: Vector of transformed parameter values on optimization scale
#'   \item `backtransformfunc`: Function to convert optimized values back to original scale
#'   \item `d_psi_d_psitrans_long`: Function for computing long-form parameter derivatives
#'   \item `d_psi_d_psitrans_short`: Function for computing short-form parameter derivatives
#'   \item `d_bi_d_omega`: Derivatives of random effects with respect to Omega elements
#'   \item `d_omega_d_Omega`: Derivatives of transformed Omega with respect to untransformed
#' }
#'
#' @details
#' Parameter Transformations:
#' \itemize{
#'   \item Population Parameters (beta):
#'     - Log transformation for positive parameters
#'     - Identity transformation for unrestricted parameters
#'     - Logit transformation for parameters bounded between 0 and 1
#'
#'   \item Variance Components (Omega diagonal):
#'     - Log transformation to ensure positivity
#'     - Typically represents between-subject variability
#'
#'   \item Correlation Components (Omega off-diagonal):
#'     - Inverse hyperbolic tangent (atanh) transformation
#'     - Ensures correlations remain between -1 and 1
#'
#'   \item Residual Error (Sigma):
#'     - Log transformation for variance parameters
#'     - Handles both proportional and additive error structures
#' }
#'
#' Derivative Computations:
#' \itemize{
#'   \item First-order derivatives for optimization algorithms
#'   \item Chain rule applied for composed transformations
#'   \item Separate handling of variance and correlation parameters
#'   \item Support for both dense and sparse matrices
#' }
#'
#' Special Features:
#' \itemize{
#'   \item Handles fixed parameters (specified as character strings)
#'   \item Preserves parameter names throughout transformations
#'   \item Automatic conversion of exponential error to proportional
#'   \item Validates parameter values and transformations
#' }
#'
#' @examples
#' # Define a two-compartment model parameters
#' p <- list(
#'   # Population parameters (fixed effects)
#'   beta = c(cl = 5,    # Clearance (L/h)
#'           v1 = 10,    # Central volume (L)
#'           v2 = 30,    # Peripheral volume (L)
#'           q = 10,     # Inter-compartmental clearance (L/h)
#'           ka = 1),    # Absorption rate (1/h)
#'   
#'   # Between-subject variability (30% CV on all parameters)
#'   Omega = matrix(c(0.09, 0, 0, 0, 0,
#'                   0, 0.09, 0, 0, 0,
#'                   0, 0, 0.09, 0, 0,
#'                   0, 0, 0, 0.09, 0,
#'                   0, 0, 0, 0, 0.09), nrow = 5, ncol = 5),
#'   
#'   # Residual error (20% CV)
#'   Sigma_prop = 0.04
#' )
#'
#' # Convert to optimization scale
#' p_optim <- p_to_optim(p)
#'
#' # Back-transform to original scale
#' p_orig <- p_optim$backtransformfunc(p_optim$values)
#'
#' @export
p_to_optim <- function(p) {
  # Handle legacy exponential error specification
  if (any(names(p) == "Sigma_exp")) {
    message("Message from p_to_optim:\n Sigma_exp will be transformed to Sigma_prop!")
    if (any(names(p) == "Sigma_prop")) {
      # Combine exponential and proportional error if both present
      p$Sigma_prop <- p$Sigma_prop + p$Sigma_exp
    } else {
      # Convert exponential to proportional error
      p$Sigma_prop <- p$Sigma_exp
    }
    p$Sigma_exp <- NULL
  }

  # Transform parameters while preserving structure
  # For matrices (Omega) and vectors (beta), apply appropriate transformations
  p2 <- map(p, ~if (is.list(.) || is.matrix(.)) {
    # Transform matrices (e.g., Omega) using matrix-specific transformations
    p_transform_mat(.)
  } else {
    # Transform vectors (e.g., beta) using vector-specific transformations
    transformed <- p_transform_vec(.)
    if (is.vector(.)) names(transformed[[1]]) <- names(.)  # Preserve parameter names
    transformed
  })

  # Extract transformed values while maintaining parameter names
  p3 <- map(p2, ~{
    res <- .[[1]]
    if (!is.null(names(res))) names(res) <- names(.[[1]])
    res
  })

  # Flatten parameter vector and clean up names
  valuess <- do.call(c, p3)
  names(valuess) <- gsub("^beta\\.", "", names(valuess))

  # Create back-transformation function
  backtransformfunc <- function(pp) {
    # Validate input parameters
    if (any(is.na(pp))) stop("Parameters contain NA values!!")
    if (missing(pp)) pp <- valuess
    
    # Restructure and back-transform parameters
    p5 <- relist(pp, p3)  # Restore parameter structure
    p5 <- map(seq_along(p5), ~p2[[.]][[2]](p5[[.]]))  # Apply back-transformations
    names(p5) <- names(p)  # Restore parameter names
    p5
  }

  # Create function for computing long-form derivatives
  d_psi_d_psitrans_long <- function(pp) {
    if (missing(pp)) pp <- valuess
    # Compute derivatives for each parameter type
    res <- do.call(c, map(names(p3), function(n) {
      r <- p2[[n]]$d_psi_d_psitrans_long(pp[grepl(n, names(valuess))])
      r
    }))
    # Label derivatives with parameter types
    names(res) <- rep(names(p3), map_dbl(p3, length))
    res
  }

  # Combine transformation functions
  tottransfunclist <- do.call(c, map(p2, ~.$transfunclist))
  
  # Create function for computing short-form derivatives
  d_psi_d_psitrans_short <- function(pp) {
    if (missing(pp)) pp <- valuess
    # Split parameters by type and compute derivatives
    pp <- split(unlist(pp), 
                factor(str_extract(names(tottransfunclist), ".*[^[:digit:]]")))
    res <- map2(p2, pp, ~.x$d_psi_d_psitrans_short(.y)) %>% unlist()
    names(res) <- names(tottransfunclist)
    res
  }

  # Return list of transformed parameters and associated functions
  list(
    values = valuess,                    # Transformed parameter values
    backtransformfunc = backtransformfunc,  # Function to restore original scale
    d_psi_d_psitrans_long = d_psi_d_psitrans_long,   # Long-form derivatives
    d_psi_d_psitrans_short = d_psi_d_psitrans_short, # Short-form derivatives
    d_bi_d_omega = p2$Omega$d_bi_d_omega,           # Random effects derivatives
    d_omega_d_Omega = p2$Omega$d_omega_d_Omega      # Omega derivatives
  )
}
