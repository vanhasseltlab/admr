#' Create a covariance matrix with specified diagonal and off-diagonal values
#'
#' @description
#' `omegas` creates a covariance matrix with specified diagonal and off-diagonal values.
#' This function is useful for creating initial or fixed covariance matrices for
#' random effects in pharmacometric models.
#'
#' @param diag Value for the diagonal elements of the matrix. This represents the
#'             variance of each random effect. For log-normal distributions, this is
#'             typically the squared coefficient of variation (CV²) on the log scale.
#' @param offdiag Value for the off-diagonal elements of the matrix. This represents
#'                the covariance between random effects. A value of 0 indicates
#'                independence between random effects.
#' @param n_om Size of the matrix (number of random effects). This should match the
#'            number of random effects in your model (e.g., 2 for CL and V in a
#'            one-compartment model).
#'
#' @returns A symmetric matrix of size `n_om` x `n_om` with:
#' \itemize{
#'   \item Diagonal elements equal to `diag` (variances)
#'   \item Off-diagonal elements equal to `offdiag` (covariances)
#' }
#'
#' @details
#' The function creates a symmetric covariance matrix for random effects where:
#' \itemize{
#'   \item Diagonal elements (ω²) represent the between-subject variability
#'   \item Off-diagonal elements (ω_ij) represent correlations between parameters
#'   \item The resulting matrix must be positive definite for valid computations
#' }
#'
#' Common use cases include:
#' \itemize{
#'   \item Initial estimates for model fitting:
#'     - Setting diagonal elements to expected variability (e.g., 0.09 for 30% CV)
#'     - Starting with zero correlations (offdiag = 0)
#'   \item Simulation studies:
#'     - Specifying known parameter variability
#'     - Testing impact of parameter correlations
#'   \item Sensitivity analysis:
#'     - Evaluating model behavior under different variability assumptions
#'     - Assessing impact of parameter correlations
#' }
#'
#' Mathematical details:
#' \itemize{
#'   \item For log-normal distributions: CV ≈ sqrt(exp(ω²) - 1)
#'   \item Correlation ρ_ij = ω_ij / sqrt(ω_ii * ω_jj)
#'   \item Matrix must be positive definite: all eigenvalues > 0
#' }
#'
#' @examples
#' # Create a diagonal matrix for a one-compartment model
#' # 30% CV on CL and V (ω² = 0.09 for each)
#' omega1 <- omegas(0.09, 0, 2)
#'
#' # Create a matrix with correlations for a two-compartment model
#' # 30% CV on all parameters (CL, V1, Q, V2)
#' # Correlation of 0.3 between parameters
#' omega2 <- omegas(0.09, 0.03, 4)
#'
#' # Create a matrix for testing parameter correlations
#' # High variability (50% CV, ω² = 0.25) and strong correlations (0.1)
#' omega3 <- omegas(0.25, 0.1, 3)
#'
#' @export
omegas <- function(diag, offdiag, n_om) {
  # Input validation (could be added)
  # if (diag < 0) stop("Diagonal elements must be non-negative")
  # if (abs(offdiag) >= abs(diag)) warning("Off-diagonal elements larger than diagonal may result in non-positive definite matrix")
  
  # Create matrix filled with off-diagonal elements
  # This ensures symmetry in the covariance matrix
  m <- matrix(offdiag, n_om, n_om)
  
  # Set diagonal elements
  # diag() extracts/replaces the diagonal of a matrix
  diag(m) <- diag
  
  # Return the covariance matrix
  # The matrix will be:
  # [diag   offdiag offdiag ...]
  # [offdiag diag   offdiag ...]
  # [offdiag offdiag diag   ...]
  # [   ...    ...    ...   ...]
  return(m)
}
