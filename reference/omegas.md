# Create a covariance matrix with specified diagonal and off-diagonal values

`omegas` creates a covariance matrix with specified diagonal and
off-diagonal values. This function is useful for creating initial or
fixed covariance matrices for random effects in pharmacometric models.

## Usage

``` r
omegas(diag, offdiag, n_om)
```

## Arguments

- diag:

  Value for the diagonal elements of the matrix. This represents the
  variance of each random effect. For log-normal distributions, this is
  typically the squared coefficient of variation (CV²) on the log scale.

- offdiag:

  Value for the off-diagonal elements of the matrix. This represents the
  covariance between random effects. A value of 0 indicates independence
  between random effects.

- n_om:

  Size of the matrix (number of random effects). This should match the
  number of random effects in your model (e.g., 2 for CL and V in a
  one-compartment model).

## Value

A symmetric matrix of size `n_om` x `n_om` with:

- Diagonal elements equal to `diag` (variances)

- Off-diagonal elements equal to `offdiag` (covariances)

## Details

The function creates a symmetric covariance matrix for random effects
where:

- Diagonal elements (ω²) represent the between-subject variability

- Off-diagonal elements (ω_ij) represent correlations between parameters

- The resulting matrix must be positive definite for valid computations

Common use cases include:

- Initial estimates for model fitting: - Setting diagonal elements to
  expected variability (e.g., 0.09 for 30% CV) - Starting with zero
  correlations (offdiag = 0)

- Simulation studies: - Specifying known parameter variability - Testing
  impact of parameter correlations

- Sensitivity analysis: - Evaluating model behavior under different
  variability assumptions - Assessing impact of parameter correlations

Mathematical details:

- For log-normal distributions: CV ≈ sqrt(exp(ω²) - 1)

- Correlation ρ_ij = ω_ij / sqrt(ω_ii \* ω_jj)

- Matrix must be positive definite: all eigenvalues \> 0

## Examples

``` r
# Load required libraries
library(admr)
library(rxode2)

# Create a diagonal matrix for a one-compartment model
# 30% CV on CL and V (ω² = 0.09 for each)
omega1 <- omegas(0.09, 0, 2)

# Create a matrix with correlations for a two-compartment model
# 30% CV on all parameters (CL, V1, Q, V2)
# Correlation of 0.3 between parameters
omega2 <- omegas(0.09, 0.03, 4)

# Create a matrix for testing parameter correlations
# High variability (50% CV, ω² = 0.25) and strong correlations (0.1)
omega3 <- omegas(0.25, 0.1, 3)
```
