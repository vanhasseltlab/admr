# Convert parameters to optimizable form

`p_to_optim` converts a parameter list into a form suitable for
optimization by transforming the parameters and computing their
derivatives. This function handles the transformation of both fixed and
random effect parameters.

## Usage

``` r
p_to_optim(p)
```

## Arguments

- p:

  A list containing the parameter structure:

  - `beta`: Vector of population parameters (fixed effects) e.g.,
    clearance (CL), volume (V), absorption rate (ka)

  - `Omega`: Covariance matrix for random effects (between-subject
    variability) diagonal elements are variances, off-diagonal are
    covariances

  - `Sigma_prop`: Proportional error variance (optional) represents CV²
    of residual error

  - `Sigma_add`: Additive error variance (optional) represents constant
    error magnitude

## Value

A list containing:

- `values`: Vector of transformed parameter values on optimization scale

- `backtransformfunc`: Function to convert optimized values back to
  original scale

- `d_psi_d_psitrans_long`: Function for computing long-form parameter
  derivatives

- `d_psi_d_psitrans_short`: Function for computing short-form parameter
  derivatives

- `d_bi_d_omega`: Derivatives of random effects with respect to Omega
  elements

- `d_omega_d_Omega`: Derivatives of transformed Omega with respect to
  untransformed

## Details

Parameter Transformations:

- Population Parameters (beta): - Log transformation for positive
  parameters - Identity transformation for unrestricted parameters -
  Logit transformation for parameters bounded between 0 and 1

- Variance Components (Omega diagonal): - Log transformation to ensure
  positivity - Typically represents between-subject variability

- Correlation Components (Omega off-diagonal): - Inverse hyperbolic
  tangent (atanh) transformation - Ensures correlations remain between
  -1 and 1

- Residual Error (Sigma): - Log transformation for variance parameters -
  Handles both proportional and additive error structures

Derivative Computations:

- First-order derivatives for optimization algorithms

- Chain rule applied for composed transformations

- Separate handling of variance and correlation parameters

- Support for both dense and sparse matrices

Special Features:

- Handles fixed parameters (specified as character strings)

- Preserves parameter names throughout transformations

- Automatic conversion of exponential error to proportional

- Validates parameter values and transformations

## Examples

``` r
# Load required libraries
library(admr)
library(rxode2)

# Define a two-compartment model parameters
p <- list(
  # Population parameters (fixed effects)
  beta = c(cl = 5,    # Clearance (L/h)
          v1 = 10,    # Central volume (L)
          v2 = 30,    # Peripheral volume (L)
          q = 10,     # Inter-compartmental clearance (L/h)
          ka = 1),    # Absorption rate (1/h)

  # Between-subject variability (30% CV on all parameters)
  Omega = matrix(c(0.09, 0, 0, 0, 0,
                  0, 0.09, 0, 0, 0,
                  0, 0, 0.09, 0, 0,
                  0, 0, 0, 0.09, 0,
                  0, 0, 0, 0, 0.09), nrow = 5, ncol = 5),

  # Residual error (20% CV)
  Sigma_prop = 0.04
)

# Convert to optimization scale
p_optim <- p_to_optim(p)

# Back-transform to original scale
p_orig <- p_optim$backtransformfunc(p_optim$values)
```
