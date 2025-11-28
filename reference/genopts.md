# Generate options for aggregate data modeling

`genopts` initializes and generates core options and settings for
aggregate data modeling and optimization. It creates a comprehensive
options object that contains all necessary information for model
fitting, including random effects, simulation settings, and likelihood
approximations. This function is the main entry point for setting up
aggregate data modeling in the package.

## Usage

``` r
genopts(
  f,
  time,
  p,
  h,
  nsim = 1,
  n = 30,
  adist = NULL,
  interact = TRUE,
  fo_appr = (nsim < 10),
  biseq = NA,
  omega_expansion = 1,
  single_betas = NA,
  p_thetai = function(p, origbeta, bi) {
     dmnorm(bi, mean = log(p$beta/origbeta),
    sigma = p$Omega, log = TRUE)$den
 },
  g = function(beta, bi = rep(0, length(beta)), ai) {
beta * exp(bi)
 },
  no_cov = F
)
```

## Arguments

- f:

  The prediction function that simulates the model output given
  parameters and time points. This function should have the signature
  `function(time, theta_i, ...)` where:

  - `time`: Vector of time points

  - `theta_i`: Matrix of individual parameters

  - Returns: Matrix of predictions

- time:

  Vector of time points at which to evaluate the model predictions.

- p:

  List containing initial parameter values and structure:

  - `beta`: Vector of population parameters (fixed effects)

  - `Omega`: Covariance matrix for random effects (between-subject
    variability)

  - `Sigma_prop`: Proportional error variance (optional)

  - `Sigma_add`: Additive error variance (optional)

- h:

  The error function that computes the variance of the predictions. If
  not provided, a default function is used that adds proportional and
  additive error components.

- nsim:

  Number of Monte Carlo samples per iteration. Default is 1.

- n:

  Number of individuals in the dataset. Used for OFV, AIC, and BIC
  calculation. Default is 30.

- adist:

  Distribution of random effects. Default is NULL (normal distribution).

- interact:

  Logical indicating whether to use FOCEI interaction. Default is TRUE.

- fo_appr:

  Logical indicating whether to use first-order approximation. Default
  is TRUE if nsim \< 10, FALSE otherwise.

- biseq:

  Sequence of random effects. Default is NA (generated internally).

- omega_expansion:

  Factor by which to expand the covariance matrix during estimation.
  Default is 1.

- single_betas:

  Matrix of beta parameters for multiple models. Default is NA.

- p_thetai:

  Function to compute the log-density of random effects. Default is a
  multivariate normal density.

- g:

  Function to transform population parameters to individual parameters.
  Default is exponential transformation.

- no_cov:

  Logical indicating whether to ignore covariance in the error model.
  Default is FALSE.

## Value

A list containing:

- `f`: The prediction function

- `time`: Time points for evaluation

- `p`: Parameter structure and initial values

- `h`: The error function

- `nsim`: Number of Monte Carlo samples

- `n`: Number of individuals

- `adist`: Distribution of random effects

- `interact`: FOCEI interaction flag

- `fo_appr`: First-order approximation flag

- `biseq`: Random effects sequence

- `omega_expansion`: Covariance expansion factor

- `single_betas`: Beta parameters for multiple models

- `p_thetai`: Random effects density function

- `g`: Parameter transformation function

- `pt`: Transformed initial parameters

- `ptrans`: Function to back-transform parameters

- `pderiv`: Function to compute parameter derivatives

- `d_g_d_beta`: Derivative of g with respect to beta

- `d_g_d_bi`: Derivative of g with respect to random effects

- `d_bi_d_omega`: Derivative of random effects with respect to Omega

- `d_omega_d_Omega`: Derivative of transformed Omega with respect to
  untransformed

- `no_cov`: Logical indicating whether to ignore covariance

## Details

The function performs several key operations:

- Parameter Transformation: - Converts parameters to optimization
  scale - Computes derivatives for optimization - Handles fixed
  parameters

- Random Effects Generation: - Uses Sobol sequences for quasi-random
  sampling - Applies normal quantile transformation - Supports custom
  distributions

- Error Function Setup: - Handles proportional error: y = f(t,θ)(1 +
  ε) - Handles additive error: y = f(t,θ) + ε - Combines both error
  types

Key features:

- Automatic derivative computation for parameter transformations

- Support for multiple models and parameter structures

- Flexible error model specification

- Efficient random effects generation

## Examples

``` r
# Load required libraries
library(admr)
library(rxode2)
library(nlmixr2)
library(dplyr)
library(tidyr)
library(mnorm)


# Define prediction function for a two-compartment model
predder <- function(time, theta_i, dose = 100) {
  n_individuals <- nrow(theta_i)
  if (is.null(n_individuals)) n_individuals <- 1

  # Create event table for dosing and sampling
  ev <- eventTable(amount.units="mg", time.units="hours")
  ev$add.dosing(dose = dose, nbr.doses = 1, start.time = 0)
  ev$add.sampling(time)

  # Solve ODE system
  out <- rxSolve(rxModel, params = theta_i, events = ev, cores = 0)

  # Return matrix of predictions
  matrix(out$cp, nrow = n_individuals, ncol = length(time), byrow = TRUE)
}

# Create options for a two-compartment model
opts <- genopts(
  f = predder,
  time = c(.1, .25, .5, 1, 2, 3, 5, 8, 12),
  p = list(
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
  ),
  nsim = 2500,  # Number of Monte Carlo samples
  n = 500,      # Number of individuals
  fo_appr = FALSE,  # Use Monte Carlo approximation
  omega_expansion = 1.2  # Expand covariance during estimation
)
```
