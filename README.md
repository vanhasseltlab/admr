admr: Aggregate Data Modeling in R
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# admr: Aggregate Data Modeling in R <img src="man/figures/logo.png" width="120" align="right" />

<!-- badges: start -->
<!-- badges: end -->

admr (Aggregate Data Modeling in R) is an open-source R package designed
to facilitate pharmacometric modeling using summary-level data. It
enables users to work with aggregate data, such as mean observations and
variance-covariance matrices, to fit pharmacokinetic and pharmacodynamic
(PK/PD) models efficiently. This package implements a newly developed
Expectation-Maximization (EM) algorithm to enhance computational
performance and provides tools for advanced modeling applications.

## Overview

The `admr` package provides a comprehensive framework for aggregate data
modeling in pharmacometrics, offering several key advantages:

- **Efficient Parameter Estimation**: Uses the Iterative Reweighting
  Monte Carlo (IRMC) algorithm for robust and fast parameter estimation.
- **Flexible Data Integration**: Works with both individual-level and
  aggregate data, making it ideal for meta-analyses and literature-based
  modeling.
- **Advanced Modeling Features**: Supports complex PK/PD models with
  various error structures and parameter transformations.
- **Comprehensive Diagnostics**: Built-in tools for model assessment,
  convergence checking, and parameter stability analysis.

## Key Features

- **Efficient Model Fitting**: Uses iterative reweighted Monte Carlo
  (IRMC) for robust parameter estimation, improving speed and
  scalability compared to traditional Monte Carlo methods.

- **Flexible Data Formats**: Supports both raw and aggregate data
  formats, allowing for the integration of summary-level data from
  diverse sources, including published literature and simulated models.

- **Comprehensive Diagnostics**: Built-in tools for model assessment,
  convergence checking, and parameter stability analysis.

- **Meta-Analysis Support**: Facilitates model-based meta-analyses by
  enabling the combination of summary data across studies.

- **R Integration**: Fully compatible with R, leveraging popular
  pharmacometric modeling libraries like rxode2.

- **Open-Source**: Developed for accessibility and ease of use by the
  pharmacometric community.

## Installation

This is an *R* package. [*R*](https://www.r-project.org/) is required,
[*RStudio*](https://posit.co/downloads/) is recommended.

You can install the development version of admr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hiddevandebeek/admr")
```

## Quick Start

Below is a complete example of how to use admr to fit a pharmacokinetic
model to aggregate data:

``` r
# Load required libraries
library(admr)
library(rxode2)
#> Warning: package 'rxode2' was built under R version 4.4.2
#> rxode2 3.0.4 using 7 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(nlmixr2)
#> Loading required package: nlmixr2data
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(mnorm)

# Load and prepare data
data(examplomycin)
examplomycin_wide <- examplomycin %>%
  filter(EVID != 101) %>%
  dplyr::select(ID, TIME, DV) %>%
  pivot_wider(names_from = TIME, values_from = DV) %>%
  dplyr::select(-c(1))

# Create aggregated data
examplomycin_aggregated <- examplomycin_wide %>%
  admr::meancov()

# Define RxODE model
rxModel <- RxODE({
  cp = linCmt(
    cl,           # Clearance
    v1,           # Volume of the central compartment
    v2,           # Volume of the peripheral compartment
    q,            # Inter-compartmental clearance
    ka            # Absorption rate constant
  )
})

# Define prediction function
predder <- function(time, theta_i, dose = 100) {
  n_individuals <- nrow(theta_i)
  if (is.null(n_individuals)) n_individuals <- 1
  
  ev <- eventTable(amount.units="mg", time.units="hours")
  ev$add.dosing(dose = dose, nbr.doses = 1, start.time = 0)
  ev$add.sampling(time)
  
  out <- rxSolve(rxModel, params = theta_i, events = ev, cores = 0)
  cp_matrix <- matrix(out$cp, nrow = n_individuals, ncol = length(time),
                      byrow = TRUE)
  return(cp_matrix)
}

# Create options
opts <- genopts(
  time = c(.1, .25, .5, 1, 2, 3, 5, 8, 12),
  p = list(
    beta = c(cl = 5, v1 = 10, v2 = 30, q = 10, ka = 1),
    Omega = matrix(c(0.09, 0, 0, 0, 0,
                    0, 0.09, 0, 0, 0,
                    0, 0, 0.09, 0, 0,
                    0, 0, 0, 0.09, 0,
                    0, 0, 0, 0, 0.09), nrow = 5, ncol = 5),
    Sigma_prop = 0.04
  ),
  nsim = 2500,
  n = 500,
  fo_appr = FALSE,
  omega_expansion = 1.2,
  f = predder
)

# Fit model to data
result <- fitIRMC(opts = opts, obs = examplomycin_aggregated)
#> Chain 1:
#> Iter | NLL and Parameters (11 values)
#> --------------------------------------------------------------------------------
#>    1: -1839.577    1.609    2.303    3.401    2.303    0.000   -2.408   -2.408   -2.408   -2.408   -2.408   -3.219
#> 
#> ### Wide Search Phase ###
#>    2: -1845.238    1.601    2.309    3.404    2.284    0.020   -2.279   -2.167   -2.333   -2.245   -2.439   -3.235
#>    3: -1845.277    1.600    2.307    3.406    2.285    0.017   -2.287   -2.195   -2.342   -2.267   -2.416   -3.235
#>    4: -1845.277    1.600    2.306    3.406    2.285    0.017   -2.287   -2.195   -2.342   -2.266   -2.416   -3.235
#>    5: -1845.277    1.600    2.306    3.406    2.285    0.017   -2.287   -2.195   -2.342   -2.266   -2.416   -3.235
#> Phase Wide Search Phase converged at iteration 5.
#> 
#> ### Focussed Search Phase ###
#>    6: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#>    7: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#>    8: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#>    9: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#>   10: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#>   11: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#> Phase Focussed Search Phase converged at iteration 11.
#> 
#> ### Fine-Tuning Phase ###
#>   12: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.200   -2.343   -2.264   -2.411   -3.235
#> Phase Fine-Tuning Phase converged at iteration 12.
#> 
#> ### Precision Phase ###
#>   13: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.201   -2.343   -2.265   -2.410   -3.235
#>   14: -1845.280    1.601    2.307    3.405    2.285    0.018   -2.286   -2.201   -2.343   -2.265   -2.410   -3.235
#> Phase Precision Phase converged at iteration 14.
#> 
#> Chain 1 Complete: Final NLL = -1845.280, Time Elapsed = 6.28 seconds
#> 
print(result)
#> -- FitIRMC Summary --
#> 
#> -- Objective Function and Information Criteria --
#>  Log-likelihood: -1845.2803
#>             AIC: 3701.56
#>             BIC: 3758.92
#> Condition#(Cov): 143.71
#> Condition#(Cor): 202.97
#> 
#> -- Timing Information --
#>      Best Chain: 6.2810 seconds
#>      All Chains: 6.2834 seconds
#>      Covariance: 28.8993 seconds
#>         Elapsed: 35.18 seconds
#> 
#> -- Population Parameters --
#> # A tibble: 6 × 6
#>   Parameter        Est.      SE  `%RSE` `Back-transformed(95%CI)` `BSV(CV%)`
#>   <chr>           <dbl>   <dbl>   <dbl> <chr>                          <dbl>
#> 1 cl             1.60    0.0153   0.955 4.96 (4.81, 5.11)               31.9
#> 2 v1             2.31    0.0839   3.64  10.05 (8.52, 11.84)             33.3
#> 3 v2             3.41    0.0393   1.16  30.12 (27.88, 32.53)            31.0
#> 4 q              2.28    0.0213   0.931 9.82 (9.42, 10.24)              32.2
#> 5 ka             0.0181  0.0791 437.    1.02 (0.87, 1.19)               30.0
#> 6 Residual Error 0.0394 NA       NA     0.0394                          NA  
#> 
#> -- Iteration Diagnostics --
#>  Iter | NLL and Parameters
#> --------------------------------------------------------------------------------
#>    1: -1839.577 1.609 2.303 3.401 2.303 0.000 -2.408 -2.408 -2.408 -2.408 -2.408 -3.219
#>    2: -1845.238 1.601 2.309 3.404 2.284 0.020 -2.279 -2.167 -2.333 -2.245 -2.439 -3.235
#>    3: -1845.277 1.600 2.307 3.406 2.285 0.017 -2.287 -2.195 -2.342 -2.267 -2.416 -3.235
#>    4: -1845.277 1.600 2.306 3.406 2.285 0.017 -2.287 -2.195 -2.342 -2.266 -2.416 -3.235
#>    5: -1845.277 1.600 2.306 3.406 2.285 0.017 -2.287 -2.195 -2.342 -2.266 -2.416 -3.235
#>    ... (omitted iterations) ...
#>   10: -1845.280 1.601 2.307 3.405 2.285 0.018 -2.286 -2.200 -2.343 -2.264 -2.411 -3.235
#>   11: -1845.280 1.601 2.307 3.405 2.285 0.018 -2.286 -2.200 -2.343 -2.264 -2.411 -3.235
#>   12: -1845.280 1.601 2.307 3.405 2.285 0.018 -2.286 -2.200 -2.343 -2.264 -2.411 -3.235
#>   13: -1845.280 1.601 2.307 3.405 2.285 0.018 -2.286 -2.201 -2.343 -2.265 -2.410 -3.235
#>   14: -1845.280 1.601 2.307 3.405 2.285 0.018 -2.286 -2.201 -2.343 -2.265 -2.410 -3.235
```

## Documentation

The package documentation is available at
<https://hiddevandebeek.github.io/admr/>. Key documentation sections
include:

- [Getting
  Started](https://hiddevandebeek.github.io/admr/articles/getting_started.html)
- [Advanced
  Topics](https://hiddevandebeek.github.io/admr/articles/advanced_modeling.html)
- [Function
  Reference](https://hiddevandebeek.github.io/admr/reference/index.html)
- [Examples](https://hiddevandebeek.github.io/admr/articles/examples.html)

## Use Cases

The `admr` package is particularly useful for:

- **Meta-Analysis**: Combining data from multiple studies where
  individual-level data is not available
- **Literature-Based Modeling**: Fitting models to published summary
  statistics
- **Simulation Studies**: Evaluating model performance with aggregate
  data
- **Population PK/PD**: Fitting complex models to summary-level data

## Getting Help

- Check the [documentation](https://hiddevandebeek.github.io/admr/)
- Browse [GitHub issues](https://github.com/hiddevandebeek/admr/issues)
- Create a new issue with a reproducible example

## Citation

If you use `admr` in your research, please cite it as:

``` r
citation("admr")
```

## License

This project is licensed under the GPL-2 License - see the
[LICENSE](LICENSE) file for details.
