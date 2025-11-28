# Compute mean and covariance of a matrix

`meancov` computes the mean and covariance of a matrix, optionally with
weights. This function is used to convert raw data into aggregate form
(mean and covariance) for use in aggregate data modeling.

## Usage

``` r
meancov(m, wt)
```

## Arguments

- m:

  A numeric matrix or data frame containing the observations. Each row
  represents an individual, and each column represents a time point. For
  pharmacometric data, columns typically represent concentration
  measurements at different time points.

- wt:

  Optional vector of weights for each observation. If not provided, all
  observations are weighted equally. Weights can be used to account for
  different sample sizes or reliability of different data sources.

## Value

A list containing:

- `E`: Vector of means for each time point (population typical values)

- `V`: Covariance matrix representing the variability between
  individuals

## Details

The function computes:

- The mean of each column (time point) using `colMeans` for unweighted
  data or weighted means for weighted data

- The covariance matrix using `cov.wt` with maximum likelihood
  estimation, which provides unbiased estimates of the population
  covariance

The maximum likelihood estimation method is used because:

- It provides unbiased estimates of the covariance matrix

- It is appropriate for aggregate data modeling where we want to
  estimate population parameters

- It handles both balanced and unbalanced designs through optional
  weights

Key features:

- Handles missing data automatically through the underlying `cov.wt`
  function

- Provides numerically stable computations

- Can be used with both raw PK data and simulated data

- Supports weighted calculations for meta-analysis or combined analysis

## Examples

``` r
# Load required libraries
library(admr)
library(rxode2)

# Create a matrix of concentration measurements
# 10 subjects measured at 10 time points
m <- matrix(rnorm(100), nrow = 10, ncol = 10)

# Compute unweighted mean and covariance
# Useful for single-study analysis
result <- meancov(m)

# Compute weighted mean and covariance
# Useful for meta-analysis or when combining studies
weights <- runif(10)  # weights could represent study sizes
result_weighted <- meancov(m, weights)
```
