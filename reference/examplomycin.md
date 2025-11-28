# Examplomycin Dataset

A simulated dataset for the fictional drug examplomycin. The dataset
contains 500 subjects, each with 9 timepoints. It was generated using a
two-compartment pharmacokinetic model with first-order absorption and
elimination. Random effects and residual errors are included to simulate
variability and noise.

## Usage

``` r
examplomycin
```

## Format

A data frame with 4500 rows and 6 variables:

- `ID`: Subject ID.

- `TIME`: Observation time (hours).

- `DV`: Observed drug concentration (mg/L).

- `AMT`: Amount of drug administered (mg).

- `EVID`: Event type indicator (0 for observation, 101 for dosing).

- `CMT`: Compartment number (1 for depot, 2 for central).

## Source

Generated using the `generate_data` function.

## Examples

``` r
#' # Load required libraries
library(admr)

# Load the dataset
data("examplomycin")

# View the first few rows
head(examplomycin)
#>    ID TIME    DV AMT EVID CMT
#> 1 460 0.00 0.000 100  101   1
#> 2 460 0.10 0.752   0    0   2
#> 3 460 0.25 1.932   0    0   2
#> 4 460 0.50 3.694   0    0   2
#> 5 460 1.00 3.479   0    0   2
#> 6 460 2.00 4.003   0    0   2
```
