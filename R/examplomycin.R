#' Examplomycin Dataset
#'
#' A simulated dataset for the fictional drug examplomycin. The dataset contains 500 subjects,
#' each with 9 timepoints. It was generated using a two-compartment pharmacokinetic model with
#' first-order absorption and elimination. Random effects and residual errors are included to
#' simulate variability and noise.
#'
#' @format A data frame with 4500 rows and 6 variables:
#'   - `ID`: Subject ID.
#'   - `TIME`: Observation time (hours).
#'   - `DV`: Observed drug concentration (mg/L).
#'   - `AMT`: Amount of drug administered (mg).
#'   - `EVID`: Event type indicator (0 for observation, 101 for dosing).
#'   - `CMT`: Compartment number (1 for depot, 2 for central).
#'
#' @examples
#' # Load the dataset
#' data("examplomycin")
#'
#' # View the first few rows
#' head(examplomycin)
#'
#' @source Generated using the `generate_data` function.
"examplomycin"
