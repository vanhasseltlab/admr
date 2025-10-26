library(testthat)
library(admr)
library(rxode2)
library(dplyr)
library(tidyr)

# Helper function to create test data from examplomycin
create_test_data <- function() {
  # Use examplomycin dataset
  data(examplomycin)

  # Prepare data in wide format
  examplomycin_wide <- examplomycin %>%
    filter(EVID != 101) %>%
    dplyr::select(ID, TIME, DV) %>%
    pivot_wider(names_from = TIME, values_from = DV) %>%
    dplyr::select(-c(1))

  # Create aggregated data
  examplomycin_aggregated <- examplomycin_wide %>%
    admr::meancov()

  # Define RxODE model
  rxModel <- function(){
    model({
      cp = linCmt(
        cl,           # Clearance
        v1,           # Volume of central compartment
        v2,           # Volume of peripheral compartment
        q,            # Inter-compartmental clearance
        ka            # Absorption rate constant
      )})
  }

  rxModel <- rxode2(rxModel)
  rxModel <- rxModel$simulationModel

  # Define prediction function
  predder <- function(time, theta_i, dose = 100) {
    n_individuals <- nrow(theta_i)

    if (is.null(n_individuals)) {
      n_individuals <- 1
    }

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

  return(list(opts = opts, obs = examplomycin_aggregated))
}

test_that("fitIRMC basic functionality works with examplomycin data", {
  test_data <- create_test_data()

  # Run fitIRMC with basic settings
  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 10,
    chains = 1
  )

  # Check basic structure
  expect_s3_class(result, "fit_admr_result")
  expect_named(result, c("final_params", "transformed_params", "param_df",
                        "covariance_matrix", "convergence_info", "diagnostics",
                        "chain_results", "iteration_history", "data"))

  # Check parameter estimates are reasonable
  expect_true(all(!is.na(result$final_params)))
  expect_true(all(!is.na(result$param_df$Est.)))

  # Check that we have the correct number of parameters
  expect_equal(length(result$final_params), 11)  # 5 beta parameters + 5 Omega parameters + 1 residual error
})

test_that("fitIRMC handles multiple chains correctly with examplomycin data", {
  test_data <- create_test_data()

  # Run fitIRMC with multiple chains
  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 10,
    chains = 3,
    pertubation = 0.1
  )

  # Check chain results
  expect_equal(length(result$chain_results), 3)
  expect_true(all(sapply(result$chain_results, function(x) !is.null(x$best_nll))))

  # Check best chain selection
  expect_true(result$diagnostics$best_chain %in% 1:3)
})

test_that("fitIRMC handles convergence criteria correctly with examplomycin data", {
  test_data <- create_test_data()

  # Run with strict convergence criteria
  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 100,
    convcrit_nll = 1e-6,
    chains = 1
  )

  # Check convergence information
  expect_true(!is.null(result$convergence_info$converged))
  expect_true(!is.null(result$convergence_info$total_iterations))
})

test_that("fitIRMC handles phase fractions correctly with examplomycin data", {
  test_data <- create_test_data()

  # Test with custom phase fractions
  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 10,
    phase_fractions = c(0.3, 0.3, 0.2, 0.2),
    chains = 1
  )

  # Check that the function runs without errors
  expect_s3_class(result, "fit_admr_result")
})

test_that("fitIRMC handles errors gracefully with examplomycin data", {
  test_data <- create_test_data()

  # Test with invalid phase fractions
  expect_error(
    fitIRMC(
      opts = test_data$opts,
      obs = test_data$obs,
      phase_fractions = c(0.3, 0.3, 0.2, 0.3)  # Sums to 1.1
    ),
    "The sum of phase_fractions must be 1"
  )
})

test_that("fitIRMC produces valid covariance matrix with examplomycin data", {
  test_data <- create_test_data()

  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 10,
    chains = 1
  )

  # Check covariance matrix properties
  expect_true(is.matrix(result$covariance_matrix))
  expect_true(all(dim(result$covariance_matrix) == c(5, 5)))  # For all parameters
  expect_true(all(!is.na(result$covariance_matrix)))
})

test_that("fitIRMC produces valid parameter transformations with examplomycin data", {
  test_data <- create_test_data()

  result <- fitIRMC(
    opts = test_data$opts,
    obs = test_data$obs,
    maxiter = 10,
    chains = 1
  )

  # Check transformed parameters
  expect_true(!is.null(result$transformed_params))
  expect_true(!is.null(result$transformed_params$beta))
  expect_true(all(result$transformed_params$beta > 0))  # Since we're using exp transformation

  # Check that we have the correct parameter names
  expected_params <- c("cl", "v1", "v2", "q", "ka")
  expect_true(all(expected_params %in% names(result$transformed_params$beta)))
})
