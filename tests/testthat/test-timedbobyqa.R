library(testthat)
library(admr)
library(rxode2)
library(dplyr)
library(tidyr)

# Helper function to create test data
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

test_that("timedbobyqa works with examplomycin data", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Initial parameters
  init <- opts$pt

  # Run timedbobyqa
  result <- timedbobyqa(init, opts, obs)

  # Basic checks
  expect_true(!is.null(result))
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true(all(!is.na(result$nll)))
  expect_true(all(!is.na(result$time)))
  expect_true(all(!is.na(result$iter)))

  # Check parameter estimates
  expect_true(length(result$p[[1]]) == 11)  # 5 beta + 5 Omega + sigma parameters
  expect_true(all(sapply(result$p, length) == 11))

  # Check convergence
  expect_true(min(result$nll) < result$nll[1])

  # Check iteration history
  expect_true(all(diff(result$iter) == 1))
  expect_true(all(diff(result$time) >= 0))
})

test_that("timedbobyqa handles multiple models correctly", {
  test_data <- create_test_data()
  opts <- list(test_data$opts, test_data$opts)  # Create two identical models
  obs <- list(test_data$obs, test_data$obs)     # Create two identical observations

  # Initial parameters
  init <- opts[[1]]$pt

  # Run timedbobyqa with multiple models
  result <- timedbobyqa(init, opts, obs, nomap = FALSE)

  # Basic checks
  expect_true(!is.null(result))
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true(all(!is.na(result$nll)))

  # Check parameter estimates
  expect_true(length(result$p[[1]]) == 11)  # 5 beta + 5 Omega + sigma parameters
  expect_true(all(sapply(result$p, length) == 11))

  # Check convergence
  expect_true(min(result$nll) < result$nll[1])
})

test_that("timedbobyqa handles boundary constraints", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Initial parameters
  init <- opts$pt

  # Run timedbobyqa
  result <- timedbobyqa(init, opts, obs)

  # Check that parameters stay within bounds
  for (i in seq_along(result$p)) {
    params <- result$p[[i]]
    expect_true(all(params >= init - 2))
    expect_true(all(params <= init + 2))
  }
})

test_that("timedbobyqa handles iteration history", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Initial parameters
  init <- opts$pt

  # Run timedbobyqa
  result <- timedbobyqa(init, opts, obs)

  # Check iteration history
  expect_true(all(diff(result$iter) == 1))
  expect_true(all(diff(result$time) >= 0))
  expect_true(all(!is.na(result$nll)))
  expect_true(all(!is.na(result$time)))
})

test_that("timedbobyqa handles convergence", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Initial parameters
  init <- opts$pt

  # Run timedbobyqa
  result <- timedbobyqa(init, opts, obs)

  # Check convergence
  expect_true(min(result$nll) < result$nll[1])
  expect_true(all(diff(result$nll[1:3]) <= 0))  # First 3 iterations should show improvement
})
