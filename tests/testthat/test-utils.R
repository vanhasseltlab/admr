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

test_that("meancov computes correct statistics", {
  # Create test matrix
  set.seed(123)
  test_matrix <- matrix(rnorm(100), nrow = 10)

  # Test without weights
  result <- meancov(test_matrix)
  expect_named(result, c("E", "V"))
  expect_equal(result$E, colMeans(test_matrix))
  expect_equal(dim(result$V), c(ncol(test_matrix), ncol(test_matrix)))

  # Test with weights
  weights <- runif(10)
  weighted_result <- meancov(test_matrix, weights)
  expect_named(weighted_result, c("E", "V"))
  expect_equal(length(weighted_result$E), ncol(test_matrix))

  # Test with examplomycin data
  test_data <- create_test_data()
  expect_named(test_data$obs, c("E", "V"))
  expect_true(is.matrix(test_data$obs$V))
  expect_true(is.numeric(test_data$obs$E))
})

test_that("omegas creates correct covariance matrices", {
  # Test diagonal matrix
  diag_matrix <- omegas(0.09, 0, 3)
  expect_equal(dim(diag_matrix), c(3, 3))
  expect_equal(diag(diag_matrix), rep(0.09, 3))
  expect_equal(diag_matrix[upper.tri(diag_matrix)], rep(0, 3))

  # Test matrix with off-diagonal elements
  corr_matrix <- omegas(0.09, 0.01, 4)
  expect_equal(dim(corr_matrix), c(4, 4))
  expect_equal(diag(corr_matrix), rep(0.09, 4))
  expect_true(all(corr_matrix[upper.tri(corr_matrix)] == 0.01))

  # Test with examplomycin data
  test_data <- create_test_data()
  omega_matrix <- test_data$opts$p$Omega
  expect_equal(dim(omega_matrix), c(5, 5))
  expect_equal(diag(omega_matrix), rep(0.09, 5))
  expect_true(all(omega_matrix[upper.tri(omega_matrix)] == 0))
})

test_that("genfitfunc handles different input types", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Test with aggregate data
  fitfun <- genfitfunc(opts, obs)
  expect_type(fitfun, "closure")
  nll <- fitfun(opts$pt)
  expect_type(nll, "double")

  # Test with raw data matrix
  raw_data <- gendataset(opts, seed = 123)
  fitfun_raw <- genfitfunc(opts, meancov(raw_data))
  nll_raw <- fitfun_raw(opts$pt)
  expect_type(nll_raw, "double")

  # Test with expected data from model
  expected_data <- gen_pop_EV(opts)[1:2]
  fitfun_expected <- genfitfunc(opts, expected_data)
  nll_expected <- fitfun_expected(opts$pt)
  expect_type(nll_expected, "double")
})

test_that("parameter transformations work correctly", {
  test_data <- create_test_data()
  p <- test_data$opts$p

  # Test transformation
  p_transformed <- p_to_optim(p)
  expect_named(p_transformed, c("values", "backtransformfunc",
                               "d_psi_d_psitrans_long", "d_psi_d_psitrans_short",
                               "d_bi_d_omega", "d_omega_d_Omega"))

  # Test back-transformation
  p_back <- p_transformed$backtransformfunc(p_transformed$values)
  expect_equal(p_back$beta, p$beta)
  expect_equal(p_back$Omega, p$Omega)

  # Test derivatives
  #expect_true(is.matrix(p_transformed$d_psi_d_psitrans_long))
  #expect_true(is.matrix(p_transformed$d_psi_d_psitrans_short))
  #expect_true(is.matrix(p_transformed$d_bi_d_omega))
  #expect_true(is.matrix(p_transformed$d_omega_d_Omega))
})

test_that("gendataset generates correct data formats", {
  test_data <- create_test_data()
  opts <- test_data$opts

  # Test raw format
  raw_data <- gendataset(opts, seed = 123)
  expect_true(is.matrix(raw_data))
  expect_equal(dim(raw_data), c(opts$nsim, length(opts$time)))

  # Test nlmixr format
  nlmixr_data <- gendataset(opts, seed = 123, nlmixrform = TRUE)
  expect_true(is.data.frame(nlmixr_data))
  expect_named(nlmixr_data, c("dv", "time", "id", "amt", "evid", "cmt"))

  # Test without residual error
  no_error_data <- gendataset(opts, seed = 123, reserr = FALSE)
  expect_true(is.matrix(no_error_data))
  expect_equal(dim(no_error_data), c(opts$nsim, length(opts$time)))
})

test_that("error handling works correctly", {
  # To be done
})

test_that("convergence diagnostics work correctly", {
  test_data <- create_test_data()
  opts <- test_data$opts
  obs <- test_data$obs

  # Test IRMC convergence
  result <- fitIRMC(
    opts = opts,
    obs = obs,
    maxiter = 50,
    chains = 1
  )

  # Check convergence information
  expect_true(!is.null(result$convergence_info$converged))
  expect_true(!is.null(result$convergence_info$total_iterations))
  expect_true(result$convergence_info$total_iterations <= 50)

  # Check parameter stability
  params_history <- result$iteration_history$parameters
  final_params <- tail(params_history, 2)
  expect_true(all(abs(diff(do.call(rbind, final_params))) < 0.1))
})
