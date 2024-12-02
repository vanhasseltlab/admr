#' Fitting aggregate data
#'
#'[fitEM()] implements the Expectation-Maximization(EM) algorithm for parameter
#'estimation of the given aggregate data model, iterating over maximum
#'likelihood updates with weighted MC updates. This version of the function uses nloptr instead of optimx. TOL = 1e-10
#'
#' @param p0 initial parameter values
#' @param opts options
#' @param obs observed data
#' @param maxiter maximum number of iterations
#' @param convcrit_nll convergence criterion for the negative log-likelihood
#' @param nomap single model or multiple models
#' @param phase_fractions vector of phase fractions
#' @param max_worse_iterations maximum number of consecutive worse iterations before skipping a phase
#' @returns A fitted model
#' @export
#' @examples
#' #test
#'

fitEM <- function(p0, opts, obs, maxiter = 100, convcrit_nll = 0.001, nomap = TRUE,
                  phase_fractions = c(0.1, 0.5, 0.4), max_worse_iterations = 10) {

  # Ensure phase_fractions sum to 1
  if (abs(sum(phase_fractions) - 1) > .Machine$double.eps^0.5) {
    stop("The sum of phase_fractions must be 1.")
  }

  start_time <- Sys.time()  # Initialize start time for total elapsed time

  if (nomap) {
    opts <- opts %>% p2opts(p0) %>% obs2opts(obs)
  } else {
    opts <- map(seq_along(opts), function(i) opts[[i]] %>% p2opts(p0) %>% obs2opts(obs[[i]]))
  }

  res <- tibble(
    iter = 1:maxiter,
    parameters = vector("list", maxiter),
    nll = NA_real_,
    approx_nll = NA_real_,
    iteration_time = 0
  )

  res$parameters[[1]] <- p0

  if (nomap) {
    res$nll[1] <- maxfunc(opts)(p0)
  } else {
    res$nll[1] <- Reduce('+', map(opts, ~ maxfunc(.)(p0)))
  }

  res$approx_nll[1] <- res$nll[1]

  # Display header for live output
  pvals <- rep(0, length(p0) + 1)
  param_count <- length(p0)
  cat(sprintf("Iter | %-12s\n", sprintf("NLL and Parameters (%d values)", param_count)))
  cat(sprintf("%s\n", strrep("-", 80)))
  cat(sprintf("%4d: %s\n", 1, paste(sprintf("%8.3f", c(res$nll[1], p0)), collapse = " ")))

  # Define phases with bounds, tolerance, and max evaluations
  phase_params <- list(
    list(bounds = 2, tol = 1e-10, maxeval = 2000, phase_name = "Exploration Phase"),
    list(bounds = 0.05, tol = 1e-10, maxeval = 2000, phase_name = "Refinement Phase"),
    list(bounds = 0.005, tol = 1e-10, maxeval = 2000, phase_name = "Final Phase")
  )

  # Validate the number of phase_fractions matches the number of phases
  if (length(phase_fractions) != length(phase_params)) {stop("Number of phase_fractions must match the number of phases.")}

  # Calculate start and end iterations for each phase
  cumulative_iters <- cumsum(floor(phase_fractions * maxiter))
  phase_start_iters <- c(2, head(cumulative_iters, -1) + 1)
  phase_end_iters <- cumulative_iters

  # Iterate through each phase
  best_nll <- res$nll[1]
  best_params <- res$parameters[[1]]
  current_iter <- 2

  for (phase_idx in seq_along(phase_params)) {
    current_phase <- phase_params[[phase_idx]]
    cat(sprintf("\n### %s ###\n", current_phase$phase_name))

    start_iter <- current_iter
    end_iter <- phase_end_iters[phase_idx]

    # Initialize starting point for the phase
    if (phase_idx > 1) {p0 <- best_params}  # Use best parameters from the previous phase

    phase_converged <- FALSE  # Flag to track convergence within the phase
    worse_counter <- 0  # Counter for worse iterations

    # Loop over iterations within each phase
    for (i in start_iter:end_iter) {
      iter_start_time <- Sys.time()  # Start timing this iteration

      # Define maxfunc for current parameters
      if (nomap) {
        ff <- maxfunc(p2opts(opts, p0))
      } else {
        ffs <- map(opts, ~ maxfunc(p2opts(., p0)))
        ff <- function(p) Reduce('+', map(ffs, ~ .(p)))
      }

      ff_nloptr <- function(params) {ff(params)}

      # Perform optimization
      m0 <- nloptr::nloptr(
        x0 = p0,
        eval_f = ff_nloptr,
        lb = p0 - current_phase$bounds,
        ub = p0 + current_phase$bounds,
        opts = list(
          algorithm = "NLOPT_LN_BOBYQA",
          xtol_rel = current_phase$tol,
          maxeval = current_phase$maxeval
        )
      )

      # Update parameters and results
      p0 <- m0$solution
      res$parameters[[i]] <- p0

      if (nomap) {
        res$nll[i] <- maxfunc(p2opts(opts, p0))(p0)
      } else {
        res$nll[i] <- Reduce('+', map(opts, ~ maxfunc(p2opts(., p0))(p0)))
      }

      res$approx_nll[i] <- m0$objective
      res$iteration_time[i] <- as.numeric(difftime(Sys.time(), iter_start_time, units = "secs"))

      # Update the best parameters and NLL
      if (i > 1 && res$nll[i] > res$nll[i - 1]) {
        worse_counter <- worse_counter + 1  # Increment if the current NLL is worse than the previous
      } else {
        worse_counter <- 0  # Reset if the current NLL improves or stays the same
      }

      if (res$nll[i] < best_nll) {
        best_nll <- res$nll[i]
        best_params <- p0
      }

      cat(sprintf("%4d: %s\n", i, paste(sprintf("%8.3f", c(res$nll[i], p0)), collapse = " ")))

      # Convergence check
      if (abs(res$nll[i] - res$approx_nll[i]) < convcrit_nll) {
        cat(sprintf("Converged in %d iterations during %s.\n", i - start_iter + 1, current_phase$phase_name))
        current_iter <- i + 1  # Update the next iteration start point
        phase_converged <- TRUE
        break
      }

      # Check if the phase should be skipped due to worsening NLL
      if (worse_counter >= max_worse_iterations) {
        cat(sprintf("Skipping %s due to %d consecutive worsening iterations.\n", current_phase$phase_name, max_worse_iterations))
        break  # Do not advance `current_iter` prematurely
      }

      # Update the next iteration start point
      current_iter <- i + 1
    }

    # Print message if the phase did not converge
    if (!phase_converged) {cat(sprintf("Did not converge during %s.\n", current_phase$phase_name))}
  }

  # After completing all phases, calculate the final results
  final_nll <- best_nll
  final_params <- best_params

  # Remove empty iterations
  res <- res[!is.na(res$nll), ]

  # Covariance computation and diagnostics
  cov_start_time <- Sys.time()
  func_fixed <- function(beta) {
    p_full <- final_params
    p_full[1:length(beta)] <- beta
    opts_list <- if (nomap) list(opts) else opts
    obs_list <- if (nomap) list(obs) else obs
    sum(sapply(seq_along(opts_list), function(i) genfitfunc(opts_list[[i]], obs_list[[i]])(p_full)))
  }

  cov_matrix_fixed <- tryCatch(
    solve(numDeriv::hessian(func_fixed, final_params[1:length(if (nomap) opts$p$beta else opts[[1]]$p$beta)],
                            method = "Richardson", method.args = list(r=2))),
    error = function(e) NA
  )

  cov_time <- as.numeric(difftime(Sys.time(), cov_start_time, units = "secs"))

  se_fixed <- if (is.matrix(cov_matrix_fixed)) {
    sqrt(diag(cov_matrix_fixed))
  } else {
    rep(NA, length(if (nomap) opts$p$beta else opts[[1]]$p$beta))
  }
  back_transformed_params <- tryCatch(
    if (nomap) opts$ptrans(final_params) else opts[[1]]$ptrans(final_params),
    error = function(e) NA
  )
  bic <- if (nomap) {
    -2 * final_nll + log(opts$n) * length(final_params)
  } else {
    -2 * final_nll + (log(opts[[1]]$n) + log(opts[[2]]$n)) * length(final_params)
  }

  # Generate output
  output <- list(
    final_params = final_params,
    transformed_params = back_transformed_params,
    BSV = tryCatch(sqrt(diag(back_transformed_params$Omega)) * 100, error = function(e) NA),
    standard_errors = se_fixed,
    covariance_matrix = cov_matrix_fixed,
    convergence_info = list(
      converged = nrow(res) < maxiter,
      total_iterations = nrow(res),
      final_nll = final_nll,
      total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),  # Use start_time for total elapsed time
      covariance_time = cov_time,
      aic = -2 * final_nll + length(final_params),
      bic = bic
    ),
    diagnostics = list(
      nll_trace = res$nll,
      approx_nll_trace = res$approx_nll,
      iteration_times = res$iteration_time,
      time_per_iteration_summary = summary(res$iteration_time, na.rm = TRUE)
    ),
    iteration_history = res
  )

  class(output) <- "fitEM_result"
  return(output)
}



#' @export
print.fitEM_result <- function(x, ...) {
  cat("-- FitEM Summary --\n\n")

  # Objective function and information criteria
  cat("-- Objective Function and Information Criteria --\n")
  conv <- x$convergence_info
  cat(sprintf(" Log-likelihood: %.4f\n", conv$final_nll))
  cat(sprintf("            AIC: %.2f\n", conv$aic))
  cat(sprintf("            BIC: %.2f\n", conv$bic))

  if (!is.null(x$covariance_matrix)) {
    cat(sprintf("Condition#(Cov): %.2f\n", kappa(x$covariance_matrix)))
    cat(sprintf("Condition#(Cor): %.2f\n", kappa(cov2cor(x$covariance_matrix))))
  }
  cat("\n")

  # Timing information
  cat("-- Timing Information --\n")
  diag <- x$diagnostics
  cat(sprintf("      Iteration: %.4f seconds\n", conv$total_time))
  cat(sprintf("     Covariance: %.4f seconds\n", conv$covariance_time))
  cat(sprintf("        Elapsed: %.2f seconds\n", conv$total_time + conv$covariance_time))
  cat("\n")

  # Population parameters and residual error
  cat("-- Population Parameters --\n")
  beta_params <- x$final_params[1:length(x$standard_errors)]
  beta_se <- x$standard_errors
  residual_error <- exp(x$final_params[length(x$final_params)])

  # Transformations and confidence intervals
  beta_estimates_original <- x$transformed_params$beta[1:length(beta_params)]
  beta_se_original <- beta_estimates_original * beta_se
  confint_low_original <- exp(beta_params - 1.96 * beta_se)
  confint_high_original <- exp(beta_params + 1.96 * beta_se)

  # Create the parameter table
  param_df <- tibble(
    Parameter = c(paste0("Beta", seq_along(beta_params)), "Residual Error"),
    `Est.` = c(beta_params, residual_error),
    `SE` = c(beta_se_original, NA),
    `%RSE` = c(100 * beta_se / abs(beta_params), NA),
    `Back-transformed(95%CI)` = c(
      paste0(
        sprintf("%.2f", beta_estimates_original), " (",
        sprintf("%.2f", confint_low_original), ", ",
        sprintf("%.2f", confint_high_original), ")"
      ),
      sprintf("%.4f", residual_error)
    ),
    `BSV(CV%)` = c(x$BSV, NA)
  )

  print(param_df)
  cat("\n")

  # Iteration diagnostics
  cat("-- Iteration Diagnostics --\n")
  cat(sprintf(" Iter | %-12s\n", "NLL and Parameters"))
  cat(sprintf("%s\n", strrep("-", 80)))
  for (i in seq_along(diag$nll_trace)) {
    nll <- diag$nll_trace[i]
    params <- x$iteration_history$parameters[[i]]
    cat(sprintf("%4d: %s\n", i, paste(sprintf("%.3f", c(nll, params)), collapse = " ")))
  }
  cat("\n")
}
