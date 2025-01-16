#' Fitting aggregate data
#'
#'[fitEM()] implements the Expectation-Maximization(EM) algorithm for parameter
#'estimation of the given aggregate data model, iterating over maximum
#'likelihood updates with weighted MC updates. This version of the function uses nloptr instead of optimx. TOL = 1e-10
#'
#' @param opts options
#' @param obs observed data
#' @param maxiter maximum number of iterations
#' @param convcrit_nll convergence criterion for the negative log-likelihood
#' @param nomap single model or multiple models
#' @param phase_fractions vector of phase fractions
#' @param max_worse_iterations maximum number of consecutive worse iterations before skipping a phase
#' @param chains number of chains
#' @param pertubation pertubation factor for the initial parameter values of each chain
#' @returns A fitted model
#' @export
#' @examples
#' #test
#'

fitEM <- function(opts, obs, maxiter = 100, convcrit_nll = 0.00001, nomap = TRUE,
                  phase_fractions = c(0.2, 0.4, 0.2, 0.2), max_worse_iterations = 10, chains = 1,
                  pertubation = 0.1) {
  init <- opts$pt

  # Ensure phase_fractions sum to 1
  if (abs(sum(phase_fractions) - 1) > .Machine$double.eps^0.5) {
    stop("The sum of phase_fractions must be 1.")
  }

  start_time <- Sys.time()  # Start time for the entire process

  perturb_init <- function(init) {
    init + rnorm(length(init), mean = 0, sd = pertubation * abs(init))  # Slight perturbation
  }

  ff_nloptr <- function(params) {
    ff(params)
  }

  grad_ff_nloptr <- function(params) {
    if (nomap == TRUE) {
      numDeriv::grad(maxfunc(opts), params)
    } else {
      Reduce('+', map(opts, ~ numDeriv::grad(maxfunc(.), params)))
      return(print("Gradient not yet available for nomap = FALSE"))
    }
  }

  chain_results <- vector("list", chains)  # Store results for each chain

  for (chain in seq_len(chains)) {
    chain_start_time <- Sys.time()  # Start time for the chain
    chain_init <- if (chain == 1) init else perturb_init(init)

    if (nomap) {
      opts <- opts %>% p2opts(chain_init) %>% obs2opts(obs)
    } else {
      opts <- map(seq_along(opts), function(i) opts[[i]] %>% p2opts(chain_init) %>% obs2opts(obs[[i]]))
    }

    res <- tibble(
      iter = 1,
      parameters = vector("list", maxiter),
      nll = NA_real_,
      approx_nll = NA_real_,
      iteration_time = 0
    )

    res$parameters[[1]] <- chain_init

    if (nomap) {
      res$nll[1] <- maxfunc(opts)(chain_init)
    } else {
      res$nll[1] <- Reduce('+', map(opts, ~ maxfunc(.)(chain_init)))
    }

    res$approx_nll[1] <- res$nll[1]

    if (chain == 1) {  # Print live output for single-chain mode
      cat(sprintf("Chain %d:\nIter | %-12s\n", chain, sprintf("NLL and Parameters (%d values)", length(init))))
      cat(sprintf("%s\n", strrep("-", 80)))
      cat(sprintf("%4d: %s\n", 1, paste(sprintf("%8.3f", c(res$nll[1], chain_init)), collapse = " ")))
    }

    phase_params <- list(
      list(bounds = 2, maxeval = 5000, phase_name = "Exploration Phase"),
      list(bounds = 1, maxeval = 5000, phase_name = "Refinement Phase"),
      list(bounds = 0.5, maxeval = 5000, phase_name = "Final Phase"),
      list(bounds = 0.01, maxeval = 5000, phase_name = "Final Phase 2")
    )

    if (length(phase_fractions) != length(phase_params)) {
      stop("Number of phase_fractions must match the number of phases.")
    }

    cumulative_iters <- cumsum(c(floor(phase_fractions[-length(phase_fractions)] * maxiter), maxiter))
    phase_start_iters <- c(2, head(cumulative_iters, -1) + 1)
    phase_end_iters <- cumulative_iters

    best_nll <- res$nll[1]
    best_params <- res$parameters[[1]]
    current_iter <- 2

    for (phase_idx in seq_along(phase_params)) {
      current_phase <- phase_params[[phase_idx]]
      if (chain == 1) {
        cat(sprintf("\n### %s ###\n", current_phase$phase_name))
      }

      # Start each phase with the best parameters found so far
      chain_init <- best_params
      phase_converged <- FALSE
      worse_counter <- 0

      # Start this phase from the current global iteration number

      for (i in current_iter:phase_end_iters[phase_idx]) {
        if (current_iter > maxiter) {
          cat("Maximum iterations reached. Terminating optimization.\n")
          break
        }

        if (phase_idx == 1) {
          algorithm <- "NLOPT_LN_BOBYQA"
          gradient <- FALSE
        } else if (phase_idx == 2) {
          algorithm <- "NLOPT_LN_BOBYQA"
          gradient <- FALSE
        } else {
          algorithm <- "NLOPT_LN_BOBYQA"
          gradient <- FALSE
        }

        iter_start_time <- Sys.time()

        if (nomap) {
          ff <- maxfunc(p2opts(opts, chain_init))
        } else {
          ffs <- map(opts, ~ maxfunc(p2opts(., chain_init)))
          ff <- function(p) Reduce('+', map(ffs, ~ .(p)))
        }

        if (gradient == TRUE) {
          m0 <- nloptr::nloptr(
            x0 = chain_init,
            eval_f = ff_nloptr,
            eval_grad_f = grad_ff_nloptr,
            lb = chain_init - current_phase$bounds,
            ub = chain_init + current_phase$bounds,
            opts = list(
              algorithm = algorithm,
              ftol_rel = 1e-10,
              maxeval = current_phase$maxeval
            )
          )
        } else {
          m0 <- nloptr::nloptr(
            x0 = chain_init,
            eval_f = ff_nloptr,
            lb = chain_init - current_phase$bounds,
            ub = chain_init + current_phase$bounds,
            opts = list(
              algorithm = algorithm,
              ftol_rel = 1e-10,
              maxeval = current_phase$maxeval
            )
          )
        }

        # Optimization step
        chain_init <- m0$solution
        res$parameters[[current_iter]] <- chain_init

        if (nomap) {
          res$nll[current_iter] <- maxfunc(p2opts(opts, chain_init))(chain_init)
        } else {
          res$nll[current_iter] <- Reduce('+', map(opts, ~ maxfunc(p2opts(., chain_init))(chain_init)))
        }

        res$approx_nll[current_iter] <- m0$objective
        res$iteration_time[current_iter] <- as.numeric(difftime(Sys.time(), iter_start_time, units = "secs"))

        # Update best parameters if a better NLL is found
        if (res$nll[current_iter] < best_nll) {
          best_nll <- res$nll[current_iter]
          best_params <- chain_init
        }

        if (chain == 1) {
          cat(sprintf("%4d: %s\n", current_iter, paste(sprintf("%8.3f", c(res$nll[current_iter], chain_init)), collapse = " ")))
        }

        # Convergence check
        if (abs(res$nll[current_iter] - res$approx_nll[current_iter]) < convcrit_nll) {
          phase_converged <- TRUE
          current_iter <- current_iter + 1
          break  # Skip the rest of this phase, continue to the next phase
        }

        # Early phase termination due to worse counter
        worse_counter <- if (res$nll[current_iter] > best_nll) worse_counter + 1 else 0
        if (worse_counter >= max_worse_iterations) {
          current_iter <- current_iter + 1
          break  # Skip the rest of this phase, continue to the next phase
        }
        current_iter <- current_iter + 1  # Increment after every iteration
      }

      # If the phase converged, note it but allow the next phase to start
      if (phase_converged) {
        cat(sprintf("Phase %s converged at iteration %d. \n",
                    current_phase$phase_name, current_iter - 1))
      }
    }

    res <- res[!is.na(res$nll), ]
    chain_time <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))
    cat(sprintf("\nChain %d Complete: Final NLL = %.3f, Time Elapsed = %.2f seconds\n \n",
                chain, best_nll, chain_time))
    chain_results[[chain]] <- list(res = res, best_nll = best_nll, best_params = best_params, time = chain_time)
  }

  # Select the best chain based on NLL
  best_chain <- which.min(sapply(chain_results, function(x) x$best_nll))
  final_result <- chain_results[[best_chain]]
  final_params <- final_result$best_params
  final_nll <- final_result$best_nll
  final_res <- final_result$res  # Final iteration history

  # Covariance computation
  cov_start_time <- Sys.time()
  func_fixed <- function(beta) {
    p_full <- final_params
    p_full[1:length(beta)] <- beta
    opts_list <- if (nomap) list(opts) else opts
    obs_list <- if (nomap) list(obs) else obs
    sum(sapply(seq_along(opts_list), function(i) genfitfunc(opts_list[[i]], obs_list[[i]])(p_full)))
  }

  cov_matrix_fixed <- tryCatch(
    solve(numDeriv::hessian(func_fixed, final_params[1:length(if (nomap) opts$p$beta else opts[[2]]$p$beta)],
                            method = "Richardson", method.args = list(r=6, v=2))),
    error = function(e) NA
  )

  cov_time <- as.numeric(difftime(Sys.time(), cov_start_time, units = "secs"))

  tryCatch(
    se_fixed <- if (is.matrix(cov_matrix_fixed)) {
      sqrt(diag(cov_matrix_fixed))
    } else {
      rep(NA, length(if (nomap) opts$p$beta else opts[[1]]$p$beta))
    }, error = function(e) NULL)

  back_transformed_params <- tryCatch(
    if (nomap) opts$ptrans(final_params) else opts[[1]]$ptrans(final_params),
    error = function(e) NA
  )

  bic <- if (nomap) {
    -2 * final_nll + log(opts$n) * length(final_params)
  } else {
    -2 * final_nll + (log(opts[[1]]$n) + log(opts[[2]]$n)) * length(final_params)
  }

  param_names <- names(back_transformed_params$beta)
  beta_params <- final_params[1:length(se_fixed)]
  beta_se <- se_fixed
  residual_error <- exp(final_params[length(final_params)])

  # Transformations and confidence intervals
  beta_estimates_original <- back_transformed_params$beta[1:length(beta_params)]
  confint_low_original <- exp(beta_params - 1.96 * beta_se)
  confint_high_original <- exp(beta_params + 1.96 * beta_se)

  # Generate output
  output <- list(
    final_params = final_params,
    transformed_params = back_transformed_params,
    param_df = tibble(
      Parameter = c(param_names, "Residual Error"),
      `Est.` = c(beta_params, residual_error),
      `SE` = c(beta_se, NA),
      `%RSE` = c(100 * beta_se / abs(beta_params), NA),
      `Back-transformed(95%CI)` = c(
        paste0(
          sprintf("%.2f", beta_estimates_original), " (",
          sprintf("%.2f", confint_low_original), ", ",
          sprintf("%.2f", confint_high_original), ")"
        ),
        sprintf("%.4f", residual_error)
      ),
      `BSV(CV%)` = c(tryCatch(sqrt(diag(back_transformed_params$Omega)) * 100, error = function(e) NA), NA)
    ),
    covariance_matrix = cov_matrix_fixed,
    convergence_info = list(
      converged = nrow(final_res) < maxiter,
      total_iterations = nrow(final_res),
      final_nll = final_nll,
      total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      chain_time = final_result$time,
      covariance_time = cov_time,
      aic = -2 * final_nll + length(final_params),
      bic = bic
    ),
    diagnostics = list(
      best_chain = best_chain,
      nll_trace = final_res$nll,
      approx_nll_trace = final_res$approx_nll,
      iteration_times = final_res$iteration_time,
      time_per_iteration_summary = summary(final_res$iteration_time, na.rm = TRUE)
    ),
    chain_results = chain_results,  # Save results of all chains
    iteration_history = final_res,
    data = list(
      opts = opts,
      obs = obs
    )
  )

  class(output) <- "fitEM_result"
  return(output)
}

#' Print fitEM results
#'
#' @param x A fitted model object returned by fitEM
#' @param ... Additional arguments (not used)
#' @export
#'

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
  cat(sprintf("     Best Chain: %.4f seconds\n", conv$chain_time))
  cat(sprintf("     All Chains: %.4f seconds\n", conv$total_time - conv$covariance_time))
  cat(sprintf("     Covariance: %.4f seconds\n", conv$covariance_time))
  cat(sprintf("        Elapsed: %.2f seconds\n", conv$total_time))
  cat("\n")

  # Population parameters and residual error
  cat("-- Population Parameters --\n")

  # Create the parameter table
  print(x$param_df)
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


#' Plot diagnostics of fitEM results
#'
#' @param x A fitted model object returned by fitEM
#' @param ... Additional arguments (not used)
#' @export
#'

plot.fitEM_result <- function(x, ...) {
  if (!inherits(x, "fitEM_result")) {
    stop("The input must be a fitEM_result object.")
  }

  chain_results <- x$chain_results

  # Prepare data for NLL trace
  nll_data <- do.call(rbind, lapply(seq_along(chain_results), function(chain) {
    data.frame(
      Iteration = seq_len(nrow(chain_results[[chain]]$res)),
      NLL = chain_results[[chain]]$res$nll,
      Chain = factor(chain)
    )
  }))

  # Extract parameter names from the first iteration of the first chain
  param_names <- names(chain_results[[1]]$res$parameters[[1]])

  # Prepare data for parameter convergence
  param_data <- do.call(rbind, lapply(seq_along(chain_results), function(chain) {
    do.call(rbind, lapply(seq_len(nrow(chain_results[[chain]]$res)), function(i) {
      param_iter <- chain_results[[chain]]$res$parameters[[i]]

      if (is.null(param_iter) || length(param_iter) == 0) {
        # Handle missing or empty parameters for this iteration
        data.frame(
          Iteration = i,
          Parameter = param_names,  # Use names from the first iteration of the first chain
          Value = NA,  # No values for missing parameters
          Chain = factor(chain),
          stringsAsFactors = FALSE
        )
      } else {
        # Align parameter values with the reference names
        values <- unlist(param_iter)
        if (length(values) != length(param_names)) {
          stop("Parameter values do not match the length of reference names.")
        }
        data.frame(
          Iteration = i,
          Parameter = param_names,
          Value = values,
          Chain = factor(chain),
          stringsAsFactors = FALSE
        )
      }
    }))
  }))

  # Reset row names to avoid mismatched indexing
  row.names(param_data) <- NULL

  # Plot 1: NLL Convergence Trace
  p1 <- ggplot(nll_data, aes(x = .data$Iteration, y = .data$NLL, color = .data$Chain, group = .data$Chain)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = min(nll_data$NLL, na.rm = TRUE), linetype = "dashed", color = "red") +
    labs(title = "NLL Convergence Trace", x = "Iteration", y = "Negative Log-Likelihood") +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot 2: Parameter Convergence
  p2 <- ggplot(param_data, aes(x = .data$Iteration, y = .data$Value, color = .data$Chain, group = interaction(.data$Chain, .data$Parameter))) +
    geom_line() +
    geom_point() +
    facet_wrap(~Parameter, scales = "free_y") +
    labs(
      title = "Parameter Convergence by Chain",
      x = "Iteration",
      y = "Parameter Value"
    ) +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(plot.title = element_text(hjust = 0.5))

  # Plot 3: Sensitivity Convergence

  # Print the plots
  print(p1)
  print(p2)
}
