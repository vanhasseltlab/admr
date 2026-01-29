utils::globalVariables(c(
  "time", "value", "lower_95", "upper_95",
  "lower_q1", "upper_q3", "time1", "time2"
))

#' Plot diagnostics of fitIRMC results
#'
#' @param x A fitted model object returned by fitIRMC
#' @param ... Additional arguments (not used)
#' @export

plot.fit_admr_result <- function(x, ...) {
  if (!inherits(x, "fit_admr_result")) {
    stop("The input must be a fit_admr_result object.")
  }

  chain_results <- x$chain_results
  nomap <- x$settings$single_dataframe
  opts <- x$data$opts

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
  if (is.null(param_names)) {
    param_names <- names(x$data$opts$pt)
  }
  if (is.null(param_names)) {
    param_names <- names(x$data$opts[[1]]$pt)
  }


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
  print(p1)

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
  print(p2)

  if (nomap == TRUE) {
    # Extract observations from the fitIRMC result
    time_points <- x$data$opts$time
    observations <- x$data$obs
    #predictions <- MCapprEV(upd_opts(x$data$opts, list(p = x$transformed_params)))
    predictions <- gen_pop_EV(upd_opts(x$data$opts, list(p = x$transformed_params)))

    # Create plot-ready boxplot data for observations and predictions
    obs_boxplot <- prepare_boxplot_df(observations$E, observations$V, time_points)
    pred_boxplot <- prepare_boxplot_df(predictions$E, predictions$V, time_points)
    obs_boxplot$source <- "Observed"
    pred_boxplot$source <- "Predicted"
    combined_boxplot <- rbind(obs_boxplot, pred_boxplot)

    # Prepare data frames for E and V heatmaps
    obs_E_df <- prepare_E_df(observations$E, "Observed", time_points)
    pred_E_df <- prepare_E_df(predictions$E, "Predicted", time_points)
    obs_V_df <- prepare_V_df(observations$V, "Observed", time_points)
    pred_V_df <- prepare_V_df(predictions$V, "Predicted", time_points)

    # --- Add residuals and standardized residuals ---
    res_E_df <- obs_E_df
    res_E_df$mean <- obs_E_df$mean - pred_E_df$mean
    res_E_df$source <- "Residual"

    std_res_E_df <- obs_E_df
    std_res_E_df$mean <- (obs_E_df$mean - pred_E_df$mean) / sqrt(pred_E_df$mean)
    std_res_E_df$source <- "Standardized Residual"

    combined_E_df <- rbind(obs_E_df, pred_E_df, res_E_df, std_res_E_df)

    res_V_df <- obs_V_df
    res_V_df$value <- obs_V_df$value - pred_V_df$value
    res_V_df$source <- "Residual"

    std_res_V_df <- pred_V_df
    std_res_V_df$value <- (obs_V_df$value - pred_V_df$value) / sqrt(abs(pred_V_df$value))
    std_res_V_df$source <- "Standardized Residual"

    combined_V_df <- rbind(obs_V_df, pred_V_df, res_V_df, std_res_V_df)

    # --- Plots ---
    # Boxplot
    p3 <- ggplot(combined_boxplot, aes(x = time)) +
      geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = source), alpha = 0.2, show.legend = FALSE) +
      geom_ribbon(aes(ymin = lower_q1, ymax = upper_q3, fill = source), alpha = 0.35) +
      geom_line(aes(y = mean, color = source, linetype = source), linewidth = 1) +
      scale_fill_manual(values = c("Observed" = "#3B8AC4", "Predicted" = "#D1495B")) +
      scale_color_manual(values = c("Observed" = "#3B8AC4", "Predicted" = "#D1495B")) +
      scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
      labs(title = "Observed vs Predicted: Mean and Variance Bands", x = "Time (h)", y = "Concentration") +
      theme_minimal(base_size = 14)
    print(p3)

    # E heatmaps with residuals
    p_E <- ggplot(combined_E_df, aes(x = time, y = 1, fill = mean)) +
      geom_tile() +
      facet_wrap(~source, ncol = 1) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = median(combined_E_df$mean, na.rm = TRUE)
      ) +
      labs(title = "Mean Vector (E) with Residuals", x = "Time", y = "", fill = "Value") +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    print(p_E)

    # V heatmaps with residuals
    combined_V_df$time1 <- factor(combined_V_df$time1, levels = sort(unique(combined_V_df$time1)))
    combined_V_df$time2 <- factor(combined_V_df$time2, levels = sort(unique(combined_V_df$time2)))

    p_V <- ggplot(combined_V_df, aes(x = time1, y = time2, fill = value)) +
      geom_tile(color = NA) +
      facet_wrap(~source, ncol = 2) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(title = "Variance-Covariance Matrix (V) with Residuals", x = "Time 1", y = "Time 2", fill = "Value") +
      coord_fixed() +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.ticks = element_blank())
    print(p_V)
  } else {
    for (i in 1:length(opts)){
      # Extract observations from the fitIRMC result
      time_points <- x$data$opts[[i]]$time
      observations <- x$data$obs[[i]]
      predictions <- gen_pop_EV(upd_opts(x$data$opts[[i]], list(p = x$transformed_params)))
      #predictions <- MCapprEV(upd_opts(x$data$opts[[i]], list(p = x$transformed_params)))

      # Create plot-ready boxplot data for observations and predictions
      obs_boxplot <- prepare_boxplot_df(observations$E, observations$V, time_points)
      pred_boxplot <- prepare_boxplot_df(predictions$E, predictions$V, time_points)
      obs_boxplot$source <- "Observed"
      pred_boxplot$source <- "Predicted"
      combined_boxplot <- rbind(obs_boxplot, pred_boxplot)

      # Prepare data frames for E and V heatmaps
      obs_E_df <- prepare_E_df(observations$E, "Observed", time_points)
      pred_E_df <- prepare_E_df(predictions$E, "Predicted", time_points)
      obs_V_df <- prepare_V_df(observations$V, "Observed", time_points)
      pred_V_df <- prepare_V_df(predictions$V, "Predicted", time_points)


      res_E_df <- obs_E_df
      res_E_df$mean <- obs_E_df$mean - pred_E_df$mean
      res_E_df$source <- "Residual"

      std_res_E_df <- obs_E_df
      std_res_E_df$mean <- (obs_E_df$mean - pred_E_df$mean) / sqrt(pred_E_df$mean)
      std_res_E_df$source <- "Standardized Residual"

      combined_E_df <- rbind(obs_E_df, pred_E_df, res_E_df, std_res_E_df)


      res_V_df <- obs_V_df
      res_V_df$value <- obs_V_df$value - pred_V_df$value
      res_V_df$source <- "Residual"

      std_res_V_df <- pred_V_df
      std_res_V_df$value <- (obs_V_df$value - pred_V_df$value) / sqrt(abs(pred_V_df$value))
      std_res_V_df$source <- "Standardized Residual"

      combined_V_df <- rbind(obs_V_df, pred_V_df, res_V_df, std_res_V_df)


      # Plot 3: Combined Boxplot with Mean and Variance Bands
      p3 <- ggplot(combined_boxplot, aes(x = time)) +
        # 95% interval ribbons
        geom_ribbon(
          aes(ymin = lower_95, ymax = upper_95, fill = source),
          alpha = 0.2,
          show.legend = FALSE
        ) +
        # IQR ribbons
        geom_ribbon(
          aes(ymin = lower_q1, ymax = upper_q3, fill = source),
          alpha = 0.35
        ) +
        # Mean lines
        geom_line(
          aes(y = mean, color = source, linetype = source),
          linewidth = 1
        ) +
        scale_fill_manual(values = c("Observed" = "#3B8AC4", "Predicted" = "#D1495B")) +
        scale_color_manual(values = c("Observed" = "#3B8AC4", "Predicted" = "#D1495B")) +
        scale_linetype_manual(values = c("Observed" = "solid", "Predicted" = "dashed")) +
        labs(
          title = paste("Observed vs Predicted: Mean and Variance Bands (Dataframe", i, ")"),
          x = "Time (h)",
          y = "Concentration",
          fill = "",
          color = "",
          linetype = ""
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size = 13)
        )
      print(p3)

      p_E <- ggplot(combined_E_df, aes(x = time, y = 1, fill = mean)) +
        geom_tile() +
        facet_wrap(~source, ncol = 1) +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = median(combined_E_df$mean, na.rm = TRUE)
        ) +
        labs(title = "Mean Vector (E)", x = "Time", y = "", fill = "Value") +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

      print(p_E)

      combined_V_df$time1 <- factor(combined_V_df$time1, levels = sort(unique(combined_V_df$time1)))
      combined_V_df$time2 <- factor(combined_V_df$time2, levels = sort(unique(combined_V_df$time2)))

      p_V <- ggplot(combined_V_df, aes(x = time1, y = time2, fill = value)) +
        geom_tile(color = NA) +
        facet_wrap(~source, ncol = 2) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        labs(
          title = "Variance-Covariance Matrix (V)",
          x = "Time 1", y = "Time 2", fill = "Value"
        ) +
        coord_fixed() +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_blank()
        )

      print(p_V)
    }
  }
}
