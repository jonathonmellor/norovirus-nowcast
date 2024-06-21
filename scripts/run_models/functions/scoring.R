# Score the last 7 days of predictions

score <- function(data, add_log = FALSE, output_path, model_name) {

  predictions_of_interest <- data |>
    # select only the predictions for the last 7 days
    dplyr::filter(specimen_date <= prediction_end_date,
      specimen_date > prediction_end_date - 7) |>
    tidyr::pivot_longer(cols = dplyr::starts_with("pi_"), names_to = "quantile", values_to = "prediction") |>
    # keep only PIs we care about
    dplyr::filter(quantile %in% c("pi_50", "pi_25", "pi_75", "pi_95", "pi_5")) |>
    dplyr::mutate(quantile = as.numeric(stringr::str_remove(quantile, "pi_")) / 100)

  predictions_of_interest_score <- predictions_of_interest |>
    dplyr::rename(true_value = target_value) |>
    scoringutils::transform_forecasts(fun = scoringutils::log_shift, offset = 1) |>
    scoringutils::score() |>
    scoringutils::add_coverage(by = c("model", "prediction_end_date", "scale"), ranges = c(50, 90)) |>
    scoringutils::summarise_scores(
      by = c("model", "prediction_end_date", "scale"),
      na.rm = TRUE
    ) |>
    scoringutils::summarise_scores(fun = round, digits = 3)

  # overall score aggregated over multiple weeks
  predictions_of_interest_score_aggregated <- predictions_of_interest |>
    dplyr::rename(true_value = target_value) |>
    scoringutils::transform_forecasts(fun = scoringutils::log_shift, offset = 1) |>
    scoringutils::score() |>
    scoringutils::add_coverage(by = c("model", "scale"), ranges = c(50, 90)) |>
    scoringutils::summarise_scores(by = c("model", "scale"), na.rm = TRUE) |>
    scoringutils::summarise_scores(fun = round, digits = 3) |>
    dplyr::mutate(time_period = "Overall") |>
    dplyr::select(model,
      scale,
      time_period,
      interval_score,
      underprediction,
      overprediction,
      ae_median,
      coverage_deviation)

  # Plot the 50% and 90% coverage, with dotted horizontal lines indicating the
  # optimal values of these scores, i.e. 0.5 and 0.9. For example, a 50%
  # coverage score of 0.5 indicates 50% of all observed values were actually
  # covered by all 50% prediction intervals. 50% coverage is plotted with a
  # solid line; 90% coverage is plotted with a dashed line.
  coverage <- predictions_of_interest_score |>
    dplyr::filter(scale == "natural") |>
    dplyr::rename(deviation = coverage_deviation, "50%" = coverage_50, "90%" = coverage_90) |>
    tidyr::pivot_longer(cols = ends_with("%"), names_to = "coverage_perc", values_to = "coverage_value") |>
    ggplot2::ggplot(ggplot2::aes(x = prediction_end_date)) +
    ggplot2::geom_line(ggplot2::aes(y = coverage_value, color = model, linetype = coverage_perc)) +
    ggplot2::geom_point(ggplot2::aes(y = coverage_value, color = model)) +
    ggplot2::geom_hline(yintercept = 0.5, color = "darkgrey", alpha = 0.7, linetype = "dotted") +
    ggplot2::geom_hline(yintercept = 0.9, color = "darkgrey", alpha = 0.7, linetype = "dotted") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_date(breaks = "1 week", date_labels = "%d-%b") +
    ggplot2::labs(
      x = "Prediction start date",
      y = "Interval coverage",
      color = "Model",
      linetype = "Coverage"
    )

  # Plot the bias, with a dotted horizontal line indicating the optimal score,
  # i.e. 0.
  bias <- predictions_of_interest_score |>
    dplyr::filter(scale == "natural") |>
    ggplot2::ggplot(ggplot2::aes(x = prediction_end_date)) +
    ggplot2::geom_line(ggplot2::aes(y = bias, color = model)) +
    ggplot2::geom_point(ggplot2::aes(y = bias, color = model)) +
    ggplot2::ylim(-1, 1) +
    ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7, linetype = "dotted") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_date(breaks = "1 week", date_labels = "%d-%b") +
    ggplot2::labs(
      x = "Prediction start date",
      y = "Bias",
      color = "Model"
    )

  # Plot the weighted interval score: the closer this is to 0, the better the
  # calibration of the model.
  is <- predictions_of_interest_score |>
    dplyr::filter(scale == "natural") |>
    ggplot2::ggplot(ggplot2::aes(x = prediction_end_date)) +
    ggplot2::geom_line(ggplot2::aes(y = interval_score, color = model)) +
    ggplot2::geom_point(ggplot2::aes(y = interval_score, color = model)) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_date(breaks = "1 week", date_labels = "%d-%b") +
    ggplot2::labs(
      x = "Prediction start date",
      y = "Interval Score",
      color = "Model"
    )

  # Plot the relative weighted interval score: the closer this is to 0, the better the
  # calibration of the model. This is the weighted interval score of log transformed targets and predictions
  # and can be interpreted as a probabalistic version of relative (as opposed to absolute) error
  # It also can be interpreted as scoring how well the model predicts growth rates

  ris <- predictions_of_interest_score |>
    dplyr::filter(scale == "log") |>
    ggplot2::ggplot(ggplot2::aes(x = prediction_end_date)) +
    ggplot2::geom_line(ggplot2::aes(y = interval_score, color = model)) +
    ggplot2::geom_point(ggplot2::aes(y = interval_score, color = model)) +
    # ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_date(breaks = "1 week", date_labels = "%d-%b") +
    ggplot2::labs(
      x = "Prediction start date",
      y = "Relative Interval Score",
      color = "Model"
    )

  # Plot the median absolute error: the closer this is to 0, the better the
  # calibration of the model.
  median_ae <- predictions_of_interest_score |>
    dplyr::filter(scale == "natural") |>
    ggplot2::ggplot(ggplot2::aes(x = prediction_end_date)) +
    ggplot2::geom_line(ggplot2::aes(y = ae_median, color = model)) +
    ggplot2::geom_point(ggplot2::aes(y = ae_median, color = model)) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_date(breaks = "1 week", date_labels = "%d-%b") +
    ggplot2::labs(
      x = "Prediction start date",
      y = "Median Absolute Error",
      color = "Model"
    )

  # Include ris plot or not depending on add_log argument
  if (add_log) {
    ggpubr::ggarrange(coverage, bias, is, ris, median_ae,
      ncol = 1,
      nrow = 5,
      common.legend = TRUE,
      legend = "bottom") |>
      ggplot2::ggsave(filename = paste0(output_path, "/", model_name, "_scoring_plots.png"), width = 16, height = 12)
  } else {
    ggpubr::ggarrange(coverage, bias, is, median_ae,
      ncol = 1,
      nrow = 4,
      common.legend = TRUE,
      legend = "bottom") |>
      ggplot2::ggsave(filename = paste0(output_path, "/", model_name, "_scoring_plots.png"), width = 16, height = 12)
  }

  score_table <- scoringutils::plot_score_table(
    predictions_of_interest_score_aggregated |> dplyr::filter(scale == "natural"),
    y = "model",
    metrics = c("interval_score", "underprediction", "overprediction", "ae_median", "coverage_deviation")) +
    ggplot2::facet_wrap(~time_period)

  ggplot2::ggsave(filename = paste0(output_path, "/", model_name, "_score_table.png"),
    plot = score_table, width = 8, height = 12)

  if (add_log) {
    score_table_log <- scoringutils::plot_score_table(
      predictions_of_interest_score_aggregated |> dplyr::filter(scale == "log"),
      y = "model",
      metrics = c("interval_score", "overprediction", "underprediction")) +
      ggplot2::facet_wrap(~time_period)

    ggplot2::ggsave(filename = paste0(output_path, "/", model_name, "_log_score_table.png"),
      plot = score_table_log, width = 8, height = 12)
  }

  return(predictions_of_interest_score_aggregated)
}
