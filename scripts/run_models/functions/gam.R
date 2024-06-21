run_gam <- function(.data,
                    model_formula,
                    prediction_end_date,
                    n_pi_samples,
                    output_columns,
                    model_hyperparams,
                    model_description

) {
  prediction_end_date <- as.Date(prediction_end_date)

  known_data <- .data |>
    dplyr::filter(specimen_date + days_to_reported <= prediction_end_date)

  # format data for model
  reporting_triangle <- .data |>
    # filter to training length and max delay
    dplyr::filter(specimen_date <= prediction_end_date,
      specimen_date > prediction_end_date - model_hyperparams$training_length) |>
    dplyr::filter(days_to_reported <= model_hyperparams$max_delay) |>
    dplyr::mutate(
      origin = as.numeric(specimen_date - min(specimen_date)) + 1, # +1 to make model happy
      # add NAs for specimen date + days reported after the max date, i.e. in the future
      target = ifelse(specimen_date + days_to_reported > prediction_end_date, NA, target),
      dow_specimen_date = lubridate::wday(specimen_date, week_start = 1),
      dow_report_date = lubridate::wday(specimen_date + days_to_reported, week_start = 1),
      weekend_reporting = ifelse((dow_report_date == 6 | dow_report_date == 7), 1, 0),
      dow_specimen_date_factor = as.factor(dow_specimen_date),
      dow_report_date_factor = as.factor(dow_report_date),
      weekend_reporting_factor = as.factor(weekend_reporting),
    )

  # known data to fit model to
  reporting_triangle_fit <- reporting_triangle |>
    dplyr::filter(!is.na(target))

  # number of knots
  k_specimen <- round(model_hyperparams$training_length / model_hyperparams$denom_specimen)
  k_report <- round(model_hyperparams$training_length / model_hyperparams$denom_report)

  ## Fit the GAM
  gam_fit <- mgcv::gam(
    as.formula(glue::glue(model_formula)),
    data = reporting_triangle_fit,
    family = "nb"
  )

  # Generate samples from model fit coefficients for unknown data
  reporting_triangle_out <- reporting_triangle |>
    dplyr::filter(is.na(target))

  fits_out <- generate_samples(
    .model = gam_fit,
    .data = reporting_triangle_out,
    .n_pi_samples = n_pi_samples,
    .method = "mh"
  )

  # tidy samples
  output_data_samples <- fits_out |>
    dplyr::mutate(model = "GAM") |>
    dplyr::select(dplyr::any_of(output_columns)) |>
    dplyr::mutate(week_starting = lubridate::floor_date(specimen_date, unit = "week", week_start = 1))

  # sum into daily total positive tests for each sample and add target
  daily_sample_predictions <- output_data_samples |>
    dplyr::group_by(specimen_date, week_starting, model, .sample) |>
    dplyr::summarise(.value = sum(.value), .groups = "drop") |>
    dplyr::left_join(known_data |>
      dplyr::summarise(target = sum(target, na.rm = TRUE), .by = specimen_date),
    by = c("specimen_date"))

  # get quantiles
  daily_preds <- samples_to_quantiles(
    .sample_predictions = daily_sample_predictions,
    remove_identifiers = c("week_starting")) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("pi_"), ~ .x + target
      )
    ) |>
    # add on additional training data within training length
    dplyr::bind_rows(known_data |>
      dplyr::filter(specimen_date > prediction_end_date - model_hyperparams$training_length,
        specimen_date <= prediction_end_date - model_hyperparams$max_delay) |>
      dplyr::summarise(target = sum(target, na.rm = TRUE), .by = specimen_date)) |>
    dplyr::mutate(t_aggregation = "daily")

  # sum into weekly total positive tests for each sample and add target
  weekly_sample_predictions <- daily_sample_predictions |>
    dplyr::group_by(week_starting, model, .sample) |>
    dplyr::summarise(.value = sum(.value),
      target = sum(target),
      .groups = "drop") |>
    dplyr::rename(specimen_date = week_starting)

  weekly_preds <- samples_to_quantiles(
    .sample_predictions = weekly_sample_predictions,
    remove_identifiers = c()) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("pi_"), ~ .x + target
      )
    ) |>
    # add on additional training data within training length
    dplyr::bind_rows(known_data |>
      dplyr::mutate(week_starting = lubridate::floor_date(specimen_date, unit = "week", week_start = 1)) |>
      dplyr::filter(specimen_date > prediction_end_date - model_hyperparams$training_length,
        specimen_date <= prediction_end_date - model_hyperparams$max_delay) |>
      dplyr::summarise(target = sum(target, na.rm = TRUE), .by = week_starting) |>
      dplyr::rename(specimen_date = week_starting)) |>
    dplyr::mutate(t_aggregation = "weekly")

  # aggregate samples by reporting delay

  # sum into daily total positive tests for each sample and add target
  reporting_delay_sample_predictions <- output_data_samples |>
    dplyr::left_join(known_data,
      by = c("specimen_date", "days_to_reported"))

  reporting_delay_preds <- samples_to_quantiles(
    .sample_predictions = reporting_delay_sample_predictions,
    remove_identifiers = c("specimen_date", "week_starting"))

  training_data <- .data |>
    # filter to training length and max delay
    dplyr::filter(specimen_date <= prediction_end_date,
      specimen_date > prediction_end_date - model_hyperparams$training_length) |>
    dplyr::filter(days_to_reported <= model_hyperparams$max_delay)

  known_data_quantiles <- training_data |>
    # filter to known
    dplyr::filter(specimen_date + days_to_reported <= prediction_end_date) |>
    dplyr::summarise("pi" = list(
      c(0.5, 0.05, 0.95, 0.025, 0.975, 0.25, 0.75, 0.17, 0.83) |>
        stats::quantile(target, probs = _, na.rm = TRUE) |>
        purrr::set_names(\(x) paste0("pi_", x) |> sub("%$", "", x = _))
    ),
    .by = c("days_to_reported")) |>
    tidyr::unnest_wider(pi)

  unknown_data_quantiles <- training_data |>
    # filter to unknown
    dplyr::filter(specimen_date + days_to_reported > prediction_end_date) |>
    dplyr::summarise("pi" = list(
      c(0.5, 0.05, 0.95, 0.025, 0.975, 0.25, 0.75, 0.17, 0.83) |>
        stats::quantile(target, probs = _, na.rm = TRUE) |>
        purrr::set_names(\(x) paste0("pi_", x) |> sub("%$", "", x = _))
    ),
    .by = c("days_to_reported")) |>
    tidyr::unnest_wider(pi)

  reporting_delay_preds |>
    ggplot2::ggplot(ggplot2::aes(x = days_to_reported)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pi_2.5, ymax = pi_97.5, fill = "unknown fit"),
      alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = pi_50, color = "unknown fit")) +
    ggplot2::geom_ribbon(data = known_data_quantiles,
      ggplot2::aes(ymin = pi_5, ymax = pi_95, fill = "known data"),
      alpha = 0.3) +
    ggplot2::geom_line(data = known_data_quantiles,
      ggplot2::aes(y = pi_50, color = "known data")) +
    ggplot2::geom_ribbon(data = unknown_data_quantiles,
      ggplot2::aes(ymin = pi_5, ymax = pi_95, fill = "unknown data"),
      alpha = 0.3) +
    ggplot2::geom_line(data = unknown_data_quantiles,
      ggplot2::aes(y = pi_50, color = "unknown data")) +
    ggplot2::scale_colour_manual(values = c("known data" = "blue", "unknown data" = "red", "unknown fit" = "black")) +
    ggplot2::scale_fill_manual(values = c("known data" = "blue", "unknown data" = "red", "unknown fit" = "black")) +
    ggplot2::coord_cartesian(ylim = c(0, 10)) +
    ggplot2::labs(y = "count of tests",
      color = "data",
      fill = "data")

  ggplot2::ggsave(path = "./outputs/plots/", filename = "days_to_reported_fitted_dist.png")


  gam_preds <- dplyr::bind_rows(
    daily_preds,
    weekly_preds
  ) |>
    dplyr::mutate(model = paste0("GAM", model_description),
      prediction_end_date = prediction_end_date) |>
    dplyr::rename(target_value = target)

  return(
    list(
      quantile_predictions = gam_preds,
      model = gam_fit
    )
  )

}
