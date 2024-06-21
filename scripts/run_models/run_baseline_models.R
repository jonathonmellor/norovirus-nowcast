# Script to run baseline (compare against) models.

# SET UP ####
wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source(paste0(wd, "/scripts/depends.R"))
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))

## CONFIG ####

# reproducibility
set.seed(8675309)

config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")

# depending on if tuning or not, set dates later
tuning <- FALSE

combined_path <- "./outputs/data/cases_with_noise.csv"

output_path <- "./outputs"
data_output_path <- glue::glue("{output_path}/data")


if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}

# Load Data ####

training_data <- aws.s3::s3read_using(
  vroom::vroom,
  object = training_data_path,
  show_col_types = FALSE
) |>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date))

latest_data <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")


# PARTIAL DATA MODEL ####
# model using only the reported to date data to generate prediction

partial_baseline_daily <-
  # generate a data set per max report date
  tidyr::expand_grid(training_data,
    prediction_end_date = max_reporting_dates) |>
  # can't have tests after end date
  dplyr::filter(specimen_date <= prediction_end_date) |>
  # create data triangle
  dplyr::mutate(target =
    dplyr::if_else(specimen_date + days_to_reported > prediction_end_date,
      NA_integer_,
      target)) |>
  # aggregate to data that would be available at prediction time
  dplyr::summarise(target = sum(target, na.rm = TRUE), .by = c("prediction_end_date", "specimen_date")) |>
  # we only need predictions for the past 7 days
  dplyr::filter(specimen_date > prediction_end_date - 7) |>
  # set prediction intervals as the collected data
  dplyr::mutate(
    "pi_50" = target,
    "pi_5" = target,
    "pi_95" = target,
    "pi_2.5" = target,
    "pi_97.5" = target,
    "pi_25" = target,
    "pi_75" = target,
    "pi_17" = target,
    "pi_83" = target,
    "model" = "Baseline: partial",
    "t_aggregation" = "daily"
  )

# use the fact our predictions are just our target to aggregate up
partial_baseline_weekly <- training_data |>
  # generate a data set per max report date
  tidyr::expand_grid(prediction_end_date = max_reporting_dates) |>
  # can't have tests after end date
  dplyr::filter(specimen_date <= prediction_end_date) |>
  # create data triangle
  dplyr::mutate(target =
    dplyr::if_else(specimen_date + days_to_reported > prediction_end_date,
      NA_integer_,
      target)) |>
  dplyr::filter(specimen_date > prediction_end_date - 7) |>
  # convert specimen date to first day of week
  dplyr::mutate(specimen_date = min(specimen_date), .by = prediction_end_date) |>
  # aggregate over the week
  dplyr::summarise(
    "pi_50" = sum(target, na.rm = TRUE),
    "pi_5" = sum(target, na.rm = TRUE),
    "pi_95" = sum(target, na.rm = TRUE),
    "pi_2.5" = sum(target, na.rm = TRUE),
    "pi_97.5" = sum(target, na.rm = TRUE),
    "pi_25" = sum(target, na.rm = TRUE),
    "pi_75" = sum(target, na.rm = TRUE),
    "pi_17" = sum(target, na.rm = TRUE),
    "pi_83" = sum(target, na.rm = TRUE),
    .by = c("specimen_date", "prediction_end_date")
  ) |>
  dplyr::mutate("t_aggregation" = "weekly",
    "model" = "Baseline: partial")

partial_baseline <- dplyr::bind_rows(
  partial_baseline_daily,
  partial_baseline_weekly
) |>
  dplyr::select(-target)



# PREVIOUS WEEK MODEL ####

previous_week_baseline_daily <-
  # generate a data set per max report date
  tidyr::expand_grid(training_data,
    prediction_end_date = max_reporting_dates) |>
  # can't have tests after end date
  dplyr::filter(specimen_date <= prediction_end_date) |>
  # create data triangle
  dplyr::mutate(target =
    dplyr::if_else(specimen_date + days_to_reported > prediction_end_date,
      NA_integer_,
      target)) |>
  # filter only to the previous weeks time series
  dplyr::filter(specimen_date <= prediction_end_date - 7 &
    specimen_date > prediction_end_date - 7 * 2) |>
  # shift previous week forward to this week to generate prediction
  dplyr::mutate(specimen_date = specimen_date + 7) |>
  # aggregate to data that would be available at prediction time
  dplyr::summarise(target = sum(target, na.rm = TRUE), .by = c("prediction_end_date", "specimen_date")) |>
  # set prediction intervals as the collected data
  dplyr::mutate(
    "pi_50" = target,
    "pi_5" = target,
    "pi_95" = target,
    "pi_2.5" = target,
    "pi_97.5" = target,
    "pi_25" = target,
    "pi_75" = target,
    "pi_17" = target,
    "pi_83" = target,
    "model" = "Baseline: previous week",
    "t_aggregation" = "daily"
  )


previous_week_baseline_weekly <-
  # generate a data set per max report date
  tidyr::expand_grid(training_data,
    prediction_end_date = max_reporting_dates) |>
  # can't have tests after end date
  dplyr::filter(specimen_date <= prediction_end_date) |>
  # create data triangle
  dplyr::mutate(target =
    dplyr::if_else(specimen_date + days_to_reported > prediction_end_date,
      NA_integer_,
      target)) |>
  # filter only to the previous weeks time series
  dplyr::filter(specimen_date <= prediction_end_date - 7 &
    specimen_date > prediction_end_date - 7 * 2) |>
  # shift previous week forward to this week to generate prediction
  dplyr::mutate(specimen_date = specimen_date + 7) |>
  # convert specimen date to first day of week
  dplyr::mutate(specimen_date = min(specimen_date),
    .by = "prediction_end_date") |>
  # aggregate over the week
  dplyr::summarise(
    "pi_50" = sum(target, na.rm = TRUE),
    "pi_5" = sum(target, na.rm = TRUE),
    "pi_95" = sum(target, na.rm = TRUE),
    "pi_2.5" = sum(target, na.rm = TRUE),
    "pi_97.5" = sum(target, na.rm = TRUE),
    "pi_25" = sum(target, na.rm = TRUE),
    "pi_75" = sum(target, na.rm = TRUE),
    "pi_17" = sum(target, na.rm = TRUE),
    "pi_83" = sum(target, na.rm = TRUE),
    .by = c("specimen_date", "prediction_end_date")
  ) |>
  dplyr::mutate("t_aggregation" = "weekly",
    "model" = "Baseline: previous week")

previous_week_baseline <-
  dplyr::bind_rows(
    previous_week_baseline_daily,
    previous_week_baseline_weekly
  ) |>
  dplyr::select(-target)

# SAVE OUTPUTS ####

dir.create(data_output_path, recursive = TRUE)
write.csv(partial_baseline, glue::glue("{data_output_path}/baseline_partial_predictions_summary.csv"),
  row.names = FALSE)
write.csv(previous_week_baseline, glue::glue("{data_output_path}/baseline_prevweek_predictions_summary.csv"),
  row.names = FALSE)
