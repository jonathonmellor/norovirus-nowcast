# uni-variate bsts model

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
set.seed(8675309)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))

output_path <- paste0(wd, "/outputs")
datapath <- "./outputs/data/cases_with_noise.csv"
config <- yaml::read_yaml(paste0(wd, "/scripts/run_models/norovirus_nowcast_config.yaml"))
hyperparams <- config$hyperparams$bsts

tuning <- FALSE
unit_root_test <- FALSE

# selecting dates
if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}

# load data
latest_data <- aws.s3::s3read_using(
  vroom::vroom,
  object = datapath)


# testing for a unit root process with a augmented dicky fuller test.
if (unit_root_test) {
  latest_daily_sums <-  latest_data |>
    dplyr::filter(specimen_date + days_to_reported <= max(max_reporting_dates)) |>
    dplyr::summarise(target = sum(target), .by = specimen_date)

  # test on un-transformed data
  ggplot2::ggplot(latest_daily_sums) +
    ggplot2::geom_line(ggplot2::aes(specimen_date, target))
  adf_natural <- fUnitRoots::adfTest(latest_daily_sums$target, lags = 10, type = "ct")

  # test on log data
  ggplot2::ggplot(latest_daily_sums) +
    ggplot2::geom_line(ggplot2::aes(specimen_date, log(target)))
  adf_log <- fUnitRoots::adfTest(log(latest_daily_sums$target), lags = 10, type = "ct")
}


#### Run model ####

bsts_combined_outputs <- c()

for (i in as.character(max_reporting_dates)) {
  data <- latest_data |>
    dplyr::filter(specimen_date + days_to_reported <= i) |>
    dplyr::summarise(target = sum(target), .by = specimen_date) |>
    dplyr::filter(specimen_date <= as.Date(i) - hyperparams$horizon,
      specimen_date > max(specimen_date) - hyperparams$training_length)

  model_components <- list()

  model_components <- bsts::AddLocalLinearTrend(
    model_components,
    # set the priors
    # fixed = false to allow the model to choose params
    slope.sigma.prior = Boom::SdPrior(
      sigma.guess = hyperparams$slope$sigma.guess,
      fixed = FALSE,
      upper.limit = hyperparams$slope$upper.limit),
    level.sigma.prior = Boom::SdPrior(
      sigma.guess = hyperparams$level$sigma.guess,
      fixed = FALSE,
      upper.limit = hyperparams$level$upper.limit),
    y = data$target)

  model_components <- bsts::AddSeasonal(
    model_components,
    y = data$target,
    nseasons  = hyperparams$seasonality)

  fit <- bsts::bsts(formula = data$target,
    family = "poisson",
    model_components,
    data = data,
    niter = hyperparams$niter)

  daily_sample_predictions <- bsts_generate_samples(.model = fit,
    .newdata = NULL,
    horizon = hyperparams$horizon,
    burnin = hyperparams$burnin)

  components <- bsts_generate_components(.model = fit,
    .newdata = data)
  components |>
    dplyr::group_by(component, .t) |>
    generate_intervals(method = "quantile") |>
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x = .t, ymax = pi_95, ymin = pi_5),
      fill = "steelblue", alpha = 0.9) +
    ggplot2::geom_line(ggplot2::aes(x = .t, y = pi_50)) +
    ggplot2::facet_wrap(~component, scales = "free")

  daily_predictions <- daily_sample_predictions |>
    dplyr::group_by(.h) |>
    generate_intervals(method = "quantile") |>
    dplyr::mutate(prediction_end_date = as.Date(i),
      specimen_date = as.Date(prediction_end_date) - hyperparams$horizon + .h,
      t_aggregation = "daily") |>
    dplyr::select(-.h)

  weekly_predictions <- daily_sample_predictions |>
    dplyr::summarise(.value = sum(.value), .by = .sample) |>
    generate_intervals(method = "quantile") |>
    dplyr::mutate(prediction_end_date = as.Date(i),
      specimen_date = as.Date(i) - 6,
      t_aggregation = "weekly")

  # bind on predictions
  bsts_combined_outputs <- bsts_combined_outputs |>
    dplyr::bind_rows(daily_predictions,
      weekly_predictions)

}

bsts_combined_outputs <- bsts_combined_outputs |>
  dplyr::mutate(model = "BSTS")

#### Plotting ####

plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))

plot_nowcast(
  data = bsts_combined_outputs,
  training_data = latest_data,
  model_name = "bsts",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6)

plot_nowcast(
  data = bsts_combined_outputs,
  training_data = latest_data,
  model_name = "bsts",
  plot_type = "lookbacks_weekly",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6)


#### Save outputs ####
# in format for comparison with other models

data_output_path <- paste0(wd, "/outputs/data/")
dir.create(data_output_path, recursive = TRUE)

bsts_formatted_summary <- bsts_combined_outputs |>
  dplyr::filter(!is.na(pi_50))

readr::write_csv(bsts_formatted_summary, file =  file.path(paste0(data_output_path, "bsts_predictions_summary.csv")))

#### Scoring ####

scoring_output_path <- fs::dir_create(fs::path(output_path, "scoring"))

# full latest data by day and week
latest_data_grouped <- latest_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(t_aggregation = "daily") |>
  dplyr::bind_rows(latest_data_grouped |>
    dplyr::mutate(
      specimen_date = lubridate::floor_date(specimen_date,
        unit = "week",
        week_start = 1)
    ) |>
    dplyr::group_by(specimen_date) |>
    dplyr::summarise(target_value = sum(target_value, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(t_aggregation = "weekly"))

bsts_formatted_scoring <- bsts_formatted_summary |>
  dplyr::left_join(latest_data_grouped,
    by = c("specimen_date", "t_aggregation"))


score(data = bsts_formatted_scoring |>
  dplyr::filter(t_aggregation == "daily"),
add_log = FALSE,
output_path = scoring_output_path,
model_name = "bsts")

score(data = bsts_formatted_scoring |>
  dplyr::filter(t_aggregation == "weekly"),
add_log = FALSE,
output_path = scoring_output_path,
model_name = "bsts_weekly")
