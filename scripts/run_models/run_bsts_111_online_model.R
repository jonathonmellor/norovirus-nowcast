### bsts with indicators

# requires 111 online data to be saved locally - produced by scripts/data_processing/save_oneoneone_online.R

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
set.seed(8675309)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))
output_path <- paste0(wd, "/outputs")
config <- yaml::read_yaml(paste0(wd, "/scripts/run_models/norovirus_nowcast_config.yaml"))
hyperparams <- config$hyperparams$bsts_111

# loading in data

datapath <- "./outputs/data/cases_with_noise.csv" # nolint: line_length_linter.
# locally saved 111 online data - produced by scripts/data_processing/save_oneoneone_online.R
nhs111_datapath <- paste0(wd, "/outputs/data/oneoneone_online_noro_symptoms.csv")

latest_data <- vroom::vroom(file = datapath)

nhs111_data <- vroom::vroom(file = nhs111_datapath) |>
  dplyr::mutate(specimen_date = as.Date(specimen_date))

regressors <- c("online_gi_symptoms",
  "online_fever",
  "online_headache",
  # "online_all_pain",
  "online_limb_pain",
  "online_stomach_pain",
  "online_not_limb_not_stomach_pain")


### setting up forecast intervals

# selecting forecast period
tuning <- FALSE
if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}


#### model running ####

bsts_combined_outputs <- c()

for (max_date in max_reporting_dates) {
  data <- latest_data |>
    dplyr::filter(specimen_date + days_to_reported <= max_date) |>
    dplyr::summarise(target = sum(target), .by = specimen_date) |>
    dplyr::filter(specimen_date > max(specimen_date) - hyperparams$training_length) |>
    # bind the explanatory variables on
    dplyr::left_join(nhs111_data,
      by = "specimen_date")

  # setting up data in train sets and test sets.
  training_data <- data |>
    dplyr::filter(specimen_date <= max_date - hyperparams$horizon)

  test_data <- data |>
    dplyr::filter(specimen_date > max_date - hyperparams$horizon) |>
    dplyr::select(-target)

  # create model
  model_components <- list()

  model_components <- bsts::AddLocalLinearTrend(
    model_components,
    # set the priors
    # fixed = false to allow the model to choose params
    slope.sigma.prior = Boom::SdPrior(sigma.guess = hyperparams$slope$sigma.guess,
      fixed = FALSE,
      upper.limit = hyperparams$slope$upper.limit),
    level.sigma.prior = Boom::SdPrior(sigma.guess = hyperparams$level$sigma.guess,
      fixed = FALSE,
      upper.limit = hyperparams$level$upper.limit),
    y = training_data$target)

  model_components <- bsts::AddSeasonal(
    model_components,
    y = training_data$target,
    nseasons  = hyperparams$seasonality)

  fit <- bsts::bsts(
    # convert character vectors into formula
    formula = str2str::v2frm(c("target", regressors)),
    model_components,
    family = "poisson",
    data = training_data,
    niter = hyperparams$niter,
    # how many regressors to be included
    expected.model.size = hyperparams$model.size)

  # plot which explanatory variables are chosen
  plot(fit, "coef", main = max_date)

  # plot model components
  components <- bsts_generate_components(.model = fit,
    .newdata = test_data)
  components |>
    dplyr::group_by(component, .t) |>
    generate_intervals(method = "quantile") |>
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x = .t, ymax = pi_95, ymin = pi_5),
      fill = "steelblue", alpha = 0.9) +
    ggplot2::geom_line(ggplot2::aes(x = .t, y = pi_50)) +
    ggplot2::facet_wrap(~component, scales = "free")

  # generate sample predictions
  daily_sample_predictions <- bsts_generate_samples(.model = fit,
    .newdata = test_data, ## if not including regressors, set to null
    horizon = hyperparams$horizon,
    burnin = hyperparams$burnin)

  daily_predictions <- daily_sample_predictions |>
    dplyr::group_by(.h) |>
    generate_intervals(method = "quantile") |>
    dplyr::mutate(prediction_end_date = as.Date(max_date),
      specimen_date = prediction_end_date - hyperparams$horizon + .h,
      t_aggregation = "daily") |>
    dplyr::select(-.h)

  weekly_predictions <- daily_sample_predictions |>
    dplyr::summarise(.value = sum(.value), .by = .sample) |>
    generate_intervals(method = "quantile") |>
    dplyr::mutate(prediction_end_date = as.Date(max_date),
      specimen_date = prediction_end_date - 6,
      t_aggregation = "weekly")

  # bind on predictions
  bsts_combined_outputs <- bsts_combined_outputs |>
    dplyr::bind_rows(daily_predictions,
      weekly_predictions)

}

bsts_combined_outputs <- bsts_combined_outputs |>
  dplyr::mutate(model = "BSTS + NHS 111 online")

#### Plotting ####

plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))

plot_nowcast(
  data = bsts_combined_outputs,
  training_data = latest_data,
  model_name = "bsts_111_online",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6)

plot_nowcast(
  data = bsts_combined_outputs,
  training_data = latest_data,
  model_name = "bsts_111_online",
  plot_type = "lookbacks_weekly",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6)


#### Save outputs ####
# in format for comparison with other models

data_output_path <- paste0(wd, "norovirus/papers/nowcast/outputs/data/")
dir.create(data_output_path, recursive = TRUE)

bsts_formatted_summary <- bsts_combined_outputs |>
  dplyr::filter(!is.na(pi_50))

readr::write_csv(
  bsts_formatted_summary,
  file =  file.path(
    paste0(data_output_path, "bsts_111_online_predictions_summary.csv")))

#### Scoring ####

scoring_output_path <- fs::dir_create(fs::path(output_path, "scoring"))

# full latest data by day and week
latest_data_grouped <- latest_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")

latest_data_grouped <- latest_data_grouped |>
  dplyr::mutate(t_aggregation = "daily") |>
  dplyr::bind_rows(latest_data_grouped |>
    dplyr::mutate(specimen_date = lubridate::floor_date(specimen_date, unit = "week", week_start = 1)) |>
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
model_name = "bsts_111_online")

score(data = bsts_formatted_scoring |>
  dplyr::filter(t_aggregation == "weekly"),
add_log = FALSE,
output_path = scoring_output_path,
model_name = "bsts_111_online")
