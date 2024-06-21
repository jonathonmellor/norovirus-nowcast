# tune univariate bsts model

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
set.seed(8675309)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "./scripts/run_models/functions/plotting.R"))
source(paste0(wd, "./scripts/run_models/functions/scoring.R"))

output_path <- fs::dir_create(fs::path(paste0(wd, "./outputs/tuning/bsts")))
datapath <- "./outputs/data/cases_with_noise.csv"
config <- yaml::read_yaml(paste0(wd, "/scripts/run_models/norovirus_nowcast_config.yaml"))
hyperparams <- config$hyperparams$bsts

# selecting forecast period
tuning <- TRUE
if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}

# load training data
latest_data <- readr::read_csv(datapath)

# tune, looping over s.d. values for slope and level

# initialise
scores_all <- c()

for (training_length in c(60, 120)) {
  for (sigma_slope_value in c(1.01, 1.1, 2)) {
    for (sigma_level_value in c(1.01, 1.1, 2)) {

      hyperparams$training_length <- training_length
      hyperparams$slope$upper.limit <- log(sigma_slope_value)
      hyperparams$level$upper.limit <- log(sigma_level_value)

      bsts_combined_outputs <- c()

      for (max_date in as.character(max_reporting_dates)) {
        data <- latest_data |>
          dplyr::filter(specimen_date + days_to_reported <= max_date) |>
          dplyr::summarise(target = sum(target), .by = specimen_date) |>
          dplyr::filter(specimen_date <= as.Date(max_date) - hyperparams$horizon,
            specimen_date > max(specimen_date) - hyperparams$training_length)

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

        daily_predictions <- daily_sample_predictions |>
          dplyr::group_by(.h) |>
          generate_intervals(method = "quantile") |>
          dplyr::mutate(prediction_end_date = as.Date(max_date),
            specimen_date = as.Date(prediction_end_date) - hyperparams$horizon + .h,
            t_aggregation = "daily") |>
          dplyr::select(-.h)

        weekly_predictions <- daily_sample_predictions |>
          dplyr::summarise(.value = sum(.value), .by = .sample) |>
          generate_intervals(method = "quantile") |>
          dplyr::mutate(prediction_end_date = as.Date(max_date),
            specimen_date = as.Date(max_date) - 6,
            t_aggregation = "weekly")

        # bind on predictions
        bsts_combined_outputs <- bsts_combined_outputs |>
          dplyr::bind_rows(daily_predictions,
            weekly_predictions)

      }

      bsts_combined_outputs <- bsts_combined_outputs |>
        dplyr::mutate(model = "BSTS")

      #### Plotting ####

      plot_nowcast(
        data = bsts_combined_outputs,
        training_data = latest_data,
        model_name = glue::glue("bsts_{training_length}_{sigma_slope_value}_{sigma_level_value}"),
        plot_type = "lookbacks",
        output_path = output_path,
        y_limit = NA,
        x_limit_upper = NA,
        x_limit_lower = as.Date(config$dates$start_date) - 6)


      #### Scoring ####

      bsts_formatted_scoring <- bsts_combined_outputs |>
        dplyr::filter(!is.na(pi_50))

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

      bsts_formatted_scoring <- bsts_formatted_scoring |>
        dplyr::left_join(latest_data_grouped,
          by = c("specimen_date", "t_aggregation"))


      # calculate scores
      scores_daily <- score(data = bsts_formatted_scoring |>
        dplyr::filter(t_aggregation == "daily"),
      add_log = FALSE,
      output_path = output_path,
      model_name = "bsts") |>
        dplyr::mutate(training_length = training_length,
          sigma_slope = sigma_slope_value,
          sigma_level = sigma_level_value,
          t_aggregation = "daily") |>
        dplyr::filter(scale == "natural")

      scores_all <- scores_all |>
        dplyr::bind_rows(scores_daily)

      gc()
    }
  }
}


# save all scores
write.csv(scores_all,
  file = paste0(output_path, "/scores.csv"),
  row.names = FALSE)
