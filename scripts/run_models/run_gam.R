# Script to run GAM norovirus nowcast

# # # # # # # # # # # #
####    SETUP     ####
# # # # # # # # # # # #

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/gam.R"))
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))

# SET GLOBAL SEED for reproducibility
set.seed(8675309)


# # # # # # # # # # # #
#### CONFIGURATION ####
# # # # # # # # # # # #

config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")

# depending on if tuning or not, set dates later
tuning <- FALSE

training_data_path <- "./outputs/data/cases_with_noise.csv"
output_path <- "./outputs"

if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}


# # # # # # # # # #
#### LOAD DATA ####
# # # # # # # # # #

training_data <- aws.s3::s3read_using(
  vroom::vroom,
  object = training_data_path,
  show_col_types = FALSE
) |>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date)) |>
  dplyr::select(-filename)


# # # # # # # # # # # # # # #
# ####  MODEL: GAM       ####
# # # # # # # # # # # # # # #

model_formula <-
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "cr") +
  s(dow_report_date_factor, bs = "re")'

gam_outputs <- run_scripted_model(wd = wd,
  model_name = "gam",
  training_data = training_data,
  prediction_end_dates = max_reporting_dates,
  model_formula = model_formula,
  output_columns = config$output_columns,
  model_hyperparams = config$hyperparams$gam,
  n_pi_samples = 1000)

gam_models <- extract_from_list(gam_outputs)$models
gam_formatted <- extract_from_list(gam_outputs)$quantile_predictions


# # # # # # # # # # # # #
####  SAVE OUTPUTS  ####
# # # # # # # # # # # # #

data_output_path <- glue::glue("{output_path}/data")
dir.create(data_output_path, recursive = TRUE)

gam_formatted_scoring <- gam_formatted |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50))
write.csv(
  x = gam_formatted_scoring,
  file = glue::glue("{data_output_path}/gam_predictions_summary.csv"),
  row.names = FALSE)


# # # # # # # # # #
#### PLOTTING ####
# # # # # # # # #

message("Now plotting models")
plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))

plot_nowcast(
  data = gam_formatted,
  training_data = training_data,
  model_name = "gam",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = "2023-10-02")

plot_nowcast(
  data = gam_formatted,
  training_data = training_data,
  model_name = "gam",
  plot_type = "lookbacks_weekly",
  output_path = plotting_output_path,
  y_limit = NA,
  x_limit_upper = NA,
  x_limit_lower = "2023-10-02")


# # # # # # # # # #
#### SCORING ####
# # # # # # # # #
message("Now scoring models")

scoring_output_path <- fs::dir_create(fs::path(output_path, "scoring"))

# full latest data by day and week
latest_data <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")

latest_data <- latest_data |>
  dplyr::mutate(t_aggregation = "daily") |>
  dplyr::bind_rows(latest_data |>
    dplyr::mutate(
      specimen_date = lubridate::floor_date(
        specimen_date, unit = "week", week_start = 1)) |>
    dplyr::group_by(specimen_date) |>
    dplyr::summarise(target_value = sum(target_value, na.rm = TRUE),
      .groups = "drop") |>
    dplyr::mutate(t_aggregation = "weekly"))

gam_formatted_scoring <- gam_formatted_scoring |>
  dplyr::left_join(
    latest_data,
    by = c("specimen_date", "t_aggregation"))


score(
  data = gam_formatted_scoring |>
    dplyr::filter(t_aggregation == "daily"),
  add_log = FALSE,
  output_path = scoring_output_path,
  model_name = "gam")

score(
  data = gam_formatted_scoring |>
    dplyr::filter(t_aggregation == "weekly"),
  add_log = FALSE,
  output_path = scoring_output_path,
  model_name = "gam_weekly")
