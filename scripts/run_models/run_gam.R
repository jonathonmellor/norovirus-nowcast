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



training_data <- vroom::vroom(training_data_path)


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
fs::dir_create(data_output_path)

gam_formatted_scoring <- gam_formatted |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50))
readr::write_csv(
  x = gam_formatted_scoring,
  file = glue::glue("{data_output_path}/gam_predictions_summary.csv"))


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

# NOTE: scoring has been removed.
