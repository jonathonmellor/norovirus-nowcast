# script to develop implimentation of `baselinenowcast` model on norovirus case
# study data.

# # # # # # # # # # # #
####    SETUP     ####
# # # # # # # # # # # #

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))
source(paste0(wd, "/scripts/run_models/functions/baselinenowcast.R"))


# we want to use the most recent versions of the packages on GitHub
remotes::install_github(repo = "epinowcast/baselinenowcast")
remotes::install_github(repo = "epinowcast/epinowcast")

library(ggplot2)


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



# Run model



baselinenowcast_model1_results <- run_scripted_model(wd = wd,
                                              model_name = "baselinenowcast",
                                              training_data = training_data,
                                              prediction_end_dates = max_reporting_dates,
                                              model_formula = "",
                                              output_columns = config$output_columns,
                                              model_hyperparams = config$hyperparams$baselinenowcast,
                                              n_pi_samples = 1000) |>
  purrr::list_rbind()|>
  dplyr::mutate(model = "baselinenowcast_model1")

baselinenowcast_formatted <- baselinenowcast_model1_results

plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))


plot_nowcast(
  data = baselinenowcast_formatted,
  training_data = training_data,
  model_name = "baselinenowcast",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = 150,
  x_limit_upper = NA,
  x_limit_lower = "2023-10-02")


# # # # # # # # # # # # #
####  SAVE OUTPUTS  ####
# # # # # # # # # # # # #

data_output_path <- glue::glue("{output_path}/data")
fs::dir_create(data_output_path)


readr::write_csv(
  x = baselinenowcast_formatted,
  file = glue::glue("{data_output_path}/baselinenowcast_model1_predictions_summary.csv"))
