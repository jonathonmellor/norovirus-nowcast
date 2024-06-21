# Script to run and score different variations on the GAM norovirus nowcast

# # # # # # # # # # # #
####    SETUP     ####
# # # # # # # # # # # #

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")
source("./run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/gam.R"))
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))

# SET GLOBAL SEED for reproducibility
set.seed(8675309)


# # # # # # # # # # # #
#### CONFIGURATION ####
# # # # # # # # # # # #

config <- yaml::read_yaml(paste0(wd, "/scripts/run_models/norovirus_nowcast_config.yaml"))

# depending on if tuning or not, set dates later
tuning <- TRUE

training_data_path <- "./outputs/data/cases_with_noise.csv"

output_path <- fs::dir_create(fs::path(paste0(wd, "norovirus/papers/nowcast/outputs/tuning/gam")))

if (tuning) {
  max_reporting_dates <- seq(from = as.Date("2023-10-08"),
    to = as.Date("2023-10-29"),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date("2023-11-05"),
    to = as.Date("2024-03-10"),
    by = 7)
}


# # # # # # # # # #
#### LOAD DATA ####
# # # # # # # # # #

training_data <- readr::read_csv(training_data_path) |>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date))

latest_data <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")


# # # # # # # # # # # # # # # # # # # #
#### TRY DIFFERENT MODEL STRUCTURES  ####
# # # # # # # # # # # # # # # # # # # #
# for testing different model structures in the lists formula_names and formula_values below

# Set parameter values

config$hyperparams$gam$training_length <- 56
config$hyperparams$gam$max_delay <- 56
# knots:
# higher denom -> fewer knots -> more uncertain at the end
config$hyperparams$gam$denom_specimen <- 7
config$hyperparams$gam$denom_report <- 14


k_specimen <- round(config$hyperparams$gam$training_length / config$hyperparams$gam$denom_specimen)
k_report <- round(config$hyperparams$gam$max_delay / config$hyperparams$gam$denom_report)

# try different dow effect structures
# try different spline types: tp or cr or gp
# m penalisation parameter, can try more knots and higher penalisation at 2 or 3


### Model structures to try:

formula_names <- c(
  "specimen-date_num-re_report-date_num-re_splines_cr",
  "specimen-date_num-re_report-date_num-re_splines_tp_m",
  "specimen-date_num-re_report-date_num-re_splines_tp",
  "specimen-date_num-re_report-date_num-re_splines_gp"
  # "specimen-date_num-re_report-date_num by weekend binary",
  # "specimen-date_num-re_report-date_num-cyclic",
  # "specimen-date_num-cyclic_report-date_num-re",
  # "specimen-date_num-cyclic_report-date_num-cyclic"
)

formula_values <- c(
  # specimen date: num + re, report date: num + re
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "cr") +
  s(dow_report_date_factor, bs = "re")',

  'target ~
  s(as.numeric(origin), k = 7, bs = "tp") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = 5, bs = "tp", m = 3) +
  s(dow_report_date_factor, bs = "re")',

  # specimen date: num + re, report date: num + re
  'target ~
  s(as.numeric(origin), k = 7, bs = "tp") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = 2, bs = "tp") +
  s(dow_report_date_factor, bs = "re")',

  # specimen date: num + re, report date: num + re; gp splines
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "gp") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "gp") +
  s(dow_report_date_factor, bs = "re")',

  # specimen date: num + re, report date: num by weekend binary
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), by = weekend_reporting_factor, k = {k_report}, bs = "cr")',

  # specimen date: num + re, report date: num + cyclic
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "cr") +
  s(dow_report_date, k = 7, bs="cc")',

  # specimen date: num + cyclic, report date: num + re
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date, k = 7, bs="cc") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "cr") +
  s(dow_report_date_factor, bs="re")',

  # specimen date: num + cyclic, report date: num + cyclic
  'target ~
  s(as.numeric(origin), k = {k_specimen}, bs = "cr") +
  s(dow_specimen_date, k = 7, bs="cc") +
  s(as.numeric(days_to_reported), k = {k_report}, bs = "cr") +
  s(dow_report_date, k = 7, bs="cc")'
)


### Run and plot models

# for each formula and associated model name, saves into a tibble:
# model_output: dataframe of quantile predictions for all prediction weeks
# plot: single plot of quantile predictions for all prediction weeks
#
# note that purrr::map2 runs this for each formula sequentially, but within
# the .f function, run_scripted_model runs each prediction week (as given by
# max_reporting_dates) in parallel

model_variations <- tibble::tibble(
  formula = formula_values,
  model_name = formula_names
) |>
  dplyr::mutate(
    model_output = purrr::map2(
      .x = formula,
      .y = model_name,
      .f = ~ extract_from_list(
        run_scripted_model(wd = wd,
          model_name = "gam",
          training_data = training_data,
          prediction_end_dates = max_reporting_dates,
          model_formula = glue::glue(.x),
          output_columns = config$output_columns,
          model_hyperparams = config$hyperparams$gam,
          n_pi_samples = 1000,
          model_description = .y)
      )$quantile_predictions,
      .progress = "model_output"),
    plot = purrr::map2(
      .x = model_output,
      .y = model_name,
      .f = ~ plot_nowcast(
        data = .x,
        training_data = training_data,
        model_name = .y,
        plot_type = "lookbacks",
        output_path = output_path,
        y_limit = 45,
        x_limit_upper = NA,
        x_limit_lower = "2023-10-02") +
        ggplot2::ggtitle(.y)
    )
  )


### Score models
print("Now scoring models")

model_variations_scoring <- model_variations$model_output |>
  dplyr::bind_rows() |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50)) |>
  dplyr::left_join(latest_data,
    by = "specimen_date")

scores_daily <- score(data = model_variations_scoring |>
  dplyr::filter(t_aggregation == "daily"),
add_log = FALSE,
output_path = output_path,
model_name = formula_names[choose_formula])

scores_weekly <- score(data = model_variations_scoring |>
  dplyr::filter(t_aggregation == "weekly"),
add_log = FALSE,
output_path = output_path,
model_name = formula_names[choose_formula])

scores_all <- scores_daily |>
  dplyr::mutate(t_aggregation = "daily") |>
  dplyr::bind_rows(
    scores_weekly |>
      dplyr::mutate(t_aggregation = "weekly")
  ) |>
  dplyr::mutate(
    training_length = stringr::str_match(model, "tl_(.*?)_")[, 2],
    max_delay = stringr::str_match(model, "delay_(.*?)_")[, 2],
    denom_specimen = stringr::str_match(model, "denom_spec_(.*?)_")[, 2],
    denom_report = stringr::str_match(model, "denom_rep_(.*)")[, 2]
  )


write.csv(scores_all,
  file = paste0(output_path, "/scores.csv"),
  row.names = FALSE)




# # # # # # # # # # # # # # # # # # # #
#### TRY DIFFERENT PARAMETER VALUES  ####
# # # # # # # # # # # # # # # # # # # #

# for testing parameter values for one model, model is chosen by choose_formula parameter
choose_formula <- 1

# parameter values to test
training_length_values <- c(28, 42, 56, 70)
max_delay_values <- c(14, 21, 28)
denom_specimen_values <- c(7, 14)
denom_report_values <- c(3, 7, 14)


### Choice of model formulas

# formula name to name output folder
formula_names <- c(
  "specimen-date_num-re_report-date_num-re",
  "specimen-date_num-re_report-date_num by weekend binary",
  "specimen-date_num-re_report-date_num-cyclic",
  "specimen-date_num-cyclic_report-date_num-re",
  "specimen-date_num-cyclic_report-date_num-cyclic"
)

# model description for labelling models with different parameters in plots/scoring
model_descriptions <- rep(
  "tl_{training_length}_delay_{max_delay}_denom_spec_{denom_specimen}_denom_rep_{denom_report}",
  5
)

formula_values <- c(
  # specimen date: num + re, report date: num + re
  'target ~
  s(as.numeric(origin), k = {round(training_length/denom_specimen)}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {round(max_delay/denom_report)}, bs = "cr") +
  s(dow_report_date_factor, bs = "re")',

  # specimen date: num + re, report date: num by weekend binary
  'target ~
  s(as.numeric(origin), k = {round(training_length/denom_specimen)}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), by = weekend_reporting_factor, k = {round(max_delay/denom_report)}, bs = "cr")',

  # specimen date: num + re, report date: num + cyclic
  'target ~
  s(as.numeric(origin), k = {round(training_length/denom_specimen)}, bs = "cr") +
  s(dow_specimen_date_factor, bs = "re") +
  s(as.numeric(days_to_reported), k = {round(max_delay/denom_report)}, bs = "cr") +
  s(dow_report_date, k = 7, bs="cc")',

  # specimen date: num + cyclic, report date: num + re
  'target ~
  s(as.numeric(origin), k = {round(training_length/denom_specimen)}, bs = "cr") +
  s(dow_specimen_date, k = 7, bs="cc") +
  s(as.numeric(days_to_reported), k = {round(max_delay/denom_report)}, bs = "cr") +
  s(dow_report_date_factor, bs="re")',

  # specimen date: num + cyclic, report date: num + cyclic
  'target ~
  s(as.numeric(origin), k = {round(training_length/denom_specimen)}, bs = "cr") +
  s(dow_specimen_date, k = 7, bs="cc") +
  s(as.numeric(days_to_reported), k = 4, bs = "cr") +
  s(dow_report_date, k = 7, bs="cc")'
)


### Run and plot models

output_path_parameter_tune <- fs::dir_create(fs::path(output_path, formula_names[choose_formula]))

# create tibble of model parameters
model_parameters <- expand.grid(training_length = training_length_values,
  max_delay = max_delay_values,
  denom_specimen = denom_specimen_values,
  denom_report = denom_report_values) |>
  tibble::tibble() |>
  dplyr::mutate(
    model_path = "norovirus/papers/nowcast/scripts/run_models/functions/gam.R",
    formula = glue::glue(formula_values[choose_formula]),
    model_description = glue::glue(model_descriptions[choose_formula]),
    model_name = formula_names[choose_formula]
  ) |>
  # this takes parameter values and nests them within a list called hyperparams for each row
  # this is the format that the model function used later within purrr::pmap takes parameters as input
  tidyr::nest(hyperparams = -c(formula, model_description, model_name)) |>
  dplyr::select(hyperparams, formula, model_description, model_name)


# given model formula, for each set of hyperparameters, saves into a tibble:
# model_output: dataframe of quantile predictions for all prediction weeks
# plot: single plot of quantile predictions for all prediction weeks
#
# note that purrr::pmap runs this for each formula sequentially, but within
# the .f function, run_scripted_model runs each prediction week (as given by
# max_reporting_dates) in parallel

parameter_variations <- model_parameters |>
  dplyr::mutate(
    model_output = purrr::pmap(
      .l = list(
        formula,
        hyperparams,
        model_description
      ),
      .f = \(x, y, z) {
        extract_from_list(
          run_scripted_model(wd = wd,
            model_name = "gam",
            training_data = training_data,
            prediction_end_dates = max_reporting_dates,
            model_formula = x,
            output_columns = config$output_columns,
            model_hyperparams = y,
            n_pi_samples = 1000,
            model_description = z)
        )$quantile_predictions
      }
    ),

    plot = purrr::map2(
      .x = model_output,
      .y = model_description,
      .f = ~ plot_nowcast(
        data = .x,
        training_data = training_data,
        model_name = .y,
        plot_type = "lookbacks",
        output_path = output_path_parameter_tune,
        y_limit = NA,
        x_limit_upper = NA,
        x_limit_lower = "2023-10-02") +
        ggplot2::ggtitle(.y)
    )
  )

print("Now scoring models")

parameter_variations_scoring <- parameter_variations$model_output |>
  dplyr::bind_rows() |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50)) |>
  dplyr::left_join(latest_data,
    by = "specimen_date")

scores_daily <- score(data = parameter_variations_scoring |>
  dplyr::filter(t_aggregation == "daily"),
add_log = FALSE,
output_path = output_path_parameter_tune,
model_name = formula_names[choose_formula])

scores_weekly <- score(data = parameter_variations_scoring |>
  dplyr::filter(t_aggregation == "weekly"),
add_log = FALSE,
output_path = output_path_parameter_tune,
model_name = formula_names[choose_formula])

scores_all <- scores_daily |>
  dplyr::mutate(t_aggregation = "daily") |>
  dplyr::bind_rows(
    scores_weekly |>
      dplyr::mutate(t_aggregation = "weekly")
  ) |>
  dplyr::mutate(
    training_length = stringr::str_match(model, "tl_(.*?)_")[, 2],
    max_delay = stringr::str_match(model, "delay_(.*?)_")[, 2],
    denom_specimen = stringr::str_match(model, "denom_spec_(.*?)_")[, 2],
    denom_report = stringr::str_match(model, "denom_rep_(.*)")[, 2]
  )


write.csv(scores_all,
  file = paste0(output_path_parameter_tune, "/scores.csv"),
  row.names = FALSE)
