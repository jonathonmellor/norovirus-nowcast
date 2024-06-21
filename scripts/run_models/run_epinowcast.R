# Script to run epinowcast norovirus nowcast

# SETUP -------------------------------------------------------------------------------

# NOTE: this script is highly computationally intensive and parallalised, optomised for a cloud computing environment
purrr::walk(
  c(
    here::here("/scripts/depends.R"),
    here::here("./scripts/run_models/functions", c(
      "epinowcast.R",
      "model_running_functions.R",
      "plotting.R",
      "scoring.R"
    ))
  ),
  source
)


# SET GLOBAL SEED for reproducibility
set.seed(8675309)


# CONFIGURATION ----------------------------------------------------------------------


config <- yaml::read_yaml(here::here("./scripts/run_models/norovirus_nowcast_config.yaml"))

# depending on if tuning or not, set dates later
tuning <- FALSE

training_data_path <-"./outputs/data/cases_with_noise.csv"

output_path <- here::here("outputs")

if (tuning) {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$tune_end_date),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date(config$dates$start_date),
    to = as.Date(config$dates$evaluate_end_date),
    by = 7)
}

# LOAD DATA --------------------------------------------------------------------------
box::use(box / s3)

training_data <- s3$read_using(
  training_data_path,
  \(.) vroom::vroom(., show_col_types = FALSE)
) |>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date)) |>
  dplyr::select(!filename)

# Prepare data for epinowcast

model_training_data <- training_data |>
  dplyr::mutate(report_date = specimen_date + days_to_reported) |>
  dplyr::group_by(report_date, specimen_date) |>
  dplyr::summarise(new_confirm = sum(target)) |>
  dplyr::ungroup() |>
  dplyr::mutate(dow_date = as.factor(lubridate::wday(specimen_date, week_start = 1)),
    dow_report_date = as.factor(lubridate::wday(report_date, week_start = 1)),
    weekend_reporting = ifelse((dow_report_date == 6 | dow_report_date == 7), 1, 0)) |>
  dplyr::rename(reference_date = specimen_date) |>
  # Note that epinowcast uses cumulative confirmed cases for a given reference (specimen) date by date of report
  # So need to run this function to get this from our nowcast prod
  # Also note this is different from the other nowcasting model
  epinowcast::enw_add_cumulative() |>
  # de-data.table-ify (data.table's in-place-modification paradigm might clash with our parallelisation...)
  as.data.frame() |>
  dplyr::mutate(dplyr::across(dplyr::where(\(.) inherits(., "IDate")), as.Date))


# MODEL: Epinowcast ------------------------------------------------------------------

# for a short run we want to use as many fits as we have reported dates,
# for longer runs we may not have that many cores, so select a smaller number.
# length(max_reporting_dates) # nolint commented_code_linter
n_parallel_fits <- 6

# Must have AT LEAST this many cores available
n_cores_required <- n_parallel_fits * config$hyperparams$epinowcast$parallel_chains + 1

if (n_cores_required > parallel::detectCores()) {
  cli::cli_abort(c(
    "x" = "At least {.val {n_cores_required}} CPU cores are required! Try reducing {.arg n_parallel_fits}?",
    " " = "{.emph This instance has {.val {parallel::detectCores()}} cores.}"
  ))
}

future::plan(future::multisession, workers = n_parallel_fits)

start_time <- Sys.time()
epinowcast_outputs <- run_scripted_model(
  wd = here::here(),
  model_name = "epinowcast",
  training_data = model_training_data,
  prediction_end_dates = max_reporting_dates,
  model_formula = "", # Not used for epinowcast
  output_columns = config$required_covariates,
  model_hyperparams = config$hyperparams$epinowcast,
  n_pi_samples = NA
) # Not used for epinowcast

print(Sys.time() - start_time)
epinowcast_models <- extract_from_list(epinowcast_outputs)$models
epinowcast_formatted <- extract_from_list(epinowcast_outputs)$quantile_predictions |>
  dplyr::rename(target_value = target)
epinowcast_delay_cdf_pdf <- purrr::map(epinowcast_outputs, "delay_cdf_pdf")


#  SAVE OUTPUTS --------------------------------------------------------------------

data_output_path <- glue::glue("{output_path}/data")
dir.create(data_output_path, recursive = TRUE)
epinowcast_formatted_scoring <- epinowcast_formatted |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50)) |>
  dplyr::mutate(
    specimen_date = as.Date(specimen_date, format = "%d/%m/%Y"),
    prediction_end_date = as.Date(prediction_end_date, format = "%d/%m/%Y"))

write.csv(x = epinowcast_formatted_scoring,
  file = glue::glue("{data_output_path}/epinowcast_predictions_summary.csv"),
  row.names = FALSE)


# PLOTTING -------------------------------------------------------------------------


plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))

plot_nowcast(
  data = epinowcast_formatted |> dplyr::filter(t_aggregation == "daily"),
  training_data = training_data,
  model_name = "epinowcast",
  plot_type = "daily",
  output_path = plotting_output_path,
  y_limit = 300,
  x_limit_upper = NA,
  x_limit_lower = NA)

plot_nowcast(
  data = epinowcast_formatted |> dplyr::filter(t_aggregation == "daily"),
  training_data = training_data,
  model_name = "epinowcast",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = 300,
  x_limit_upper = NA,
  x_limit_lower = NA)


# SCORING --------------------------------------------------------------------------

print("Now scoring models")

scoring_output_path <- fs::dir_create(fs::path(output_path, "scoring"))

# Aggregate data by specimen_date for scoring
latest_data <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")

epinowcast_formatted_scoring <- merge(epinowcast_formatted_scoring, latest_data, by = "specimen_date", all.x = TRUE)

score(data = epinowcast_formatted_scoring |> dplyr::filter(t_aggregation == "daily"),
  add_log = FALSE,
  output_path = scoring_output_path,
  model_name = "epinowcast")
