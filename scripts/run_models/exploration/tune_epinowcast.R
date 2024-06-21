# Script to run epinowcast norovirus nowcast

# SETUP -------------------------------------------------------------------------------

# This script is highly compute intensive doing parallel runs.

# Recommend running this on high compute instance

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source(paste0(wd, "/scripts/depends.R"))
source(paste0(wd, "/scripts/run_models/functions/epinowcast.R"))
source(paste0(wd, "/scripts/run_models/functions/model_running_functions.R"))
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))

# SET GLOBAL SEED for reproducibility
set.seed(8675309)


# CONFIGURATION ----------------------------------------------------------------------

config <- yaml::read_yaml(paste0(wd, "/scripts/run_models/norovirus_nowcast_config.yaml"))

# depending on if tuning or not, set dates later
tuning <- TRUE

training_data_path <- "./outputs/data/cases_with_noise.csv"
output_path <- "./outputs"

scoring_output_path <- fs::dir_create(fs::path(
  output_path, "tuning", "epinowcast", "scoring")
)

data_output_path <- fs::dir_create(fs::path(
  output_path, "tuning", "epinowcast", "data")
)



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

training_data <- readr::read_csv(training_data_path)|>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date)) |>
  dplyr::select(-filename)

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

n_parallel_fits <- length(max_reporting_dates)

# Must have AT LEAST this many cores available
n_cores_required <- n_parallel_fits * config$hyperparams$epinowcast$parallel_chains + 1

if (n_cores_required > parallel::detectCores()) {
  cli::cli_abort(c(
    "x" = "At least {.val {n_cores_required}} CPU cores are required! Try reducing {.arg n_parallel_fits}?",
    " " = "{.emph This instance has {.val {parallel::detectCores()}} cores.}"
  ))
}

future::plan(future::multisession, workers = n_parallel_fits)


# Set tuning parameters

training_lengths <- c(21, 28, 35, 42, 49)
max_delays <- c(7, 14, 21)

epinowcast_models <- list()

iter <- 1
for (tl in training_lengths) {
  for (md in max_delays) {
    config$hyperparams$epinowcast$training_length <- tl
    config$hyperparams$epinowcast$max_delay <- md
    print(config$hyperparams$epinowcast)

    if (tl >= md) {
      start <- Sys.time()
      epinowcast_outputs <- run_scripted_model(wd = wd,
        model_name = "epinowcast",
        training_data = model_training_data,
        prediction_end_dates = max_reporting_dates,
        model_formula = "", # Not used for epinowcast
        output_columns = config$required_covariates,
        model_hyperparams = config$hyperparams$epinowcast,
        n_pi_samples = NA) # Not used for epinowcast

      print(Sys.time() - start)

      epinowcast_models[[iter]] <- extract_from_list(epinowcast_outputs)$models
      epinowcast_formatted <- extract_from_list(epinowcast_outputs)$quantile_predictions |>
        dplyr::rename(target_value = target) |>
        dplyr::mutate(tl = tl, md = md)

      path_name <- "{data_output_path}/epinowcast_predictions_tuning_{tl}_{md}.rds"
      saveRDS(epinowcast_formatted, glue::glue(path_name)) # save as we go
      print("Output saved to RDS")



      iter <- iter + 1
    }
  }
}


# load in that combined data
# this will break if not all are created
combined_results <- c()

for (tl in training_lengths) {
  for (md in max_delays) {

    path_name <- "{data_output_path}/epinowcast_predictions_tuning_{tl}_{md}.rds"
    combined_results <- dplyr::bind_rows(combined_results,
      readRDS(glue::glue(path_name)))

  }
}



# SCORING --------------------------------------------------------------------------


print("Now scoring models")


epinowcast_formatted_scoring <- combined_results |>
  mutate(model = paste0("tl_", tl, "_md_", md)) |>
  # score based on latest data
  dplyr::select(-target_value) |>
  dplyr::filter(!is.na(pi_50)) |>
  dplyr::distinct()

# Aggregate data by specimen_date for scoring
latest_data <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop")

epinowcast_formatted_scoring <- merge(epinowcast_formatted_scoring, latest_data, by = "specimen_date", all.x = TRUE)

scoring_epinowcast <- score(data = epinowcast_formatted_scoring |> dplyr::filter(t_aggregation == "daily"),
  add_log = FALSE,
  output_path = scoring_output_path,
  model_name = "epinowcast")

scoring_epinowcast_final <- scoring_epinowcast |>
  dplyr::mutate(params = stringr::str_remove_all(model, pattern = "(tl_)|(md_)")) |>
  tidyr::separate(params, into = c("training_length", "max_delay"))

readr::write_csv(x = scoring_epinowcast_final,
  file = glue::glue("{output_path}/epinowcast/scoring/scores.csv"))
