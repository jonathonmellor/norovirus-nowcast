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
nowcast_date <- max_reporting_dates[[1]]

test_results <- run_baselinenowcast(.data = training_data,
                    prediction_end_date = nowcast_date,
                    n_pi_samples = 100,
                    model_hyperparams = config$hyperparams$baselinenowcast)
test_results

test_results |>
  ggplot() +
  geom_point(aes(x=specimen_date, y=target)) +
  geom_ribbon(aes(x=specimen_date, ymax=pi_95, ymin=pi_5, alpha="90%")) +
  geom_ribbon(aes(x=specimen_date, ymax=pi_75, ymin=pi_25, alpha="50%")) +
  scale_alpha_manual(values = c("90%"=0.3,
                                "50%" = 0.5)) +
  scale_y_continuous(trans = "sqrt") +
  xlab("Reference date") +
  ylab("cases") +
  theme(legend.position = "bottom")


baselinenowcast_results <- run_scripted_model(wd = wd,
                                              model_name = "baselinenowcast",
                                              training_data = training_data,
                                              prediction_end_dates = max_reporting_dates,
                                              model_formula = "",
                                              output_columns = config$output_columns,
                                              model_hyperparams = config$hyperparams$baselinenowcast,
                                              n_pi_samples = 1000) |>
  purrr::list_rbind()

baselinenowcast_formatted <- baselinenowcast_results

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

final_observed <- training_data_raw |>
  dplyr::summarise(target = sum(new_confirm, na.rm=TRUE),
                   .by="reference_date")

baselinenowcast_results |>
  dplyr::filter(specimen_date > prediction_end_date - 7) |>
  dplyr::select(-target) |>
  dplyr::right_join(final_observed, by=c("specimen_date" = "reference_date")) |>
  dplyr::mutate(prediction_end_date = factor(prediction_end_date)) |>
  dplyr::filter(specimen_date > min(max_reporting_dates) - 14) |>
  ggplot() +
  geom_point(aes(x=specimen_date, y=target)) +
  geom_line(aes(x=specimen_date, y=pi_50, group=prediction_end_date, color=prediction_end_date)) +
  geom_ribbon(aes(x=specimen_date,ymax=pi_95, ymin=pi_5, group=prediction_end_date, fill=prediction_end_date), alpha=0.5) +
  coord_cartesian(ylim=c(0,125))
