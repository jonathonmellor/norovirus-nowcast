# script to develop implementation of day-of-week version of
# `baselinenowcast` model on norovirus case study data.

# # # # # # # # # # # #
####    SETUP     ####
# # # # # # # # # # # #

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))
source(paste0(wd, "/scripts/run_models/functions/baselinenowcast_dow.R"))


# we want to use the most recent versions of the packages on GitHub
remotes::install_github(repo = "epinowcast/baselinenowcast")
remotes::install_github(repo = "epinowcast/epinowcast")

# SET GLOBAL SEED for reproducibility
set.seed(8675309)

# # # # # # # # # # # #
#### CONFIGURATION ####
# # # # # # # # # # # #

config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")
training_data_path <- "./outputs/data/cases_with_noise.csv"
output_path <- "./outputs"

# will only use evaluation dates
max_reporting_dates <- seq(from = as.Date(config$dates$tune_end_date) + lubridate::days(7),
                             to = as.Date(config$dates$evaluate_end_date),
                             by = 7)

# # # # # # # # # #
#### LOAD DATA ####
# # # # # # # # # #

training_data <- vroom::vroom(training_data_path)


# Specifications
nowcast_date <- max_reporting_dates[3]

test_results <- run_baselinenowcast_dow(.data = training_data,
                                    prediction_end_date = nowcast_date,
                                    n_pi_samples = 100,
                                    model_hyperparams = config$hyperparams$baselinenowcast_dow)
test_results

test_results |>
  ggplot() +
  geom_point(aes(x=specimen_date, y=target)) +
  geom_ribbon(aes(x=specimen_date, ymax=pi_95, ymin=pi_5, alpha="90%")) +
  geom_ribbon(aes(x=specimen_date, ymax=pi_75, ymin=pi_25, alpha="50%")) +
  scale_alpha_manual(values = c("90%"=0.3,
                                "50%" = 0.5)) +
  scale_y_continuous(trans = "sqrt") +
  xlab("Specimen date") +
  ylab("cases") +
  theme(legend.position = "bottom")

# fit all dates

baselinenowcast_dow_results <- run_scripted_model(wd = wd,
                                              model_name = "baselinenowcast_dow",
                                              training_data = training_data,
                                              # TODO issues fitting first two weeks
                                              prediction_end_dates = max_reporting_dates[3:23],
                                              model_formula = "",
                                              output_columns = config$output_columns,
                                              model_hyperparams = config$hyperparams$baselinenowcast_dow,
                                              n_pi_samples = 1000) |>
  purrr::list_rbind()

plotting_output_path <- fs::dir_create(fs::path(output_path, "plots"))

plot_nowcast(
  data = baselinenowcast_dow_results,
  training_data = training_data,
  model_name = "baselinenowcast_dow",
  plot_type = "lookbacks",
  output_path = plotting_output_path,
  y_limit = 150,
  x_limit_upper = NA,
  x_limit_lower = "2023-10-02")
