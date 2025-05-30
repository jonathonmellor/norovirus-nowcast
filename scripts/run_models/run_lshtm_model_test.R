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



training_data_raw <- vroom::vroom(training_data_path) |>
  # convert to epinowcast naming conventions
  dplyr::rename(reference_date = specimen_date,
                confirm = target) |>
  dplyr::mutate(report_date = reference_date + days_to_reported) |>
  dplyr::select(-days_to_reported) |>
  # adding this to try handle the problem calculating pobs
  dplyr::filter(!is.na(confirm))



# Run model

# Approach: transform synthetic data to format needed for package and run getting
# started page code in order.

# select a test date to get the code working in line with the package getting started page
nowcast_date <- max_reporting_dates[[1]]

# apply required filtering
target_data <- training_data_raw |>
  dplyr::rename(new_confirm = confirm) |>
  epinowcast::enw_filter_report_dates(latest_date = nowcast_date + 30) |>
  epinowcast::enw_filter_reference_dates(
    latest_date = nowcast_date
  ) |>
  epinowcast::enw_add_cumulative()

latest_data <- epinowcast::enw_latest_data(target_data)

observed_data <- epinowcast::enw_filter_report_dates(
  target_data,
  latest_date = nowcast_date
)

obs_data_by_reference_date <- epinowcast::enw_latest_data(observed_data)

ggplot() +
  geom_line(
    data = obs_data_by_reference_date,
    aes(x = reference_date, y = confirm), color = "darkred"
  ) +
  geom_line(
    data = latest_data,
    aes(x = reference_date, y = confirm), color = "black"
  ) +
  theme_bw() +
  xlab("Reference date") +
  ylab("Confirmed admissions") +
  #scale_y_continuous(trans = "log10") +
  ggtitle("Comparing real-time and later observed cases")

# Specify the maximum delay, which will determine the length of your delay
# distribution. Empirical data outside this delay window will not be used for
# training.
max_delay <- config$hyperparams$gam$max_delay
n_training_volume <- config$hyperparams$gam$training_length

# Specify the number of reference times to use to estimate the delay
# distribution. Note this assumes you want the most recent observations.
# NOTE: check which to chose for this??
n_history_delay <- 0.5 * n_training_volume

# Specify the number of retrospective nowcast datasets
# to use for uncertainty estimation.
# NOTE: check which to chose for this??
n_retrospective_nowcasts <- 0.5 * n_training_volume

training_data <- epinowcast::enw_filter_reference_dates(
  observed_data,
  include_days = n_training_volume - 1
)

latest_training_data <- epinowcast::enw_latest_data(training_data)

target_data <- epinowcast::enw_filter_reference_dates(
  latest_data,
  include_days = n_training_volume - 1
)


# Get the reporting triangle, adding an additional day because epinowcast
# we want the max_delay + 1 entries since 0 is a valid delay.
# This also validates that the data is in the correct format and
# runs preprocessing see ?enw_preprocess_data for more details
pobs <- epinowcast::enw_preprocess_data(
  obs = training_data,
  max_delay = max_delay + 1
)

# as we only have one group here we only need reference_date, delay,
# and new_confirm
reporting_triangle_df <- select(
  pobs$new_confirm[[1]],
  reference_date,
  delay,
  new_confirm
)

# we now pivot to wide format, dropping the reference_date column, and
# convert to a matrix
# this is the format that baselinenowcast expects
reporting_triangle <- reporting_triangle_df |>
  pivot_wider(names_from = delay, values_from = new_confirm) |>
  select(-reference_date) |>
  as.matrix()

triangle_df <- as.data.frame(reporting_triangle) |>
  mutate(time = row_number()) |>
  pivot_longer(!time,
               values_to = "count",
               names_prefix = "V",
               names_to = "delay"
  ) |>
  mutate(delay = as.numeric(delay))

ggplot(
  triangle_df,
  aes(x = delay, y = time, fill = count)
) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Reporting triangle", x = "Delay", y = "Time") +
  theme_bw() +
  scale_y_reverse()

delay_pmf <- baselinenowcast::get_delay_estimate(
  reporting_triangle = reporting_triangle,
  max_delay = max_delay,
  n = n_history_delay
)

delay_df <- data.frame(
  delay = 0:(length(delay_pmf) - 1),
  pmf = delay_pmf
)

ggplot(delay_df) +
  geom_line(aes(x = delay, y = cumsum(pmf))) +
  xlab("Delay") +
  ylab("Cumulative proportion reported") +
  ggtitle("Empirical point estimate of cumulative proportion reported by delay") + # nolint
  theme_bw()

ggplot(delay_df) +
  geom_point(aes(x = delay, y = pmf)) +
  xlab("Delay") +
  ylab("Proportion reported") +
  ggtitle("Empirical point estimate of proportion reported by delay") +
  theme_bw()

point_nowcast_matrix <- baselinenowcast::apply_delay(
  rep_tri_to_nowcast = reporting_triangle,
  delay_pmf = delay_pmf
)

point_nowcast_df <- target_data |>
  mutate(nowcast = rowSums(point_nowcast_matrix))

prep_latest_data <- latest_training_data |>
  mutate(type = "Real-time data") |>
  select(type, reference_date, count = confirm)

# Combine data into a single dataframe for plotting
plot_data <- point_nowcast_df |>
  pivot_longer(
    cols = c(confirm, nowcast),
    names_to = "type",
    values_to = "count"
  ) |>
  mutate(type = case_when(
    type == "confirm" ~ "Final observed data",
    type == "nowcast" ~ "Point nowcast",
    TRUE ~ type
  )) |>
  bind_rows(prep_latest_data)

# Create plot with data type as a variable
ggplot(plot_data, aes(x = reference_date, y = count, color = type)) +
  geom_line() +
  scale_color_manual(values = c(
    "Real-time data" = "darkred",
    "Final observed data" = "black",
    "Point nowcast" = "darkblue"
  )) +
  theme_bw() +
  xlab("Reference date") +
  ylab("Confirmed admissions") +
  scale_y_continuous(trans = "log10") +
  ggtitle("Comparing real-time, nowcasted, and later observed cases") +
  theme(legend.position = "bottom") +
  labs(color = "Type")


trunc_rep_tri_list <- baselinenowcast::truncate_triangles(reporting_triangle,
                                         n = n_retrospective_nowcasts
)
retro_rep_tri_list <- baselinenowcast::generate_triangles(trunc_rep_tri_list)

retro_pt_nowcast_mat_list <- baselinenowcast::generate_pt_nowcast_mat_list(
  reporting_triangle_list = retro_rep_tri_list,
  n = n_history_delay
)

disp_params <- baselinenowcast::estimate_dispersion(
  pt_nowcast_mat_list = retro_pt_nowcast_mat_list,
  trunc_rep_tri_list = trunc_rep_tri_list,
  reporting_triangle_list = retro_rep_tri_list,
  n = n_retrospective_nowcasts
)

nowcast_draws_df <- baselinenowcast::get_nowcast_draws(
  point_nowcast_matrix, reporting_triangle,
  dispersion = disp_params,
  draws = 100
)

head(nowcast_draws_df)

# Prepare the datasets for joining
latest_data_prepped <- latest_training_data |>
  mutate(time = row_number()) |>
  rename(obs_confirm = confirm) |>
  mutate(reference_date = as.Date(reference_date))

final_data_prepped <- target_data |>
  select(reference_date, final_confirm = confirm) |>
  mutate(reference_date = as.Date(reference_date))

# Join the datasets
obs_with_nowcast_draws_df <- nowcast_draws_df |>
  left_join(latest_data_prepped, by = "time") |>
  left_join(final_data_prepped, by = "reference_date")
head(obs_with_nowcast_draws_df)

# Create a separate dataframe for observed and final data
combined_data <- obs_with_nowcast_draws_df |>
  select(reference_date, obs_confirm, final_confirm) |>
  distinct() |>
  pivot_longer(
    cols = c(obs_confirm, final_confirm),
    names_to = "type",
    values_to = "count"
  ) |>
  mutate(type = case_when(
    type == "obs_confirm" ~ "Observed data",
    type == "final_confirm" ~ "Final observed data"
  ))

# Plot with draws for nowcast only
ggplot() +
  # Add nowcast draws as thin gray lines
  geom_line(
    data = obs_with_nowcast_draws_df,
    aes(
      x = reference_date, y = pred_count, group = draw,
      color = "Nowcast draw", linewidth = "Nowcast draw"
    )
  ) +
  # Add observed data and final data once
  geom_line(
    data = combined_data,
    aes(
      x = reference_date,
      y = count,
      color = type,
      linewidth = type
    )
  ) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "Nowcast draw" = "gray",
      "Observed data" = "darkred",
      "Final observed data" = "black"
    ),
    name = ""
  ) +
  scale_linewidth_manual(
    values = c(
      "Nowcast draw" = 0.2,
      "Observed data" = 1,
      "Final observed data" = 1
    ),
    name = ""
  ) +
  scale_y_continuous(trans = "log10") +
  xlab("Reference date") +
  ylab("cases") +
  theme(legend.position = "bottom") +
  ggtitle("Comparison of cases as of the nowcast date, later observed counts, \n and probabilistic nowcasted counts") # nolint
