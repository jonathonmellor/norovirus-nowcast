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



training_data <- vroom::vroom(training_data_path) |>
  # convert to epinowcast naming conventions
  dplyr::rename(reference_date = specimen_date,
                confirm = target) |>
  dplyr::mutate(report_date = reference_date + days_to_reported) |>
  dplyr::select(-days_to_reported)



# Run model

# Approach: transform synthetic data to format needed for package and run getting
# started page code in order.

# select a test date to get the code working in line with the package getting started page
nowcast_date <- max_reporting_dates[[1]]

observed_long <- training_data |>
  epinowcast::enw_filter_report_dates(latest_date = nowcast_date) |>
  epinowcast::enw_filter_reference_dates(include_days = config$hyperparams$gam$training_length - 1)

pobs <- epinowcast::enw_preprocess_data(
  obs = observed_long,
  max_delay = config$hyperparams$gam$max_delay + 1
)

reporting_triangle <- pobs$reporting_triangle[[1]] |>
  dplyr::select(-c(".group", "reference_date")) |>
  as.matrix() |>
  unname()

cleaned_reporting_triangle <- baselinenowcast::replace_lower_right_with_NA(reporting_triangle)

delay_pmf <- baselinenowcast::get_delay_estimate(
  reporting_triangle = cleaned_reporting_triangle,
  max_delay = config$hyperparams$gam$max_delay,
  n = config$hyperparams$gam$training_length
)

point_nowcast_matrix <- baselinenowcast::apply_delay(
  rep_tri_to_nowcast = cleaned_reporting_triangle,
  delay_pmf = delay_pmf
)

trunc_rep_mat_list <- baselinenowcast::truncate_triangles(
  reporting_triangle = cleaned_reporting_triangle
)

retro_rep_tri_list <- baselinenowcast::generate_triangles(
  reporting_triangle_list = trunc_rep_mat_list
)

# problem here?
retro_pt_nowcast_mat_list <- baselinenowcast::generate_pt_nowcast_mat_list(
  reporting_triangle_list = retro_rep_tri_list
)


disp_params <- baselinenowcast::estimate_dispersion(
  pt_nowcast_mat_list = retro_pt_nowcast_mat_list,
  trunc_rep_mat_list = trunc_rep_mat_list
)
n_for_delay_estimate <- min(
  sapply(retro_rep_tri_list, nrow))
for (i in 1:length(retro_rep_tri_list)) {
  # Get number of rows in the retrospective reporting triangle 
  nr0 <- nrow(retro_rep_tri_list[[i]])
  
  # For each reporting triangle, it will use the bottom `n_for_delay_estimate`
  # rows in the current function defaults (which are not the same as the default 
  # configuration, and should probably be considered more carefully).
  # Therefore, we need to check that any of these triangles contain 0s in the bottom 16 rows. 
  if (all(retro_rep_tri_list[[i]][(nr0-n_for_delay_estimate + 1):nr0,1] == 0)) print(i)
  
  # We see that 1, 2, 3, and 4 all contain the bottom 16 rows as NA. 
}

# Try again with `n_for_delay_estimate=28`-------------------------------------
# Running again being specific about how many triangles to make based on 
# specification of 56 as the total training data. 

n_for_delay_estimate <- config$hyperparams$gam$training_length/2
n_for_uncertainty_estimate <- config$hyperparams$gam$training_length - n_for_delay_estimate
delay_pmf <- baselinenowcast::get_delay_estimate(
  reporting_triangle = cleaned_reporting_triangle,
  max_delay = config$hyperparams$gam$max_delay,
  n = n_for_delay_estimate 
)

point_nowcast_matrix <- baselinenowcast::apply_delay(
  rep_tri_to_nowcast = cleaned_reporting_triangle,
  delay_pmf = delay_pmf
)

trunc_rep_mat_list <- baselinenowcast::truncate_triangles(
  reporting_triangle = cleaned_reporting_triangle,
  n = n_for_uncertainty_estimate # This n is the number of triangles to make, not the number
  # to use for delay estimation, though they happen to be the same 
)

retro_rep_tri_list <- baselinenowcast::generate_triangles(
  reporting_triangle_list = trunc_rep_mat_list
)
# Generate 28 retrospective reporting triangles, of which the delay is estimated
# from the last 28 rows of data. Since n
retro_pt_nowcast_mat_list <- baselinenowcast::generate_pt_nowcast_mat_list(
  reporting_triangle_list = retro_rep_tri_list
)

for (i in 1:length(retro_rep_tri_list)) {
  # Get number of rows in the retrospective reporting triangle 
  nr0 <- nrow(retro_rep_tri_list[[i]])
  
  # For each reporting triangle, it will use the bottom `n_for_delay_estimate`
  # Therefore, we need to check that any of these triangles contain 0s in the bottom 28 rows. 
  if (all(retro_rep_tri_list[[i]][(nr0-n_for_delay_estimate + 1):nr0,1] == 0)) print(i)
  
}
# No zeros in bottom 28 rows of any of the retrospective reporting triangles,
# therefore we can generate the retrospective point nowcasts 

disp_params <- baselinenowcast::estimate_dispersion(
  pt_nowcast_mat_list = retro_pt_nowcast_mat_list,
  trunc_rep_mat_list = trunc_rep_mat_list
)

