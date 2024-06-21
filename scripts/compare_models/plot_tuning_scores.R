# Load tuning scores and plot as a table

#### SETUP ####

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source(paste0(wd, "./scripts/depends.R"))
# webshot::install_phantomjs() # nolint
output_path <- paste0(wd, "./outputs/tuning/")

#### GAM tuning score table ####

scores_gam <- vroom::vroom(
  file = paste0(wd, "./outputs/tuning/gam/specimen-date_num-re_report-date_num-re/scores.csv")
) |>
  dplyr::filter(scale == "natural",
    t_aggregation == "daily",
    denom_report < max_delay) |>
  dplyr::arrange(training_length,
    max_delay,
    denom_specimen,
    denom_report) |>
  dplyr::select(training_length,
    max_delay,
    denom_specimen,
    denom_report,
    interval_score
  )

# plot table
scores_gam |>
  gt::gt() |>
  gt::data_color(columns = interval_score,
    method = "numeric",
    palette = c("white", "red")
  ) |>
  gt::data_color(columns = c("training_length",
    "max_delay",
    "denom_specimen",
    "denom_report"),
  method = "factor",
  palette = "Blues") |>
  gt::cols_label(
    training_length = "Training length (days)",
    max_delay = "Maximum reporting delay (days)",
    denom_specimen = "Specimen date spline knot frequency (days)",
    denom_report = "Reporting delay spline knot frequency (days)",
    interval_score = "Weighted interval score"
  ) |>
  gt::cols_width(everything() ~ px(100)) |>
  gt::gtsave(paste0(output_path, "gam_tuning_scores.html"))

webshot::webshot(url = paste0(output_path, "gam_tuning_scores.html"),
  file = paste0(output_path, "gam_tuning_scores.png"))


#### BSTS tuning score table ####

scores_bsts <- vroom::vroom(
  file = paste0(wd, "./outputs/tuning/bsts/scores.csv")
) |>
  dplyr::filter(scale == "natural",
    t_aggregation == "daily") |>
  dplyr::select(
    training_length,
    sigma_level,
    sigma_slope,
    interval_score
  ) |>
  dplyr::arrange(
    training_length,
    sigma_level,
    sigma_slope
  )

# plot table
scores_bsts |>
  gt::gt() |>
  gt::data_color(columns = interval_score,
    method = "numeric",
    palette = c("white", "red")
  ) |>
  gt::data_color(columns = c(
    "training_length",
    "sigma_level",
    "sigma_slope"),
  method = "factor",
  palette = "Blues") |>
  gt::cols_label(
    training_length = "Training length (days)",
    sigma_level = gt::html("exp(&sigma;<sub>&mu;</sub>)"),
    sigma_slope = gt::html("exp(&sigma;<sub>&delta;</sub>)"),
    interval_score = "Weighted interval score"
  ) |>
  gt::cols_width(everything() ~ px(100)) |>
  gt::gtsave(paste0(output_path, "bsts_tuning_scores.html"))

webshot::webshot(url = paste0(output_path, "bsts_tuning_scores.html"),
  file = paste0(output_path, "bsts_tuning_scores.png"))



#### BSTS + 111 online tuning score table ####

scores_bsts_111 <- vroom::vroom(
  file = paste0(wd, "./outputs/tuning/bsts_111/scores.csv")
) |>
  dplyr::filter(scale == "natural",
    t_aggregation == "daily",
    training_length > 60,
    model_size > 1) |>
  dplyr::select(
    training_length,
    model_size,
    sigma_level,
    sigma_slope,
    interval_score
  ) |>
  dplyr::arrange(
    training_length,
    model_size,
    sigma_level,
    sigma_slope
  )

# plot table
scores_bsts_111 |>
  gt::gt() |>
  gt::data_color(columns = interval_score,
    method = "numeric",
    palette = c("white", "red")
  ) |>
  gt::data_color(columns = c(
    "training_length",
    "model_size",
    "sigma_level",
    "sigma_slope"),
  method = "factor",
  palette = "Blues") |>
  gt::cols_label(
    training_length = "Training length (days)",
    model_size = "Expected model size",
    sigma_level = gt::html("exp(&sigma;<sub>&mu;</sub>)"),
    sigma_slope = gt::html("exp(&sigma;<sub>&delta;</sub>)"),
    interval_score = "Weighted interval score"
  ) |>
  gt::cols_width(everything() ~ px(100)) |>
  gt::gtsave(paste0(output_path, "bsts_111_tuning_scores.html"))

webshot::webshot(url = paste0(output_path, "bsts_111_tuning_scores.html"),
  file = paste0(output_path, "bsts_111_tuning_scores.png"))


# EPINOWCAST SCORING ####
scores_epinowcast <- vroom::vroom(
  file = glue::glue("{output_path}/epinowcast/scoring/scores.csv")
) |>
  dplyr::filter(scale == "natural") |>
  dplyr::select(
    "training_length",
    "max_delay",
    "interval_score"
  ) |>
  dplyr::arrange(
    training_length,
    max_delay
  )

scores_epinowcast |>
  gt::gt() |>
  gt::data_color(columns = interval_score,
    method = "numeric",
    palette = c("white", "red")
  ) |>
  gt::data_color(columns = c(
    "training_length",
    "max_delay"),
  method = "factor",
  palette = "Blues") |>
  gt::cols_label(
    training_length = "Training length (days)",
    max_delay = "Maximum reporting delay (days)",
    interval_score = "Weighted interval score"
  ) |>
  gt::cols_width(everything() ~ px(100)) |>
  gt::gtsave(paste0(output_path, "epinowcast_tuning_scores.html"))

webshot::webshot(url = paste0(output_path, "epinowcast_tuning_scores.html"),
  file = paste0(output_path, "epinowcast_tuning_scores.png"))
