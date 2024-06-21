# Load data from each model and score / plot together

# IMPORTANT: Gets prediction dataframes from S3, so must be uploaded to model_data_path

#### SETUP ####

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source(paste0(wd, "./scripts/depends.R"))
source(paste0(wd, "./scripts/run_models/functions/plotting.R"))
output_path <- paste0(wd, "./outputs/combined/")
config <- yaml::read_yaml(paste0(wd, "./scripts/run_models/norovirus_nowcast_config.yaml"))
model_data_path <- "./outputs/data/"

training_data_path <- "./outputs/data/cases_with_noise.csv"

prediction_filenames <- c(
  "gam_predictions_summary.csv",
  "bsts_predictions_summary.csv",
  "bsts_111_online_predictions_summary.csv",
  "epinowcast_predictions_summary.csv",
  "baseline_prevweek_predictions_summary.csv"
)


#### LOAD DATA ####

# load prediction data
predictions <- tibble::tibble(
  prediction_filename = prediction_filenames
) |>
  dplyr::mutate(
    prediction_path = paste0(model_data_path, prediction_filename),
    predictions = purrr::map(
      .x = prediction_path,
      .f = ~ aws.s3::s3read_using(
        vroom::vroom,
        object = .x)
    )
  ) |>
  dplyr::pull(predictions) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    specimen_date = as.Date(specimen_date),
    prediction_end_date = as.Date(prediction_end_date),
    # when using a single baseline
    model = ifelse(model == "Baseline: previous week", "Baseline", model)
  ) |>
  dplyr::select(-dplyr::any_of(c("target_value")))


# load training data
training_data <- aws.s3::s3read_using(
  vroom::vroom,
  object = training_data_path,
  show_col_types = FALSE
) |>
  dplyr::mutate(specimen_date = lubridate::ymd(specimen_date))

training_data_daily <- training_data |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(t_aggregation = "daily")

training_data_weekly <- training_data_daily |>
  dplyr::mutate(specimen_date = lubridate::floor_date(specimen_date, unit = "week", week_start = 1)) |>
  dplyr::group_by(specimen_date) |>
  dplyr::summarise(target_value = sum(target_value, na.rm = TRUE),
    .groups = "drop") |>
  dplyr::mutate(t_aggregation = "weekly")

# join latest values to predictions dataframe
predictions <- predictions |>
  dplyr::full_join(dplyr::bind_rows(training_data_daily,
    training_data_weekly),
  by = c("specimen_date", "t_aggregation"))

#### PLOTTING ####

plot_nowcast(
  data = predictions |>
    dplyr::mutate(model = dplyr::case_when(
      model == "BSTS + NHS 111 online" ~ "BSTS + \nNHS 111 online",
      TRUE ~ model)) |>
    dplyr::filter(t_aggregation == "daily"),
  training_data = training_data,
  model_name = "all",
  plot_type = "lookbacks",
  output_path = output_path,
  y_limit = 110,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6,
  save_plot = FALSE)

ggplot2::ggsave(
  filename = paste0(output_path, "/all_lookbacks.png"),
  width = 12,
  height = 12,
  dpi = 500
)

plot_nowcast(
  data = predictions |>
    dplyr::mutate(model = dplyr::case_when(
      model == "BSTS + NHS 111 online" ~ "BSTS + \nNHS 111 online",
      TRUE ~ model)) |>
    dplyr::filter(t_aggregation == "weekly"),
  training_data = training_data,
  model_name = "all",
  plot_type = "lookbacks_weekly",
  output_path = output_path,
  y_limit = 750,
  x_limit_upper = NA,
  x_limit_lower = as.Date(config$dates$start_date) - 6,
  save_plot = FALSE)

ggplot2::ggsave(
  filename = paste0(output_path, "/all_lookbacks_weekly.png"),
  width = 12,
  height = 12,
  dpi = 500
)



#### SCORING ####

##### Calculate scores #####

# full set of scoring metrics used
# names formatted for plot
scoring_metrics <- c(
  "interval_score" = "WIS",
  "interval_score_log" = "Log WIS",
  "relative_is_baseline" = "Weighted interval \nskill score",
  "ae_median" = "Median absolute error",
  "bias" = "Bias",
  "interval_score_mean" = "WIS \n(daily mean)",
  "interval_score_log_mean" = "Log WIS \n(daily mean)",
  "relative_is_baseline_mean" = "Weighted interval \nskill score \n(daily mean)",
  "ae_median_mean" = "Median absolute error \n(daily mean)",
  "bias_mean" = "Bias \n(daily mean)",
  "coverage_90_mean" = "90% PI coverage \n(daily mean)",
  "coverage_50_mean" = "50% PI coverage \n(daily mean)",
  "coverage_deviation_mean" = "Coverage \ndeviation \n(daily mean)"
)

# score predictions
predictions_scored <- predictions |>
  # select only the predictions for the last 7 days
  dplyr::filter(specimen_date <= prediction_end_date,
    specimen_date > prediction_end_date - 7) |>
  # transform specimen date for weekly predictions to be end of week
  dplyr::mutate(
    specimen_date = ifelse(t_aggregation == "weekly", specimen_date + 6, specimen_date),
    specimen_date = as.Date(specimen_date, format = "%Y-%m-%d")) |>
  tidyr::pivot_longer(cols = dplyr::starts_with("pi_"), names_to = "quantile", values_to = "prediction") |>
  # keep only PIs we care about
  dplyr::filter(quantile %in% c("pi_50", "pi_25", "pi_75", "pi_95", "pi_5")) |>
  dplyr::mutate(quantile = as.numeric(stringr::str_remove(quantile, "pi_")) / 100) |>
  dplyr::rename(true_value = target_value) |>
  # score
  scoringutils::transform_forecasts(fun = scoringutils::log_shift, offset = 1) |>
  scoringutils::score() |>
  scoringutils::add_coverage(
    by = c("model", "specimen_date", "prediction_end_date", "scale", "t_aggregation"),
    ranges = c(50, 90)) |>
  scoringutils::summarise_scores(
    by = c("model", "specimen_date", "prediction_end_date", "scale", "t_aggregation"),
    na.rm = TRUE
  )

# calculate mean daily scores per week
scores_mean_daily <- predictions_scored |>
  dplyr::filter(scale == "natural",
    t_aggregation == "daily") |>
  dplyr::select(-c("specimen_date", "scale", "t_aggregation")) |>
  scoringutils::summarise_scores(
    by = c("model", "prediction_end_date"),
    na.rm = TRUE
  ) |>
  # set specimen date as end of week
  dplyr::rename(specimen_date = prediction_end_date) |>
  tidyr::pivot_longer(cols = -c(specimen_date, model),
    names_to = "metric",
    values_to = "value") |>
  dplyr::mutate(t_aggregation = "daily",
    metric = paste0(metric, "_mean"))

# long format for plotting
predictions_scored_long <- predictions_scored |>
  # calculate interval score relative to baselines
  dplyr::group_by(specimen_date, prediction_end_date, scale, t_aggregation) |>
  # if baseline WIS is 0, define differently
  dplyr::mutate(
    relative_is_baseline = dplyr::case_when(
      interval_score[model == "Baseline"] != 0 ~ 1 - interval_score / interval_score[model == "Baseline"],
      TRUE ~ 0)) |>
  dplyr::ungroup() |>
  # add scale into scoring metric names
  tidyr::pivot_wider(
    id_cols = c(model, specimen_date, t_aggregation, prediction_end_date),
    names_from = c(scale),
    values_from = -c(model, specimen_date, t_aggregation, scale, prediction_end_date)
  ) |>
  dplyr::rename_with(~ str_replace(., "_natural", "")) |>
  dplyr::select(-c("prediction_end_date")) |>
  tidyr::pivot_longer(cols = -c(specimen_date, model, t_aggregation),
    names_to = "metric",
    values_to = "value") |>
  # add daily averaged coverage
  dplyr::bind_rows(scores_mean_daily) |>
  dplyr::mutate(model = ifelse(model == "BSTS + NHS 111 online", "BSTS + \nNHS 111 online", model)) |>
  # select metrics
  dplyr::filter(metric %in% names(scoring_metrics))


# add mean daily interval score relative to baseline
# (unable to do this as part of scores_mean_daily as not a recognised scoring metric)
predictions_scored_long <- predictions_scored_long |>
  dplyr::bind_rows(
    predictions_scored_long |>
      dplyr::filter(t_aggregation == "daily",
        metric == "relative_is_baseline") |>
      dplyr::mutate(prediction_end_date = lubridate::floor_date(specimen_date,
        unit = "week",
        week_start = 1) + 6) |>
      dplyr::summarise(value = mean(value),
        .by = c("prediction_end_date", "t_aggregation", "model")) |>
      dplyr::mutate(metric = "relative_is_baseline_mean") |>
      dplyr::rename(specimen_date = prediction_end_date)
  )

# set prediction end dates for x axis ticks
prediction_end_dates <- seq(min(predictions$prediction_end_date, na.rm = TRUE) - 7,
  max(predictions$prediction_end_date, na.rm = TRUE), by = "1 week")

##### Overall scores #####

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_table <- c(
  "interval_score" = "WIS",
  "ae_median" = "Median absolute error",
  "bias" = "Bias",
  "coverage_deviation" = "Coverage deviation"
)

# overall score - daily and weekly
overall_scores_daily_weekly <- predictions_scored |>
  dplyr::filter(scale == "natural") |>
  scoringutils::summarise_scores(by = c("model", "scale", "t_aggregation"), na.rm = TRUE) |>
  scoringutils::summarise_scores(fun = round, digits = 3) |>
  dplyr::select(model, t_aggregation, !!names(scoring_metrics_table)) |>
  # calculate interval score relative to baselines
  # dplyr::group_by(t_aggregation) |>
  # dplyr::mutate(
  #   `Weighted interval \nskill score` = interval_score / interval_score[model == "Baseline"]) |>
  # dplyr::ungroup() |>
  dplyr::arrange(t_aggregation, model) |>
  # format names
  dplyr::rename(Model = model,
    `Temporal granularity` = t_aggregation) |>
  dplyr::rename_at(
    vars(names(scoring_metrics_table)), ~ scoring_metrics_table[.]
  )

# save table
write.csv(overall_scores_daily_weekly,
  file = paste0(output_path, "scoring_overall.csv"),
  row.names = FALSE)


##### Score - by day of prediction 1 to 7 #####

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_dow <- c(
  "interval_score" = "(a) WIS",
  "relative_is_baseline" = "(b) Weighted interval \nskill score",
  "bias" = "(c) Bias",
  "coverage_deviation" = "(d) Coverage \ndeviation"
)
# other options: interval_score_log, ae_median

overall_scores_dow_plot <- predictions_scored |>
  dplyr::filter(t_aggregation == "daily") |>
  dplyr::filter(scale == "natural") |>
  # score by dow
  dplyr::mutate(index = as.numeric(specimen_date + 7 - prediction_end_date)) |>
  scoringutils::summarise_scores(by = c("model", "index"), na.rm = TRUE) |>
  scoringutils::summarise_scores(fun = round, digits = 3) |>
  # calculate interval score relative to baselines
  dplyr::group_by(index) |>
  dplyr::mutate(
    relative_is_baseline = 1 - interval_score / interval_score[model == "Baseline"]) |>
  dplyr::ungroup() |>
  dplyr::select(model, index, !!names(scoring_metrics_dow)) |>

  # longer format for plotting
  tidyr::pivot_longer(cols = -c(index, model),
    names_to = "metric",
    values_to = "value") |>
  dplyr::mutate(metric = scoring_metrics_dow[metric],
    metric = factor(metric, levels = unname(scoring_metrics_dow))) |>
  dplyr::mutate(model = dplyr::case_when(
    model == "BSTS + NHS 111 online" ~ "BSTS + \nNHS 111 online",
    TRUE ~ model
  )) |>

  # plot
  ggplot2::ggplot(aes(x = index, y = value, color = model, linetype = model)) +
  ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_grid(metric ~ ., switch = "y", scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  ggplot2::scale_x_continuous(breaks = seq(1, 7), labels = c("Mon", "Tues", "Wed", "Thurs", "Fri", "Sat", "Sun")) +
  theme_ham() +
  theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA)) +
  ggplot2::labs(
    x = "Day of prediction",
    y = "",
    color = "Model",
    linetype = "Model"
  ) +
  theme(strip.text.y.left = element_text(angle = 0),
    legend.position = "bottom")
overall_scores_dow_plot

ggplot2::ggsave(filename = paste0(output_path, "/", "scoring_dow.png"), overall_scores_dow_plot, width = 10, height = 9)


##### Plot noro tests for scoring #####

# calculate partial data up to end of prediction weeks
partial_data <- training_data |>
  # convert dates to last Sunday of week
  dplyr::mutate(prediction_end_date = lubridate::floor_date(specimen_date,
    unit = "week",
    week_start = 1) + 6) |>
  dplyr::filter(prediction_end_date <= max(predictions$prediction_end_date, na.rm = TRUE),
    prediction_end_date >= min(predictions$prediction_end_date, na.rm = TRUE),
    # plot only last 7 days of predictions
    specimen_date + days_to_reported <= prediction_end_date,
    specimen_date + days_to_reported > prediction_end_date - 7) |>
  dplyr::summarise(n_tests_initial = sum(target, na.rm = TRUE), .by = c("specimen_date", "prediction_end_date")) |>
  dplyr::mutate(metric = "(a) Positive test count")

partial_data_weekly <- partial_data |>
  dplyr::summarise(n_tests_initial = sum(n_tests_initial),
    .by = c("prediction_end_date", "metric")) |>
  dplyr::rename(specimen_date = prediction_end_date)

# plot
noro_tests_daily <- ggplot2::ggplot(data = partial_data,
  ggplot2::aes(x = specimen_date)) +
  # plot partial data
  ggplot2::geom_line(data = partial_data,
    ggplot2::aes(
      x = specimen_date,
      y = n_tests_initial,
      col = "initial",
      group = as.factor(prediction_end_date)),
    alpha = 1) +
  ggplot2::geom_point(data = partial_data,
    ggplot2::aes(
      x = specimen_date,
      y = n_tests_initial,
      col = "initial")) +
  # plot latest data
  ggplot2::geom_point(data = training_data_daily |>
    dplyr::filter(specimen_date <= max(prediction_end_dates),
      specimen_date >= min(prediction_end_dates)),
  ggplot2::aes(
    x = specimen_date,
    y = target_value,
    col = "final")) +
  ggplot2::scale_color_manual(values = c("initial" = "#B3589A", "final" = "black")) +
  ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
  ggplot2::labs(x = "",
    y = "",
    color = "Reported data") +
  theme_ham() +
  ggplot2::theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "top") +
  ggplot2::scale_x_date(breaks = prediction_end_dates,
    labels = NULL) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  # add strip label for metric name
  ggplot2::facet_grid(metric ~ ., switch = "y") +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))

noro_tests_weekly <- ggplot2::ggplot(data = partial_data_weekly,
  ggplot2::aes(x = specimen_date)) +
  # plot partial data
  ggplot2::geom_point(ggplot2::aes(
    x = specimen_date,
    y = n_tests_initial,
    col = "initial")) +
  # plot latest data
  ggplot2::geom_point(data = training_data_weekly |>
    dplyr::filter(specimen_date <= max(prediction_end_dates),
      specimen_date >= min(prediction_end_dates)),
  ggplot2::aes(
    x = specimen_date + 6,
    y = target_value,
    col = "final")) +
  ggplot2::geom_line(data = training_data_weekly |>
    dplyr::filter(specimen_date <= max(prediction_end_dates),
      specimen_date >= min(prediction_end_dates)),
  ggplot2::aes(
    x = specimen_date + 6,
    y = target_value,
    col = "final")) +
  ggplot2::scale_color_manual(values = c("initial" = "#B3589A", "final" = "black")) +
  ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
  ggplot2::labs(x = "",
    y = "",
    color = "Reported data") +
  theme_ham() +
  ggplot2::theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "top") +
  ggplot2::scale_x_date(breaks = prediction_end_dates,
    labels = NULL) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  # add strip label for metric name
  ggplot2::facet_grid(metric ~ ., switch = "y") +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))

##### Weekly scoring plot #####

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_weekly <- c(
  "interval_score" = "(b) WIS",
  "relative_is_baseline" = "(c) Weighted interval \n skill score",
  "bias" = "(d) Bias"
)
# other options: interval_score_log, ae_median

# scores based on weekly PIs
scores_weekly_plot <- predictions_scored_long |>
  dplyr::filter(t_aggregation == "weekly",
    metric %in% names(scoring_metrics_weekly)) |>
  dplyr::mutate(metric = scoring_metrics_weekly[metric],
    metric = factor(metric, levels = unname(scoring_metrics_weekly))) |>
  # prediction end date on x axis
  ggplot2::ggplot(ggplot2::aes(x = specimen_date)) +
  ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = value, color = model)) +
  ggplot2::geom_line(ggplot2::aes(y = value, color = model, linetype = model)) +
  ggplot2::facet_grid(metric ~ ., switch = "y", scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  theme_ham() +
  theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA)) +
  ggplot2::scale_x_date(breaks = prediction_end_dates, labels = scales::label_date_short()) +
  ggplot2::labs(
    x = "End of prediction week",
    y = "",
    color = "Model",
    linetype = "Model"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))


# combine with noro tests plot and save plot
(noro_tests_weekly +
  ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "pt"),
    axis.text.y = element_text(margin = margin(r = 0)))) /
  (scores_weekly_plot +
    ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom")) +
  plot_layout(heights = c(1, 3))

ggplot2::ggsave(filename = paste0(output_path, "/", "scoring_weekly.png"), width = 12, height = 9)


##### Daily average scoring plot #####

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_daily_avg <- c(
  "interval_score_mean" = "(b) WIS \n(daily mean)",
  "relative_is_baseline_mean" = "(c) Weighted interval \nskill score \n(daily mean)",
  "bias_mean" = "(d) Bias \n(daily mean)",
  "coverage_deviation_mean" = "(e) Coverage deviation \n(daily mean)"
)
# other options: interval_score, interval_score_log_mean, ae_median_mean,
# coverage_90_mean, coverage_50_mean

# only keep specified scoring metrics
scores_daily_avg <- predictions_scored_long |>
  dplyr::filter(t_aggregation == "daily",
    metric %in% names(scoring_metrics_daily_avg)) |>
  dplyr::mutate(metric = scoring_metrics_daily_avg[metric],
    metric = factor(metric, levels = unname(scoring_metrics_daily_avg)))

scores_daily_avg_plot <- scores_daily_avg |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_date)) +
  ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7) +
  ggplot2::geom_point(ggplot2::aes(y = value, color = model)) +
  ggplot2::geom_line(ggplot2::aes(y = value, color = model, linetype = model)) +
  ggplot2::facet_grid(metric ~ ., switch = "y", scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  theme_ham() +
  theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA)) +
  ggplot2::scale_x_date(breaks = prediction_end_dates, labels = scales::label_date_short()) +
  ggplot2::labs(
    x = "End of prediction week",
    y = "",
    color = "Model",
    linetype = "Model"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))

# combine with noro tests plot and save plot
(noro_tests_daily +
  ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "pt"),
    axis.text.y = element_text(margin = margin(r = 0)))) /
  (scores_daily_avg_plot +
    ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom")) +
  plot_layout(heights = c(1, 4))

ggplot2::ggsave(filename = paste0(output_path, "/", "scoring_daily_avg.png"), width = 12, height = 12)


##### Daily scoring plot #####

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_daily <- c(
  "interval_score" = "(b) WIS \n(daily)",
  "relative_is_baseline" = "(c) Weighted interval \nskill score \n(daily)"
)
# other options: bias, ae_median

# only keep specified scoring metrics
scores_daily <- predictions_scored_long |>
  dplyr::filter(t_aggregation == "daily",
    metric %in% names(scoring_metrics_daily)) |>
  dplyr::mutate(metric = scoring_metrics_daily[metric],
    metric = factor(metric, levels = unname(scoring_metrics_daily)))

scores_daily_plot <- scores_daily |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_date)) +
  ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = value, color = model)) +
  ggplot2::facet_grid(metric ~ ., switch = "y", scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  theme_ham() +
  theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA)) +
  ggplot2::scale_x_date(breaks = prediction_end_dates, labels = scales::label_date_short()) +
  ggplot2::labs(
    x = "Prediction date",
    y = "",
    color = "Model",
    linetype = "Model"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))

# combine with noro tests plot and save plot
(noro_tests_daily +
  ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "pt"),
    axis.text.y = element_text(margin = margin(r = 0)))) /
  (scores_daily_plot +
    ggplot2::theme(plot.margin = margin(0, 0, 0, 0, "pt"),
      legend.position = "bottom")) +
  plot_layout(heights = c(1, 2))

ggplot2::ggsave(filename = paste0(output_path, "/", "scoring_daily.png"), width = 12, height = 7)



##### Combined daily and weekly scoring plot #####

# Not included in paper

# choose scoring metrics to include
# and formatted metric names to display
scoring_metrics_daily_weekly <- c(
  "interval_score_daily" = "(a) WIS \n(daily)",
  "interval_score_log_mean_daily" = "(b) Log WIS \n(daily mean)",
  "ae_median_mean_daily" = "(c) Median absolute error \n(daily mean)",
  "interval_score_log_weekly" = "(d) Log WIS \n(weekly)",
  "ae_median_weekly" = "(e) Median absolute error \n(weekly)"
)
# other options: relative_is_baseline_mean_daily, bias_mean_daily,
# coverage_90_mean_daily, coverage_50_mean_daily,
# interval_score_weekly, relative_is_baseline_weekly, bias_weekly


scores_daily_weekly <- predictions_scored_long |>
  # format scoring metric names
  dplyr::mutate(
    metric = paste0(metric, "_", t_aggregation)) |>
  # only the scoring metrics we care about
  dplyr::filter(metric %in% names(scoring_metrics_daily_weekly)) |>
  dplyr::mutate(metric = scoring_metrics_daily_weekly[metric],
    metric = factor(metric, levels = unname(scoring_metrics_daily_weekly)))

scores_daily_weekly_plot <- scores_daily_weekly |>
  ggplot2::ggplot(ggplot2::aes(x = specimen_date)) +
  ggplot2::geom_hline(yintercept = 0, color = "darkgrey", alpha = 0.7) +
  ggplot2::geom_line(ggplot2::aes(y = value, color = model, linetype = model)) +
  ggplot2::geom_point(data = scores_daily_weekly |>
    dplyr::filter(metric != "(a) WIS \n(daily)"),
  ggplot2::aes(y = value, color = model)) +
  ggplot2::facet_grid(metric ~ ., switch = "y", scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2") +
  theme_ham() +
  theme(strip.background = element_blank(),
    strip.placement  = "outside",
    panel.spacing    = unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA)) +
  ggplot2::scale_x_date(breaks = prediction_end_dates, labels = scales::label_date_short()) +
  ggplot2::labs(
    x = "Prediction end date",
    y = "",
    color = "Model",
    linetype = "Model"
  ) +
  ggplot2::coord_cartesian(
    xlim = c(min(prediction_end_dates), max(prediction_end_dates))
  ) +
  ggplot2::theme(strip.text.y.left = element_text(angle = 0))

scores_daily_weekly_plot

ggplot2::ggsave(filename = paste0(output_path, "/", "scoring_weekly_daily.png"), width = 12, height = 12)
