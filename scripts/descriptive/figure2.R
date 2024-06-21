# script for plots on reporting delays




source("./scripts/depends.R")
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")

# config
combined_path <- "./outputs/data/cases_with_noise.csv"
study_start_date_visual <- as.Date(config$dates$start_date) - 6 # shift to beginning of week
study_end_date_visual <- as.Date(config$dates$evaluate_end_date)

ggplot2::theme_set(ggplot2::theme_bw())

# setup
loaded <- aws.s3::s3read_using(
  vroom::vroom,
  object = combined_path
)

# note: removing incomplete data
# preferably would truncate time series by max(specimen_date) - max(days_to_reported)
# to not show partial distributions


max_days_to_report_shown <- 10

# FIGURE 2 #######
# show the time delay from specimen date to SGSS receive date
prepared_delay <- loaded |>
  # remove incomplete specimen dates from most recent data
  dplyr::filter(specimen_date < max(specimen_date) - max(days_to_reported)) |>
  # use only data from start of study period?
  dplyr::filter(specimen_date >= study_start_date_visual) |>
  # painful factor manipulation to truncate x axis
  dplyr::mutate(
    delay_trunc =
      as.factor(
        dplyr::if_else(days_to_reported >= max_days_to_report_shown,
          paste0(max_days_to_report_shown, "+"),
          as.character(days_to_reported)
        )
      )
  ) |>
  dplyr::mutate(delay = as.integer(stringr::str_remove(delay_trunc, "\\+"))) |>
  # add in this variable to split by later on
  dplyr::mutate(period =
    dplyr::if_else(specimen_date <= config$dates$tune_end_date,
      "tune", "evaluate"))


overall_delay <- prepared_delay |>
  # calculate raw proportion
  # na.rm=T because of most recent missing data
  dplyr::summarise(
    target_count = sum(target, na.rm = TRUE),
    .by = c("delay", "delay_trunc")
  ) |>
  dplyr::mutate(target_prop = target_count / sum(target_count)) |>
  dplyr::arrange(delay)

time_delay_plot <- overall_delay |>
  ggplot() +
  geom_col(aes(
    x = forcats::fct_reorder(delay_trunc, delay),
    y = 100 * target_prop
  ), fill = "royalblue", color = "black") +
  guides(fill = "none") +
  ylab("Proportion of tests (%)") +
  xlab("Days to report")

time_delay_plot

ggplot2::ggsave(
  glue::glue("./outputs/figure2.png"),
  time_delay_plot,
  width = 8,
  height = 6
)

# exploration of a breakdown across time periods

split_delay <- prepared_delay |>
  # calculate raw proportion
  # na.rm=T because of most recent missing data
  dplyr::summarise(
    target_count = sum(target, na.rm = TRUE),
    .by = c("delay", "delay_trunc", "period")
  ) |>
  dplyr::mutate(target_prop = target_count / sum(target_count), .by = "period") |>
  dplyr::arrange(delay)

time_delay_plot_tune_evaluate <- split_delay |>
  dplyr::mutate(period = factor(period, levels = c("tune", "evaluate"))) |>
  ggplot() +
  geom_col(aes(
    x = forcats::fct_reorder(delay_trunc, delay),
    y = 100 * target_prop
  ), fill = "royalblue", color = "black") +
  guides(fill = "none") +
  ylab("Proportion of tests (%)") +
  xlab("Days to report") +
  facet_wrap(~period)

time_delay_plot_tune_evaluate_count <- split_delay |>
  dplyr::mutate(period = factor(period, levels = c("tune", "evaluate"))) |>
  ggplot() +
  geom_col(aes(
    x = forcats::fct_reorder(delay_trunc, delay),
    y = target_count
  ), fill = "royalblue", color = "black") +
  guides(fill = "none") +
  ylab("Count of tests") +
  xlab("Days to report") +
  facet_wrap(~period)

combined_time_delay_periods <- time_delay_plot_tune_evaluate /
  time_delay_plot_tune_evaluate_count

ggplot2::ggsave(
  glue::glue("./outputs/time_delay_supplement.png"),
  combined_time_delay_periods,
  width = 8,
  height = 6
)
