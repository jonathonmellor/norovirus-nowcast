# script for descriptive time series plots




source("./scripts/depends.R")

config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")


# config
combined_path <- "./outputs/data/cases_with_noise.csv"
# shift to beginning of week
study_start_date_visual <- as.Date(config$dates$start_date) - 6
study_end_date_visual <- as.Date(config$dates$evaluate_end_date)


ggplot2::theme_set(ggplot2::theme_bw())


# setup
raw <- aws.s3::s3read_using(
  vroom::vroom,
  object = combined_path
)

spine <- data.frame(specimen_date = seq(min(raw$specimen_date),
  max(raw$specimen_date), by = "day"))

loaded <- spine |>
  dplyr::left_join(raw, by = "specimen_date") |>
  dplyr::mutate(week = lubridate::isoweek(specimen_date),
    year = lubridate::isoyear(specimen_date)) |>
  # max brings it to Sunday... add 2
  dplyr::mutate(prediction_end_date = max(specimen_date),
    .by = c(year, week)) |>
  dplyr::filter(prediction_end_date != max(prediction_end_date))


# Figure 1B, weekly ########

# explore day of week

# we aren't ingesting data consistently...
loaded$prediction_end_date |> unique() |> lubridate::wday(label = TRUE)

# data within the week of ingest
initial <- loaded |>
  # dplyr::mutate(week_start = prediction_end_date - 6) |>
  dplyr::filter(specimen_date + days_to_reported <= prediction_end_date) |>
  dplyr::summarise(n_tests_initial = sum(target, na.rm = TRUE),
    .by = c("specimen_date", "prediction_end_date"))

# finalised data tests
revisions <- loaded |>
  # dplyr::mutate(week_start = prediction_end_date - 6) |>
  dplyr::summarise(n_tests_final = sum(target, na.rm = TRUE),
    .by = c("specimen_date", "prediction_end_date"))

# both initial and revised data, can calculate the change
combined <- revisions |>
  dplyr::left_join(initial, by = c("specimen_date", "prediction_end_date")) |>
  dplyr::mutate(n_tests_initial = dplyr::coalesce(n_tests_initial, 0)) |>
  dplyr::mutate(n_revisions = n_tests_final - n_tests_initial) |>
  tidyr::pivot_longer(cols = c("n_tests_final",
    "n_revisions",
    "n_tests_initial"),
  names_to = "state",
  values_to = "n_tests")

# show initial with revisions on top
weekly_plot <- combined |>
  dplyr::filter(state != "n_tests_initial") |>
  dplyr::mutate(week_start = prediction_end_date - 6) |>
  dplyr::summarise(n_tests = sum(n_tests), .by = c("week_start", "state")) |>
  dplyr::filter(week_start >= study_start_date_visual &
    week_start <= study_end_date_visual) |>
  ggplot() +
  geom_col(aes(x = week_start, y = n_tests, group = state, fill = state),
    color = "grey10") +
  ylab("Positive test count") +
  xlab("Week of ingest") +
  theme(legend.position = "bottom") +
  ggtitle("b.") +
  scale_fill_manual("Test reporting", values = c("maroon", "royalblue"),
    labels = c("final", "initial"))

weekly_plot



# Figure 1A, daily #######

# show how the data tails off nearer max_reporting_date
# need a different data structure (wide not long) for this kind of ggplot
line_plot <- combined |>
  tidyr::pivot_wider(names_from = state, values_from = n_tests) |>
  dplyr::mutate(prediction_end_date_fct = as.factor(prediction_end_date)) |>
  dplyr::filter(specimen_date >= study_start_date_visual &
    specimen_date <= study_end_date_visual) |>
  ggplot() +
  geom_line(aes(x = specimen_date, y = n_tests_final, linetype = "final")) +
  geom_line(aes(x = specimen_date, y = n_tests_initial,
    color = prediction_end_date_fct, linetype = "initial"),
  linewidth = 0.8) +
  geom_point(aes(x = specimen_date, y = n_tests_initial,
    color = prediction_end_date_fct), size = 2) +
  xlab("Specimen date") +
  ylab("Positive test count") +
  ggtitle("a.") +
  scale_color_discrete(name = "Ingest date", na.translate = FALSE) +
  guides(color = "none") +
  theme(legend.position = "bottom") +
  scale_linetype_manual(name = "Report type", values = c(8, 1))

line_plot

figure1 <- line_plot / weekly_plot

figure1

ggplot2::ggsave(
  glue::glue("./outputs/figure1.png"),
  figure1,
  width = 10,
  height = 11
)
