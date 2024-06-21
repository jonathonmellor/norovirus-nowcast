# Plot to show 111 trends compared with noro cases



source("./scripts/depends.R")
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")


ggplot2::theme_set(ggplot2::theme_bw())


# config

combined_path <- "./outputs/data/cases_with_noise.csv"
nhs111_datapath <- "./outputs/data/oneoneone_online_noro_symptoms.csv"

study_start_date_visual <- as.Date(config$dates$start_date) - 6 # shift to beginning of week
study_end_date_visual <- as.Date(config$dates$evaluate_end_date)

nhs111_data <- vroom::vroom(file = nhs111_datapath) |>
  dplyr::mutate(specimen_date = as.Date(specimen_date)) |>
  dplyr::select(-c("online_all_pain"))

# PROCESS ####

# collabse tests
cases <- aws.s3::s3read_using(
  vroom::vroom,
  object = combined_path
) |>
  dplyr::summarise(target = sum(target, na.rm = TRUE), .by = "specimen_date")

case_111 <- nhs111_data |>
  dplyr::filter(specimen_date >= study_start_date_visual &
    specimen_date <= study_end_date_visual) |>
  # long for easier facet when plotting
  tidyr::pivot_longer(cols = dplyr::starts_with("online"),
    names_to = "indicator", values_to = "metric_value") |>
  # clean up names, may want to revise depending on paper choices
  dplyr::mutate(indicator = stringr::str_remove(indicator, "online_") |>
    stringr::str_replace_all("_", " ")) |>
  # as each facet will have cases signal, we want to join rather than bind rows
  dplyr::left_join(cases, by = "specimen_date") |>
  # create a rolling mean version, a scaled version, and both
  dplyr::mutate(metric_value_smooth =
    zoo::rollmean(metric_value, k = 7,
      na.pad = TRUE,
      align = "right"),
  target_smooth = zoo::rollmean(target,
    k = 7,
    na.pad = TRUE,
    align = "right"),
  metric_value_smooth_scaled = scales::rescale(metric_value_smooth),
  target_smooth_scaled = scales::rescale(target_smooth),
  metric_value_scaled = scales::rescale(metric_value),
  target_scaled = scales::rescale(target),
  .by = "indicator")


unsmoothed_plot <- case_111 |>
  ggplot() +
  geom_line(aes(x = specimen_date, y = metric_value_scaled,
    color = indicator)) +
  geom_line(aes(x = specimen_date, y = target_scaled), color = "black") +
  facet_wrap(~indicator, nrow = 2) +
  ggthemes::scale_color_few() +
  guides(color = "none") +
  xlab("Specimen date") +
  ggtitle("b.") +
  ylab("Scaled metric")

unsmoothed_plot


smoothed_plot <- case_111 |>
  ggplot() +
  geom_line(aes(x = specimen_date, y = metric_value_smooth_scaled, color = indicator)) +
  geom_line(aes(x = specimen_date, y = target_smooth_scaled), color = "black") +
  facet_wrap(~indicator, nrow = 2) +
  ggthemes::scale_color_few() +
  guides(color = "none") +
  xlab("") +
  ggtitle("a.") +
  ylab("7-day average scaled metric")

smoothed_plot

nhs111_plot <- smoothed_plot / unsmoothed_plot

nhs111_plot


ggplot2::ggsave(
  glue::glue("./outputs/111_figure.png"),
  nhs111_plot,
  width = 10,
  height = 10)

# Make another graph to show relative scale better.
# keep separate as each y axis will mess up the patchwork x alignment

unscaled_plot <- case_111 |>
  ggplot() +
  geom_line(aes(x = specimen_date, y = metric_value / 1000, color = indicator)) +
  facet_wrap(~indicator, nrow = 2, scales = "free_y") +
  ggthemes::scale_color_few() +
  guides(color = "none") +
  xlab("Date") +
  ylab("Metric/1000")


ggplot2::ggsave(
  glue::glue("./111_figure_unscaled.png"),
  unscaled_plot,
  width = 10,
  height = 6)
