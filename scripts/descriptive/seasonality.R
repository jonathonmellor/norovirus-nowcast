# script to show day of week seasonality




source("./scripts/depends.R")
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")


# config
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")
study_start_date_visual <- as.Date(config$dates$start_date) - 6 # shift to beginning of week
study_end_date_visual <- as.Date(config$dates$evaluate_end_date)


ggplot2::theme_set(ggplot2::theme_bw())

cases <- aws.s3::s3read_using(
  vroom::vroom,
  object = combined_path
) |>
  dplyr::filter(specimen_date >= study_start_date_visual &
    specimen_date <= study_end_date_visual) |>
  dplyr::summarise(target = sum(target, na.rm = TRUE),
    .by = "specimen_date") |>
  # get the week to group by in plots
  # (if we have more than 1 year we should do this better)
  dplyr::mutate(specimen_week = lubridate::epiweek(specimen_date),
    # shift dow week to monday to align with nowcasts
    dow = lubridate::wday(specimen_date,
      label = TRUE,
      abbr = FALSE,
      week_start = getOption("lubridate.week.start", 1))) |>
  dplyr::mutate(target_scaled = target - mean(target), .by = specimen_week) |>
  dplyr::mutate(mean_target_scaled = mean(target_scaled), .by = dow)

# calculate some acf's to show autocorrelation of the time series (and hence seasonality)
acf_result <- data.frame(
  acf = acf(cases$target)$acf,
  lag = acf(cases$target)$lag
)


# seasonality plot over day of week
dow_plot <- cases |>
  ggplot() +
  geom_line(aes(x = dow, y = target_scaled, group = specimen_week), alpha = 0.5, color = "brown") +
  geom_line(aes(x = dow, y = mean_target_scaled, group = specimen_week), color = "black", linewidth = 1.2) +
  ylab("Difference from week average count") +
  ggtitle("a.") +
  xlab("Specimen date day of week")

dow_plot

# autocorrelation plot, highlighting 7 visually
auto_corr_plot <- acf_result |>
  dplyr::mutate(is_seven = lag %% 7 == 0) |>
  ggplot() +
  geom_col(aes(x = lag, y = acf, fill = is_seven)) +
  xlab("Lag (days)") +
  ylab("Correlation") +
  guides(fill = "none") +
  ggtitle("b.") +
  ggthemes::scale_fill_hc()

auto_corr_plot

seasonality_plot <- dow_plot / auto_corr_plot

seasonality_plot

ggplot2::ggsave(
  glue::glue("{wd}/norovirus/papers/nowcast/outputs/seasonality_plot.png"),
  seasonality_plot,
  width = 8,
  height = 8)
