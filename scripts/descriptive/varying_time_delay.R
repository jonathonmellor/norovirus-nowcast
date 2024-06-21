# script to show varying time delay distribution
# NOTE: the `norovirus_tests` object uses individual patient data and has thus been removed from this repo.
# therefore this specific result cannot be reproduced directly while complying with patient confidentiallity.


source("./scripts/depends.R")
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")


# config
config <- yaml::read_yaml("./scripts/run_models/norovirus_nowcast_config.yaml")
study_start_date_visual <- as.Date(config$dates$start_date) - 6 # shift to beginning of week
study_end_date_visual <- as.Date(config$dates$evaluate_end_date)



ggplot2::theme_set(ggplot2::theme_bw())

tuning <- FALSE

# generate week ends
if (tuning) {
  max_reporting_dates <- seq(from = as.Date("2023-10-08"),
    to = as.Date("2023-10-29"),
    by = 7)
} else {
  max_reporting_dates <- seq(from = as.Date("2023-10-08"),
    to = as.Date("2024-03-10"),
    by = 7)
}




# line level approach to distribution ####
path <- "EXCLUDED"



noro_delay <- norovirus_tests |>
  tidyr::expand_grid(prediction_end_date = max_reporting_dates) |>
  # can't have tests after end date
  dplyr::filter(specimen_date <= prediction_end_date) |>
  dplyr::filter(specimen_date > prediction_end_date - 7) |>
  # calculate time delay
  dplyr::mutate(
    time_delay = sgss_received_date - specimen_date
  ) |>
  # segmet to time period of interest
  dplyr::filter(specimen_date >= study_start_date_visual & specimen_date <= study_end_date_visual) |>
  # generate distribution quantiles
  dplyr::summarise(
    ci_50 = quantile(time_delay, 0.5),
    ci_95 = quantile(time_delay, 0.95),
    ci_05 = quantile(time_delay, 0.05),
    ci_75 = quantile(time_delay, 0.75),
    ci_25 = quantile(time_delay, 0.25),
    mean = mean(time_delay),
    .by = c("prediction_end_date")
  )

noro_delay |>
  ggplot() +
  geom_ribbon(aes(x = prediction_end_date,
    ymax = ci_95,
    ymin = ci_05,
    alpha = "95% interval"),
  fill = "royalblue") +
  geom_ribbon(aes(x = prediction_end_date, ymax = ci_75, ymin = ci_25, alpha = "50% interval"),
    fill = "royalblue") +
  geom_line(aes(x = prediction_end_date, y = ci_50, linetype = "median")) +
  geom_line(aes(x = prediction_end_date, y = mean, linetype = "mean")) +
  scale_alpha_manual(name = "", values = c("50% interval" = 0.5, "95% interval" = 0.2)) +
  scale_linetype_manual(name = "", values = c("median" = 3, "mean" = 2)) +
  scale_y_continuous(name = "reporting delay (days)", breaks = seq(0, 14, 1), limits = c(0, 14)) +
  theme(legend.position = "bottom") +
  scale_x_date(name = "end of prediction week",
    breaks = scales::date_breaks(width = "1 week"),
    labels = scales::label_date_short()) -> time_varying_delay_plot

time_varying_delay_plot


ggplot2::ggsave(
  glue::glue("./outputs/time_varying_delay.png"),
  time_varying_delay_plot,
  width = 9,
  height = 7
)
