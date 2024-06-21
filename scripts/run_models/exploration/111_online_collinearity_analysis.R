# this script check the correlation of variables for noro now casting both visual and with ccf tests

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")

# loading in data

datapath <- "./outputs/data/cases_with_noise.csv"
nhs111_datapath <- "./outputs/data/oneoneone_online_noro_symptoms.csv"

max_date <- c("2024-03-10") ## setting the max date of the analysis

latest_data <- readr::read_csv(datapath) |>
  dplyr::filter(specimen_date <= max(as.Date(max_date))) |>
  dplyr::summarise(target = sum(target, na.rm = T), .by = specimen_date)


# getting the 111 data

regressors <- c("online_gi_symptoms",
  "online_fever",
  "online_headache",
  "online_all_pain",
  "online_limb_pain",
  "online_stomach_pain",
  "online_not_limb_not_stomach_pain")

nhs111_data <- vroom::vroom(file = nhs111_datapath) |>
  dplyr::mutate(specimen_date = as.Date(specimen_date)) |>
  dplyr::select(specimen_date, dplyr::all_of(regressors))

## join data
noro_111_combo <- latest_data |>
  dplyr::left_join(nhs111_data, by = "specimen_date") |>
  tidyr::pivot_longer(!specimen_date, names_to = "symptom", values_to = "count") |>
  dplyr::filter(specimen_date >= "2023-10-01") |>
  dplyr::group_by(symptom) |>
  # rescale between 0 and 1
  dplyr::mutate(count = scales::rescale(count, to = c(0, 1), from = range(count, na.rm = TRUE, finite = TRUE), ))

# basic plot of all symptoms
all_symptom_and_target_plot <- ggplot2::ggplot(noro_111_combo, ggplot2::aes(specimen_date, y = count)) +
  ggplot2::geom_line(group = 1) +
  ggplot2::labs(x = "prediction end date", y = "") +
  ggplot2::facet_wrap(. ~ symptom)
all_symptom_and_target_plot


# getting a variable for just target for easier plotting

target <- noro_111_combo |> dplyr::filter(symptom == "target")

# plot loop
for (x in regressors) {
  one_var <- noro_111_combo |>  dplyr::filter(symptom == x)

  symptom <- as.character(x)
  plot <- ggplot2::ggplot(one_var) +
    ggplot2::geom_line(mapping = ggplot2::aes(target$specimen_date, y = target$count, color = "red")) + # creating plots to allow for a visual
    ggplot2::geom_line(ggplot2::aes(specimen_date, y = count)) +
    ggplot2::ggtitle(glue::glue("{symptom}"))
  ggplot2::ggsave(glue::glue("target_&_{symptom}.png"), plot, path = "./outputs/colinearity_plots/")

  for (y in regressors) {
    one_var2 <- noro_111_combo |>  dplyr::filter(symptom == y)
    symptom2 <- as.character(y)
    ccf(one_var$count, one_var2$count, lag = 8, main = glue::glue("target & {symptom2}")) ### this creates ccf plots for all variables against each other,
    # change one_var to target to see correlation between explanatory and output variables


  }

}
