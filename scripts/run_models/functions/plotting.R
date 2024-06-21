
#' Create nowcast plots for all lookbacks
#'
#' note: plots only last 7 days of predictions
#'
#' @param data Dataframe that need to be in the formatted model output and must include the columns:
#'  - specimen_date
#'  - prediction_end_date
#'  - model
#'  - target_value
#'  - pi_5
#'  - pi_90
#'  - pi_50
#' @param training_data Dataframe of the training data that must include the columns:
#'  - specimen_date
#'  - days_to_reported
#'  - target
#' @param model_name Name of model to be plotted
#' @param plot_type Type of nowcast plot:
#'  - "lookbacks": plot all daily lookbacks on a single plot
#'  - "lookbacks_weekly": plot all weekly lookbacks on a single plot
#' @param output_path Local output folder path to save plots into
#' @param y_limit numeric value for the maximum y-axis value to show on plots.
#' If you're happy to let the script pick this, leave it as NA.
#' @param x_limit_upper Date to provide upper limit on plot's x-axis.
#' If you're happy to let the script pick this, leave it as NA.
#' @param x_limit_lower Date to provide lower limit on plot's x-axis.
#' If you're happy to let the script pick this, leave it as NA.
#' @param save_plot logical: whether to save plot as png
#'
plot_nowcast <- function( # nolint: cyclocomp_linter.
    # TODO reduce cyclomatic complexity
    data,
    training_data,
    model_name,
    plot_type,
    output_path,
    y_limit = NA,
    x_limit_upper = NA,
    x_limit_lower = "2023-10-02",
    save_plot = TRUE
    ) {

  latest_data <- training_data |>
    dplyr::group_by(specimen_date) |>
    dplyr::summarise(target_value = sum(target, na.rm = TRUE),
      .groups = "drop")

  # define tune and evaluate weeks
  train_weeks <- seq(from = as.Date("2023-10-08"),
    to = as.Date("2023-10-29"),
    by = 7)
  test_weeks <- seq(from = as.Date("2023-11-05"),
    to = as.Date("2024-03-10"),
    by = 7)
  prediction_end_dates <- c(train_weeks, test_weeks)

  # format prediction data for plotting
  data <- data |>
    dplyr::filter(!is.na(model),
      # plot only last 7 days of predictions
      specimen_date > prediction_end_date - 7) |>
    dplyr::mutate(
      across(starts_with("pi_"), ~ ifelse(specimen_date > prediction_end_date - 7, ., NA)),
      train_or_test = dplyr::case_when(
        prediction_end_date %in% train_weeks ~ "tune",
        prediction_end_date %in% test_weeks ~ "evaluate",
        TRUE ~ NA),
      train_or_test = factor(train_or_test, levels = c("tune", "evaluate")),
      # only plot pi_50 for baseline models
      pi_50 = ifelse(stringr::str_detect(model, "aseline"), pi_50, NA))


  # calculate partial data up to end of prediction weeks
  partial_data <- training_data |>
    # convert dates to last Sunday of week
    dplyr::mutate(prediction_end_date = lubridate::floor_date(specimen_date,
      unit = "week",
      week_start = 1) + 6) |>
    dplyr::filter(prediction_end_date <= max(data$prediction_end_date, na.rm = TRUE),
      prediction_end_date >= min(data$prediction_end_date, na.rm = TRUE),
      # plot only last 7 days of predictions
      specimen_date + days_to_reported <= prediction_end_date,
      specimen_date + days_to_reported > prediction_end_date - 7) |>
    dplyr::summarise(n_tests_initial = sum(target, na.rm = TRUE), .by = c("specimen_date", "prediction_end_date"))

  # fprmat axis limits
  if (is.null(x_limit_lower) || is.na(x_limit_lower)) {
    x_limit_lower <- data$specimen_date |>
      min() |>
      as.Date()
  } else {
    x_limit_lower <- as.Date(x_limit_lower)
    data <- data |>
      dplyr::filter(specimen_date >= x_limit_lower)
    latest_data <- latest_data |>
      dplyr::filter(specimen_date >= x_limit_lower)
  }


  # Daily #####
  if (plot_type == "lookbacks") {

    data <- data |>
      dplyr::filter(t_aggregation == "daily")

    plt <- ggplot2::ggplot(
      data = data,
      ggplot2::aes(x = specimen_date)) +
      # plot PIs
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = pi_95,
          ymax = pi_5,
          fill = as.factor(train_or_test),
          group = as.factor(prediction_end_date),
          alpha = "90%")
      ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = pi_75,
          ymax = pi_25,
          fill = as.factor(train_or_test),
          group = as.factor(prediction_end_date),
          alpha = "50%")
      ) +
      ggplot2::scale_alpha_manual(values = c("90%" = 0.4, "50%" = 0.7)) +
      ggplot2::scale_fill_manual(values = c("tune" = "grey60", "evaluate" = "#9BBF85"))

    if (any(grepl("Baseline", data$model))) {
      # plot pi_50 for baseline models if any
      plt <- plt +
        ggplot2::geom_line(
          ggplot2::aes(y = pi_50, group = as.factor(prediction_end_date)),
          color = "#9BBF85",
          linewidth = 1) +
        ggplot2::guides(color = "none") +
        ggnewscale::new_scale_color()
    }

    for (partial_data_model in c("GAM", "epinowcast")) {
      if (partial_data_model %in% data$model) {
        # plot partial data
        partial_data <- partial_data |>
          dplyr::mutate(model = partial_data_model)
        plt <- plt +
          ggplot2::geom_line(
            data = partial_data,
            ggplot2::aes(
              x = specimen_date,
              y = n_tests_initial,
              col = "initial",
              group = as.factor(prediction_end_date)
            ),
            alpha = 1) +
          ggplot2::geom_point(
            data = partial_data,
            ggplot2::aes(x = specimen_date, y = n_tests_initial, col = "initial"))
      }
    }

    plt <- plt +
      # plot latest data
      ggplot2::geom_point(
        data = latest_data |> dplyr::filter(specimen_date <= max(data$specimen_date)),
        ggplot2::aes(x = specimen_date, y = target_value, col = "final")) +
      ggplot2::scale_color_manual(values = c("initial" = "#B3589A", "final" = "black")) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(x = "Specimen date",
        y = "Positive test count",
        color = "Reported data",
        fill = "Data split",
        alpha = "PI") +
      theme_ham() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(colour = "black", linewidth = 0.05, linetype = "dashed"),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        text = ggplot2::element_text(size = 12),
        legend.position = "bottom") +
      # set axis limits
      ggplot2::coord_cartesian(
        xlim = c(x_limit_lower, x_limit_upper),
        ylim = c(0, y_limit)
      ) +
      ggplot2::facet_wrap(~model,
        ncol = 1,
        strip.position = "right") +
      ggplot2::scale_x_date(breaks = prediction_end_dates,
        labels = scales::label_date_short())

    if (save_plot) {
      ggplot2::ggsave(filename = paste0(output_path,
        "/",
        model_name,
        "_lookbacks.png"),
      plt,
      width = 12,
      height = 10,
      dpi = 500
      )
    }
  }


  # Weekly #####
  if (plot_type == "lookbacks_weekly") {

    data <- data |>
      dplyr::filter(t_aggregation == "weekly")

    partial_data <- partial_data |>
      dplyr::summarise(n_tests_initial = sum(n_tests_initial), .by = c("prediction_end_date")) |>
      dplyr::rename(specimen_date = prediction_end_date)

    latest_data <- latest_data |>
      # group by end of week
      dplyr::mutate(specimen_date = lubridate::floor_date(specimen_date, unit = "week", week_start = 1) + 6) |>
      dplyr::summarise(target_value = sum(target_value), .by = specimen_date)

    plt <- ggplot2::ggplot(data = data) +
      # plot PIs
      # note that no middle line is displayed so input for middle is irrelevant but
      # cannot use pi_50 as that may be all NAs due to line 144
      ggplot2::geom_boxplot(
        ggplot2::aes(
          x = specimen_date + 6,
          ymin = pi_5,
          lower = pi_25,
          middle = pi_25,
          upper = pi_75,
          ymax = pi_95,
          fill = as.factor(train_or_test),
          group = as.factor(prediction_end_date)),
        stat = "identity",
        fatten = NULL) +
      ggplot2::scale_fill_manual(values = c("tune" = "grey60", "evaluate" = "#9BBF85"))

    for (partial_data_model in c("GAM", "epinowcast")) {
      if (partial_data_model %in% data$model) {
        # plot partial data
        partial_data <- partial_data |>
          dplyr::mutate(model = partial_data_model)
        plt <- plt +
          ggplot2::geom_point(
            data = partial_data,
            ggplot2::aes(x = specimen_date, y = n_tests_initial, col = "initial"))
      }
    }

    plt <- plt +
      # plot latest data
      ggplot2::geom_point(
        data = latest_data |> dplyr::filter(specimen_date <= max(data$prediction_end_date)),
        ggplot2::aes(x = specimen_date, y = target_value, col = "final")) +
      ggplot2::geom_line(
        data = latest_data |> dplyr::filter(specimen_date <= max(data$prediction_end_date)),
        ggplot2::aes(x = specimen_date, y = target_value, col = "final")) +
      ggplot2::scale_color_manual(values = c("initial" = "#B3589A", "final" = "black")) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::labs(x = "End of prediction week",
        y = "Positive test count",
        color = "Reported data",
        fill = "Data split",
        alpha = "PI") +
      theme_ham() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_line(colour = "black", linewidth = 0.05, linetype = "dashed"),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
        text = ggplot2::element_text(size = 12),
        legend.position = "bottom") +
      # set axis limits
      ggplot2::coord_cartesian(
        xlim = c(x_limit_lower, x_limit_upper),
        ylim = c(0, y_limit)
      ) +
      ggplot2::facet_wrap(~model,
        ncol = 1,
        strip.position = "right") +
      ggplot2::scale_x_date(breaks = prediction_end_dates,
        labels = scales::label_date_short())

    if (save_plot) {
      ggplot2::ggsave(filename = paste0(output_path,
        "/",
        model_name,
        "_lookbacks_weekly.png"),
      plt,
      width = 12,
      height = 10,
      dpi = 500
      )
    }
  }

  return(plt)

}


#' Create single nowcast plot for one lookbacks
#'
#' @param data Dataframe that need to be in the formatted model output and must include the columns:
#'  - specimen_date
#'  - prediction_end_date
#'  - model
#'  - target_value
#'  - pi_5
#'  - pi_90
#'  - pi_50
#' @param training_data Dataframe of the training data that must include the columns:
#'  - specimen_date
#'  - days_to_reported
#'  - target
#' @param model_name Name of model to be plotted
#' @param plot_type Type of nowcast plot:
#'  - "daily": plot separate daily nowcast for each lookback
#' @param output_path Local output folder path to save plots into
#' @param y_limit numeric value for the maximum y-axis value to show on plots.
#' If you're happy to let the script pick this, leave it as NA.
#' @param x_limit_upper Date to provide upper limit on plot's x-axis.
#' If you're happy to let the script pick this, leave it as NA.
#' @param x_limit_lower Date to provide lower limit on plot's x-axis.
#' If you're happy to let the script pick this, leave it as NA.
#'
plot_nowcast_single <- function(
    data,
    training_data,
    model_name,
    plot_type,
    output_path,
    y_limit = NA,
    x_limit_upper = NA,
    x_limit_lower = "2023-10-02"
    ) {

  latest_data <- training_data |>
    dplyr::group_by(specimen_date) |>
    dplyr::summarise(target_value = sum(target, na.rm = TRUE),
      .groups = "drop")

  plot_colors <- c("Data at time of nowcast" = "black", "Data since nowcast" = "red", "Nowcast" = "blue")

  for (nowcast_date in unique(data$prediction_end_date)) {

    data_filtered <- data |>
      dplyr::filter(prediction_end_date == as.Date(nowcast_date))

    min_date <- min(data_filtered$specimen_date, na.rm = TRUE)
    max_date <- max(data_filtered$specimen_date, na.rm = TRUE)

    if (is.null(x_limit_lower) || is.na(x_limit_lower)) {
      x_limit_lower_nowcast <- min_date |>
        as.Date()
    } else {
      x_limit_lower_nowcast <- as.Date(x_limit_lower)
    }

    latest_data_filtered <- latest_data |>
      dplyr::filter(specimen_date <= max_date,
        specimen_date >= min_date)

    if (plot_type == "daily") {

      plt <- ggplot2::ggplot(data_filtered, ggplot2::aes(x = specimen_date)) +
        ggplot2::geom_point(ggplot2::aes(y = target_value, color = "Data at time of nowcast"),
          size = 2) +
        ggplot2::geom_line(ggplot2::aes(y = target_value, color = "Data at time of nowcast")) +
        ggplot2::geom_point(data = latest_data_filtered,
          ggplot2::aes(y = target_value, color = "Data since nowcast")) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = pi_95, ymax = pi_5, fill = "Nowcast"),
          alpha = 0.2) +
        ggplot2::geom_line(ggplot2::aes(y = pi_50),
          color = "black",
          linewidth = 0.7,
          linetype = "dashed",
          alpha = 0.6) +
        ggplot2::labs(x = "Specimen date",
          y = glue::glue("Daily norovirus {gsub('_', ' ', 'positive_tests')}"),
          fill = paste0("Nowcast up to ", max(data_filtered$prediction_end_date), ":"),
          color = "Data",
          title = glue::glue(
            "Nowcast for norovirus {gsub('_', ' ', 'positive_tests')} ",
            "up to {max(data_filtered$prediction_end_date)}"
          )
        ) +
        ggplot2::scale_color_manual(values = plot_colors) +
        ggplot2::scale_fill_manual(values = plot_colors) +
        ggplot2::scale_x_date(label = scales::label_date_short(),
          breaks = "1 week") +
        theme_ham() +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_line(colour = "black", linewidth = 0.05, linetype = "dashed"),
          panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
          text = ggplot2::element_text(size = 12),
          legend.position = "bottom")
    }

    # set axis limits
    plt <- plt +
      ggplot2::coord_cartesian(
        xlim = c(x_limit_lower_nowcast, x_limit_upper),
        ylim = c(0, y_limit)
      )

    ggplot2::ggsave(filename = paste0(output_path,
      "/",
      model_name,
      "_",
      plot_type,
      "_",
      as.Date(nowcast_date),
      ".png"),
    plt,
    width = 9,
    height = 7,
    dpi = 500
    )
  }
}



#' Plotting theme for Health Analysis Modelling team
#'
#' @inheritParams ggthemes::theme_fivethirtyeight
#'
#' @export

theme_ham <- function(base_size = 14, base_family = "sans") {
  ggplot2::theme_bw(base_size, base_family) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(1.5, "lines"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.5, linetype = 1),
      plot.background = ggplot2::element_rect(fill = "white", colour = "white"),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white"),
      panel.grid = ggplot2::element_line(colour = "#D9D9D9"),
      strip.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_rect(colour = "black", fill = "white"),
      legend.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linetype = 1, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(linewidth = 0),
      plot.caption = ggplot2::element_text(size = 12, lineheight = 0.5, color = "black", hjust = 0),
      plot.caption.position = "plot",
      plot.title = ggtext::element_textbox(
        width = ggplot2::unit(1, "npc"), size = 16, padding = ggplot2::margin(t = 10)),
      plot.title.position = "plot"
    )
}
