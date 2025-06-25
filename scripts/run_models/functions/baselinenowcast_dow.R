run_baselinenowcast_dow <- function(.data,
                                prediction_end_date,
                                n_pi_samples = 500,
                                model_hyperparams,
                                timestep = "day",
                                by = NULL,
                                ...) {

  # parameters
  max_delay <- model_hyperparams$max_delay
  n_training_volume <- model_hyperparams$training_length

  n_history_delay <- model_hyperparams$n_history_delay
  n_retrospective_nowcasts <- model_hyperparams$n_retrospective_nowcasts

  eval_timeframe <- model_hyperparams$eval_timeframe

  # prepare data
  target_data <- .data |>
    dplyr::rename(reference_date = specimen_date,
                  count = target) |>
    dplyr::mutate(report_date = reference_date + days_to_reported)

  # model each dow separately
  all_nowcasts <- data.frame()
  for( i in 1:7){

    # apply filters to simulate real-time truncation
    train_df_i <- epinowcast::enw_filter_report_dates(
      # filter to one dow
      obs = target_data[lubridate::wday(target_data$reference_date) == i,],
      latest_date = prediction_end_date
    ) |>
      epinowcast::enw_filter_reference_dates(
        latest_date = prediction_end_date
      ) |>
      dplyr::mutate(reference_date = as.Date(reference_date))

    rep_tri <- train_df_i |>
      dplyr::mutate(delay = as.integer(difftime(report_date,
                                         reference_date,
                                         units = "days"
      ))
      ) |>
      dplyr::select(reference_date, delay, count) |>
      dplyr::filter(delay <= max_delay, delay >= 0) |>
      tidyr::pivot_wider(
        names_from = delay,
        values_from = count
      ) |>
      dplyr::select(-reference_date) |>
      as.matrix()

    delay_pmf <- baselinenowcast::get_delay_estimate(
      reporting_triangle = rep_tri,
      max_delay = max_delay,
      n = n_history_delay
    )
    pt_nowcast_mat <- baselinenowcast::apply_delay(
      rep_tri_to_nowcast = rep_tri,
      delay_pmf = delay_pmf
    )

    # This will throw a warning because not all triangles can be nowcasted, in
    # practince n_retrospective_nowcasts will be less than 28
    trunc_rts <- baselinenowcast::truncate_triangles(
      reporting_triangle = rep_tri,
      n = n_retrospective_nowcasts
    )

    retro_rts <- baselinenowcast::generate_triangles(
      trunc_rep_tri_list = trunc_rts,
      structure = c(1, 7)
    )
    # This is going to throw a bunch of warnings because not all
    # triangles may be nowcastable depending on prediction_end_date.
    # It will only use those that are nowcastable.
    retro_nowcasts <- baselinenowcast::generate_pt_nowcast_mat_list(
      reporting_triangle_list = retro_rts,
      n = n_history_delay
    )
    disp_params <- baselinenowcast::estimate_dispersion(
      pt_nowcast_mat_list = retro_nowcasts,
      trunc_rep_tri_list = trunc_rts,
      reporting_triangle_list = retro_rts
    )
    nowcast_draws_df <- baselinenowcast::get_nowcast_draws(
      point_nowcast_matrix = pt_nowcast_mat,
      reporting_triangle = rep_tri,
      dispersion = disp_params,
      draws = n_pi_samples
    )


    date_df <- data.frame(
      reference_date = unique(as.Date(train_df_i$reference_date), 'desc')) |>
      dplyr::arrange(reference_date) |>
      dplyr::mutate(time = row_number())
    # Join with the original filtered data
    obs_data <- train_df_i |>
      dplyr::group_by(reference_date) |>
      dplyr::summarise(data_as_of = sum(count, na.rm = TRUE))

    nowcast_w_data <- nowcast_draws_df |>
      dplyr::left_join(date_df, by = "time") |>
      dplyr::left_join(obs_data, by = "reference_date")

    # Bind together the nowcasts for each weekday
    all_nowcasts <- dplyr::bind_rows(all_nowcasts, nowcast_w_data)

  }



  # Order by reference dates
  all_nowcasts <- all_nowcasts |>
    dplyr::arrange(reference_date, "desc") |>
    # only need last 14 days
    dplyr::filter(reference_date >= as.Date(prediction_end_date) - days(14))

  target_data_summarised <- .data |>
    dplyr::rename(
      reference_date = specimen_date
    ) |>
    dplyr::mutate(
      report_date = reference_date + days_to_reported
    ) |>
    epinowcast::enw_filter_report_dates(
      latest_date = as.Date(prediction_end_date) + days(eval_timeframe)
    ) |>
    dplyr::group_by(reference_date)|>
    dplyr::summarise(target = sum(target, na.rm = TRUE)) |># This is the actual target
    dplyr::filter(reference_date <= prediction_end_date) |>
    dplyr::arrange(reference_date,'desc') |>
    dplyr::mutate(reference_date = as.Date(reference_date))

  obs_with_nowcast_draws_df <- all_nowcasts |>
    dplyr::left_join(target_data_summarised, by = "reference_date") |>
    dplyr::rename(.value = pred_count,
                  specimen_date = reference_date,
                  .sample = draw) |>
    dplyr::select(specimen_date, .sample, target, .value, data_as_of) |>
    dplyr::mutate(model = "baselinenowcast_dow")

  nowcast_quantiles <- samples_to_quantiles(
    .sample_predictions = obs_with_nowcast_draws_df,
    remove_identifiers = c()) |>
    dplyr::mutate(t_aggregation = "daily",
                  prediction_end_date = prediction_end_date,
                  model = "baselinenowcast_dow")

}
