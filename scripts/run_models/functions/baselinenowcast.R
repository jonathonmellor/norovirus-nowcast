
run_baselinenowcast <- function(.data,
                                prediction_end_date,
                                n_pi_samples = 500,
                                model_hyperparams,
                                ...) {

  # parameters
  max_delay <- model_hyperparams$max_delay
  n_training_volume <- model_hyperparams$training_length

  n_history_delay <- model_hyperparams$n_history_delay
  n_retrospective_nowcasts <- model_hyperparams$n_retrospective_nowcasts
  eval_timeframe <- model_hyperparams$eval_timeframe
  
  # Ignoring the vignette and epinowcast preprocessing for the most part
  training_data <- .data |>
    rename(
      reference_date = specimen_date
    ) |>
    mutate(
      report_date = reference_date + days_to_reported
      ) |>
    epinowcast::enw_filter_report_dates(
    latest_date = prediction_end_date
  ) 
  reporting_triangle <- training_data |>
    # Make reporting triangle by hand
    select(reference_date, days_to_reported, target) |>
    filter(days_to_reported <= max_delay, days_to_reported >= 0) |>
    pivot_wider(
      names_from = days_to_reported,
      values_from = target
    ) |>
    select(-reference_date) |>
    as.matrix()

 
  delay_pmf <- baselinenowcast::get_delay_estimate(
    reporting_triangle = reporting_triangle,
    max_delay = max_delay,
    n = n_history_delay
  )

  point_nowcast_matrix <- baselinenowcast::apply_delay(
    rep_tri_to_nowcast = reporting_triangle,
    delay_pmf = delay_pmf
  )

  trunc_rep_tri_list <- baselinenowcast::truncate_triangles(reporting_triangle,
                                                            n = n_retrospective_nowcasts
  )
  retro_rep_tri_list <- baselinenowcast::generate_triangles(trunc_rep_tri_list)

  retro_pt_nowcast_mat_list <- baselinenowcast::generate_pt_nowcast_mat_list(
    reporting_triangle_list = retro_rep_tri_list,
    n = n_history_delay
  )

  disp_params <- baselinenowcast::estimate_dispersion(
    pt_nowcast_mat_list = retro_pt_nowcast_mat_list,
    trunc_rep_tri_list = trunc_rep_tri_list,
    reporting_triangle_list = retro_rep_tri_list,
    n = n_retrospective_nowcasts
  )

  nowcast_draws_df <- baselinenowcast::get_nowcast_draws(
    point_nowcast_matrix, reporting_triangle,
    dispersion = disp_params,
    draws = n_pi_samples
  )

  training_data_summarised<- training_data |>
    dplyr::mutate(reference_date = as.Date(reference_date)) |>
    group_by(reference_date) |>
    summarise(data_as_of = sum(target, na.rm = TRUE)) 
  # This isn't a target
    # its the data as of the nowcast date summed across the delays that have
    # been observed
  
  target_data_summarised <- .data |>
    rename(
      reference_date = specimen_date
    ) |>
    mutate(
      report_date = reference_date + days_to_reported
    ) |>
    epinowcast::enw_filter_report_dates(
      latest_date = prediction_end_date + days(eval_timeframe)
    ) |> group_by(reference_date)|>
    summarise(target = sum(target, na.rm = TRUE)) |># This is the actual target
   filter(reference_date <= prediction_end_date) |>
    arrange(reference_date,'desc') |>
    mutate(
      time = row_number()
    )
    

  obs_with_nowcast_draws_df <- nowcast_draws_df |>
    dplyr::left_join(target_data_summarised, by = "time") |>
    dplyr::mutate(reference_date = as.Date(reference_date))|>
    dplyr::left_join(training_data_summarised, by = "reference_date") |>
    dplyr::rename(.value = pred_count,
                  specimen_date = reference_date,
                  .sample = draw) |>
    dplyr::select(specimen_date, .sample, target, .value, data_as_of) |>
    dplyr::mutate(model = "baselinenowcast")

  nowcast_quantiles <- samples_to_quantiles(
    .sample_predictions = obs_with_nowcast_draws_df,
    remove_identifiers = c()) |>
    # Remove this step because observations + predicted nowcast draws already
    # happened in `get_nowcast_draws()`
    # dplyr::mutate(
    #   dplyr::across(
    #     dplyr::starts_with("pi_"), ~ .x + target
    #   )
    # )|>
    dplyr::mutate(t_aggregation = "daily",
                  prediction_end_date = prediction_end_date,
                  model = "baselinenowcast")

}
