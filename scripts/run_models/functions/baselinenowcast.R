
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

  # prepare data
  target_data <- .data |>
    # convert to epinowcast naming conventions
    dplyr::rename(reference_date = specimen_date,
                  new_confirm = target) |>
    dplyr::mutate(report_date = reference_date + days_to_reported) |>
    dplyr::select(-days_to_reported) |>
    # adding this to try handle the problem calculating pobs
    dplyr::filter(!is.na(new_confirm)) |>
    ## apply filters to simulate real-time truncation
    epinowcast::enw_filter_report_dates(latest_date = prediction_end_date) |>
    epinowcast::enw_filter_reference_dates(
      latest_date = prediction_end_date
    ) |>
    # convert to cumulative for `epinowcast` required delay structure
    epinowcast::enw_add_cumulative()

  # again, not really doing anything procesing-wise, keeping for consistent naming with vignette
  observed_data <- epinowcast::enw_filter_report_dates(
    target_data,
    latest_date = prediction_end_date
  )

  # enforce training length
  training_data <- epinowcast::enw_filter_reference_dates(
    observed_data,
    include_days = n_training_volume - 1
  )

  latest_training_data <- epinowcast::enw_latest_data(training_data)

  pobs <- epinowcast::enw_preprocess_data(
    obs = training_data,
    max_delay = max_delay + 1
  )

  reporting_triangle <- dplyr::select(
    pobs$new_confirm[[1]],
    reference_date,
    delay,
    new_confirm
  ) |>
    pivot_wider(names_from = delay, values_from = new_confirm) |>
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
    draws = 100
  )

  latest_data_prepped <- latest_training_data |>
    dplyr::mutate(time = row_number()) |>
    dplyr::rename(obs_confirm = confirm) |>
    dplyr::mutate(reference_date = as.Date(reference_date))

  obs_with_nowcast_draws_df <- nowcast_draws_df |>
    dplyr::left_join(latest_data_prepped, by = "time") |>
    dplyr::rename(.value = pred_count,
                  specimen_date = reference_date,
                  target = obs_confirm,
                  .sample = draw) |>
    dplyr::select(specimen_date, .sample, target, .value) |>
    dplyr::mutate(model = "baselinenowcast")

  nowcast_quantiles <- samples_to_quantiles(
    .sample_predictions = obs_with_nowcast_draws_df,
    remove_identifiers = c()) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("pi_"), ~ .x + target
      )
    )|>
    dplyr::mutate(t_aggregation = "daily",
                  prediction_end_date = prediction_end_date,
                  model = "baselinenowcast")

}
