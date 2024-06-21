# Epinowcast function

run_epinowcast <- function(.data,
                           model_formula,
                           n_pi_samples, # NB this isn't used here but needed to run with run_scripted_models
                           prediction_end_date,
                           output_columns,
                           model_hyperparams,
                           model_description) {
  # Load samples_to_quantiles() within the context of this function
  # (to avoid relying on something which is hopefully, maybe, in the global environment)
  samples_to_quantiles <- local({
    # Source within this local() environment
    source(
      here::here("./scripts/run_models/functions/model_running_functions.R"),
      local = TRUE
    )
    # Return the object we are interested in
    samples_to_quantiles
  })

  # Preprocess data ----------------------------------------------------------------------------

  model_data <- .data |>
    dplyr::filter(
      report_date >= prediction_end_date - model_hyperparams$training_length,
      report_date <= prediction_end_date
    ) |>
    epinowcast::enw_complete_dates()

  triangle_data <- epinowcast::enw_preprocess_data(
    model_data,
    max_delay = model_hyperparams$max_delay)

  # SET PRIORS ####

  # capture the defaults
  default_priors <- epinowcast::enw_reference(data = triangle_data)

  # save to update them
  new_priors <- default_priors$priors

  # reduce the sd on all parameters
  new_priors[, sd := 0.2]
  # move normal from 1 to log(2) to best represent our prior knowledge of delay
  new_priors[variable == "refp_mean_int", mean := log(3)]
  new_priors[variable == "refp_sd_int", mean := 0.3 * sqrt(log(5) - log(1))]



  # Set up epinowcast model --------------------------------------------------------------------

  # Compile in a new temporary directory - otherwise compiled object is saved in tempdir(), meaning
  # that if we have multiple parallel processes we end up overwriting the same object multiple
  # times in tempdir()!
  # TODO can we compile once, rather than in all parallel processes?
  tmpd <- withr::local_tempdir()
  model <- epinowcast::enw_model(threads = TRUE, target_dir = tmpd)

  fit <- epinowcast::enw_fit_opts(
    save_warmup = FALSE, output_loglik = TRUE, pp = FALSE,
    chains = model_hyperparams$chains,
    parallel_chains = model_hyperparams$parallel_chains,
    threads_per_chain = 2, # our EC2 instances have max 2 threads per core
    iter_sampling = model_hyperparams$n_iterations,
    iter_warmup = model_hyperparams$n_warmup,
    show_messages = FALSE, refresh = 0,
    adapt_delta = 0.98, max_treedepth = 15
  ) # Can set pp = TRUE to generate posterior predictions for the observed data (useful to evaluate performance)

  # Defines the generative process for the expected incidence

  expectation_module <- epinowcast::enw_expectation(
    # Default here is a daily random effect, also adding a day of the week effect due to periodicity in data
    r = ~ 0 + (1 | day) + (1 | day_of_week),
    data = triangle_data,
    # weighting to apply to current and past expected observations (default is 1)
    latent_reporting_delay = 1
  )

  # Formula describing the non-parametric logit hazard model

  # NOTE: I have removed the day of week random effect from the hazard model
  # in some tuning instances this caused intractable runtimes (particularly for long)
  # delays and training lengths. I believe this is because of having two day of week
  # effects, which the model would struggle to fit. Jon
  report_module <- epinowcast::enw_report(~1, data = triangle_data)

  # Not included at the moment due to convergence issues (also expect less day of the week effect in specimen date)
  # reference_module <- enw_reference(~ (1 | day_of_week), data = triangle_data) # nolint: commented_code_linter

  # Run nowcast ------------------------------------------------------------------

  # Note default observation model in the likelihood is negative binomial
  nowcast <- epinowcast::epinowcast(triangle_data,
    expectation = expectation_module,
    report = report_module,
    fit = fit,
    priors = new_priors,
    # obs = epinowcast::enw_obs(family = "poisson", data = triangle_data), # nolint: commented_code_linter
    model = model)

  cat("Finished nowcasting. Compilation and fitting took:", nowcast$run_time, "\n")
  cat("Max Rhat is:", nowcast$max_rhat, "\n")
  cat("Proportion divergent transitions:", nowcast$per_divergent_transitions, "\n")

  # Plot to check
  # plot(nowcast, type = "nowcast", latest_obs = .data |> filter(reference_date >= as.Date("2023-08-01"))) # nolint: commented_code_linter

  # Post-processing --------------------------------------------------------------

  epinowcast_samples <- epinowcast::enw_nowcast_samples(nowcast$fit[[1]], nowcast$latest[[1]])
  epinowcast_posteriors <- summary(nowcast, type = "fit")
  delay_cdf_pdf <- extract_epinowcast_cdf(nowcast, model_hyperparams$max_delay)

  epinowcast_daily_predictions <- epinowcast::enw_nowcast_summary(nowcast$fit[[1]], nowcast$latest[[1]],
    probs = c(0.025, 0.05, 0.17, 0.25, 0.5, 0.75, 0.83, 0.95, 0.975)) |>
    dplyr::rename(specimen_date = reference_date, target = confirm) |>
    dplyr::rename_with(~ tolower(gsub("q", "pi_", .x))) |>
    # NB epinowcast will predict up to the maximum report date wheras other model predicts up to maximum specimen date
    # so filtering here to harmonize outputs
    dplyr::filter(specimen_date <= prediction_end_date) |>
    dplyr::mutate(prediction_end_date = prediction_end_date) |>
    dplyr::select(specimen_date, prediction_end_date, dplyr::starts_with("pi"), target) |>
    dplyr::mutate(t_aggregation = "daily")

  # Aggregate samples to weekly level
  epinowcast_weekly_samples <- epinowcast_samples |>
    dplyr::mutate(week_starting = lubridate::floor_date(reference_date, "weeks", week_start = 1)) |>
    dplyr::group_by(week_starting, .draw) |>
    dplyr::summarise(sample = sum(sample), confirm = sum(confirm)) |>
    dplyr::ungroup() |>
    dplyr::rename(.sample = .draw, .value = sample, target = confirm) |>
    dplyr::mutate(model = "epinowcast")

  # Calculate quantiles for weekly predictions
  epinowcast_weekly_predictions <- samples_to_quantiles(
    .sample_predictions = epinowcast_weekly_samples
  ) |>
    dplyr::rename(specimen_date = week_starting) |>
    dplyr::mutate(
      specimen_date = data.table::as.IDate(specimen_date),
      prediction_end_date = prediction_end_date,
      t_aggregation = "weekly"
    )

  # Bind together daily and weekly predictions
  epinowcast_predictions <- dplyr::bind_rows(epinowcast_daily_predictions, epinowcast_weekly_predictions) |>
    dplyr::mutate(model = "epinowcast")

  return(
    list(
      quantile_predictions = epinowcast_predictions,
      model = nowcast$fit[[1]],
      delay_cdf_pdf = delay_cdf_pdf
  ))
}

# Function to extract the posterior distribution for the reporting delay from specimen date to report date
# Adapted from https://package.epinowcast.org/articles/single-timeseries-rt-estimation.html

extract_epinowcast_cdf <- function(nowcast, max_delay) {
  # Extract samples for mean and sd of delay distribution
  draws <- nowcast |>
    (\(x) {
      x$fit[[1]]$draws(variables = c("refp_mean", "refp_sd"), format = "df")
    })() |>
    data.table::as.data.table()

  # Get cdf for delay distribution (for reporting delay from 1 to max_delay)
  draws[
    ,
    cdf := purrr::map2(
      `refp_mean[1]`, `refp_sd[1]`,
      ~ data.table::data.table(
        delay = 1:max_delay, cdf = stats::plnorm(1:max_delay, .x, .y) / stats::plnorm(max_delay, .x, .y)
      )
    )
  ]
  draws <- draws[, data.table::rbindlist(cdf, idcol = "row_id")]
  # Calculate pdf from cdf
  draws <- draws |>
    dplyr::group_by(row_id) |>
    dplyr::mutate(pdf = cdf - dplyr::lag(cdf)) |>
    dplyr::mutate(pdf = dplyr::case_when(is.na(pdf) ~ cdf,
      TRUE ~ pdf)) |>
    dplyr::ungroup() |>
    data.table::as.data.table()
  # Calculate quantiles from samples for cdf
  cdf_draws <- draws[,
    .(
      mean = mean(cdf),
      pi_2.5 = stats::quantile(cdf, probs = 0.025),
      pi_97.5 = stats::quantile(cdf, probs = 0.975)
    ),
    by = "delay"
  ]
  # Calculate quantiles from samples for pdf
  pdf_draws <- draws[,
    .(
      mean = mean(pdf),
      pi_2.5 = stats::quantile(pdf, probs = 0.025),
      pi_97.5 = stats::quantile(pdf, probs = 0.975)
    ),
    by = "delay"
  ]

  list(
    pdf = as.data.frame(pdf_draws),
    cdf = as.data.frame(cdf_draws)
  )
}
