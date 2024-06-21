
## this function is used to generate predictions from the fit and also out put them in more useable format
#' @param .model the model fit
#' @param .newdata the new data that the model is going to base the predictions off, null If uni-variate
#' @param horizon the forecast horizon a numeric value which shows how many observations ahead will be forecast
#' @param burnin a numeric value of the number of iteration which will be ignored when running the model

bsts_generate_samples <- function(.model, .newdata, horizon, burnin = 0) {

  prediction_object <- predict(
    .model,
    newdata = .newdata,
    horizon = horizon,
    burnin = burnin
  )

  prediction_samples <- prediction_object[["distribution"]] |>
    data.frame() |>
    dplyr::mutate(.sample = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("X"),
      names_to = ".h", # .h and H is an index for the forecast day
      values_to = ".value"
    ) |>
    dplyr::mutate(.h = as.integer(stringr::str_remove(.h, "X")))


}

## used in creating plots that look at the fit when creating bsts models.
#' @param .model the model fit
#' @param .newdata the new data that the model is going to base the predictions off, null If uni-variate
#' @param burnin a numeric value of the number of iteration which will be ignored when running the model
bsts_generate_components <- function(.model, .newdata, burnin = 0) {

  state_component_samples <- .model[["state.contributions"]] |>
    as.data.frame.table()  |>
    dplyr::mutate(
      .sample = as.integer(mcmc.iteration),
      .t = as.integer(time)
    ) |>
    dplyr::filter(.sample > burnin) |>
    dplyr::mutate(.sample = .sample - min(.sample) + 1) |>
    dplyr::rename(.value = Freq) |>
    dplyr::select(-c("time", "mcmc.iteration"))


}

# function for aggregating symptoms by name
#' @param symptom_name the symptom name that is going to aggregated

aggregate_111 <- function(symptom_name) {

  symptom_name <- dplyr::across(starts_with(symptom_name)) |>  rowSums()
}
#' Wrapper function for running the selected model from script
#'
#' Takes data, parameters, and model selection to run over each lookback.
#' This wrapper can also run multiple models in parallel on different lookbacks.
#'
#' @param wd string of working directory, defaults to git root folder.
#' @param model_name String of the function name of the model to be run
#' @param training_data Data frame of the input data.
#' @param prediction_end_dates A sequence of dates defining the start of each look back.
#' @param required_covariates A list of strings of data's column names to keep,
#' generally already listed in the config.
#' @param model_hyperparams A named list of numerical priors for the model,
#' also generally found in the config.
#' @param n_pi_samples Number of posterior samples to calculate
#'
#'
#' @returns Combined model outputs, as a data frame.

run_scripted_model <- function(
    wd = system("echo $(git rev-parse --show-toplevel)", intern = TRUE),
    model_name,
    training_data,
    prediction_end_dates,
    model_formula,
    output_columns,
    model_hyperparams,
    n_pi_samples,
    model_description = ""
    ) {
  model_path <- fs::path(wd, model_hyperparams$model_path)

  if (file.exists(model_path)) {
    suppressWarnings(source(model_path, local = TRUE))
  } else {
    stop(paste(
      "No such model file found. Please check this path for mistakes:",
      model_path
    ))
  }

  print(glue::glue("Fitting {model_name} model: started"))


  model_outputs <- furrr::future_map(
    prediction_end_dates |>
      purrr::set_names(),
    \(prediction_end_date) {
      do.call(
        paste0("run_", stringr::str_remove(model_name, ".R$")),
        args = list(
          .data = training_data,
          model_formula = model_formula,
          n_pi_samples = n_pi_samples,
          prediction_end_date = prediction_end_date,
          output_columns = output_columns,
          model_hyperparams = model_hyperparams,
          model_description = model_description
        )
      )
    },
    .options = furrr::furrr_options(
      seed = NULL,
      packages = c(
        "data.table",
        "dplyr",
        "epinowcast",
        "here",
        "lubridate",
        "purrr",
        "stats"
      )
    )
  )

  print(glue::glue("Fitting {model_name} model: done"))

  return(model_outputs)

}




#' Function to extract formatted dataframes and model objects
#'
#' The function binds the list elements by row.
#'
#' @param model_outputs Individual model output from run_model
#'
#' @returns List with dataframe with formatted predictions and model objects

extract_from_list <- function(model_outputs) {
  # we want to bind all prediction dates, and access model by prediction date
  quantile_predictions <- model_outputs |>
    purrr::map("quantile_predictions") |>
    purrr::list_rbind()

  # Yes okay, we have to iterate twice, but that's a small price to pay
  models <- purrr::map(model_outputs, "model")

  return(list(
    quantile_predictions = quantile_predictions,
    models = models
  ))
}



#' Generate prediction samples from fit `mgcv` and new data using the `gratia` package
#'
#' Thin wrapper for the `gratia::posterior_samples()` function, which gives us
#' model prediction samples.
#'
#' [gratia package version](https://gavinsimpson.github.io/gratia/reference/predicted_samples.html)
#'
#' @param .data Dataframe containing the historic (train) and future (test) data.
#' @param .model Fitted `mgcv` object.
#' @param .n_pi_samples Integer number of samples to generate from model.
#'
#' @returns Matrix with one row per `.data` input, and one column per `.n_pi_samples`.
#'
generate_samples <- function(.data, .model, .n_pi_samples = 1000, .method = "gaussian") {

  results <- gratia::posterior_samples(
    model = .model,
    data = .data,
    n = .n_pi_samples,
    method = .method
  ) |>
    dplyr::rename(
      .sample = .draw,
      .value = .response
    )

  formatted_results <- .data |>
    # generate the row number in the raw data
    dplyr::mutate(row = dplyr::row_number()) |>
    dplyr::left_join(results, by = c("row" = ".row")) |>
    dplyr::select(-row)

  formatted_results

}



#' Take prediction samples and convert to prediction intervals.
#'
#' Essentially a wrapper function of:
#'  - `aggregate_samples()`
#'  - `generate_intervals()`
#'  Which summarizes our prediction samples into communicable formats.
#'
#' @param .sample_predictions Dataframe of one row per sample per set of covariates.
#' @param remove_identifiers Vector of column names which appear in the model predictions, but we
#' do not want to be aggregated by. All covariates not specified in this arguement, or aggregated,
#' will be grouped by. In the case where the data is at the correct aggregation, or no variables
#' need to be removed, a NULL should be provided.
#'
#' @returns Dataframe of one row per covariate group, with quantiles and probabilities of trends.
#'
samples_to_quantiles <- function(
    .sample_predictions,
    remove_identifiers = NULL,
    method = "quantile"
    ) {

  unavailable_identifiers <- setdiff(
    remove_identifiers,
    names(.sample_predictions)
  )

  if (length(unavailable_identifiers) != 0) {
    stop(
      "Columns named in `remove_identifiers` must be present in `.data`!\n",
      "These columns were not available: ",
      paste(unavailable_identifiers, collapse = ", "))
  }

  .sample_predictions |>
    dtplyr::lazy_dt() |>
    aggregate_samples(remove_identifiers = remove_identifiers) |>
    # group by the unique covariates
    dplyr::group_by(dplyr::across(!c(.sample, .value, model))) |>
    # find the summary stats we want
    generate_intervals(method = method) |>
    dplyr::ungroup()
}






#' Standardised set of quantiles for creating prediction intervals.
#'
#' Takes sample data and produces prediction intervals for pre-defined
#' quantiles. Note that quantiles of this form are an assumption. Data must be
#' grouped by whatever covariates *before* calling
#' this function.
#'
#' Alternative quantile approaches are available [see this
#' reference](https://easystats.github.io/bayestestR/articles/credible_interval.html)
#'
#' @param .data Input data frame of processed prediction samples (not directly
#'   from [intervals$generate_samples()]) which contains the column `.value`,
#'   the quantity to be summarized. One row per sample per covariate group.
#' @param method Whether to calculate using standard quantile approach or HDI
#'   (`"quantile"` or `"hdi"`).
#'
#' @returns Data frame with one row per set of groupings, rather than one row
#'   per sample (input)
#'
generate_intervals <- function(.data, method = c("quantile", "hdi")) {

  method <- match.arg(method)

  nested <- switch(
    method,

    # allowing NA's through because this will be applied to empty data where there are
    # no predictions
    "quantile" = dplyr::summarise(
      .data,
      "pi" = list(
        c(0.5, 0.05, 0.95, 0.025, 0.975, 0.25, 0.75, 0.17, 0.83) |>
          stats::quantile(.value, probs = _, na.rm = TRUE) |>
          purrr::set_names(\(x) paste0("pi_", x) |> sub("%$", "", x = _))
      ),
      .groups = "drop"
    ) |>
      dplyr::collect(),


    # unclear whether NA can be passed through
    # unclear what the appropriate central value should be (pi_50) so have gone
    # with median
    "hdi" = dplyr::summarise(
      .data,
      "pi_50" = stats::quantile(.value, probs = 0.5, na.rm = TRUE), # median
      "pi" = list(
        c(0.9, 0.95, 0.5, 0.66) |>
          bayestestR::hdi(.value, ci = _, na.rm = TRUE) |>
          tidyr::pivot_wider(names_from = CI, values_from = !CI, names_vary = "slowest") |>
          purrr::set_names(paste0("pi_", 100 * c(0.05, 0.95, 0.025, 0.975, 0.25, 0.75, 0.17, 0.83)))
      ),
      .groups = "drop"
    ) |>
      dplyr::collect()
  )

  tidyr::unnest_wider(nested, pi)
}




#' Take low-level prediction samples and aggregate them to coarser covariates.
#'
#' Prediction samples are produced at the lowest level possible, which gives the
#' flexibility to aggregate predictions so e.g. higher geographies (Trust -> region).
#' This is prefered to aggregating intervals, as it preserves the uncertainty of the
#' modelled predictions more robustly.
#'
#' This is quite opinionated with what columns it expects within the dataframe, beware.
#' These include: .value, target
#'
#' @param .sample_predictions Dataframe of one row per sample per set of covariates.
#' @param remove_identifiers Vector of column names which appear in the model predictions, but we
#' do not want to be aggregated by. All covariates not specified in this arguement, or aggregated,
#' will be grouped by. In the case where the data is at the correct aggregation, or no variables
#' need to be removed, a NULL should be provided.
#'
#' @returns Dataframe of one row per sample per set of non-removed covariates.
#'
aggregate_samples <- function(
    .sample_predictions,
    remove_identifiers = NULL
    ) {

  .sample_predictions |>
    # calculate the summed up predictions to the level we want
    dplyr::summarise(
      .value = sum(.value, na.rm = TRUE),
      target = if (all(is.na(target))) NA_real_ else sum(target, na.rm = TRUE),
      .by = !dplyr::all_of(c(remove_identifiers, ".value", "target"))
    )
}
