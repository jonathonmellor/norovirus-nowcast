# # # # # # # # # # # #
####    SETUP     ####
# # # # # # # # # # # #

wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)
source("./scripts/depends.R")
source("./scripts/run_models/functions/model_running_functions.R")
source(paste0(wd, "/scripts/run_models/functions/plotting.R"))
source(paste0(wd, "/scripts/run_models/functions/scoring.R"))
source(paste0(wd, "/scripts/run_models/functions/baselinenowcast.R"))


# we want to use the most recent versions of the packages on GitHub
remotes::install_github(repo = "epinowcast/baselinenowcast")
remotes::install_github(repo = "epinowcast/epinowcast")

library(ggplot2)
library(baselinenowcast)
library(epinowcast)
library(dplyr)
library(lubridate)
library(glue)



training_data_path <- "./outputs/data/cases_with_noise.csv"
training_data_raw <- vroom::vroom(training_data_path)

training_data <- training_data_raw |>
  mutate(
    report_date = specimen_date + days(days_to_reported)
  ) |>
  rename(count = target,
         reference_date = specimen_date) |>
  select(reference_date, report_date, count)

# Specifications
n_draws <- 100
nowcast_date <- "2024-01-21"
max_delay <- 14
n_history_delay <- 28 # We are going to use in this case, 28 x 7 total
# reference dates
# we don't actually have enough here
n_retrospective_nowcasts <- 28
filter_ref_dates <- TRUE

all_nowcasts <- data.frame()
for( i in 1:7){
  train_df_i <- enw_filter_report_dates(
    obs = training_data[wday(training_data$reference_date) == i,],
    latest_date = nowcast_date
  ) |>
    mutate(reference_date = as.Date(reference_date))
  rep_tri <- train_df_i |>
    mutate(delay = as.integer(difftime(report_date,
                                       reference_date,
                                       units = "days"
    ))
    ) |>
    select(reference_date, delay, count) |>
    filter(delay <= max_delay, delay >= 0) |>
    pivot_wider(
      names_from = delay,
      values_from = count
    ) |> select(-reference_date) |> as.matrix()

  delay_pmf <- get_delay_estimate(
    reporting_triangle = rep_tri,
    max_delay = max_delay,
    n = n_history_delay
  )
  pt_nowcast_mat <- apply_delay(
    rep_tri_to_nowcast = rep_tri,
    delay_pmf = delay_pmf
  )

  # This will throw a warning because not all triangles can be nowcasted, in
  # practince n_retrospective_nowcasts will be less than 28
 trunc_rts <- truncate_triangles(
    reporting_triangle = rep_tri,
    n = n_retrospective_nowcasts
  )

 retro_rts <- generate_triangles(
   trunc_rep_tri_list = trunc_rts,
   structure = c(1, 7)
 )
 # This is going to throw a bunch of warnings because only the first 15
 # triangles are actually nowcastable. It will use those 15.
 retro_nowcasts <- generate_pt_nowcast_mat_list(
   reporting_triangle_list = retro_rts,
   n = n_history_delay
 )
 disp_params <- estimate_dispersion(
   pt_nowcast_mat_list = retro_nowcasts,
   trunc_rep_tri_list = trunc_rts,
   reporting_triangle_list = retro_rts
 )
 nowcast_draws_df <- get_nowcast_draws(
   point_nowcast_matrix = pt_nowcast_mat,
   reporting_triangle = rep_tri,
   dispersion = disp_params,
   draws = n_draws
 )


 date_df <- data.frame(
   reference_date = unique(as.Date(train_df_i$reference_date), 'desc')) |>
   arrange(reference_date) |>
   mutate(time = row_number())
 # Join with the original filtered data
 obs_data <- train_df_i |>
   group_by(reference_date) |>
   summarise(data_as_of = sum(count, na.rm = TRUE))

 nowcast_w_data <- nowcast_draws_df |>
   left_join(date_df, by = "time") |>
   left_join(obs_data, by = "reference_date")

 # Bind together the nowcasts for each weekday
  all_nowcasts <- bind_rows(all_nowcasts, nowcast_w_data)

}

# Order by reference dates
all_nowcasts <- all_nowcasts |>
  arrange(reference_date, "desc")

# Quick plot of nowcasts compared to data as of the nowcast date
ggplot(all_nowcasts |> filter(reference_date >= as.Date(nowcast_date) - days(14))) +
  geom_line(aes(x = reference_date, y = pred_count, group = draw),
            linewidth = 0.2, alpha = 0.2) +
  geom_line(aes(x = reference_date, y = data_as_of), color = "magenta4") +
  xlab("Specimen date") +
  ylab("Positive test count")+
  theme_bw() +
  ggtitle(glue("Nowcast vs data observed as of the nowcast date for {nowcast_date}"))
