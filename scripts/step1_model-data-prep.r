###project: Snook Activity Patterns
###author(s): MW
###goal(s): set up activity data for analyses
###date(s): May 2026
###note(s): 

# Housekeeping ------------------------------------------------------------

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, scales, ggplot2, suncalc, 
                 dataRetrieval, lunar)

### define custom functions ----
fill_next_nonzero <- function(x) {
      for (i in seq_along(x)) {
            if (!is.na(x[i]) && x[i] == 0) {  # Only process non-NA values
                  # Find the next non-zero, non-NA value
                  next_nonzero <- x[(i+1):length(x)][which(!is.na(x[(i+1):length(x)]) & x[(i+1):length(x)] != 0)[1]]
                  
                  if (!is.na(next_nonzero)) {
                        x[i] <- next_nonzero
                  }
            }
      }
      return(x)
}

### load necessary data ----
### acceleration data and metadata ---
dat <- read_csv('local-data/snook-accelerometer-data.csv') |> 
      filter(distance >= 21)
tags <- dat |> dplyr::select(id, tl, sl, weight, date_tagged) |> distinct()
sites <- dat |> dplyr::select(station, distance, latitude, longitude) |> distinct()

### read in high-res temperature data from bottle creek ---
temp <- read_csv("local-data/botcreek_temp_15mins.csv") |> 
      mutate(datetime_est = as.POSIXct(datetime, format = "%m/%d/%y %H:%M"))

### read in marsh stage data ---
hydro <- readxl::read_xlsx('../1. archive/MAP/data/hydrology/mo215_current.xlsx') |> 
      mutate(date = as.Date(Date)) |> 
      filter(date >= as.Date("2024-01-13") & date <= as.Date("2024-12-14")) |> 
      rename(stage = `Stage (cm)`) |> 
      select(date, stage) |> 
      arrange(date) |> 
      mutate(delta_stage = stage - lag(stage)) |> 
      mutate(delta_stage = case_when(
            is.na(delta_stage) ~ -0.30480000,
            TRUE ~ delta_stage
      ))

### bin the acceleration and detection data hourly to account for serial correlation ----
acc <- dat |> 
      mutate(hour_block = floor_date(datetime, unit = "hour")) |> 
      group_by(id, hour_block) |>                                     
      summarize(mean_acceleration = mean(acceleration, na.rm = TRUE),
                most_visited_station = names(which.max(table(station))),
                n = n(),
                .groups = "drop") |> 
      rename(datetime = hour_block,
             station = most_visited_station) |> 
      filter(n>=3) |> 
      mutate(datetime_eastern = with_tz(datetime, tzone = "America/New_York"),
             date = as.Date(datetime_eastern),
             year = year(datetime_eastern),
             month = month(datetime_eastern),
             day = day(datetime_eastern),
             time = format(datetime_eastern, format = "%H:%M:%S")) |> 
      select(-datetime) |> 
      rename(datetime_est = datetime_eastern) |> 
      select(datetime_est, date, year, month, day, time, everything()) |> 
      left_join(temp) |> 
      left_join(hydro) |> 
      left_join(tags) |> 
      left_join(sites)

# estimate length following tagging --------------------------------------

### set parameters from Taylor et al 2000 paper ---
### Von Bertalanffy growth params ---
L_asymp <- 947.3
K <- 0.175
t0 <- -1.352
### length to weight conv ---
a <- 10^(-5.3564)
b <- 3.1117

bz <- df <- acc |> 
      select(date, id, date_tagged, weight, sl, tl) |> 
      distinct() |> 
      mutate(sl_mm = sl*10,
             sl_cm = sl,
             tl_cm = tl,
             tl_mm = tl*10,
             weight_g = weight,
             weight_kg = weight*0.001) |> 
      select(-sl, -tl, -weight, -date_tagged) |> 
      filter(!is.na(date)) |> 
      mutate(fl_mm = -14.8606 + 0.9512*tl_mm) |> 
      mutate(age_years = t0 -1/K*log(1-fl_mm/L_asymp),
             age_days = age_years*365.25) |> 
      arrange(id, date) |>
      group_by(id) |>
      mutate(first_date = min(date),
             last_date = max(date)) |> 
      ungroup() |>
      select(-sl_cm, -tl_cm, -weight_kg) |> 
      rename(sl_obs = sl_mm, 
             tl_obs = tl_mm, 
             weight_obs = weight_g, 
             fl_est = fl_mm,
             age_at_tagging_days = age_days, 
             age_at_tagging_years = age_years)
glimpse(bz)

bz1 <- bz |> 
      group_by(id) |> 
      summarise(first_date = min(first_date), last_date = max(last_date), .groups = "drop") |> 
      rowwise() |> 
      mutate(date = list(seq.Date(from = first_date, to = last_date, by = "day"))) |> 
      unnest(date)

bz2 <- bz1 |> 
      left_join(bz, by = c('id', 'date', 'first_date', 'last_date')) |> 
      mutate(date_filled = case_when(
            is.na(age_at_tagging_years) ~ "filled",
            TRUE ~ "observed"
      )) |> 
      group_by(id) |> 
      mutate(across(sl_obs:age_at_tagging_days, 
                    ~ifelse(is.na(.),
                            mean(., na.rm = TRUE),
                            .))) |> 
      ungroup() |> 
      group_by(id) |> 
      mutate(days_since_tagged = as.numeric(
            difftime(date, first_date, units = "days")
      )) |> 
      mutate(years_since_tagged = days_since_tagged/365.25) |> 
      ungroup() |> 
      group_by(id) |> 
      mutate(age_days = age_at_tagging_days + days_since_tagged) |> 
      ungroup() |> 
      mutate(age_years = age_at_tagging_years + years_since_tagged) |> 
      mutate(predicted_fl_mm = L_asymp * (1 - exp(-K * (age_years - t0)))) |> 
      mutate(predicted_weight_g = a * (predicted_fl_mm^b)) |> 
      select(date, id, age_years, predicted_fl_mm, predicted_weight_g, date_filled) |> 
      rename(fl_vbf = predicted_fl_mm,
             weight_vbf = predicted_weight_g) |> 
      arrange(id, date) |> 
      group_by(id) |> 
      mutate(daily_growth_g_vbf = weight_vbf - lag(weight_vbf, default = first(weight_vbf))) |> 
      ungroup() |> 
      arrange(id, date) |> 
      filter(date_filled == "observed") |> 
      select(date, id, age_years, fl_vbf, weight_vbf, daily_growth_g_vbf)
glimpse(bz2)

acc1 <- acc |> 
      left_join(bz2, by = c('date', 'id')) |> 
      rename(weight_obs = weight, 
             tl_obs = tl,
             sl_obs = sl) |> 
      select(-date_tagged) |> 
      arrange(id, date) |>
      # handle zeros for growth rate on first day of observation
      # not ultimately used in MS, but retained for posterity
      group_by(id) |> 
      mutate(daily_growth_g_vbf = fill_next_nonzero(daily_growth_g_vbf)) |> 
      ungroup()
glimpse(acc1)

write_csv(acc1, 'local-data/snook-acc-model-data.csv')
