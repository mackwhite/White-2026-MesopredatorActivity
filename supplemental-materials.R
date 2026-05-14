###project: EctothermPredatorActivityDrivers
###author(s): MW, ROS
###goal(s): 
###date(s): February 2025
###note(s): 
# Housekeeping ------------------------------------------------------------

### define custom functions ----
## none

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, suncalc,
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot)

### load necessary data ----
all <- read_csv('local-data/snook-acc-model-data.csv') |>
      mutate(time_og = time,
             time = hour(hms(time_og)),
             y = mean_acceleration) |>
      rename(weight = weight_vbf,
             temp = temp_c) |> 
      mutate(id = as.factor(id), station = as.factor(station)) |>
      dplyr::select(y, time, weight, temp, stage, id, station)
glimpse(all)

### rms conversions ---
slope <- 0.01922
intercept <- 0

og <- read_rds('local-data/archive/accelerometer-model-data-012025.RDS') |> 
      filter(species == "common snook") |> 
      rename(acceleration_raw = acceleration) |> 
      mutate(acceleration = slope*acceleration_raw + intercept) |> 
      rename(stage = marsh_stage,
             y = acceleration) |> 
      mutate(id = as.factor(id),
             station = as.factor(station),
             month = as.factor(month))
      

### supplemental figure one panel a ---
summ <- og |> group_by(id) |> summarize(median = median(y)) |> arrange(median)
print(summ)
# order <- c('5362', '5039', '5043', '5363', '5040', '5361', '5367', '5042')
order <- c('5362', '5039', '5363', '5361', '5042', '5043', '5040', '5367')

b <- og |> 
      mutate(x = id) |> 
      mutate(x = factor(x, levels = order)) |> 
      ggplot(aes(x = x, y = y)) + 
      geom_boxplot(fill = "white", outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(y = expression(bold("Activity (m/s"^2*")"))) +
      # scale_y_continuous(breaks = c(0,50,100,150,200,250), limits = c(0,260)) +
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold')) +
      geom_text(aes(x = max(as.numeric(factor(og$id))), 
                    y = max(og$y, na.rm = TRUE), 
                    label = 'Raw Data'),
                vjust = 0.25, hjust = 1.0, size = 5.0, fontface = "bold", inherit.aes = FALSE)

b

### supplemental figure one panel b ---
summ <- all |> group_by(id) |> summarize(median = median(y)) |> arrange(median)
print(summ)
order <- c('5362', '5039', '5363', '5361', '5042', '5043', '5040', '5367')

a <- all |> 
      mutate(x = id) |> 
      mutate(x = factor(x, levels = order)) |> 
      ggplot(aes(x = x, y = y)) + 
      geom_boxplot(fill = "white", outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      labs(y = expression(bold("Activity (m/s"^2*")"))) +
      # scale_y_continuous(breaks = c(10,30,50,70,90), limits = c(0,90)) +
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold')) +
      geom_text(aes(x = max(as.numeric(factor(all$id))), 
                    y = max(all$y, na.rm = TRUE), 
                    label = 'Data Binned Hourly'),
                vjust = 0.25, hjust = 1.0, size = 5.0, fontface = "bold", inherit.aes = FALSE)
a

### combine a and b plots
s1 <- ggarrange(a, b, ncol = 1, align = "v", labels = c("a", "b"))
### add figure axis titles
s1 |> 
      annotate_figure(
            left = text_grob(label = expression(bold("Activity (m/s"^2*")")),
                             just = 'centre',
                             color = "black",
                             face = "bold",
                             size = 18,
                             rot = 90)
      )

ggsave('output/figs/supplemental-figure-one.png',
       dpi = 600, units= 'in', height = 8, width = 8)

### create supplemental figure 2 with time series of all individual fishes ---
### set up breaks for each month in dataset
month_breaks <- og |> 
      mutate(first_of_month = floor_date(date, "month")) |>   # Get first day of month
      pull(first_of_month) |> 
      unique() |> 
      sort()

### plot it out
s2 <- og |> 
      group_by(id) |> 
      mutate(obs = n()) |> 
      ungroup() |> 
      ggplot(aes(x = date, y = y)) +
      geom_point() +
      facet_wrap(~id, ncol = 1) +
      geom_text(aes(x = max(date), y = 4.5, label = paste("obs =", obs)), 
                hjust = 0.5, size = 3, fontface = "bold") +
      labs(y = expression(bold("Activity (m/s"^2*")")),
           x = "Month") +
      # scale_y_continuous(breaks = c(50,150,250), limits = c(0,260)) +
      scale_x_date(
            breaks = month_breaks,  
            labels = format(month_breaks, "%b")
      ) +
      theme(axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
            axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
            axis.title.y = element_text(size = 16, face = "bold", colour = "black"),
            axis.title.x = element_text(size = 16, face = "bold", colour = "black"),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_line(colour = "grey70", linetype = "dashed"),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(fill = 'lightgrey'),
            strip.text = element_text(size = 14, face = "bold", colour = "black", hjust = 0.5))
s2

ggsave('output/figs/supplemental-figure-two.png',
       dpi = 600, units= 'in', height = 8, width = 8)

### create summary table for taggers ---
total_obs <- tibble(
      id = c(5039,5040,5042,5043,5361,5362,5363,5367),
      obs = c(5313,18422,3565,9355,13696,2756,8876,6177)
      ) |> 
      mutate(id = as.character(id))

all |> 
      mutate(id = as.character(id)) |> 
      group_by(id) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup() |> 
      left_join(total_obs, by = c("id")) |> 
      capture.output(file = "output/tables/t1-tagger-info.csv")

all |> 
      mutate(id = as.character(id)) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup() |> 
      capture.output(file = "output/tables/t2-alltagger-summary-info.csv")

### local summaries

summ_ind <- all |> 
      mutate(id = as.character(id)) |> 
      group_by(id) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            min_acc = min(y),
            max_acc = max(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup() |> 
      left_join(total_obs, by = c("id"))

summ_all <- all |> 
      mutate(id = as.character(id)) |> 
      summarize(
            weight = mean(weight*.001),
            fork_length = mean(fl_vbf/10),
            median_acc = median(y),
            mean_acc = mean(y),
            min_acc = min(y),
            max_acc = max(y),
            sd_acc = sd(y),
            cv_acc = sd(y)/mean(y),
            days_obs = n_distinct(date)) |> 
      ungroup()

glimpse(all)

station_data <- all |> 
      group_by(station, latitude, longitude) |>
      summarize(mean_acc = mean(y, na.rm = TRUE),
                median_acc = median(y, na.rm = TRUE),
                sd_acc = sd(y, na.rm = TRUE),
                n = n(),
                .groups = 'drop')

writexl::write_xlsx(station_data, 'local-data/station-data-for-map06132025.xlsx')

m <- lm(mean_acc ~ log(n), data = station_data)
summary(m)

station_data |> 
      ggplot(aes(x= log(n), y= mean_acc)) +
      geom_point() +
      geom_abline()

station_order <- all |> 
      distinct(station, longitude) |> 
      arrange(longitude) |> 
      pull(station)

all |> 
      mutate(station = factor(station, levels = station_order)) |> 
      ggplot(aes(x = station, y = y)) +
      geom_boxplot(aes(color = station), alpha = 0.3, outliers = FALSE) +
      geom_jitter(aes(color = station), alpha = 0.7)

# supplemental temp and stage figure --------------------------------------

all <- read_csv('local-data/snook-model-data-march.csv') |> 
      mutate(time_og = time,
             time = time_num,
             y = mean_hr_acc) |> 
      rename(weight = weight_vbf, 
             stage = stage_d,
             temp = temp_hr) |> 
      mutate(id = as.factor(id),
             station = as.factor(station),
             month = as.factor(month))
summary(all)

temp_hr <- read_csv("local-data/botcreek_temp_15mins.csv")
temphr1 <- temp_hr |> 
      mutate(datetime_est = as.POSIXct(datetime, format = "%m/%d/%y %H:%M")) |> 
      filter(datetime_est >= as.Date("2024-01-14") & datetime_est <= as.Date("2024-12-14"))

temp <- temphr1 |> 
      mutate(date = as.Date(datetime_est)) |> 
      group_by(date) |> 
      summarise(
            # min_temp = min(temp_c, na.rm = TRUE),
            # max_temp = max(temp_c, na.rm = TRUE),
            temp = mean(temp_c, na.rm = TRUE),
            # sd_temp = sd(temp_c, na.rm = TRUE)
      )

stage <- readxl::read_xlsx('../1. archive/MAP/data/hydrology/mo215_current.xlsx') |> 
      mutate(date = as.Date(Date)) |> 
      filter(date >= as.Date("2024-01-14") & date <= as.Date("2024-12-14")) |> 
      rename(stage = `Stage (cm)`) |> 
      select(date, stage) |> 
      arrange(date)

env <- temp |> left_join(stage)
glimpse(env)

env_long <- env |> 
      pivot_longer(cols = c('temp', 'stage'),
                   names_to = 'metric',
                   values_to = 'value')

keep <- c("env_long")
rm(list = setdiff(ls(), keep))

env_wide <- env_long |>
      filter(metric %in% c("temp", "stage")) |>
      select(date, metric, value) |>
      pivot_wider(names_from = metric, values_from = value)

### defining range of data for transformations ----
min_temp <- min(env_wide$temp, na.rm = TRUE)
max_temp <- max(env_wide$temp, na.rm = TRUE)
min_stage <- min(env_wide$stage, na.rm = TRUE)
max_stage <- max(env_wide$stage, na.rm = TRUE)

### function to scale marsh stage to the range of temperature ----
scale_stage <- function(x) {
      (x - min_stage) / (max_stage - min_stage) * (max_temp - min_temp) + min_temp
}

### inverse of above function to rescale seconday y axis with real values ----
inv_scale_stage <- function(x) {
      (x - min_temp) / (max_temp - min_temp) * (max_stage - min_stage) + min_stage
}

### scale data for secondary y axis using function from above ----
env_wide <- env_wide |> 
      mutate(stage_scaled = scale_stage(stage))

### grab monthly breaks for below plot ----
month_breaks <- env_wide |> 
      mutate(first_of_month = floor_date(date, "month")) |>   # Get first day of month
      pull(first_of_month) |> 
      unique() |> 
      sort()

### set y axis breaks - manually because of the two axes ----
y_primary_limits <- c(16, 33)
y_primary_breaks <- seq(16, 32, by = 4)

# Plot with second axis
env_wide |> 
ggplot(aes(x = date)) +
      geom_line(aes(y = temp, color = "Temperature"), linetype = 'solid', linewidth = 1.25) +
      geom_line(aes(y = stage_scaled, color = "Marsh Stage"), linetype = 'solid', linewidth = 1.25) +
      scale_y_continuous(
            name = "Temperature (°C)",
            limits = y_primary_limits,
            breaks = y_primary_breaks,
            sec.axis = sec_axis(
                  trans = ~ inv_scale_stage(.),
                  name = "Marsh Stage (cm)"
            )
      ) +
      scale_color_manual(values = c("Temperature" = "#e41a1c", 
                                    "Marsh Stage" = "#377eb8")) +
      scale_x_date(
            breaks = month_breaks,  
            labels = format(month_breaks, "%b")
      ) +
      labs(x = "Month", color = NULL) +
      theme(
            axis.text.x = element_text(size = 14, face = "bold", colour = "black", angle = 45, vjust = 0.8), 
            axis.text.y = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
            axis.title.y.left = element_text(size = 15, face = "bold", colour = "#e41a1c"),
            axis.title.y.right = element_text(size = 15, face = "bold", colour = "#377eb8"),
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            # legend.position = c(0.18,0.90),
            legend.position = 'none',
            legend.text = element_text(size = 10, color = "black", face = 'bold'),
            legend.title = element_text(size = 10, color = "black", face = 'bold')
      )

ggsave('output/figs/supplemental-figure-3.png',
       dpi = 600, units= 'in', bg = 'white', height = 5, width = 6)
