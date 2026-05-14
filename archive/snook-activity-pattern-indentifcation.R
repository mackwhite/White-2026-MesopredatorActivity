###project: mesoconsumer behaviors
###author(s): MW
###goal(s): quantify effects of temerature on diel activity patterns of snook
###date(s): February 2025
###note(s): 
# Housekeeping ------------------------------------------------------------

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, suncalc,
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot)

all <- read_csv('local-data/snook-model-data-march.csv') |> 
      mutate(time_og = time,
             time = time_num,
             y = mean_hr_acc,
             total_dens = total_dens + 1) |> 
      rename(weight = weight_vbf, 
             stage = stage_d,
             temp = temp_hr) |> 
      mutate(id = as.factor(id),
             station = as.factor(station),
             month = as.factor(month))
glimpse(all)
hist(all$y)

# fit glmms with trigonometric terms --------------------------------------
trig_null <- glmmTMB(y ~ 1 + (1|id) + (1|station), 
                     family = gaussian(link = "log"),
                     data = all)

trig_unimodal <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + (1|id) + (1|station),
                         family = gaussian(link = "log"),
                         data = all)

trig_bimodal <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + cos(2*pi*time/12) + sin(2*pi*time/12) + (1|id) + (1|station),
                        family = gaussian(link = "log"),
                        data = all)

performance::compare_performance(trig_null, trig_unimodal, trig_bimodal) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
      capture.output(file = "output/tables/q1-trig-glmm-model-comparison.csv")

# summary(trig_bimodal) |> 
#       capture.output(file = "output/tables/q1-bestfit-glmm-summary.xlsx")
r.squaredGLMM(trig_bimodal)
performance(trig_bimodal)

new_data <- expand.grid(
      time = seq(from = 0, to = 23, by = 0.1)) |> 
      mutate(
            cos_24 = cos(2 * pi * time / 24),
            sin_24 = sin(2 * pi * time / 24),
            cos_12 = cos(2 * pi * time / 12),
            sin_12 = sin(2 * pi * time / 12),
            id = NA,  
            station = NA,
            month = NA
      )

### use best-fit model to predict new data ---
pred <- predict(trig_bimodal, newdata = new_data, se.fit = TRUE)
new_data$fit <- pred$fit
new_data$se.fit <- pred$se.fit

pred_fit <- new_data |> 
      mutate(predicted_log = pred$fit,
             predicted = exp(predicted_log),
             lower_log = predicted_log - pred$se.fit,
             upper_log = predicted_log + pred$se.fit) |>
      mutate(lower = exp(lower_log),
             upper = exp(upper_log)) |>
      rename(x = time,
             y = predicted)

glmm_fit <- pred_fit |> mutate(model = "glmm") |> 
      select(model, x, y, lower, upper)

bfm_glmm <- trig_bimodal

### run gamms now ---
snook_gamm <- mgcv::gam(y ~ s(time, bs="cc") + s(id, bs="re")+ s(station, bs="re"),
                        family = gaussian(link = 'log'),
                        data = all)
snook_gamm_vis <- visreg(snook_gamm, "time", type = "conditional", scale = "response")
gamm_fit <- snook_gamm_vis$fit |> 
      rename(predicted = visregFit,
             lower = visregLwr,
             upper = visregUpr) |> 
      mutate(model = "gamm") |> 
      select(-y) |> 
      rename(x = time,
             y = predicted) |> 
      select(model, x, y, lower, upper)
glimpse(gamm_fit)
glimpse(glmm_fit)

### compare performance of each model ---
compare_performance(bfm_glmm, snook_gamm) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
      capture.output(file = "output/tables/q1-gamm-glmm-comparison.csv")

summary(snook_gamm)
performance(snook_gamm)
# summary.gam(snook_gamm) |> 
#       capture.output(file = "output/tables/q1-gamm-summary.xlsx")

### join models together ---
snook_fit <- rbind(glmm_fit, gamm_fit) |> 
      mutate(model = as.factor(model))

### clean up environment
keep <- c("all", "snook_fit", 'gamm_fit')
rm(list = setdiff(ls(), keep))

a <- snook_fit |> 
      rename(Model = model) |> 
      mutate(Model = case_when(
            Model == "gamm" ~ "GAMM",
            Model == "glmm" ~ "GLMM"
      )) |> 
      ggplot(aes(x = x)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.4) +
      geom_line(aes(y = y, linetype = Model), linewidth = 2) +
      theme_bw() +
      labs(x = "Time of Day (h)", y = expression(bold("Predicted Activity (m/s"^2*")")),
           fill = "Model",
           color = "Model",
           linetype = "Model") +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(breaks = c(0.35,0.40,0.45,0.50,0.55,0.60), limits = c(0.32,0.6251)) +
      scale_fill_manual(values = c("GLMM" = "#1f78b4", "GAMM" = "#33a02c")) +
      scale_color_manual(values = c("GLMM" = "#1f78b4", "GAMM" = "#33a02c")) +
      scale_linetype_manual(values = c("GLMM" = "solid", "GAMM" = "dashed")) +
      theme(axis.text = element_text(size = 10, face = "bold", colour = "black"),
            axis.title = element_text(size = 12, face = "bold", colour = "black"),
            plot.title = element_text(size = 10, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            # panel.background = element_blank(),
            legend.position = c(0.13,0.85),
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 8, color = "black"))
a

# ggsave('output/figs/figure-two.png',
#        dpi = 600, units= 'in', height = 4, width = 4)

temp <- read_csv("local-data/botcreek_temp_15mins.csv") |>
      mutate(
            datetime_est = as.POSIXct(datetime, format = "%m/%d/%y %H:%M"),
            date = as.Date(datetime_est),
            time = as.numeric(hour(datetime_est))
            # y = temp_c
            # y = scale(temp_c, center = TRUE, scale = TRUE)[,1]
      ) |>
      group_by(date) |>
      mutate(
            y = temp_c - mean(temp_c, na.rm = TRUE)
      ) |>
      ungroup() |>
      select(y, time)

temp_summary <- temp |> 
      rename(x = time) |> 
      group_by(x) |> 
      summarize(
            ymean = mean(y, na.rm = TRUE),
            mean_temp = mean(y, na.rm = TRUE),
            median_temp = median(y, na.rm = TRUE),
            sd_temp = sd(y, na.rm = TRUE),
            se_temp = sd_temp / sqrt(n()),
            # upper = mean_temp + sd_temp,
            # lower = mean_temp - sd_temp,
            upper = mean_temp + se_temp,
            lower = mean_temp - se_temp,
            .groups = 'drop') |> 
      rename(y = ymean)

temp_summary_clean <- temp_summary |> 
      mutate(across(everything(), as.numeric))

new_x <- tibble(x = seq(0, 23, by = 0.01))

vars_to_interp <- c("y", "mean_temp", "median_temp", "sd_temp", "se_temp", 
                    "upper", "lower")

interp_values <- map_dfc(vars_to_interp, function(var) {
      interp_result <- approx(
            x = temp_summary_clean$x,
            y = temp_summary_clean[[var]],
            xout = new_x$x
      )
      tibble(!!var := interp_result$y)
})

interp_df <- bind_cols(tibble(x = new_x$x), interp_values)
glimpse(interp_df)
# Smooth mean temperature
temp_spline <- smooth.spline(x = interp_df$x, y = interp_df$mean_temp, spar = 0.6)
# Smooth upper and lower bounds
upper_spline <- smooth.spline(x = interp_df$x, y = interp_df$upper, spar = 0.6)
lower_spline <- smooth.spline(x = interp_df$x, y = interp_df$lower, spar = 0.6)

# Assemble into a tibble
temp_smooth_df <- tibble(
      x = temp_spline$x,
      y = temp_spline$y,
      upper = upper_spline$y,
      lower = lower_spline$y
)

temp_summary1 <- temp_smooth_df
# glimpse(crepusc_summary)
glimpse(temp_summary1)
glimpse(gamm_fit)

min_gamm <- min(gamm_fit$lower, na.rm = TRUE)
max_gamm <- max(gamm_fit$upper, na.rm = TRUE)
min_temp <- min(temp_summary1$lower, na.rm = TRUE)
max_temp <- max(temp_summary1$upper, na.rm = TRUE)

scale_temp <- function(x) {
      (x - min_temp) / (max_temp - min_temp) * (max_gamm - min_gamm) + min_gamm
}
inv_scale_temp <- function(x) {
      (x - min_gamm) / (max_gamm - min_gamm) * (max_temp - min_temp) + min_temp
}

temp_summary_scaled <- temp_summary1 |> 
      mutate(
            y_scaled = scale_temp(y),
            lower_scaled = scale_temp(lower),
            upper_scaled = scale_temp(upper)
      )

lat <- 25.458015
long <- -80.875852
dates <- seq(as.Date("2024-01-13"), as.Date("2024-12-14"), by = "day")

sun_time <- getSunlightTimes(
      date = dates,
      lat = lat,
      lon = long,
      keep = c('sunsetStart', 'night',
               'nightEnd','sunriseEnd'),
      tz = 'EST'
)

sun_time_summary <- sun_time %>%
      mutate(
            s_morning = hour(nightEnd) + minute(nightEnd)/60 + second(nightEnd)/3600,
            e_morning = hour(sunriseEnd) + minute(sunriseEnd)/60 + second(sunriseEnd)/3600,
            s_night= hour(sunsetStart) + minute(sunsetStart)/60 + second(sunsetStart)/3600,
            e_night = hour(night) + minute(night)/60 + second(night)/3600
      ) |> 
      summarise(
            dusk_min = min(s_night),
            dusk_max = max(e_night),
            dawn_min = min(s_morning),
            dawn_max = max(e_morning)
      )
print(sun_time_summary)

dawn_min <- 4.068056
dawn_max <- 7.249444
dusk_min <- 17.52556
dusk_max <- 20.81667
      
ggplot() +
      geom_rect(aes(xmin = dawn_min, xmax = dawn_max,
                    ymin = -Inf, ymax = Inf),
                fill = "grey70", alpha = 0.25) +
      geom_rect(aes(xmin = dusk_min, xmax = dusk_max,
                    ymin = -Inf, ymax = Inf),
                fill = "grey70", alpha = 0.25) +
      # geom_vline(xintercept = dawn_min, color = "grey70", size = 1, linetype = 1) + 
      # geom_vline(xintercept = dawn_max, color = "grey70", size = 1, linetype = 1) + 
      # geom_vline(xintercept = dusk_min, color = "grey70", size = 1, linetype = 1) + 
      # geom_vline(xintercept = dusk_max, color = "grey70", size = 1, linetype = 1) + 
      # geom_ribbon(data = temp_summary_scaled, aes(x = x, ymin = lower_scaled, ymax = upper_scaled), 
      #             fill = "#6a3d9a", alpha = 0.4) +
      geom_ribbon(data = temp_summary_scaled, aes(x = x, ymin = lower_scaled, ymax = upper_scaled),
                  fill = "#6a3d9a", alpha = 0.4) +
      geom_line(data = temp_summary_scaled, aes(x = x, y = y_scaled), 
                color = "#6a3d9a", linewidth = 2) +
      geom_ribbon(data = gamm_fit, aes(x = x, ymin = lower, ymax = upper), 
                  fill = "#33a02c", alpha = 0.4) +
      geom_line(data = gamm_fit, aes(x = x, y = y), 
                color = "black", linewidth = 2, linetype = 'dashed') +
      scale_y_continuous(
            name = "Predicted Activity (m/s²)",
            sec.axis = sec_axis(trans = ~inv_scale_temp(.), name = "Change in Temperature (°C)")
      ) +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      labs(x = "Time of Day (h)") +
      theme(axis.text = element_text(size = 10, face = "bold", colour = "black"),
            axis.title = element_text(size = 12, face = "bold", colour = "black"),
            axis.title.y.left = element_text(color = "#33a02c", face = "bold"),
            axis.title.y.right = element_text(color = "#6a3d9a", face = "bold"),
            plot.title = element_text(size = 10, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "right",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold'))

ggsave('output/figs/supplemental-figure-x-temptimeofday-activity.png',
       dpi = 600, units= 'in', height = 4, width = 5)

morning_act_peak <- 7.13
morning_temp_valley <- 8.08
8.08-7.13
0.95 # ~1hr before
evening_act_peak <- 17.94
evening_temp_peak <- 16.40
17.94-16.40
1.54 # ~1.5hrs after
