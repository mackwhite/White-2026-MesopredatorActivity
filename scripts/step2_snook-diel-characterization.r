###project: Snook Activity Patterns
###author(s): MW, JR, WRJ, ROS
###goal(s): 
###date(s): May 2026
###note(s): 
# Housekeeping ------------------------------------------------------------

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, suncalc, nlme, DHARMa,
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot, suncalc)

### load necessary data ----
all <- read_csv('local-data/snook-acc-model-data.csv') |>
      mutate(time_og = time,
             time    = hour(hms(time_og)),
             y       = mean_acceleration) |>
      mutate(id      = as.factor(id), station = as.factor(station)) |>
      dplyr::select(y, time, id, station, year, month, day) |> 
      arrange(id, year, month, day)
glimpse(all)

# fit glmms with trigonometric terms --------------------------------------

### models fit with Gaussian error distribution
trig_null     <- glmmTMB(y ~ 1 + (1|id) + (1|station), 
                         family = gaussian(link = "log"),
                         data = all)

trig_unimodal <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + (1|id) + (1|station),
                         family = gaussian(link = "log"),
                         data = all)

trig_bimodal  <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + cos(2*pi*time/12) + sin(2*pi*time/12) + (1|id) + (1|station),
                         family = gaussian(link = "log"),
                         data = all)

### models fit with Gammma error distributions
trig_null_gamma     <- glmmTMB(y ~ 1 + (1|id) + (1|station), 
                               family = Gamma(link = "log"),
                               data = all)

trig_unimodal_gamma <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + (1|id) + (1|station),
                               family = Gamma(link = "log"),
                               data = all)

trig_bimodal_gamma  <- glmmTMB(y ~ cos(2*pi*time/24) + sin(2*pi*time/24) + cos(2*pi*time/12) + sin(2*pi*time/12) + (1|id) + (1|station),
                               family = Gamma(link = "log"),
                               data = all)

### model comparison ----
performance::compare_performance(trig_null, trig_unimodal, trig_bimodal,
                                 trig_null_gamma, trig_unimodal_gamma, trig_bimodal_gamma) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
      capture.output(file = "output/q1-trig-glmm-model-comparison.csv")

performance::check_model(trig_bimodal_gamma)
performance::check_collinearity(trig_bimodal_gamma)
performance::check_convergence(trig_bimodal_gamma)

### model summary ---
# summary(trig_bimodal) |> 
#       capture.output(file = "output/q1-bestfit-glmm-summary.xlsx")

### generate file for model predictions/visualizations ---
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
pred <- predict(trig_bimodal_gamma, newdata = new_data, se.fit = TRUE)
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

bfm_glmm <- trig_bimodal_gamma

# fit gamm analog for comparison -----------------------------------------------
# double-check temporal autocorrelation ----------------------------------------

tac <- mgcv::gam(y ~ s(time, bs="cc", k = 10) + s(id, bs="re", k = 10)+ s(station, bs="re"),
                        family = Gamma(link = 'log'), data = all, method = "REML")

sim <- simulateResiduals(tac)
sim_agg <- recalculateResiduals(sim, group = all$time) 
testTemporalAutocorrelation(sim_agg, time = sort(unique(all$time)))

# fit gamm with no autocorrelation structure ------------------------------
snook_gamm <- mgcv::gam(y ~ s(time, bs="cc", k = 10) + s(id, bs="re", k = 10)+ s(station, bs="re"),
                        family = Gamma(link = 'log'),
                        data = all, method = "REML")
summary(snook_gamm)
set.seed(1)
k.check(snook_gamm)
set.seed(1)
gam.check(snook_gamm)

### visualize model predictions ---
snook_gamm_vis <- visreg(snook_gamm, "time", type = "conditional", scale = "response")

### generate file for model predictions/visualizations ---
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

### compare performance of gamm and glmm model ---
compare_performance(bfm_glmm, snook_gamm) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
      capture.output(file = "output/tables/q1-gamm-glmm-comparison.csv")

### model summary for gamm ---
summary(snook_gamm)
performance(snook_gamm)
summary.gam(snook_gamm) |>
      capture.output(file = "output/tables/q1-gamm-summary.xlsx")

### join model fit datasets together for visualizations ---
snook_fit <- rbind(glmm_fit, gamm_fit) |> 
      mutate(model = as.factor(model))

### clean up environment
keep <- c("all", "snook_fit", 'gamm_fit')
rm(list = setdiff(ls(), keep))


# Pull crepuscular information for Shark River ---------------------------------

### designate coordinates for reference
lat <- 25.4413
lon <- -80.9060

### provide range of dates for taking range/mean
dates <- seq(as.Date("2024-01-13"), as.Date("2024-12-14"), by = "day")

### pull information
sun_times <- getSunlightTimes(
      date = dates,
      lat = lat, lon = lon,
      keep = c("dawn", "sunrise", "sunset", "dusk"),
      tz = "America/New_York"
) |> 
      
      ### extract the decimal hour to align with data in figure
      mutate(
            dawn_hour    = hour(dawn)    + minute(dawn)    / 60,
            sunrise_hour = hour(sunrise) + minute(sunrise) / 60,
            sunset_hour  = hour(sunset)  + minute(sunset)  / 60,
            dusk_hour    = hour(dusk)    + minute(dusk)    / 60
      ) |> 
      
      ### take the summary stats for all the metrics :D
      summarize(
            ### anywhere a column ends in '_hour'
            across(ends_with('_hour'), 
            ### take the mean, min, and max and label it with the function
            list(mean = mean, min = min, max = max),
            ### but in addition to the function (i.e., 'fn'), start column with og column (i.e., col) name
            .names = '{.col}_{.fn}'
      ))

print(sun_times)

# Visualize the model predictions together w crepusc info ----------------------
a <- snook_fit |> 
      rename(Model = model) |> 
      mutate(Model = case_when(
            Model == "gamm" ~ "GAMM",
            Model == "glmm" ~ "GLMM"
      )) |> 
      ggplot(aes(x = x)) +
      ### morning crepuscular period
      annotate("rect",
               xmin = sun_times$dawn_hour_min,
               xmax = sun_times$sunrise_hour_max,
               ymin = -Inf, ymax = Inf,
               fill = "grey", alpha = 0.15) +
      ### evening crepuscular period
      annotate("rect",
               xmin = sun_times$sunset_hour_min,
               xmax = sun_times$dusk_hour_max,
               ymin = -Inf, ymax = Inf,
               fill = "grey", alpha = 0.15) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Model), alpha = 0.4) +
      geom_line(aes(y = y, linetype = Model), linewidth = 2) +
      theme_bw() +
      labs(x = "Time of Day (h)", y = expression(bold("Predicted Activity (m/s"^2*")")),
           fill = "Model",
           linetype = "Model") +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(breaks = c(0.35,0.45,0.55,0.65), limits = c(0.32,0.65)) +
      scale_fill_manual(values = c("GLMM" = "#1f78b4", "GAMM" = "#33a02c")) +
      scale_linetype_manual(values = c("GLMM" = "solid", "GAMM" = "dashed")) +
      theme(axis.text = element_text(size = 10, face = "bold", colour = "black"),
            axis.title = element_text(size = 12, face = "bold", colour = "black"),
            plot.title = element_text(size = 10, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = c(0.13,0.85),
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 8, color = "black"))
a

ggsave('output/figs/snook-diel-pattern.png',
       dpi = 600, units= 'in', height = 4, width = 4)
