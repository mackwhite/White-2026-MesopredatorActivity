###project: Snook Activity Patterns
###author(s): MW, JR, WRJ, ROS
###goal(s): 
###date(s): May 2026
###note(s): 
# Housekeeping ------------------------------------------------------------

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, suncalc, nlme, DHARMa,
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot)

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
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) #|> 
# capture.output(file = "output/q1-trig-glmm-model-comparison.csv")

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

tac_test <- all |>
      mutate(datetime = make_datetime(year, month, day, time, tz = "America/New_York")) |>
      arrange(id, datetime) |>
      group_by(id) |>
      mutate(time_cont = as.numeric(difftime(datetime, min(datetime), units = "hours"))) |>
      ungroup()

tac <- mgcv::gam(y ~ s(time, bs="cc") + s(id, bs="re")+ s(station, bs="re"),
                        family = Gamma(link = 'log'),
                        data = tac_test)
plot(tac, shade = TRUE)
sim <- simulateResiduals(tac)
sim_agg <- recalculateResiduals(sim, group = tac_test$time) 
testTemporalAutocorrelation(sim_agg, time = sort(unique(tac_test$time)))

# fit gamm with no autocorrelation structure ------------------------------
snook_gamm <- mgcv::gam(y ~ s(time, bs="cc") + s(id, bs="re")+ s(station, bs="re"),
                        family = Gamma(link = 'log'),
                        data = all)
plot(snook_gamm, shade = TRUE)

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
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) #|> 
      # capture.output(file = "output/tables/q1-gamm-glmm-comparison.csv")

### model summary for gamm ---
summary(snook_gamm)
performance(snook_gamm)
# summary.gam(snook_gamm) |> 
#       capture.output(file = "output/tables/q1-gamm-summary.xlsx")

### join model fit datasets together for visualizations ---
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
            legend.position = c(0.13,0.85),
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 8, color = "black"))
a

ggsave('output/figs/snook-diel-pattern.png',
       dpi = 600, units= 'in', height = 4, width = 4)