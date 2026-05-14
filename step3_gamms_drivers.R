###project: mesoconsumer behaviors
###author(s): MW
###goal(s): quantify effects of temerature on diel activity patterns of snook
###date(s): February 2025
###note(s): 
# Housekeeping ------------------------------------------------------------

### load necessary libraries ----
# install.packages("librarian")
librarian::shelf(tidyverse, readr, lme4, ggpubr, performance, 
                 scales, ggstats, ggeffects, visreg, mgcv, MuMIn, glmmTMB, corrplot)

### load necessary data ----
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

### evaluate collinearity ---
numeric_vars <- all |> 
      select(weight,
             stage, 
             temp,
             time,
             # lunar_ill,
             activity = y)
cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
corrplot(cor_matrix, method = 'number', type = 'upper', tl.cex = 0.8)

# ggsave('output/figs/supp-correlation-matrices.png',
#        dpi = 600, units= 'in', height = 4, width = 4)

### clean environment ---
keep <- c("all")
rm(list = setdiff(ls(), keep))

### summarize to identify median + quantiles for each predictor ---
summary <- all |> 
      summarize(
            ### time summary ---
            min_time = min(time, na.rm = TRUE),
            q25_time = quantile(time, 0.25, na.rm = TRUE),
            median_time = quantile(time, 0.50, na.rm = TRUE),
            q75_time = quantile(time, 0.75, na.rm = TRUE),
            max_time = max(time, na.rm = TRUE),
            ### temperature summary ---
            min_temp = min(temp, na.rm = TRUE),
            q25_temp = quantile(temp, 0.25, na.rm = TRUE),
            median_temp = quantile(temp, 0.50, na.rm = TRUE),
            q75_temp = quantile(temp, 0.75, na.rm = TRUE),
            max_temp = max(temp, na.rm = TRUE),
            ### weight summary ---
            min_weight = min(weight, na.rm = TRUE),
            q25_weight = quantile(weight, 0.25, na.rm = TRUE),
            median_weight = quantile(weight, 0.50, na.rm = TRUE),
            q75_weight = quantile(weight, 0.75, na.rm = TRUE),
            max_weight = max(weight, na.rm = TRUE),
            ### stage summary ---
            min_stage = min(stage, na.rm = TRUE),
            q25_stage = quantile(stage, 0.25, na.rm = TRUE),
            median_stage = quantile(stage, 0.50, na.rm = TRUE),
            q75_stage = quantile(stage, 0.75, na.rm = TRUE),
            max_stage = max(stage, na.rm = TRUE)
      )

# influence of weight on diel activity cycles -----------------------------
t_size <- gam(y ~ te(time, weight) + s(stage, k = 3) + s(temp, k = 3) + s(station, bs = "re"),
            family = gaussian(link = 'log'),
            data = all,
            method = "REML")

s_size <- gam(y ~ s(time, bs="cc") + s(weight, k = 3) + s(stage, k = 3) + s(temp, k = 3) + s(station, bs = "re"),
              family = gaussian(link = 'log'),
              data = all,
              method = "REML")

l_size <- gam(y ~ s(time, bs="cc") + weight + s(stage, k = 3) + s(temp, k = 3) + s(station, bs = "re"),
              family = gaussian(link = 'log'),
              data = all,
              method = "REML")

performance::compare_performance(t_size, s_size, l_size)|> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc)

# performance::compare_performance(t_size, s_size, l_size)|> 
#       mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
#       capture.output(file = "output/tables/size-model-comparison.csv")

concurvity(t_size)
performance(t_size)
summary.gam(t_size)
summary.gam(s_size)
size_dev_ratio <- 19.6/15.1
bfm <- t_size

### clear up some memory ---
keep <- c('all', 'bfm', 't_size', 's_size')
rm(list = setdiff(ls(), keep))

### smoothed interaction effect of weight and time ----
new_data <- expand.grid(
      time = seq(from = 0, to = 24, by = 0.1),
      weight = seq(from = 326, to = 6393, by = 10),
      temp = 28,
      stage = 55
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(bfm, newdata=new_data, se.fit=TRUE)
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
             y = weight,
             z = predicted)

a <- pred_fit |> 
      mutate(y = y/1000) |>
      ggplot(aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", direction = 1,
                           limits = c(0.071, 0.606)) + 
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(breaks = c(0,1,2,3,4,5,6)) +
      labs(x = "Time of Day (h)", y = "Weight (kg)", 
           fill = expression(bold("Predicted Activity (m/s"^2*")"))) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),            
            legend.position = "bottom",
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 10, color = "black"))
a

### smoothed main effect of size ----
new_data <- expand.grid(
      time = 12,
      weight = seq(from = 326, to = 6393, by = 1),
      temp = 28,
      stage = 55
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(s_size, newdata=new_data, se.fit=TRUE)
new_data$fit <- pred$fit
new_data$se.fit <- pred$se.fit

pred_fit <- new_data |> 
      mutate(predicted_log = pred$fit,
             predicted = exp(predicted_log),
             lower_log = predicted_log - pred$se.fit,
             upper_log = predicted_log + pred$se.fit) |>
      mutate(lower = exp(lower_log),
             upper = exp(upper_log)) |>
      select(time, weight, predicted, lower, upper) |> 
      rename(x = weight,
             y = predicted)

b <- pred_fit |> 
      mutate(x = x/1000) |> 
      ggplot(aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line(linewidth = 2, linetype = "solid") +
      theme_bw() +
      labs(x = "Weight (kg)", y = expression(bold("Predicted Activity (m/s"^2*")"))) +
      scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
      scale_y_continuous(breaks = c(0.25,0.30,0.35,0.40,0.45), limits = c(0.25,0.45)) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),            
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold'))
b

# summary.gam(bfm) |> capture.output(file = "output/tables/q2-model-summary-size.csv")
# influence of temperature on diel activity cycles -----------------------------

t_temp <- gam(y ~ te(time, temp) + s(stage, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
            family = gaussian(link = 'log'),
            data = all,
            method = "REML")

s_temp <- gam(y ~ s(time, bs="cc") + s(temp, k = 3) + s(stage, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
              family = gaussian(link = 'log'),
              data = all,
              method = "REML")

l_temp <- gam(y ~ s(time, bs="cc") + temp + s(stage, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
              family = gaussian(link = 'log'),
              data = all,
              method = "REML")

performance::compare_performance(t_temp, s_temp, l_temp) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc)

# performance::compare_performance(t_temp, s_temp, l_temp) |> 
#       mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
#       capture.output(file = "output/tables/temp-model-comparison.csv")

concurvity(t_temp)
performance(t_temp)
summary.gam(t_temp)
summary.gam(s_temp)
temp_dev_ratio <- 18.7/15.1
bfm <- t_temp

### clear up some memory ---
keep <- c('all', 'bfm', "a", "b", 't_size', 's_size', 't_temp', 's_temp')
rm(list = setdiff(ls(), keep))

### smoothed interaction effect of temp and time ----
new_data <- expand.grid(
      time = seq(from = 0, to = 24, by = 0.1),
      weight = 1396,
      temp = seq(from = 17, to = 34, by = 0.1),
      stage = 55
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(bfm, newdata=new_data, se.fit=TRUE)
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
             y = temp,
             z = predicted)

c <- pred_fit |> 
      ggplot(aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_y_continuous(breaks = c(18,20,22,24,26,28,30,32,34)) +
      scale_fill_viridis_c(option = "plasma", direction = 1,
                           limits = c(0.071, 0.606)) +  
      labs(x = "Time of Day (h)", y = "Temperature (°C)", 
           fill = expression(bold("Predicted Activity (m/s"^2*")"))) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),            
            legend.position = "bottom",
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 10, color = "black"))
c

### smoothed main effect of temperature ----
new_data <- expand.grid(
      time = 12,
      weight = 1396,
      temp = seq(from = 17, to = 34, by = 0.1),
      stage = 55
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(s_temp, newdata=new_data, se.fit=TRUE)
new_data$fit <- pred$fit
new_data$se.fit <- pred$se.fit

pred_fit <- new_data |> 
      mutate(predicted_log = pred$fit,
             predicted = exp(predicted_log),
             lower_log = predicted_log - pred$se.fit,  
             upper_log = predicted_log + pred$se.fit) |> 
      mutate(lower = exp(lower_log),
             upper = exp(upper_log)) |> 
      select(time, temp, predicted, lower, upper) |> 
      rename(x = temp,
             y = predicted)

d <- pred_fit |> 
      ggplot(aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line(linewidth = 2, linetype = "solid") +
      theme_bw() +
      labs(x = "Temperature (°C)", y = expression(bold("Predicted Activity (m/s"^2*")"))) +
      scale_x_continuous(breaks = c(18,20,22,24,26,28,30,32,34)) +
      scale_y_continuous(breaks = c(0.25,0.30,0.35,0.40,0.45), limits = c(0.25,0.45)) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),            
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold'))
d

# summary.gam(bfm) |> capture.output(file = "output/tables/q2-model-summary-temp.csv")

# influence of stage on diel activity cycles -----------------------------
t_stage <- gam(y ~ te(time, stage) + s(temp, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
             family = gaussian(link = 'log'),
             data = all,
             method = "REML")

s_stage <- gam(y ~ s(time, bs="cc") + s(stage, k = 3) + s(temp, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
               family = gaussian(link = 'log'),
               data = all,
               method = "REML")

l_stage <- gam(y ~ s(time, bs="cc") + stage + s(temp, k = 3) + s(weight, k = 3) + s(station, bs = "re"),
               family = gaussian(link = 'log'),
               data = all,
               method = "REML")

performance::compare_performance(t_stage, s_stage, l_stage) |> 
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc)

# performance::compare_performance(t_stage, s_stage, l_stage) |> 
#       mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |> 
#       capture.output(file = "output/tables/stage-model-comparison.csv")

concurvity(t_stage)
performance(t_stage)
summary.gam(t_stage)
summary(s_stage)
stage_dev_ratio <- 16.9/15.1
bfm <- t_stage

### clear up some memory ---
keep <- c('all', 'bfm', "a", "b", 'c', 'd', 't_size', 's_size', 't_temp', 's_temp', 't_stage', 's_stage')
rm(list = setdiff(ls(), keep))

### smoothed interaction effect of stage and time ----
new_data <- expand.grid(
      time = seq(from = 0, to = 24, by = 0.1),
      weight = 1396,
      temp = 28,
      stage = seq(from = 38, to = 79, by = 0.1)
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(bfm, newdata=new_data, se.fit=TRUE)
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
             y = stage,
             z = predicted)

e <- pred_fit |> 
      ggplot(aes(x = x, y = y, fill = z)) +
      geom_tile() +
      scale_y_continuous(breaks = c(40,50,60,70,80)) +
      scale_x_continuous(breaks = c(0,4,8,12,16,20,24)) +
      scale_fill_viridis_c(option = "plasma", direction = 1,
                           limits = c(0.071, 0.606)) + 
      labs(x = "Time of Day (h)", y = "Marsh Stage (cm)", 
           fill = expression(bold("Predicted Activity (m/s"^2*")"))) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),            
            legend.position = "bottom",
            legend.text = element_text(face = 'bold', size = 8, color = "black"),
            legend.title = element_text(face = 'bold', size = 10, color = "black"))
e

### smoothed main effect of stage ----
new_data <- expand.grid(
      time = 12,
      weight = 1396,
      temp = 28,
      stage = seq(from = 38, to = 79, by = 0.1)
) |>
      mutate(id = all$id[4],
             station = all$station[2])

pred <- predict(s_stage, newdata=new_data, se.fit=TRUE)
new_data$fit <- pred$fit
new_data$se.fit <- pred$se.fit

pred_fit <- new_data |> 
      mutate(predicted_log = pred$fit,
             predicted = exp(predicted_log),
             lower_log = predicted_log - pred$se.fit,  
             upper_log = predicted_log + pred$se.fit) |> 
      mutate(lower = exp(lower_log),
             upper = exp(upper_log)) |> 
      select(stage, time, predicted, lower, upper) |> 
      rename(x = stage,
             y = predicted)

f <- pred_fit |> 
      ggplot(aes(x = x, y = y)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
      geom_line(linewidth = 2, linetype = "dashed") +
      theme_bw() +
      labs(x = "Marsh Stage (cm)", y = expression(bold("Predicted Activity (m/s"^2*")"))) +
      scale_x_continuous(breaks = c(40,50,60,70,80)) +
      scale_y_continuous(breaks = c(0.25,0.30,0.35,0.40,0.45), limits = c(0.25,0.45)) +
      theme(axis.text = element_text(size = 12, face = "bold", colour = "black"),
            axis.title = element_text(size = 14, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "none",
            legend.text = element_text(face = 'bold'),
            legend.title = element_text(face = 'bold'))
f

# summary.gam(bfm) |> capture.output(file = "output/tables/q2-model-summary-stage.csv")

### check out plots ----
keep <- c('all', "a", "b", 'c', 'd', 'e', 'f', 't_size', 's_size', 't_temp', 's_temp', 't_stage', 's_stage')
rm(list = setdiff(ls(), keep))

compare_performance(t_size, t_temp, t_stage)

### arrange panels for figures 3:5 ----

### figure three ---
fig3_arranged <- ggarrange(a,b, ncol = 2, align = "h", labels = c("a", "b"), common.legend = TRUE,
                           legend = "bottom")
fig3_arranged
# ggsave('output/figs/figure-three.png',
#        dpi = 600, units= 'in', bg = 'white', height = 4, width = 8)

### figure four ---
fig5_arranged <- ggarrange(c,d, ncol = 2, align = "h", labels = c("a", "b"), common.legend = TRUE,
                           legend = "bottom")
fig5_arranged
# ggsave('output/figs/figure-five.png',
#        dpi = 600, units= 'in', bg = 'white', height = 4, width = 8)

### figure five ---
fig4_arranged <- ggarrange(e,f, ncol = 2, align = "h", labels = c("a", "b"), common.legend = TRUE,
                           legend = "bottom")
fig4_arranged
# ggsave('output/figs/figure-four.png',
#        dpi = 600, units= 'in', bg = 'white', height = 4, width = 8)

# full model --------------------------------------------------------------

full_gamm <- mgcv::gam(y ~ s(time, bs="cc") + s(temp, k = 3) + s(weight, k = 3) + s(stage, k = 3) + s(station, bs="re"),
                        family = gaussian(link = 'log'),
                        data = all)

concurvity(full_gamm)
performance(full_gamm)
# summary.gam(full_gamm) |> capture.output(file = "output/tables/q2-model-summary-full.csv")

performance::compare_performance(t_size, t_temp, t_stage, full_gamm) |>
      mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc)

# performance::compare_performance(t_size, t_temp, t_stage, full_gamm) |>
#       mutate(dAICc = AICc - min(AICc)) |> arrange(dAICc) |>
#       capture.output(file = "output/tables/q2-gamm-model-comparisons.xlsx")

### summary stats for manuscript ----
fig3_arranged
summary.gam(t_size)
size_dev_ratio <- 1.30

fig5_arranged
summary.gam(t_temp)
temp_dev_ratio <- 1.24

fig4_arranged
summary.gam(t_stage)
stage_dev_ratio <- 1.12

glimpse(all)

act <- all |> 
      summarize(
            min = min(y),
            max = max(y), 
            mean = mean(y), 
            sd = sd(y)
      )
glimpse(act)