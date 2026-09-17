# Modeling flight phenology of Gonepteryx rhamni in Germany

##Libraries:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(raster)
library(rnaturalearth)
library(MASS)
library(DHARMa)
library(forecast)
library(lubridate)
library(randomForest)
library(lme4)
library(gam)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(ggeffects)

#1. Data -------
train_y1 <- 2016:2019
train_y2 <- 2020:2024 

g_r_1_inat <- read.csv("Data_analysis/g_r_1_inat")
g_r_2_inat <- read.csv("Data_analysis/g_r_2_inat")
g_r_1_obs <- read.csv("Datasets_analysis/g_r_1_obs")
g_r_2_obs <- read.csv("Datasets_analysis/g_r_2_obs")

g_r_1_all <- rbind(g_r_1_inat, g_r_1_obs) 
g_r_2_all <- rbind(g_r_2_inat, g_r_2_obs)

g_r_all <- rbind(g_r_1_all, g_r_2_all)

#2. Plots, general statistics ------
g_r_summary <- g_r_all %>% 
  group_by(year) %>% 
  summarise(
    total_n = sum(n),
    .groups = "drop"
  )

par(mfrow=c(1,1))
plot(g_r_summary$total_n ~ g_r_summary$year, 
     pch = 19, 
     xlab = "year", 
     ylab = "Total number of observations (G. rhamni)")

#Maps with occurrence records: 
germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")
germany <- germany[6]

gr_map_1 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = g_r_1_all, col = "blue", size = 0.02, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 

gr_map_2 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = g_r_2_all, col = "blue", size = 0.02, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 

library(gridExtra)
grid.arrange(gr_map_1, gr_map_2, ncol = 2)


  #the upper north part is only filled slowly after 2021

any(duplicated(g_r_2_all["grid.id.x",]))

#for basic statistics: 
get_hotspots <- function(species_data) { #10 3 used in all, obs models
  species_data %>% 
  #  filter(gen_time > 6) %>% 
    filter(n > 2) %>% 
    dplyr::select(-first_obs,
                  -last_obs,
                  -grid_id.x, -grid_id.y, -X, -n, -lon, -lat
                  )
}

g_r_1_obs_hot <- get_hotspots(g_r_1_obs) 
g_r_1_inat_hot <- get_hotspots(g_r_1_inat)
g_r_2_obs_hot <- get_hotspots(g_r_2_obs)
g_r_2_inat_hot <- get_hotspots(g_r_2_inat)
g_r_1_all_hot <- get_hotspots(g_r_1_all)
g_r_2_all_hot <- get_hotspots(g_r_2_all)

par(mfrow=c(1,2))
hist(g_r_2_all_hot$peak, xlab = "Second flight peak of G. rhamni (day of the year)")
hist(g_r_1_all_hot$peak, xlab = "First flight peak of G. rhamni (day of the year)")


get_per_year_stats <- function(dataset) { 
  summary <- dataset %>% 
    group_by(year) %>% 
    summarise(
      peak = mean(peak), 
      gen_size = mean(gen_time), 
      .groups = "drop"
    )
  return(summary)
}

g_r_1_y_sum <- get_per_year_stats(g_r_1_all_hot)
g_r_2_y_sum <- get_per_year_stats(g_r_2_all_hot)

plot(g_r_1_y_sum$peak ~ g_r_1_y_sum$year, pch = 19, 
     xlab = "year", ylab = "G. rhamni's mean first flight peak (day of the year)")

plot(g_r_1_y_sum$gen_size ~ g_r_1_y_sum$year, pch = 19, 
     xlab = "year", ylab = "First peak length (days)")

plot(g_r_2_y_sum$gen_size ~ g_r_2_y_sum$year, pch = 19, 
     xlab = "year", ylab = "Second peak length (days)")

ggplot(g_r_1_y_sum, aes(x = peak, y = gen_size)) + 
         geom_point() +  
         geom_smooth(method = "lm") 
         
plot(g_r_1_y_sum$gen_size ~ g_r_1_y_sum$peak, pch = 19, 
     xlab = "First peak (days)", 
     ylab = "First peak length (days)")
abline(add = TRUE)

  
#3. Data processing -------
get_hotspots <- function(species_data) { #10 3 used in all, obs models
  species_data %>% 
    filter(n > 2) %>% 
    dplyr::select(-first_obs,
                  -last_obs,
                   -gen_time,
                  -grid_id.x, -grid_id.y, -X, -n, -lon, -lat, 
                  -year, 
                  -stateProvince
    )
}

#Hotspots for all years: 
g_r_1_obs_hot <- get_hotspots(g_r_1_obs) 
g_r_1_inat_hot <- get_hotspots(g_r_1_inat)
g_r_2_obs_hot <- get_hotspots(g_r_2_obs)
g_r_2_inat_hot <- get_hotspots(g_r_2_inat)

#together inat and obs.org: 
g_r_1_all_hot <- get_hotspots(g_r_1_all)
g_r_2_all_hot <- get_hotspots(g_r_2_all)

#year split
g_r_1_y1 <-  g_r_1_all %>% filter(year %in% train_y1) %>% get_hotspots() 

g_r_1_y2 <- g_r_1_all %>% filter(year %in% train_y2) %>% get_hotspots() 

g_r_2_y1 <- g_r_2_all %>% filter(year %in% train_y1) %>% get_hotspots()

g_r_2_y2 <- g_r_2_all %>% filter(year %in% train_y2) %>% get_hotspots()

mean(g_r_1_y1$peak) #111
mean(g_r_1_y2$peak) #103
mean(g_r_1_all$peak) #104
mean(g_r_2_all$peak) #202
mean(g_r_2_y1$peak) #195
mean(g_r_2_y2$peak) #201

hist(g_r_1_y1$peak)
hist(g_r_1_y2$peak)

#keeping different periods for analysis: 


#Dataset actually used in RF, Rpart, LMMs
# only obserrvation.org, filter years: 
g_r_1_obs_sub <- g_r_1_obs %>%  #117 
  filter(year %in% train_y1) %>%
  filter(n > 2) %>% 
  dplyr::select(-first_obs,
                -last_obs, -gen_time, 
                -grid_id.x, -grid_id.y, -X, -n, -lon, -lat, 
                -year, -stateProvince)

g_r_2_obs_sub <- g_r_2_obs %>%  
  filter(year %in% train_y1) %>%
  filter(n > 2) %>% 
  dplyr::select(-first_obs,
                -last_obs, -gen_time, 
                -grid_id.x, -grid_id.y, -X, -n, -lon, -lat, 
                -year, -stateProvince)

#4. Random Forest  ----------
#only obs y1: 
set.seed(12) #to get the same result every time
rf_gr1 <- randomForest(peak ~ ., #peak = response variable - all the other columns in the dataset should be the predictors
                     data = g_r_1_obs_sub,
                     ntree = 500, #can also be 1000, plot the model to check minimum error levels 
                     ntry = round(sqrt(ncol(g_r_1_obs_sub) - 1)), 
                     importance = TRUE, 
                     na.action = na.omit)

print(rf_gr1) # for g_r_1_obs_sub: Mean of squared residuals: 282.0153
#% Var explained: 47.67
#for 1_y1: 42.22, 299.74, for all:14 319  #summary to check MSE etc
varImpPlot(rf_gr1) #plot to check the most important ones 
imp_gr1 <- importance(rf_gr1)
imp_gr1 <- imp_gr1[order(imp_gr1[,1], decreasing = TRUE), , drop = FALSE]
head(imp_gr1, 15) #table of top 15 
#total mar prec, abs d tmin, total n prec, oc mean, jan tmin, may tmax

par(mfrow= c(1,1))
plot(g_r_1_obs_hot$peak ~ g_r_1_obs_hot$abs_d_tmin)
plot(g_r_1_obs_hot$peak ~ g_r_1_obs_hot$av_oc_tmean)
#negative not clear

#when run for g_r_1_y2, g_r_2_y2: #1.69 325.76, 1.35 243.22 2 y2 - excluded 


#obs and inat y1:
rf_gr1_all <- randomForest(peak ~ ., 
                           data = g_r_1_y1,
                           ntree = 500,
                           ntry = round(sqrt(ncol(g_r_1_y1) - 1)), 
                           importance = TRUE, 
                           na.action = na.omit)


print(rf_gr1_all) #43 294
varImpPlot(rf_gr1_all) 
#n mean min d min jl max n max

#main ones: n min mean, d min, oc mean, aug mean, abs n max, av aug mean, abs jl max
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_oc_tmean)
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_d_tmin)
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$abs_jl_tmax)
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_n_tmin)
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_n_tmean)
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$abs_n_tmax)
#plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_aug_tmean) unclear
#plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$total_n_prec) unclear

#poly like, 2 degree 
# 5. Rpart  --------
set.seed(12)
rpart_gr1 <- rpart(peak ~ ., data = g_r_1_obs_sub)
rpart.plot(rpart_gr1, cex = 0.9)
cp <- as.data.frame(printcp(rpart_gr1))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_gr1, 0.081231)
rpart.plot(rpart_pruned, cex = 1.2)
#d tmin, jan tmax, mar tmin


plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_oc_tmean)
#not so clear but negative
set.seed(12)
rpart_gr1 <- rpart(peak ~ ., data = g_r_1_all_hot)
rpart.plot(rpart_gr1, cex = 0.7)
cp <- as.data.frame(printcp(rpart_gr1))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_gr1, 0.02122703)
rpart.plot(rpart_pruned, cex = 0.8)
#abs d tmin, total n prec, av jun tmax

plot(g_r_1_all_hot$peak ~ g_r_1_all_hot$total_n_prec)
#not linear, more distributed left, weird shape
plot(g_r_1_all_hot$peak ~ g_r_1_all_hot$abs_d_tmin)
plot(g_r_1_all_hot$peak ~ g_r_1_all_hot$av_jun_tmax)
#linear spread out


#6. Polynomial models ----------
#Adding year back to the dataset to run mixed models with year intercept
g_r_1_obs_sub <- g_r_1_obs %>%  
  filter(year %in% train_y1) %>%
  filter(n > 2) %>% 
  dplyr::select(-first_obs, 
                -last_obs, -gen_time, 
                -grid_id.x, -grid_id.y, -X, -n, -lon, -lat)

poly_lm <- lm(peak ~ poly(av_oc_tmean, 2) + 
                poly(av_d_tmin, 2) + 
                poly(abs_jl_tmax, 2) + 
                poly(av_n_tmin, 2) + 
                poly(av_n_tmean, 2) + 
                poly(abs_n_tmax, 2), 
                data = g_r_1_obs_sub)
summary(poly_lm) #res error 17,11, rsq adj 0.46
resid <- simulateResiduals(poly_lm)
plot(resid) #good model

#how to include October although not significant? 

plot(av_oc_tmean ~ year, data = g_r_1_obs_sub)

#jl oc not sign. Try only winter? 
poly_lm2 <- lm(peak ~ poly(av_d_tmin, 3) + #better than 2
                av_jl_tmax +
                poly(av_n_tmin, 2) + 
                poly(av_n_tmean, 3) + #not good in fitting
                poly(abs_n_tmax, 2), 
              data = g_r_1_obs_sub)
summary(poly_lm2) #0.50 16.57
#poly 1 of av dmin, both of av n min, both of av nmean, poly of nmax

resid2 <- simulateResiduals(poly_lm2)
plot(resid2) #not ok slightly problematic

#simplifying: without oct and poly effect of dec.
poly_lm3 <- lm(peak ~ av_d_tmin + 
                 poly(av_n_tmin, 2) + 
                 poly(abs_n_tmax, 2), 
               data = g_r_1_obs_sub)
summary(poly_lm3) #0.42, 17.8 al sign 
simulateResiduals(poly_lm3) %>% plot()

#only clear relationships:
poly_lm4 <- lm(peak ~  #variations of n mean - similar results
                 poly(av_d_tmin, 3) + 
                 poly(av_n_tmean, 2) +
                 poly(av_n_tmin, 2) + 
                 poly(abs_n_tmax, 2), 
               data = g_r_1_obs_sub)
summary(poly_lm4) 
#50, 16.5 (all show some sign, only linear terms not all)
resid <- simulateResiduals(poly_lm4)
plot(resid) #pretty good 

AIC(poly_lm, poly_lm3, poly_lm4)
#fourth model is the best AIC = 930.9

#7. Mixed models ----------
lmm4.1 <- lmer(peak ~  #variations of n mean - similar results
                 poly(av_d_tmin, 3) + 
                 poly(av_n_tmean, 3) +
                 poly(av_n_tmin, 2) + 
                 poly(abs_n_tmax, 2) + 
                 (1 | year), 
               data = g_r_1_obs_sub)
summary(lmm4.1) 
resid <- simulateResiduals(lmm4.1) %>% plot() #ok
#sign: same, but n tmean 2 instead of linear

#Model diagnostics: 
AIC(poly_lm4, lmm4.1) 
ΔAIC = AIC(poly_lm4) - AIC(lmm4.1) # = 75.14 big difference

plot(av_d_tmin ~ year, data = g_r_1_obs_sub)
plot(abs_n_tmax ~ year, data = g_r_1_obs_sub)

lmm5 <- lmer(peak ~  
               poly(av_d_tmin, 3) + 
               poly(av_n_tmean, 2) +
               poly(av_n_tmin, 2) + 
               poly(abs_n_tmax, 2) + 
               (0 + abs_n_tmax | year),
             data = g_r_1_obs_sub)
summary(lmm5) 
resid <- simulateResiduals(lmm5) %>% plot() #ok

AIC(lmm4.1, lmm5) #similar  - 850 856, 1 better

lmm5.1 <- lmer(peak ~  
               poly(av_d_tmin, 3) + 
               poly(av_n_tmean, 2) +
               poly(av_n_tmin, 2) + 
               poly(abs_n_tmax, 2) + 
               (0 + av_d_tmin | year),
             data = g_r_1_obs_sub)
summary(lmm5.1) 

AIC(lmm4.1, lmm5, lmm5.1) #5.1. 853, 4.1. lower

resid <- simulateResiduals(lmm5.1) %>% plot() #ok

#n teman? 
lmm5.2 <- lmer(peak ~  
                 poly(av_d_tmin, 3) + 
                 poly(av_n_tmean, 3) +
                 poly(av_n_tmin, 2) + 
                 poly(abs_n_tmax, 2) + 
                 (0 + av_d_tmin | year),
               data = g_r_1_obs_sub)
summary(lmm5.2)
resid <- simulateResiduals(lmm5.2) %>% plot()  

#since n mean is problematic, try to do groups per n mean

poly_lm5 <- lm(peak ~  #variations of n mean - similar results
                 poly(av_d_tmin, 3) + 
                 poly(av_n_tmin, 2) +
                 poly(abs_n_tmax, 2) + 
                 av_n_tmin:av_n_tmean,
               data = g_r_1_obs_sub)
summary(poly_lm5)  #all sign, 16.54 0.50
resid <- simulateResiduals(poly_lm5)
plot(resid) #ok

#Model diagnostics:
AIC(lmm4.1, lmm5, 
  lmm5.1, lmm5.2) #5.2 the best

#lmm5.2 is the best model 
ranef(lmm5.2)
summary(lmm5.2
        )

g_r_1_obs_sub$fitted <- fitted(lmm5.2)

ggplot(data = g_r_1_obs_sub, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

ggplot(data = g_r_1_obs_sub, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ year) + 
  labs(
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 


years <- unique(g_r_1_obs_sub$year)
pred_n_min <- ggpredict(lmm5.2, 
                        terms = c("av_n_tmin[all]", "year"),
                        type = "random") %>% as.data.frame()
pred_n_min <- pred_n_min %>% mutate(fac_gr = group)
  
pred_n_max <- ggpredict(lmm5.2, terms = c("abs_n_tmax[all]", "year"),
                        type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)

pred_d_tmin <- ggpredict(lmm5.2, terms = c("av_d_tmin[all]", "year"),
                         type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
pred_n_tmean <- ggpredict(lmm5.2, terms = c("av_n_tmean[all]", "year"),
                          type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
g_r_1_plot <- g_r_1_obs_sub %>% mutate(fac_gr = year)


ggplot() + 
  geom_ribbon(data = pred_n_min, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_min, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = g_r_1_plot, 
             aes(x = av_n_tmin, y = peak)) + 
  facet_wrap(~ fac_gr, scales = "free") +
  labs(x = "Average minimum temperature of November (°C)",
       y = "Flight peak of G. rhamni (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = g_r_1_plot, 
             aes(x = abs_n_tmax, y = peak)) + 
  facet_wrap(~ fac_gr, scale = "free") +
  labs(x = "Average maximum temperature of November (°C)",
       y = "Flight peak of G. rhamni (day of the year)")


ggplot() + 
  geom_ribbon(data = pred_d_tmin, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_d_tmin, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = g_r_1_plot, 
             aes(x = av_d_tmin, y = peak)) + 
  facet_wrap(~ fac_gr, scale = "free") +
  labs(x = "Average minimum temperature of December (°C)",
       y = "Flight peak of G. rhamni (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_n_tmean, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_tmean, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = g_r_1_plot, 
             aes(x = av_n_tmean, y = peak)) + 
  facet_wrap(~ fac_gr, scales = "free") +
  labs(x = "Average mean temperature of November (°C)",
       y = "Flight peak of G. rhamni (days)")
#bad fit

ranef(lmm5.2)
coef(lmm5.2)

#n min is the best predictor 

#global predictions lmm5.1

pred_n_min <- ggpredict(lmm5.1, terms = "av_n_tmin[all]")
pred_n_max <- ggpredict(lmm5.1, terms = "abs_n_tmax[all]")
pred_d_tmin <- ggpredict(lmm5.1, terms = "av_d_tmin[all]")

ggplot() + 
  geom_ribbon(data = pred_n_min, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, color = "lightgrey", fill = NA)) +
  geom_point(data = g_r_1_obs_sub, 
             aes(x = av_n_tmin, y = peak)) + 
  geom_line(data = pred_n_min, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5)  +
  labs(title = "Effects of average November minimum temperature on 
       Gonepteryx rhamni peak, 2016-2019", x = "November min temperature (C)",
       y = "Peak (days)")

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, color = "lightgrey", fill = NA)) +
  geom_point(data = g_r_1_obs_sub, 
             aes(x = abs_n_tmax, y = peak)) + 
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) + 
  labs(title = "Effects of absolute November maximum temperature on 
       Gonepteryx rhamni peak, 2016-2019", x = "November maximum temperature (C)",
       y = "Peak (days)")
#this looks ok

ggplot() + 
  geom_ribbon(data = pred_d_tmin, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, color = "lightgrey", fill = NA)) +
  geom_point(data = g_r_1_obs_sub, 
             aes(x = av_d_tmin, y = peak)) + 
  geom_line(data = pred_d_tmin, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) + 
  labs(title = "Effects of average December minimum temperature on 
       Gonepteryx rhamni peak, 2016-2019", x = "December min temperature (C)",
       y = "Peak (days)")

###different intercepts, mixed per year a bit better 

ranef(lmm4)

#earlier peaks, better predictions: 2019 

#rpart 
rpart_gr_y1 <- rpart(peak ~ ., data = g_r_1_obs_sub)
rpart.plot(rpart_gr_y1, cex = 0.7) #nov and dec temp not here
cp <- as.data.frame(printcp(rpart_gr_y1))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_gr_y1, 0.0524106)
rpart.plot(rpart_pruned, cex = 0.8)
#oc tmean, abs ap tmax, abs aug tmax last
#oc has a big effect: 87 119 (32 days) !! - heatwave (abs 36 in Aug 30 days diff)

plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$abs_aug_tmax)
#nog clear 
plot(g_r_1_obs_sub$peak ~ g_r_1_obs_sub$av_oc_tmean)
#negative clear 
#if its warm they have an early peak


#only 2 of var expl. when running for generation 2

#6. Models excluded from the analysis ---------

#for second peak:
set.seed(12) #to get the same result every time
rf_gr2 <- randomForest(peak ~ ., #peak = response variable - all the other columns in the dataset should be the predictors
                       data = g_r_2_y1,
                       ntree = 500, #can also be 1000, plot the model to check minimum error levels 
                       ntry = round(sqrt(ncol(g_r_2_y1) - 1)), 
                       importance = TRUE, 
                       na.action = na.omit)

print(rf_gr2) #second data peak: 6.49 306.11
varImpPlot(rf_gr2)

# inat with smoother filtering (6-2 not 10-3)
#only filter for n > 2: 
rf_gr1_inat <- randomForest(peak ~ ., 
                            data = g_r_1_inat_hot,
                            ntree = 500,
                            ntry = round(sqrt(ncol(g_r_1_inat_hot) - 1)), 
                            importance = TRUE, 
                            na.action = na.omit)

print(rf_gr1_inat) #6
varImpPlot(rf_gr1_inat)

plot(g_r_1_inat_hot$peak ~ g_r_1_inat_hot$abs_jun_tmax)
#negative relationship
plot(g_r_1_inat_hot$peak ~ g_r_1_inat_hot$mean_nov_prec)
#negative relationship
plot(g_r_1_inat_hot$peak ~ g_r_1_inat_hot$mean_dec_prec)
#positive 

#jun tmax last, may tmax, nov dec prec last, total mean mar prec

#All second peak: 
rf_gr2all <- randomForest(peak ~ ., 
                          data = g_r_2_merged,
                          ntree = 500,
                          ntry = round(sqrt(ncol(g_r_2_merged) - 1)), 
                          importance = TRUE, 
                          na.action = na.omit)
print(rf_gr2all) #4 / 2 when adding naturgucker
varImpPlot(rf_gr2all)

rf_gr1all <- randomForest(peak ~ ., 
                          data = g_r_1_merged,
                          ntree = 500,
                          ntry = round(sqrt(ncol(g_r_1_merged) - 1)), 
                          importance = TRUE, 
                          na.action = na.omit)
print(rf_gr1all) #11 322
varImpPlot(rf_gr1all) #total n mar pres, abs d min may max, jan tmin

#smooth filtering only n > 2:
#abs d tmin, av jan tmin, mean mar prec, jan tmin, abs may tmax, 
#total mar prec, abs f min

