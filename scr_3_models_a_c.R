# Modeling flight phenology of Anthocharis cardamines in Germany

#Libraries:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(psych)
library(raster)
library(MASS)
library(DHARMa)
library(ggeffects)
library(rnaturalearth)
library(forecast)
library(lubridate)
library(randomForest)
library(lme4)
library(rpart)
library(rpart.plot)

# 1. Data -------

anth_c_inat <- read.csv("Data_analysis/anth_c_inat") 
anth_c_ob <- read.csv("Datasets_analysis/anth_c_obs")
anth_c_all <- rbind(anth_c_inat, anth_c_ob) %>% distinct()

anth_c_all$first_obs <- as.integer(yday(as.Date(anth_c_all$first_obs)))
anth_c_all$last_obs <- as.integer(yday(as.Date(anth_c_all$last_obs)))

anth_c_summary <- anth_c_all %>% 
  group_by(year) %>% 
  summarise(
    total_n = sum(n),
    .groups = "drop"
  )

#2. Plots and basic statistics -------
par(mfrow=c(1,1))
plot(anth_c_summary$total_n ~ anth_c_summary$year, 
     pch = 19, 
     xlab = "year", 
     ylab = "Total number of observations (A. cardamines)")


anth_c_hotspots <- subset(anth_c_all, n > 2) %>% filter(gen_time > 6)
#previously 1, 6, then would filter for 3 and 6

# **: filtering not so necessary since we filter for n later anyway

anth_c_peak <- anth_c_hotspots %>% 
  group_by(year, grid_id, lon, lat) %>% 
  summarise(
    peak = as.integer(mean(c(last_obs, first_obs))),
    .groups = "drop"
  )

anth_c_dset <- inner_join(anth_c_hotspots, anth_c_peak, by = c("year", "grid_id", "lon", "lat"))

get_per_year_stats <- function(dataset) { 
  summary <- dataset %>% 
    group_by(year) %>% 
    summarise(
      peak = mean(peak), 
      gen_length = mean(gen_time), 
      .groups = "drop"
    )
  return(summary)
}

anth_c_summary <- get_per_year_stats(anth_c_dset)

plot(peak ~ year, data = anth_c_summary, pch = 19)

mean(anth_c_dset$peak)

#only obs? 
anth_c_obs <- subset(anth_c_ob, n > 3) %>% filter(gen_time > 6)

anth_c_obs$first_obs <- as.integer(yday(as.Date(anth_c_obs$first_obs)))
anth_c_obs$last_obs <- as.integer(yday(as.Date(anth_c_obs$last_obs)))

# **: filtering not so necessary since we filter for n later anyway

anth_c_peak_obs <- anth_c_obs %>% 
  group_by(year, grid_id, lon, lat) %>% 
  summarise(
    peak = as.integer(mean(c(last_obs, first_obs))),
    .groups = "drop"
  )

anth_c_obs_data <- inner_join(anth_c_obs, anth_c_peak_obs, by = c("year", "grid_id", "lon", "lat"))

par(mfrow=c(1,1))
hist(anth_c_dset$peak, xlab = "Flight peak od A. cardamines (day of the year)")

anth_c_dset <- anth_c_dset %>% filter(peak >= 90 & peak <= 170) #up to 20 of June
hist(anth_c_dset$peak, xlab = "Flight peak od A. cardamines (day of the year)")
anth_c_obs_data <- anth_c_obs_data %>% filter(peak >= 90 & peak <= 170) #up to 20 of June

#Using two periods for random forest: 
train_y1 <- 2016:2019
train_y2 <- 2020:2024 

#Maps with occurrence records:  
germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")
germany <- germany[6]

ggplot() + 
  geom_sf(data = germany, fill = NA, col = "blacK") +
  geom_point(data = anth_c_all, col = "blue", size = 0.02, fill = NA, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 
#its everywhere but the north part is slowly filled after 2020

#3. Analysis
# Random forest, rpart
# 3.1. RF, rpart 2016-2019 -------- 
#2016 - 2019 
rf_data_ac2 <- anth_c_dset %>%  #99 obs
  filter(year %in% train_y1) %>%
  dplyr::select(-X, -n, -grid_id, -first_obs,
                -last_obs, -gen_time, -year, -lon, -lat, -stateProvince)
                
set.seed(12)
rfac2 <- randomForest(peak ~ ., 
                      data = rf_data_ac2,
                      ntree = 500,
                      ntry = round(sqrt(ncol(rf_data_ac2) - 1)), 
                      importance = TRUE, 
                      na.action = na.omit)
print(rfac2) #59.0 var. expl, 91.7 MSE
varImpPlot(rfac2)
imp_ac2 <- importance(rfac2)
imp_ac2 <- imp_ac2[order(imp_ac2[,1], decreasing = TRUE), , drop = FALSE]
head(imp_ac2, 15)
#d tmin, s tmax, d tmean, d tmax, jun tmean, 
#aug tmax, d tmax, oc tmean, abs aug tmax, av n tmean, abs jan tmax,
#av n tmin, abs jl tmax
#all dec more important!

plot(rf_data_ac2$peak ~ rf_data_ac2$av_d_tmin)
plot(rf_data_ac2$peak ~ rf_data_ac2$av_d_tmean)
plot(rf_data_ac2$peak ~ rf_data_ac2$av_d_tmax)
plot(rf_data_ac2$peak ~ rf_data_ac2$av_s_tmax)
plot(rf_data_ac2$peak ~ rf_data_ac2$av_aug_tmax) #poly  3
plot(rf_data_ac2$peak ~ rf_data_ac2$av_jun_tmean) #linear
plot(rf_data_ac2$peak ~ rf_data_ac2$av_n_tmean) #poly 2
plot(rf_data_ac2$peak ~ rf_data_ac2$av_jl_tmax) #poly 2? 
#very steep negative, linear

#rpart
set.seed(12)
rpart_model_y1 <- rpart(peak ~ ., data = rf_data_ac2, method = "anova")
rpart.plot(rpart_model_y1, cex = 0.8, main = "rpart model splits: Anthocharis cardamines peak")
cp_y1 <- as.data.frame(printcp(rpart_model_y1))
(cp_y1[which.min(cp_y1[,"xerror"]), "CP"])

rpart_ac_prun_y1 <- prune(rpart_model_y1, 0.01000)
par(mfrow=c(1,1))
rpart.plot(rpart_ac_prun_y1, cex = 0.8)
#s tmax, d tmax, jl prec, d tmin, aug tmax

#av s tmax, av d tmax, av d tmin
#dmax: total jl prec - total jl prec
#dmin: av aug tmax - av may tmax
#extra info: may tmax, opossite effect!!, total jl prec 
plot(rf_data_ac2$peak ~ rf_data_ac2$total_jl_prec)
#unclear to fit model
plot(rf_data_ac2$peak ~ rf_data_ac2$av_aug_tmax)
#better, can fit linear

# 3.2. LMs 2016-2019 -------

rf_data_ac2 <- anth_c_dset %>%  #99 observations
  filter(year %in% train_y1) %>%
  dplyr::select(-X, -n, -grid_id, -first_obs,
                -last_obs, -gen_time, -lon, -lat)


lm_ac2 <- lm(peak ~ av_d_tmin 
             + av_d_tmean 
             + av_d_tmax
             + av_s_tmax 
             + av_aug_tmax, 
            data = rf_data_ac2)
summary(lm_ac2) 
simulateResiduals(lm_ac2) #one problem

lm_ac2.1 <- lm(peak ~ av_d_tmin 
             + av_d_tmean 
             + av_d_tmax
             + av_s_tmax 
             + av_aug_tmax 
             + av_jun_tmean , 
             data = rf_data_ac2)
summary(lm_ac2.1)  

resid <- simulateResiduals(lm_ac2)
plot(resid) #almost okay, 1 small problem

#interaction: 
lm_ac2.2 <- lm(peak ~ av_d_tmin * 
               av_s_tmax + 
               + av_d_tmean 
               + av_d_tmax
               + av_s_tmax 
               + av_aug_tmax 
               + av_jun_tmean
               + av_n_tmax,
               data = rf_data_ac2)
summary(lm_ac2.2)  
simulateResiduals(lm_ac2.2) %>% plot() #ok!
#too complex: reject jun tmean

plot(peak ~ av_n_tmean, data = rf_data_ac2)

#chosen model: 
lm_ac_2.3 <- lm(peak ~ 
                  av_s_tmax + 
                  av_d_tmin 
                + av_d_tmean 
                + av_d_tmax
                + av_aug_tmax
                + av_n_tmax, 
                data = rf_data_ac2)

summary(lm_ac_2.3)
#0.61 9.3
resid <- simulateResiduals(lm_ac_2.3)
plot(resid) #ok, but interaction not significant

lm_ac_2.4 <- lm(peak ~  av_s_tmax * 
                  av_d_tmax + 
                + av_d_tmean 
                + av_d_tmin
                + av_aug_tmax
                + av_n_tmax, 
                data = rf_data_ac2)

summary(lm_ac_2.4)
resid <- simulateResiduals(lm_ac_2.4)
plot(resid) #ok, interaction not sign

#make nov poly ? 
#chosen model: 
lm_ac_2.5 <- lm(peak ~ 
                av_s_tmax +
                av_d_tmax + 
                + av_d_tmin
                + av_jun_tmean
                + av_aug_tmax
                + poly(av_n_tmax, 2),
                data = rf_data_ac2)

summary(lm_ac_2.5)

lm_ac_2.6 <- lm(peak ~ 
                  av_s_tmax +
                  av_d_tmax + 
                  + av_d_tmin
                + av_jun_tmean
                + av_aug_tmax
                + poly(av_n_tmean, 2),
                data = rf_data_ac2)

summary(lm_ac_2.6)

resid <- simulateResiduals(lm_ac_2.6)
plot(resid) #ok

AIC(lm_ac2.2, lm_ac_2.3, lm_ac_2.4, lm_ac_2.5, lm_ac_2.6)
# 5 or 6 are the best

lmm_ac2_1 <- lmer(peak ~ 
                  av_d_tmin 
                + av_s_tmax 
                + av_d_tmax
                + av_jun_tmean
                + av_aug_tmax
                + poly(av_n_tmean, 2) + (1 | year),
                data = rf_data_ac2)

summary(lmm_ac2_1)
ranef(lmm_ac2_1)
simulateResiduals(lmm_ac2_1) %>% plot()
#ok

#State province/year bias - first years expected to play a role: 
lmm_ac2_2 <- lmer(peak ~ #like lm ac 5
                    av_d_tmin 
                  + av_s_tmax
                  + av_d_tmax
                  + av_aug_tmax
                  + poly(av_n_tmax, 2) + (1 | year),
                  data = rf_data_ac2)

summary(lmm_ac2_2)
resid <- simulateResiduals(lmm_ac2_2) %>% plot
#not ok
ranef(lmm_ac2_2) #-2 for 2019

summary(lm_ac_2.6)

#Model diagnostics: 
AIC(lm_ac_2.5, lm_ac_2.6, lmm_ac2_1, lmm_ac2_2) 
BIC(lm_ac_2.6, lmm_ac2_1, lmm_ac2_2) 
#mixed with state province and all predictors is better
#lmm ac2.1 the best

#lmm_ac2_1 predictions: 
rf_data_ac2$fitted <- fitted(lmm_ac2_1)

ggplot(data = rf_data_ac2, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ year, scale = "free") +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines, 2016-2019",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

#pretty good
rf_data_ac2$stateProvince <- as.factor(rf_data_ac2$stateProvince)
ranef(lmm_ac2_1)
#2017 - 1.3

plot_data <- rf_data_ac2 %>% mutate(fac_gr = year)

pred_n_tmean <- ggpredict(lmm_ac2_1, terms = c("av_n_tmean[all]", "year"), 
                          type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)

pred_aug_tmax <- ggpredict(lmm_ac2_1, terms = c("av_aug_tmax[all]", "year"), 
                           type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)

ggplot() + 
  geom_ribbon(data = pred_n_tmean, 
              aes(x = x, ymin = conf.low, group = fac_gr,  
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_tmean, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = plot_data, 
             aes(x = av_n_tmean, y = peak)) + 
  facet_wrap(~ fac_gr, scales = "free") +
  labs(x = "Average mean temperature of November (°C)",
       y = "Flight peak of A. cardamines (day of the year)") 


ggplot() + 
  geom_ribbon(data = pred_aug_tmax, 
              aes(x = x, ymin = conf.low, group = fac_gr, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_aug_tmax, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) + 
  geom_point(data = plot_data, 
             aes(x = av_aug_tmax, y = peak)) + 
  facet_wrap(~ fac_gr) +
  labs(x = "Average maximum temperature of August (°C)",
       y = "Flight peak of A. cardamines (day of the year)") #not that good

# with polynomial model: 

rf_data_ac2$fitted <- fitted(lm_ac_2.5)

ggplot(data = rf_data_ac2, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines, 2016-2019",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  facet_wrap(~ year) +
  theme_minimal() 


#predictions, most significant ones:
pred_n_max <- ggpredict(lm_ac_2.5, terms = "av_n_tmax")
pred_aug_tmax <- ggpredict(lm_ac_2.5, terms = "av_aug_tmax")
pred_d_tmin <- ggpredict(lm_ac_2.5, terms = "av_d_tmin")

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_point(data = rf_data_ac2, 
             aes(x = av_n_tmax, y = peak)) + 
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  facet_wrap(~ year) +
  labs(title = "Effects of average November max temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "Average maximum temperature of November (C)",
       y = "Peak (days)")
#great but spread out

ggplot() + 
  geom_ribbon(data = pred_aug_tmax, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_point(data = rf_data_ac2, 
             aes(x = av_aug_tmax, y = peak)) + 
  geom_line(data = pred_aug_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  facet_wrap(~ year) +
  labs(title = "Effects of average August maximum temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "Average maximum temperature of August (C)",
       y = "Peak (days)")
#missing a small group, good for the rest 

ggplot() + 
  geom_ribbon(data = pred_d_tmin, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_point(data = rf_data_ac2, 
             aes(x = av_d_tmin, y = peak)) + 
  geom_line(data = pred_d_tmin, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  facet_wrap(~ year) +
  labs(title = "Effects of average December min temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "December minimum temperature (C)",
       y = "Peak (days)")
#missing small group



# 3.3. Random foret, rpart 2020-2024 ---------- 
rf_data_ac3 <- anth_c_dset %>% #414 obs.
  filter(year %in% train_y2) %>%
  filter(n > 3) %>% #here we have enough data
  dplyr::select(-X, -n, -grid_id, -first_obs,
                -last_obs, -stateProvince, -gen_time, -year, -lon, -lat)

set.seed(12)
rfac3 <- randomForest(peak ~ ., 
                      data = rf_data_ac3,
                      ntree = 500,
                      ntry = round(sqrt(ncol(rf_data_ac3) - 1)), 
                      importance = TRUE, 
                      na.action = na.omit)
print(rfac3) #46.5 Var. 89.8 MSE

imp_ac3 <- importance(rfac3)
imp_ac3 <- imp_ac3[order(imp_ac3[,1], decreasing = TRUE), , drop = FALSE]
head(imp_ac3, 15)

varImpPlot(rfac3)
#jan jun tmax, n max, abs jan tmax, av jl tmean, abs jl tmax,
# f min, n prec, may max, aug tmax, abs d tmax
#jan jul jun more important recent years


plot(rf_data_ac3$peak ~ rf_data_ac3$abs_jan_tmax)
#negative linear
plot(rf_data_ac3$peak ~ rf_data_ac3$abs_jun_tmax)
#negative linear
plot(rf_data_ac3$peak ~ rf_data_ac3$abs_n_tmax)
#spread out not clear
plot(rf_data_ac3$peak ~ rf_data_ac3$av_jan_tmax)
plot(rf_data_ac3$peak ~ rf_data_ac3$av_jl_tmax)
plot(rf_data_ac3$peak ~ rf_data_ac3$av_jl_tmean)
#negative linear, spread out but still clear
plot(rf_data_ac3$peak ~ rf_data_ac3$abs_jl_tmax)
#linear spread out


set.seed(12)
rpart_model_y2 <- rpart(peak ~ ., data = rf_data_ac3, method = "anova")
rpart.plot(rpart_model_y2, cex = 0.8, main = "rpart model splits: Anthocharis cardamines peak")
cp_y2 <- as.data.frame(printcp(rpart_model_y2))
(cp[which.min(cp_y2[,"xerror"]), "CP"])

rpart_ac_prun_y2 <- prune(rpart_model_y2, 0.0132765)
par(mfrow=c(1,1))
rpart.plot(rpart_ac_prun_y2, cex = 0.8)

#3.4. LMs 2020-2024 ---------
rf_data_ac3 <- anth_c_dset %>% #412
  filter(year %in% train_y2) %>%
  filter(n > 3) %>% #here we have enough data - can do further filtering
  dplyr::select(-X, -n, -grid_id, -first_obs,
                -last_obs, -gen_time, -lon, -lat) %>% na.omit()

#Adding interaction for basic predictors from rpart: 

lm_ac3 <- lm(peak ~ 
                 abs_jan_tmax * 
                 av_n_tmax + 
                 av_d_tmax + 
                 abs_jun_tmax +
                 av_f_tmax, #excluding abs d tmax
               data = rf_data_ac3)

summary(lm_ac3) #february, av d not sign
resid <- simulateResiduals(lm_ac3) %>% plot()

lm_ac3.1 <- lm(peak ~ 
               abs_jan_tmax *
               av_n_tmax + 
               abs_jun_tmax +
               abs_d_tmax, #excl abs d tmax
             data = rf_data_ac3)

summary(lm_ac3.1) 
resid <- simulateResiduals(lm_ac3.1) %>% plot()
#dmax not sign

#can get simpler:
lm_ac3.2 <- lm(peak ~ 
                 abs_jan_tmax *
                 av_n_tmax + 
                 abs_jun_tmax, #excl abs d tmax
               data = rf_data_ac3)

summary(lm_ac3.2) 
resid <- simulateResiduals(lm_ac3.2) %>% plot()
#All sign: abs jan -2, av n m-6, abs jun -3, inter positive 0.2

AIC(lm_ac3, lm_ac3.1, lm_ac3.2)
#last one is the  best 
par(mfrow=c(2,2))
plot(lm_ac3.2) #ok  
resid <- simulateResiduals(lm_ac3.2) %>% plot() #ok

lmm_ac3 <- lmer(peak ~ 
                 abs_jan_tmax * 
                 av_n_tmax + 
                 abs_jun_tmax + (1 | year),
               data = rf_data_ac3)

#singular fit

rf_data_ac3$fitted <- fitted(lm_ac3.2)
rf_data_ac3$residuals <- resid(lm_ac3.2)

ggplot(data = rf_data_ac3, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

ggplot(data = rf_data_ac3, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ year) + 
  labs(
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

#trying per year
rf_data_ac3$year <- as.factor(rf_data_ac3$year)


pred_jan_tmax <- ggpredict(lm_ac3.2, terms = "abs_jan_tmax[all]")

pred_jun_tmax <- ggpredict(lm_ac3.2, terms = "abs_jun_tmax[all]")

pred_n_tmax <- ggpredict(lm_ac3.2, terms = "av_n_tmax[all]") 



ggplot() + 
  geom_ribbon(data = pred_jan_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jan_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac3, 
             aes(x = abs_jan_tmax, y = peak)) + 
  facet_wrap(~ year) +
  labs(x = "Absolute maximum temperature of January (°C)",
       y = "Flight peak of A. cardamines (day of the year)")
#not good

ggplot() + 
  geom_ribbon(data = pred_jun_tmax, 
              aes(x = x, ymin = conf.low,  
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jun_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac3, 
             aes(x = abs_jun_tmax, y = peak)) + 
  facet_wrap(~ year) +
  labs( x = "Absolute maximum temperature of June (°C)",
       y = "Flight peak of A. cardamines (day of the year)")
#ok


ggplot() + 
  geom_ribbon(data = pred_n_tmax, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac3, 
             aes(x = av_n_tmax, y = peak)) + 
  facet_wrap(~ year) +
  labs(x = "Average maximum temperature of November (°C)",
       y = "Flight peak of A. cardamines (day of the year)")


# More ---------

#Using the whole period: 
# RF, rpart all years 
rf_data_ac <- anth_c_dset %>%  #490 observations
  filter(n > 3) %>% #then its mainly obs.org, with inat have lower var expl.
  dplyr::select(-X, -n, -grid_id, -first_obs, 
                -last_obs, -gen_time, -year, -stateProvince, -lon, -lat
  ) 

set.seed(12)
rfac <- randomForest(peak ~ ., 
                     data = rf_data_ac,
                     ntree = 500,
                     ntry = round(sqrt(ncol(rf_data_ac) - 1)), 
                     importance = TRUE, 
                     na.action = na.omit)

print(rfac) #50 var expl, 87 MSE

imp_ac <- importance(rfac)
imp_ac <- imp_ac[order(imp_ac[,1], decreasing = TRUE), , drop = FALSE]
head(imp_ac, 15)

varImpPlot(rfac) #2 x jan tmax, abs jun tmax last, jan tmean, n tmean last, aug tmean last, 
#abs d tmax, av n tmax, av jl tmean, av n tmeam, ab jan tmean, av may tmax tmean

#use all for rpart to use the methodology without being biased
set.seed(12)
rpart_model <- rpart(peak ~ ., data = rf_data_ac, method = "anova")
rpart.plot(rpart_model, cex = 0.8, main = "rpart model splits: Anthocharis cardamines peak")
cp <- as.data.frame(printcp(rpart_model))
(cp[which.min(cp[,"xerror"]), "CP"])

rpart_ac_prun <- prune(rpart_model, 0.01)
par(mfrow=c(1,1))
rpart.plot(rpart_ac_prun, cex = 0.65)
#Jan max 

#consider interactions:
#abs jl max, abs jun max
#av jan mean, av d min
#January explains 42% of the data (abs (13) + av (5.8) max)

#observation.org data only
rf_data_ac_obs <- anth_c_obs_data %>%  #692
  dplyr::select(-X, -n, -grid_id, -first_obs, 
                -last_obs, -gen_time, -lon, -lat, -year, -stateProvince
  ) 

set.seed(12)
rfac_obs <- randomForest(peak ~ ., 
                         data = rf_data_ac_obs,
                         ntree = 500,
                         ntry = round(sqrt(ncol(rf_data_ac_obs) - 1)), 
                         importance = TRUE, 
                         na.action = na.omit)

print(rfac_obs) #49 var expl, 87 MSE - keeping all for analysis (extra 40 observations)
importance(rfac_obs) 
varImpPlot(rfac_obs)
#similar, more clear effect of winter: abs jan d max, av jan max 
plot(rf_data_ac_obs$peak ~ rf_data_ac_obs$abs_d_tmax)
plot(rf_data_ac_obs$peak ~ rf_data_ac_obs$abs_jan_tmax)
plot(rf_data_ac_obs$peak ~ rf_data_ac_obs$av_jan_tmax)
plot(rf_data_ac_obs$peak ~ rf_data_ac_obs$abs_jun_tmax)
plot(rf_data_ac_obs$peak ~ rf_data_ac_obs$av_n_tmax)
#all: very clear negative relationships

plot(rf_data_ac$peak ~ rf_data_ac$abs_d_tmax)
plot(rf_data_ac$peak ~ rf_data_ac$av_d_tmin)
#similar
plot(rf_data_ac$peak ~ rf_data_ac$abs_jan_tmax)
plot(rf_data_ac$peak ~ rf_data_ac$av_jan_tmax)
#similar
plot(rf_data_ac$peak ~ rf_data_ac$abs_jun_tmax)
plot(rf_data_ac$peak ~ rf_data_ac$av_n_tmax)
plot(rf_data_ac$peak ~ rf_data_ac$av_jl_tmean)
plot(rf_data_ac$peak ~ rf_data_ac$av_jan_tmean)
#mar feb not so clear
#all linear negative, a bit more spread out
#use linear mixed models


# LMs for all years

rf_data_ac <- anth_c_dset %>%  #723 
  filter(n > 3) %>% #then its mainly obs.org, with inat have lower var expl.
  dplyr::select(-X, -n, -grid_id, -first_obs, 
                -last_obs, -gen_time, -lon, -lat
  ) %>% na.omit()

lm_ac <- lm(peak ~ abs_d_tmax + abs_jan_tmax + av_jan_tmax + abs_jun_tmax +
              av_n_tmax + av_jl_tmean + av_jan_tmean + abs_jl_tmax, 
            data = rf_data_ac)
summary(lm_ac) #0.48 rsq 9.6
#av n tmax, abs jun tmax, av jul tmean negative effect
# -3, -2, -1
# p value < 2.2e-16
par(mfrow=c(2,2))
plot(lm_ac)
#pretty good
resid <- simulateResiduals(lm_ac)
plot(resid) #ok 

#simplify: 
lm_ac2 <- lm(peak ~ av_jan_tmax + abs_jun_tmax +
               av_n_tmax + abs_jan_tmax, 
             data = rf_data_ac)

summary(lm_ac2) #av jan tmax not sign, abs d tmax not sign, same metrics
#main predictors: jan tmax, -3.4, n max -3
resid <- simulateResiduals(lm_ac2)
plot(resid) 
#ok


lm_ac3 <- lm(peak ~ abs_jun_tmax * av_jan_tmax +
               av_n_tmax + abs_jan_tmax,
             data = rf_data_ac)
#jan mean positive effect, others nefative: jun tmax, nov tmax most imp.
summary(lm_ac3)

resid3 <- simulateResiduals(lm_ac3)
plot(resid3) #ok, linear

lm_ac4 <- lm(peak ~ av_jan_tmax * abs_jun_tmax +
               av_n_tmax + abs_jl_tmax, 
             data = rf_data_ac)
summary(lm_ac4)

resid4 <- simulateResiduals(lm_ac4)
plot(resid4) #ok not so linear

#Model diagnsotics:
AIC(lm_ac, lm_ac2, lm_ac3, lm_ac4)
#4 is the best, keep 3 because of linearity

rf_data_ac$fitted <- fitted(lm_ac3)
rf_data_ac$residuals <- resid(lm_ac3)

ggplot(data = rf_data_ac, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

ggplot(data = rf_data_ac, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ year) + 
  labs(title = "Predicted vs actual peak values Anthocharis cardamines",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 


rf_data_ac$year <- as.factor(rf_data_ac$year)

#Mixed: 
lmm_ac <- lmer(peak ~ av_jan_tmax * abs_jun_tmax +
                 av_n_tmax + abs_jan_tmax
               + (1 | year), 
               data = rf_data_ac)
summary(lmm_ac)
ranef(lmm_ac) # 2017: -6, 2021 -2, 2019 5, 2020 3
coef(lmm_ac)

resid <- simulateResiduals(lmm_ac)
plot(resid) #not good both graphs 

#predictions, most significant ones:
pred_n_max <- ggpredict(lm_ac3, terms = "av_n_tmax")
#pred_jan_tmax <- ggpredict(lm_ac4, terms = "av_jan_tmax")
#pred_jl_tmax <- ggpredict(lm_ac4, terms = "abs_jl_tmax")
pred_jun_tmax <- ggpredict(lm_ac3, terms = "abs_jun_tmax")

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac, 
             aes(x = av_n_tmax, y = peak)) + 
  facet_wrap(~ year) +
  labs(x = "Average maximum temperature of November (°C)",
       y = "Peak (days)")
#spread out but ok

ggplot() + 
  geom_ribbon(data = pred_jun_tmax, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jun_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac, 
             aes(x = abs_jun_tmax, y = peak)) + 
  facet_wrap(~year) +
  labs(x = "Absolute maximum temperature of June (°C)",
       y = "Peak (days)")
#ok

summary(lm_ac3)
summary(lm_ac3.2)

summary(lm_ac)

#trying predictions per year: 
rf_data_ac$year <- as.factor(rf_data_ac$year)

ac_mixed_plot <- rf_data_ac %>% mutate(fac_gr = year)

pred_n_max <- ggpredict(lmm_ac, terms = c("av_n_tmax[all]", "year"), 
                        type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
pred_jan_tmax <- ggpredict(lmm_ac, terms = c("av_jan_tmax[all]", "year"),
                           type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
pred_jan_tmean <- ggpredict(lmm_ac, terms = c("av_jan_tmean[all]", "year"),
                            type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
pred_jun_tmax <- ggpredict(lmm_ac, terms = c("abs_jun_tmax[all]", "year"),
                           type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)
pred_jl_tmax <- ggpredict(lmm_ac, terms = c("abs_jl_tmax[all]", "year"),
                          type = "random") %>% as.data.frame() %>% 
  mutate(fac_gr = group)

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = ac_mixed_plot, 
             aes(x = av_n_tmax, y = peak)) + 
  facet_wrap(~ fac_gr) +
  labs(x = "November maximum temperature (C)",
       y = "Peak (days)")
#spread out


ggplot() + 
  geom_ribbon(data = pred_jan_tmax,  
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jan_tmax, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = ac_mixed_plot, 
             aes(x = av_jan_tmax, y = peak)) + 
  facet_wrap(~ fac_gr) +
  labs(title = "Effects of average January max temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "January max temperature (C)",
       y = "Peak (days)")


ggplot() + 
  geom_ribbon(data = pred_jan_tmean, 
              aes(x = x, ymin = conf.low, group = fac_gr, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jan_tmean, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = ac_mixed_plot, 
             aes(x = av_jan_tmean, y = peak)) + 
  facet_wrap(~ fac_gr) +
  labs(title = "Effects of average January temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "January mean temperature (C)",
       y = "Peak (days)")
#bad 


ggplot() + 
  geom_ribbon(data = pred_jun_tmax, 
              aes(x = x, ymin = conf.low, group = fac_gr,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jun_tmax, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = ac_mixed_plot, 
             aes(x = abs_jun_tmax, y = peak)) + 
  facet_wrap(~ fac_gr) +
  labs(title = "Effects of absolute max June temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "June max abs temperature (C)",
       y = "Peak (days)")
#bad


ggplot() + 
  geom_ribbon(data = pred_jl_tmax, 
              aes(x = x, ymin = conf.low, group = fac_gr, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jl_tmax, 
            aes(x = x, y = predicted, group = fac_gr), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = ac_mixed_plot, 
             aes(x = abs_jl_tmax, y = peak)) + 
  facet_wrap(~ fac_gr, scale = "free") +
  labs(title = "Effects of absolute max July temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "July max temperature (C)",
       y = "Peak (days)")
#really bad


#Only June jan max and November are acceptable: 

lm_alt <- lm(peak ~  av_jan_tmax + 
               + abs_jun_tmax +
               av_n_tmax , 
             data = rf_data_ac)
summary(lm_alt) #0.48 10 
resid <- simulateResiduals(lm_alt) %>% plot() #ok

lmm_alt <- lmer(peak ~  abs_jun_tmax + 
                  av_n_tmax + + av_jan_tmax + (1 | stateProvince), 
                data = rf_data_ac)
summary(lmm_alt) 
resid <- simulateResiduals(lmm_alt) %>% plot()

lmm_alt2 <- lmer(peak ~  abs_jun_tmax + 
                   av_n_tmax + + av_jan_tmax + (1 | year), 
                 data = rf_data_ac)
summary(lmm_alt2) 
resid <- simulateResiduals(lmm_alt2) %>% plot() #not ok

AIC(lm_alt, lmm_alt) #same
BIC(lm_alt, lmm_alt) #lower simple

# Using the simplest model: 


rf_data_ac$fitted <- fitted(lm_alt)

ggplot(data = rf_data_ac, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines, 2016-2019",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

ggplot(data = rf_data_ac, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ stateProvince) +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines, 2016-2019",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

ggplot(data = rf_data_ac, aes(x = fitted, y = peak)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ year) +
  labs(title = "Predicted vs actual peak values Anthocharis cardamines, 2016-2019",
       x = "Fitted Values", 
       y = "Actual peak Values") +
  theme_minimal() 

pred_n_max <- ggpredict(lm_alt, terms = "av_n_tmax")
pred_jun_tmax <- ggpredict(lm_alt, terms = "abs_jun_tmax")
pred_jan_tmax <- ggpredict(lm_alt, terms = "av_jan_tmax")

ggplot() + 
  geom_ribbon(data = pred_n_max, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_n_max, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac, 
             aes(x = av_n_tmax, y = peak)) +
  labs(title = "Effects of average November max temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "November maximum temperature (C)",
       y = "Peak (days)")

ggplot() + 
  geom_ribbon(data = pred_jun_tmax, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jun_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac, 
             aes(x = abs_jun_tmax, y = peak)) + 
  labs(title = "Effects of absolute June max temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "June max temperature (C)",
       y = "Peak (days)")

ggplot() + 
  geom_ribbon(data = pred_jan_tmax, 
              aes(x = x, ymin = conf.low, 
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jan_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = rf_data_ac, 
             aes(x = av_jan_tmax, y = peak)) + 
  labs(title = "Effects of absolute max January temperature on 
       Anthocharis cardamines peak, 2016-2019", x = "January max temperature (C)",
       y = "Peak (days)")




