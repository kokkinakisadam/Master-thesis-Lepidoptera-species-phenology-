# Modeling flight phenology of Euproctis chrysorrhoea in Germany

#Libraries:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(raster)
library(MASS)
library(DHARMa)
library(ggeffects)
library(rnaturalearth)
library(forecast)
library(lubridate)
library(randomForest)
library(mgcv)
library(rpart)
library(rpart.plot)

#1. Data -----------
eu_chr_ob <- read.csv("Datasets_analysis/eu_chr_obs")
eu_chr_inat <- read.csv("Data_analysis/eu_chr_inat") 
eu_chr_all <- rbind(eu_chr_inat, eu_chr_ob) %>% distinct() #problem diff var, why inat 63

eu_chr_inat$first_obs <- as.integer(yday(as.Date(eu_chr_inat$first_obs)))
eu_chr_inat$last_obs <- as.integer(yday(as.Date(eu_chr_inat$last_obs)))

eu_chr_all$first_obs <- as.integer(yday(as.Date(eu_chr_all$first_obs)))
eu_chr_all$last_obs <- as.integer(yday(as.Date(eu_chr_all$last_obs)))

#2. Maps and basic statistics ---------

#Maps with occurrence records:
ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = eu_chr_all, col = "blue", size = 0.9, fill = "blue", aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 
#upper right part close to berlin hotspot

hist(eu_chr_inat$first_obs, breaks = 20)
hist(eu_chr_all$first_obs, breaks = 20, xlab = "First observation of E. chrysorrhoea (day of the year)")
eu_chr_all_1 <- eu_chr_all %>% filter(first_obs <= 150)
eu_chr_all_2 <- eu_chr_all %>% filter(first_obs > 150) %>% 
  filter(first_obs <= 190)

eu_chr_data_1 <- eu_chr_all_1 %>% dplyr::select(-X, -grid_id, 
                                                -n, -last_obs, -year,
                                                -stateProvince, -lon, -lat)
eu_chr_data_2 <- eu_chr_all_2 %>% dplyr::select(-X, -grid_id, 
                                                -n, -last_obs, -year,
                                                -stateProvince, -lon, -lat)
#3. Random forest ----------
set.seed(12)
rf_eu <- randomForest(first_obs ~ ., 
                      data = eu_chr_data_1,
                      ntree = 500,
                      ntry = round(sqrt(ncol(eu_chr_data_1) - 1)), 
                      importance = TRUE, 
                      na.action = na.omit)
print(rf_eu) 
#16.7, Var expl. 135.9 MSE
imp_eu <- importance(rf_eu)
imp_eu <- imp_eu[order(imp_eu[,1], decreasing = TRUE), , drop = FALSE]
head(imp_eu, 15)
varImpPlot(rf_eu) 
#abs d tmax, av jun tmax, total mar prec, av s tmax, abs f tmin, 
#av jan tmax, abs jan tmax

plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$abs_d_tmax) 
#negative polynomial or linear
plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$total_mar_prec, 
     xlab = "Total March precipitation (mm)", 
     ylab = "First observation - first part (days)")
#positive polynomial or linear
plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$av_jun_tmax)
#negative
plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$abs_jan_tmax)
#unclear poly negative 2 groups
plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$abs_f_tmin)
#negative polynomial
plot(eu_chr_data_1$first_obs ~ eu_chr_data_1$av_s_tmax)
#unclear, could be a 2 or 3 degree polynomial, 2 groups

set.seed(12)
rf_eu_2 <- randomForest(first_obs ~ ., 
                        data = eu_chr_data_2,
                        ntree = 500,
                        ntry = round(sqrt(ncol(eu_chr_data_2) - 1)), 
                        importance = TRUE, 
                        na.action = na.omit)
print(rf_eu_2) #0 Var expl. 

#more rain - later emerg. , hot June - early emergence

#For analysis: Exclude 2018, one observation 
eu_chr_all_1 <- eu_chr_all_1 %>% 
  filter(year %in% c(2019:2024))

#4. Rpart -------
par(mfrow=c(1,1))
set.seed(12)
rpart_eu <- rpart(first_obs ~ ., data = eu_chr_data_1)
rpart.plot(rpart_eu, cex = 0.8) 
cp_eu <- as.data.frame(printcp(rpart_eu))
(cp_eu[which.min(cp_eu[,"xerror"]), "CP"])
rpart_eu_pr <- prune(rpart_eu, 0.08953911)
par(mfrow=c(1,1))
rpart.plot(rpart_eu_pr, cex = 1.1)
#jun tmax, jl tmax, d tmean
#can do interaction: jl tmax by jun tmax, d tmean by jun tmax


plot(first_obs ~ av_jun_tmax, data = eu_chr_all_1)
plot(first_obs ~ av_jl_tmax, data = eu_chr_all_1)
plot(first_obs ~ abs_f_tmin, data = eu_chr_all_1)
#unclear

#try additive model with 2 versions (data split - no data split)

#5. GAMs ---------
#Trying one version combining results from rpart and rf and one only with rpart: 
gam_eu <- mgcv::gam(first_obs ~ s(av_jun_tmax) + 
                      s(av_jl_tmax) + 
                      s(av_d_tmean), data = eu_chr_all_1)
summary(gam_eu) #dev explained: 49.8, Rsq adj = 0.396, jl jn are sign 

#interactions: 
gam_eu_2 <- mgcv::gam(first_obs ~ s(av_jun_tmax) + 
                        s(av_jl_tmax, by = av_jun_tmax) + 
                        s(av_d_tmean, by = av_jun_tmax), data = eu_chr_all_1)
summary(gam_eu_2)  #49 %, 0.394 

#adding febr, mar prec: 

gam_eu_3 <- mgcv::gam(first_obs ~ s(av_jun_tmax) + 
                        s(av_jl_tmax, by = av_jun_tmax) + 
                        s(abs_d_tmax) + 
                        s(total_mar_prec),
                      data = eu_chr_all_1)
summary(gam_eu_3)
#58.9, 0.51 rsq: almost all significant - not total mar prec

gam_eu_4 <- mgcv::gam(first_obs ~ s(av_jun_tmax) + 
                        s(av_jl_tmax, by = av_jun_tmax) + 
                        s(av_d_tmean, by = av_jun_tmax) + 
                        s(abs_d_tmax), 
                      data = eu_chr_all_1)
summary(gam_eu_4)
#rsq 0.51, dev 58.6

AIC(gam_eu, gam_eu_2, gam_eu_3, gam_eu_4)
#3 or 4 better - the same
#3 is the best: 61.7 0.524 

par(mfrow=c(2,2))
gam.check(gam_eu_3)  #acceptable fit, better predictions at the later peaks
draw(gam_eu_3, constant = coef(gam_eu_3)[1], caption = FALSE) 
#av jun tmax shoult be excluded, low k (0.96) and negative values in the plot

gam_eu_5 <- mgcv::gam(first_obs ~ s(av_jl_tmax, by = av_jun_tmax) + 
                        s(abs_d_tmax) + 
                        s(total_mar_prec),
                      data = eu_chr_all_1)
summary(gam_eu_5) #58.4%, R sq = 0.506 
gam.check(gam_eu_5)  #qqplot not th ebest
AIC(gam_eu_3, gam_eu_5) #same
draw(gam_eu_5, constant = coef(gam_eu_5)[1], caption = FALSE) 
#better predictions

pred_jl_tmax <- ggpredict(gam_eu_5, terms = "av_jl_tmax[all]")
pred_abs_d_tmax <- ggpredict(gam_eu_5, terms = "abs_d_tmax[all]")
pred_av_jun_tmax <- ggpredict(gam_eu_5, terms = "av_jun_tmax[all]")

ggplot() + 
  geom_ribbon(data = pred_jl_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jl_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = eu_chr_all_1, 
             aes(x = av_jl_tmax, y = first_obs), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Average maximum temperature of July (°C)",
       y = "First observation of E. chrysorrhoea (day of the year)")


ggplot() + 
  geom_ribbon(data = pred_abs_d_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_abs_d_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = eu_chr_all_1, 
             aes(x = abs_d_tmax, y = first_obs), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Absolute maximum temperature of December (°C)",
       y = "First observation of E. chrysorrhoea (day of the year)")


ggplot() + 
  geom_ribbon(data = pred_av_jun_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_av_jun_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = eu_chr_all_1, 
             aes(x = av_jun_tmax, y = first_obs), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Average maximum temperature of January (°C)",
       y = "First observation of E. chrysorrhoea (day of the year)")


