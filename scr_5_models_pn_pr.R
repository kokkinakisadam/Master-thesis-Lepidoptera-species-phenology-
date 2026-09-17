# Modeling flight phenology of Pieri napi and Pieris rapae in Germany

#Libraries: 

library(rnaturalearth)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(psych)
library(rnaturalearthdata)
library(raster)
library(MASS)
library(DHARMa)
library(forecast)
library(lubridate)
library(randomForest)
library(lme4)
library(mgcv)
library(rpart)
library(rpart.plot)

#1.Data -------
p_n_1_inat <- read.csv("Data_analysis/p_n_1_inat_all")
p_n_2_inat <- read.csv("Data_analysis/p_n_2_inat_all")
p_r_1_inat <- read.csv("Data_analysis/p_r_1_inat_all")
p_r_2_inat <- read.csv("Data_analysis/p_r_2_inat_all")
p_n_1_ob <- read.csv("Datasets_analysis/p_n_1_obs_all")
p_n_2_ob <- read.csv("Datasets_analysis/p_n_2_obs_all")
p_r_1_ob <- read.csv("Datasets_analysis/p_r_1_obs_all")
p_r_2_ob <- read.csv("Datasets_analysis/p_r_2_obs_all")

p_n_1 <- rbind(p_n_1_inat, p_n_1_ob) %>% distinct()
p_n_2 <- rbind(p_n_2_inat, p_n_2_ob) %>% distinct()
p_r_1 <- rbind(p_r_1_inat, p_r_1_ob) %>% distinct()
p_r_2 <- rbind(p_r_2_inat, p_r_2_ob) %>% distinct()

# 2. Plots and functions --------
germany <- ne_countries(scale = "medium", country = "Germany",
                        returnclass = "sf")
germany <- germany[6]

pr_map_1 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_r_1, col = "blue", size = 0.1, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 

 ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_r_1, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "green") 

 pr_map_2 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_r_2, col = "blue", size = 0.1, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 

ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_r_2, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "green")

grid.arrange(pr_map_1, pr_map_2, ncol = 2)
#north and south ends are filled in 2024
#mean peak? 

pn_map_1 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_1, col = "blue", size = 0.1, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 
  

 ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_1, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "green")

pn_map_2 <- ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_2, col = "blue", size = 0.06, aes(x = lon, y = lat)) + 
  labs(y = "latitude", x= "longitude") +
  theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  facet_wrap(~year) 

 ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_2, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "green")

grid.arrange(pn_map_1, pn_map_2, ncol = 2)
#from central to north and south


#histogram of 2 statistical peaks  / generations:

hist(p_n_1$peak)
hist(p_n_2$peak)

mean(p_n_1$peak)
mean(p_n_2$peak)

hist(p_r_1$peak)
hist(p_r_2$peak)
mean(p_r_1$peak)
mean(p_r_2$peak)


#Using peak as response variable
get_hotspots <- function(species_data) {
  species_data %>% 
    filter(n > 2) %>% 
    dplyr::select(-X, -grid_id.x, -grid_id.y, 
                  -lon, -lat, 
                  -year,
                   -n, -first_obs, -last_obs, -gen_time, -stateProvince)
}

p_n_1_hot <- get_hotspots(p_n_1)
p_n_2_hot <- get_hotspots(p_n_2)
p_r_1_hot <- get_hotspots(p_r_1)
p_r_2_hot <- get_hotspots(p_r_2)

par(mfrow=c(1,2))
hist(p_n_1_hot$peak, xlab = "First flight peak of P. napi (day of the year)", main = NA)
hist(p_n_2_hot$peak, xlab = "Second flight peak of P. napi (day of the year)", main = NA)


par(mfrow=c(1,2))
hist(p_r_1_hot$peak, xlab = "First flight peak of P. rapae (day of the year)", main = NA)
hist(p_r_2_hot$peak, xlab = "Second flight peak of P. rapae (day of the year)", main = NA)

p_n_1_hot_coords <- p_n_1 %>% filter(n > 2) 
p_n_2_hot_coords <- p_n_2 %>% filter(n > 2) 

ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_1_hot_coords, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "pink")

ggplot() + 
  geom_sf(data = germany, fill = NA, col = "black") +
  geom_point(data = p_n_2_hot_coords, size = 1, aes(x = lon, y = lat, color = peak)) +
  scale_color_gradient(low = "blue",  high = "pink")


#only for descriptive statistics: 

get_hotspots_d <- function(species_data){
  species_data %>% 
    filter(n > 2) 
  }
  
p_n_1_hotd <- get_hotspots_d(p_n_1)
p_n_2_hotd <- get_hotspots_d(p_n_2)
p_r_1_hotd <- get_hotspots_d(p_r_1)
p_r_2_hotd <- get_hotspots_d(p_r_2)
  
  
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

p_n_1_sum <- get_per_year_stats(p_n_1_hotd)
p_n_2_sum <- get_per_year_stats(p_n_2_hotd)
p_r_1_sum <- get_per_year_stats(p_r_1_hotd)
p_r_2_sum <- get_per_year_stats(p_r_2_hotd)

#negative relationship for the first generation
#no clear relationship for the second generation
#Earlier emergence -> longer first generation 
#Pr: a later second peak makes the second generation longer
plot(peak ~ gen_size, data = p_n_1_sum, pch = 19)
#negative rel
plot(peak ~ gen_size, data = p_n_2_sum, pch = 19)
plot(peak ~ gen_size, data = p_r_1_sum, pch = 19)
plot(peak ~ gen_size, data = p_r_2_sum, pch = 19)

#Trend: Earlier peaks - big difference 
#2023 appears very different (for other species too)
plot(peak ~ year, data = p_n_1_sum, pch = 19, ylab = "First peak (days)")
#high year follows low year
plot(peak ~ year, data = p_r_1_sum, pch = 19)

#trend: later peak of second generation
plot(peak ~ year, data = p_n_2_sum, pch = 19)
plot(peak ~ year, data = p_r_2_sum, pch = 19, ylab = "Second peak (days)")

## theres something wrong: Pr and Pn 1 same results 

#3. P. napi -----
#3.1. RF peak --------

set.seed(12)
rfpn_part1 <- randomForest(peak ~ ., #372 obs (n > 2)
                             data = p_n_1_hot, #or p_n_1_hot
                             ntree = 1000,
                             ntry = round(sqrt(ncol(p_n_1_hot) - 1)), 
                             importance = TRUE, 
                             na.action = na.omit)
print(rfpn_part1) # 21.9 102.1 with no gen time filter, acceptable / 20 76 for n > 4

imp_pn1 <- importance(rfpn_part1)
imp_pn1 <- imp_pn1[order(imp_pn1[,1], decreasing = TRUE), , drop = FALSE]
head(imp_pn1, 15)

varImpPlot(rfpn_part1)
#main: abs jan and jun max, av jun max - all negative
#more: sep max(positive), jun tmean(negative)

#june tmean tmax last, may 2x max, n mean, sep max, mar mean,
#jan mean, nov min mean jun prec last
#abs jan tmax also imp. 
plot(p_n_1_hot$peak ~ p_n_1_hot$abs_jan_tmax) #negative, curvy, 12 + promote early peaks? 
plot(p_n_1_hot$peak ~ p_n_1_hot$abs_jun_tmax) #negative, linear spread out
#plot(p_n_1_hot$peak ~ p_n_1_hot$total_jun_prec) #unclear
#plot(p_n_1_hot$peak ~ p_n_1_hot$av_jun_tmax) #negative, linear
plot(p_n_1_hot$peak ~ p_n_1_hot$abs_s_tmax) #positive, but unclear
#plot(p_n_1_hot$peak ~ p_n_1_hot$total_n_prec) #unclear
plot(p_n_1_hot$peak ~ p_n_1_hot$av_jun_tmean) #negative, a bit unclear

#try both lm with abs max of s jan jun and jun mean 

set.seed(12)
rfpn_part2 <- randomForest(peak ~ ., #1431 / 628 fr n > 4
                             data = p_n_2_hot,
                             ntree = 500,
                             ntry = round(sqrt(ncol(p_n_2_hot) - 1)), 
                             importance = TRUE, 
                             na.action = na.omit)
print(rfpn_part2)  #13.3 338.4
varImpPlot(rfpn_part2) #mar min, sep mean last, total mar prec,a v mar mean, 
#d mean, sep max last,f mean, f min, abs jl tmax last, abs mar tmax, total jun 
#prec last, mean jun prec last, abs n tmin,
#sep mean, abs f tmax
#very similar when using obs only
plot(p_n_2_hot$peak ~ p_n_2_hot$av_mar_tmin)  #positive, spread out
plot(p_n_2_hot$peak ~ p_n_2_hot$av_mar_tmean) #same
plot(p_n_2_hot$peak ~ p_n_2_hot$abs_f_tmax) #positive, 2 groups
plot(p_n_2_hot$peak ~ p_n_2_hot$total_mar_prec) #unclear
plot(p_n_2_hot$peak ~ p_n_2_hot$total_d_prec) #unclear

#positive 
#too many observations? how to filter

#3.2. Rpart -------
set.seed(12)
rpart_pn1peak <- rpart(peak ~ ., data = p_n_1_hot)
rpart.plot(rpart_pn1peak, cex = 0.8) #chaos
cp <- as.data.frame(printcp(rpart_pn1peak))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_pn1peak, 0.01987373)
rpart.plot(rpart_pruned, cex = 0.8)
#jan tmax, n tmean late peaks 
#prec: Jun, febr, jan left peak, together with abs d tmax 
#jun tmax abs main distinction

plot(peak ~ total_jun_prec, data = p_n_1_hot)
plot(peak ~ total_feb_prec, data = p_n_1_hot) #group
plot(peak ~ total_jan_prec, data = p_n_1_hot) #group
plot(peak ~ av_n_tmean, data = p_n_1_hot) #values around an optimum
plot(peak ~ av_jan_tmax, data = p_n_1_hot) #negative, closer to the right
plot(peak ~ abs_d_tmax, data = p_n_1_hot) #negative

par(mfrow = c(1,2))
plot(peak ~ abs_d_tmax, data = p_n_1_hot, 
     ylab = "First peak (days)", xlab = "Absolute maximum 
     temperature of December (°C)", cex = 0.8)
plot(peak ~ abs_jan_tmax, data = p_n_1_hot, ylab = 
     "First peak (days)", xlab = "Absolute maximum 
     temperature of January (°C)", cex = 0.8)
plot(peak ~ abs_jun_tmax, data = p_n_1_hot, 
     ylab = "First peak (days)", xlab = "Absolute maximum 
     temperature of June (°C)", cex = 0.8) #linear
plot(peak ~ total_feb_prec, data = p_n_1_hot, ylab = "First peak (days)", xlab = "Total 
     precipitation of February (mm)", cex = 0.8)



#3.3. GAMs -------

#keep year for visualization: 
p_n_1_hot <- p_n_1 %>% filter(n > 2)

library(mgcv) #use this package for general additive models

pn1_gam <- gam(peak ~ s(abs_d_tmax) +
                 s(abs_jun_tmax) + 
                 s(abs_jan_tmax) + 
                 s(total_feb_prec) + 
                 s(av_jan_tmax) + 
                 s(total_jun_prec),
               data = p_n_1_hot)
summary(pn1_gam)

#Deviance explained: 30 %, R sq adj = 0.268 (all sign but not av jan)
#all significant

pn1_gam_2 <- gam(peak ~ s(abs_d_tmax) +
                 s(abs_jun_tmax) + 
                 s(abs_jan_tmax) + 
                 s(total_feb_prec) + 
                 s(total_jun_prec),
               data = p_n_1_hot)
summary(pn1_gam_2)
#27.9 %, 0.257

AIC(pn1_gam, pn1_gam_2) #the same - 
#more complex slightly better: 0.9 difference

#can try interactions: jun max d max jun prec

pn1_gam_3 <- mgcv::gam(peak ~ s(abs_d_tmax, by = total_jun_prec) +
                   s(total_jan_prec) +
                   s(abs_jun_tmax) + 
                   s(abs_jan_tmax) + 
                   s(av_jan_tmax),
                 data = p_n_1_hot)
summary(pn1_gam_3) #28.7, 0.264, jan tmax feb not sign
#without total feb prec at all: 30.7 0.276
AIC(pn1_gam, pn1_gam_2, pn1_gam_3)
#3 is the best 3662.385

#checking fit: 
par(mfrow=c(2,2))
mgcv::gam.check(pn1_gam)
gam.check(pn1_gam_2)
gam.check(pn1_gam_3) #really good distribution and QQ
#Ks are borderline, better to include one with less borderline k values

draw(pn1_gam_3, constant = coef(pn1_gam_3)[1], caption = FALSE) #checking this one too

#d max and jan max are non linear, others negative linear
library(ggeffects)
pred_d_tmax <- ggpredict(pn1_gam_3, terms = "abs_d_tmax[all]")

pred_jun_tmax <- ggpredict(pn1_gam_3, terms = "abs_jun_tmax[all]")

pred_abs_jan_tmax <- ggpredict(pn1_gam_3, terms = "abs_jan_tmax[all]")

pred_av_jan_tmax <- ggpredict(pn1_gam_3, terms = "av_jan_tmax[all]")

pred_prec <- ggpredict(pn1_gam_3, terms = "total_jan_prec[all]")


ggplot() + 
  geom_ribbon(data = pred_d_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_d_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_1_hot, 
             aes(x = abs_d_tmax, y = peak), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Absolute maximum temperature of December (°C)",
       y = "First flight peak of P. napi (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_jun_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_jun_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_1_hot, 
             aes(x = abs_jun_tmax, y = peak), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Absolute maximum temperature of June (°C)",
       y = "First flight peak of P. napi (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_av_jan_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_av_jan_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_1_hot, 
             aes(x = av_jan_tmax, y = peak), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Average maximum temperature of January (°C)",
       y = "First flight peak of P. napi (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_abs_jan_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_abs_jan_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_1_hot, 
             aes(x = abs_jan_tmax, y = peak), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Absolute maximum temperature of January (°C)",
       y = "First flight peak of P. napi (day of the year)")


ggplot() + 
  geom_ribbon(data = pred_prec,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_prec, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_1_hot, 
             aes(x = total_jan_prec, y = peak), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Total precipitation of January (mm)",
       y = "First flight peak of P. napi (day of the year)")


#not good
#4. P. rapae ------
# 4.1. RF peak -----------
set.seed(12)
rfpr_part1 <- randomForest(peak ~ ., 
                            data = p_r_1_hot, #or p_r_1_hot
                            ntree = 500,
                            ntry = round(sqrt(ncol(p_r_1_hot) - 1)), 
                            importance = TRUE, 
                            na.action = na.omit)
print(rfpr_part1) #14.7 132.0
imp_pr1 <- importance(rfpr_part1)
imp_pr1 <- imp_pr1[order(imp_pr1[,1], decreasing = TRUE), , drop = FALSE]
head(imp_pr1, 15)

varImpPlot(rfpr_part1) #jl tmax, aug tmax, jl mean, total n prec, 
#av jun tmean, av jan tmax, abs jan tmax av n tmin
#all: total n prec, all max of summer, jun mean, early cold -4, jan prec

plot(p_r_1_hot$peak ~ p_r_1_hot$abs_jun_tmax)
#negative linear 
plot(p_r_1_hot$peak ~ p_r_1_hot$abs_jl_tmax)
#negative linear spread out
plot(p_r_1_hot$peak ~ p_r_1_hot$abs_jan_tmax) 
#relatively clear, either poly(2) or linear
plot(p_r_1_hot$peak ~ p_r_1_hot$av_jl_tmax)
#negative linear but spread out
plot(p_r_1_hot$peak ~ p_r_1_hot$av_aug_tmax)
#negative spread out 
plot(p_r_1_hot$peak ~ p_r_1_hot$av_jun_tmax)
#linear negative
plot(p_r_1_hot$peak ~ p_r_1_hot$av_jun_tmean)
#negative
plot(p_r_1_hot$peak ~ p_r_1_hot$av_n_tmin)
#grouping
plot(p_r_1_hot$peak ~ p_r_1_hot$total_n_prec)
#group

set.seed(12)
rfpr_part2 <- randomForest(peak ~ ., 
                           data = p_r_2_hot, 
                           ntree = 500,
                           ntry = round(sqrt(ncol(p_r_2_hot) - 1)), 
                           importance = TRUE, 
                           na.action = na.omit)
print(rfpr_part2)  #3.78 493.8909

#4.2. Rpart --------
set.seed(12)
rpart_pr1_peak <- rpart(peak ~ ., data = p_r_1_hot)
rpart.plot(rpart_pr1_peak, cex = 0.8)
cp <- as.data.frame(printcp(rpart_pr1_peak))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_pr1_peak, 0.1763807)
rpart.plot(rpart_pruned, cex = 1.2)

  
#5. P. napi alternative method (peak 2) ------
#Here keeping only locations with data in both flight peaks

#5.1 Data processing ---------
p_n_1 <- p_n_1 %>% 
rename(first_obs_1 = first_obs) %>%
rename(peak_1 = peak) %>%
rename(last_obs_1 = last_obs) %>%
rename(gen_time_1 = gen_time) %>% 
filter(n > 1) %>%
dplyr::select(-n, -X, -grid_id.y)
    
p_n_2 <- p_n_2 %>% 
   rename(first_obs_2 = first_obs) %>%
   rename(peak_2 = peak) %>% 
   rename(last_obs_2 = last_obs) %>%
   rename(gen_time_2 = gen_time) %>%
   filter(n > 2) %>%
   dplyr::select(-n, -X, -grid_id.y)
  
p_n_all <- inner_join(p_n_1, p_n_2, relationship = "many-to-many") 
#650 / 112 

#5.2. Random forest (2) ----------
  
p_n_rf <- p_n_all %>% 
  dplyr::select(-year, -lon, -lat, -grid_id.x, -first_obs_1, 
                -first_obs_2, -last_obs_1, -last_obs_2, -gen_time_1, 
                -gen_time_2, -stateProvince)
set.seed(12)
rfpn_part2 <- randomForest(peak_2 ~ ., 
                        data = p_n_rf, 
                        ntree = 1000,
                        ntry = round(sqrt(ncol(p_n_rf) - 1)), 
                        importance = TRUE, 
                        na.action = na.omit)
print(rfpn_part2) #17.27 269.00 #slightly better than with no filtering
imp_pn2 <- importance(rfpn_part2)
imp_pn2 <- imp_pn2[order(imp_pn2[,1], decreasing = TRUE), , drop = FALSE]
head(imp_pn2, 15)
varImpPlot(rfpn_part2) 
#other model showed: mar dec may, this model shows f, mar s
  
plot(peak_2 ~ av_f_tmax, data = p_n_rf) #positive
plot(peak_2 ~ abs_f_tmax, data = p_n_rf) #positive
plot(peak_2 ~ av_mar_tmean, data = p_n_rf) #positive
plot(peak_2 ~ av_f_tmean, data = p_n_rf) #2 groups
plot(peak_2 ~ av_jun_tmean, data = p_n_rf) #positive uncl
plot(peak_2 ~ total_feb_prec, data = p_n_rf) #pos
plot(peak_2 ~ av_s_tmax, data = p_n_rf) #positive 
plot(peak_2 ~ av_s_tmean, data = p_n_rf) #positive steep,maybe poly
  
#5.3. Rpart (2) -------
set.seed(12)
rpart_pn2_peak <- rpart(peak_2 ~ ., data = p_n_rf)
rpart.plot(rpart_pn2_peak, cex = 0.8)
cp <- as.data.frame(printcp(rpart_pn2_peak))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_pn2_peak, 0.0281312)
rpart.plot(rpart_pruned, cex = 1.1) 
#abs f tmax (11), av f tmax 9.9: cold - early! 
  
#5.4. GAM (2) ------
#not included in the thesis
    
library(mgcv)
pn2_gam <- gam(peak_2 ~ s(abs_f_tmax) +
                     s(av_f_tmax) + 
                     s(av_s_tmax) + 
                     s(av_jun_tmean),
                   data = p_n_all)
summary(pn2_gam)
 
#19.3% only s  max abs fmax sign.
  
pn2_gam2 <- mgcv::gam(peak_2 ~ s(abs_f_tmax) +
               s(av_f_tmax, by = abs_f_tmax) + 
               s(av_s_tmax) + 
               s(av_jun_tmean),
               data = p_n_all)
summary(pn2_gam2) #19.3 0.177
AIC(pn2_gam, pn2_gam2) #2 better
  
gam.check(pn2_gam2) #weird group in the left, apart from that good
#k and p-value is very low, k adjustment: 
  
pn2_gam2 <- mgcv::gam(peak_2 ~ s(abs_f_tmax) +
                s(av_f_tmax, by = abs_f_tmax) + 
                s(av_s_tmax) + 
                s(av_jun_tmean),
                data = p_n_all)
summary(pn2_gam2) 
#again k and p vaues low, terms are removed: 
  
pn2_gam3 <- mgcv::gam(peak_2 ~ s(abs_f_tmax) +
                 s(av_f_tmax, by = abs_f_tmax),
                 data = p_n_all)
summary(pn2_gam3)  #19.6, 0.177 same as before
AIC(pn2_gam, pn2_gam2, pn2_gam3) # 2/3 better
gam.check(pn2_gam3)  #ok
draw(pn2_gam3, constant = coef(pn2_gam3)[1], caption = FALSE) 
#abs is a better metric: narrower confidence interval


pred_abs_f_tmax <- ggpredict(pn2_gam3, terms = "abs_f_tmax[all]")
pred_av_f_tmax <- ggpredict(pn2_gam3, terms = "av_f_tmax[all]")

ggplot() + 
  geom_ribbon(data = pred_abs_f_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_abs_f_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_all, 
             aes(x = abs_f_tmax, y = peak_2), size = 0.9) + 
  #facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Absolute maximum temperature of February (°C)",
       y = "Second flight peak of P. napi (day of the year)")

ggplot() + 
  geom_ribbon(data = pred_av_f_tmax,  
              aes(x = x, ymin = conf.low,
                  ymax = conf.high, colour = "lightgrey", fill = NA)) +
  geom_line(data = pred_av_f_tmax, 
            aes(x = x, y = predicted), 
            color = "blue", linewidth = 1.5) +
  geom_point(data = p_n_all, 
             aes(x = av_f_tmax, y = peak_2), size = 0.9) + 
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text = element_text(size = 7), 
        panel.spacing = unit(1, "line")) +
  labs(x = "Average maximum temperature of February (°C)",
       y = "Second flight peak of P. napi (day of the year)")

#6. P. rapae alternative method-------
#6.1. Data processing -------
p_r_1 <- p_r_1 %>% 
rename(first_obs_1 = first_obs) %>%
rename(peak_1 = peak) %>%
rename(last_obs_1 = last_obs) %>%
rename(gen_time_1 = gen_time) %>% 
filter(n > 1) %>%
dplyr::select(-n, -X, -grid_id.y)
  
p_r_2 <- p_r_2 %>% 
  rename(first_obs_2 = first_obs) %>%
  rename(peak_2 = peak) %>% 
  rename(last_obs_2 = last_obs) %>%
  rename(gen_time_2 = gen_time) %>% 
  filter(n > 1) %>%
  dplyr::select(-n, -X, -grid_id.y) 
  
p_r_all <- inner_join(p_r_1, p_r_2, relationship = "many-to-many") 
  #111
  
#6.2. Random forest (2) ------
p_r_rf <- p_r_all %>% 
  dplyr::select(-year, -lon, -lat, -grid_id.x, 
                -first_obs_1, 
                -first_obs_2, 
                -last_obs_1, -last_obs_2, -gen_time_1, 
                -gen_time_2, -stateProvince)
  
rfpr_part2 <- randomForest(peak_2 ~ ., 
                      data = p_r_rf, 
                      ntree = 1000,
                      ntry = round(sqrt(ncol(p_r_rf) - 1)), 
                      importance = TRUE, 
                      na.action = na.omit)
print(rfpr_part2)  #8.5 551.9
imp_pr2 <- importance(rfpr_part2)
imp_pr2 <- imp_pr2[order(imp_pr2[,1], decreasing = TRUE), , drop = FALSE]
head(imp_pr2, 15)
varImpPlot(rfpr_part2) 

#all unclear
plot(peak_2 ~ av_d_tmin, data = p_r_rf) #curvy
plot(peak_2 ~ abs_d_tmin, data = p_r_rf)#curvy
plot(peak_2 ~ av_f_tmax, data = p_r_rf) #positive?
plot(peak_2 ~ abs_f_tmax, data = p_r_rf) #positive? 
plot(peak_2 ~ abs_may_tmax, data = p_r_rf) #curvy?
plot(peak_2 ~ av_n_tmin, data = p_r_rf) #positive curvy?
plot(peak_2 ~ total_mar_prec, data = p_r_rf) #positive
plot(peak_2 ~ av_f_tmean, data = p_r_rf) #unclear
  
  
#6.3. Rpart (2) ------ 
set.seed(12)
rpart_pr2_peak <- rpart(peak_2 ~ ., data = p_r_rf)
rpart.plot(rpart_pr2_peak, cex = 0.8)
cp <- as.data.frame(printcp(rpart_pr2_peak))
(cp[which.min(cp[,"xerror"]), "CP"])
rpart_pruned <- prune(rpart_pr2_peak, 0.118778)
rpart.plot(rpart_pruned, cex = 1.2)
#only d tmin 







