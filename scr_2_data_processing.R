#Extraction and processing of species occurrence data from GBIF for 
#five Lepidoptera species in Germany. Merging with climate data

#Libraries: 

library(dplyr)
library(sf)
library(terra)
library(rnaturalearth)
library(CoordinateCleaner) 
library(Rchelsa)
library(ggplot2)
library(psych)
library(raster)
library(ncdf4)
library(MASS)
library(DHARMa)
library(forecast)
library(lubridate)
library(randomForest)
library(lme4)


# 1. Data ---------
#1.1. Reading data ---------
# Observation.org
p_napi_obs <- read.csv("Data_sp/0050642-260226173443078_pnapi.csv", sep = "",  quote = "\"", header = T) 
p_rapae_obs <- read.csv("Data_sp/0050674-260226173443078_prapae.csv", sep = "", quote = "\"", header = T)
anth_c_obs <- read.csv("Data_sp/0050704-260226173443078_ant_c.csv", sep = "",  quote = "\"", header = T)
eu_chr_obs <- read.csv("Data_sp/0050719-260226173443078_eu_chr.csv", sep = "",  quote = "\"", header = T)
gone_rh_obs <- read.csv("Data_sp/0036668-260519110011954.csv", sep = "", quote = "\"", header = T)

# iNaturalist
eu_chr_inat <- read.csv("Data_sp/0036491-260519110011954.csv", sep = "", quote = "\"", header = T)
anth_c_inat <- read.csv("Data_sp/0036548-260519110011954.csv", sep = "", quote = "\"", header = T)
p_napi_inat <- read.csv("Data_sp/0036557-260519110011954.csv", sep = "", quote = "\"", header = T)
p_rapae_inat <- read.csv("Data_sp/0036573-260519110011954.csv",sep = "", quote = "\"", header = T)
gone_rh_inat <- read.csv("Data_sp/0036605-260519110011954.csv", sep = "", quote = "\"", header = T)

# Climate data
all_clim_res <- read.csv("Datasets_analysis/all_clim_res_new")

#checking lon and lat: 
str(p_napi_obs$decimalLatitude) 
str(p_napi_obs$decimalLongitude)
#character 

# 1.2. Cleaning data --------
#function to format coordinates before processing: 

cleaning_coords <- function(df) {
  df %>% 
    mutate(
      decimalLongitude = as.numeric(as.character(decimalLongitude)), #converting the character values
      decimalLatitude = as.numeric(as.character(decimalLatitude))
    ) %>% 
    filter( 
      !is.na(decimalLongitude), 
      !is.na(decimalLatitude)) %>% 
    cc_dupl() %>% #from coordinates cleaner package
    cc_zero() %>%
    cc_sea() %>%
    cc_outl()
}

p_napi_obs <- cleaning_coords(p_napi_obs)
p_rapae_obs <- cleaning_coords(p_rapae_obs)
anth_c_obs <- cleaning_coords(anth_c_obs)
eu_chr_obs <- cleaning_coords(eu_chr_obs)
gone_rh_obs <- cleaning_coords(gone_rh_obs)

p_napi_inat <- cleaning_coords(p_napi_inat)
p_rapae_inat <- cleaning_coords(p_rapae_inat)
anth_c_inat <- cleaning_coords(anth_c_inat)
eu_chr_inat <- cleaning_coords(eu_chr_inat)
gone_rh_inat <- cleaning_coords(gone_rh_inat)


#1.3. Climate data --------
all_tmax <- readRDS("Data_clim/results_tmax2") 
all_tmin <- readRDS("Data_clim/results_tmin2") 
all_tmean <- readRDS("Data_clim/results_tmean2") 
all_prec <- readRDS("Data_clim/results_prec2") 


# 2. Grid ----
# 2.1 Grid for Germany ------

# get Germany shape file
germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")
germany <- germany[6]
crs(germany) #4326
extent <- ext(germany)
# transform to UTM Zone 32N to have a equal area grid 
germany_utm <- st_transform(germany, 32632)

# 5 km grid for analysis 
grid_5_utm <- st_make_grid(germany_utm, cellsize = 5000, square = TRUE)
# adding unique grid IDs
grid_5_utm <- st_sf(grid_id = 1:length(grid_5_utm), geometry = grid_5_utm)
grid_5_utm <- st_intersection(grid_5_utm, germany_utm)
# extracting grid centroids
grid_5_utm$centroid <- st_centroid(grid_5_utm$geometry)

# 2.2. grid for analysis 4326 transformation: -------

grid_5_proj <- st_transform(grid_5_utm, 4326)

grid_5_proj$centroid <- st_centroid(grid_5_proj$geometry)
grid_5_proj$lon <- st_coordinates(grid_5_proj$centroid)[, 1]
grid_5_proj$lat <- st_coordinates(grid_5_proj$centroid)[, 2]

#processing butterfly and moth data to get dates and locations on grid:

#3. Editting species data -------
# 3.1. adding geometry -----

p_napi_obs$geometry <- paste(p_napi_obs$decimalLongitude, p_napi_obs$decimalLatitude, sep = ", ")
p_rapae_obs$geometry <- paste(p_rapae_obs$decimalLongitude, p_rapae_obs$decimalLatitude, sep = ", ")
anth_c_obs$geometry <- paste(anth_c_obs$decimalLongitude, anth_c_obs$decimalLatitude, sep = ", ")
eu_chr_obs$geometry <- paste(eu_chr_obs$decimalLongitude, eu_chr_obs$decimalLatitude, sep = ", ")
gone_rh_obs$geometry <- paste(gone_rh_obs$decimalLongitude, gone_rh_obs$decimalLatitude, sep = ", ")

p_napi_inat$geometry <- paste(p_napi_inat$decimalLongitude, p_napi_inat$decimalLatitude, sep = ", ")
p_rapae_inat$geometry <- paste(p_rapae_inat$decimalLongitude, p_rapae_inat$decimalLatitude, sep = ", ")
anth_c_inat$geometry <- paste(anth_c_inat$decimalLongitude, anth_c_inat$decimalLatitude, sep = ", ")
eu_chr_inat$geometry <- paste(eu_chr_inat$decimalLongitude, eu_chr_inat$decimalLatitude, sep = ", ")
gone_rh_inat$geometry <- paste(gone_rh_inat$decimalLongitude, gone_rh_inat$decimalLatitude, sep = ", ")

# 3.2. Editting with dates -------------

#Compiled: 
get_coords_dates <- function(species_data) {
  species_data$eventDate <- as.Date(species_data$eventDate)
  species_coords_dates <- species_data %>% dplyr::select(eventDate, decimalLongitude, decimalLatitude, 
                                                    stateProvince) %>%
    distinct() %>% as.data.frame() %>% 
    na.omit()
  return(species_coords_dates)
}

p_napi_coords_dates_o <- get_coords_dates(p_napi_obs)
p_rapae_coords_dates_o <- get_coords_dates(p_rapae_obs) 
anth_c_coords_dates_o <- get_coords_dates(anth_c_obs)
eu_chr_coords_dates_o <- get_coords_dates(eu_chr_obs)
gone_rh_coords_dates_o <- get_coords_dates(gone_rh_obs)

p_napi_coords_dates_i <- get_coords_dates(p_napi_inat)
p_rapae_coords_dates_i <- get_coords_dates(p_rapae_inat) 
anth_c_coords_dates_i <- get_coords_dates(anth_c_inat)
eu_chr_coords_dates_i <- get_coords_dates(eu_chr_inat)
gone_rh_coords_dates_i <- get_coords_dates(gone_rh_inat)


#3.3. spatial data transformation -------
p_napi_dist_o <- st_as_sf(p_napi_coords_dates_o, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) 
p_rapae_dist_o <- st_as_sf(p_rapae_coords_dates_o, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
anth_c_dist_o <- st_as_sf(anth_c_coords_dates_o, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
eu_chr_dist_o <- st_as_sf(eu_chr_coords_dates_o, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
gone_rh_dist_o <- st_as_sf(gone_rh_coords_dates_o, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)


p_napi_dist_i <- st_as_sf(p_napi_coords_dates_i, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) 
p_rapae_dist_i <- st_as_sf(p_rapae_coords_dates_i, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
anth_c_dist_i <- st_as_sf(anth_c_coords_dates_i, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
eu_chr_dist_i <- st_as_sf(eu_chr_coords_dates_i, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
gone_rh_dist_i <- st_as_sf(gone_rh_coords_dates_i, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)


#3.3. merging species with grid -------

p_napi_grid_o <- st_join(p_napi_dist_o, grid_5_proj, join = st_within) %>% na.omit()
p_rapae_grid_o <- st_join(p_rapae_dist_o, grid_5_proj, join = st_within) %>% na.omit()
anth_c_grid_o <- st_join(anth_c_dist_o, grid_5_proj, join = st_within) %>% na.omit()
eu_chr_grid_o <- st_join(eu_chr_dist_o, grid_5_proj, join = st_within) %>% na.omit()
gone_rh_grid_o <- st_join(gone_rh_dist_o, grid_5_proj, join = st_within) %>% na.omit()

p_napi_grid_i <- st_join(p_napi_dist_i, grid_5_proj, join = st_within) %>% na.omit()
p_rapae_grid_i <- st_join(p_rapae_dist_i, grid_5_proj, join = st_within) %>% na.omit()
anth_c_grid_i <- st_join(anth_c_dist_i, grid_5_proj, join = st_within) %>% na.omit()
eu_chr_grid_i <- st_join(eu_chr_dist_i, grid_5_proj, join = st_within) %>% na.omit()
gone_rh_grid_i <- st_join(gone_rh_dist_i, grid_5_proj, join = st_within) %>% na.omit()


# 5. Creating merged datasets ----------

#summarize per grid
#return NA for grids with no observations 
#inner_join for each species 

#yearly summary
# groupping species: 

#5.1. Yearly summary for analysis!! -------
get_sp_y_summary <- function(species_grid_5) {
  sp_y_summary <- species_grid_5 %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>%
    group_by(year, grid_id, lon, lat) %>%
    summarise(
      n = n(),  
      #stateProvince = stateProvince,
      first_obs = min(eventDate), 
      last_obs = max(eventDate), 
      gen_time = as.integer(max(last_obs) - min(first_obs)),
      .groups = "drop") %>% #corrected gen time
    arrange(year, grid_id, lon, lat) %>% 
    distinct()
  return(sp_y_summary)
}

p_napi_y_sum_o <- get_sp_y_summary(p_napi_grid_o)
p_rapae_y_sum_o <- get_sp_y_summary(p_rapae_grid_o)
anth_c_y_sum_o <- get_sp_y_summary(anth_c_grid_o)
eu_chr_y_sum_o <- get_sp_y_summary(eu_chr_grid_o)
gone_rh_y_sum_o <- get_sp_y_summary(gone_rh_grid_o)

plot(gone_rh_y_sum_o$first_obs, breaks = 30)

p_napi_y_sum_i <- get_sp_y_summary(p_napi_grid_i)
p_rapae_y_sum_i <- get_sp_y_summary(p_rapae_grid_i)
anth_c_y_sum_i <- get_sp_y_summary(anth_c_grid_i)
eu_chr_y_sum_i <- get_sp_y_summary(eu_chr_grid_i)
gone_rh_y_sum_i <- get_sp_y_summary(gone_rh_grid_i)


#5.2. Country - states summary --------
#can do that with the dataset 2 (with counts) too

get_sp_de_summary <- function(sp_y_summary) {
  sp_de_summary <- sp_y_summary %>%
    filter(year %in% c(2016:2025)) %>%
    group_by(year) %>%
    summarise(
      n = n(),
      mean_gen_time = as.integer(mean(as.numeric(last_obs - first_obs)), na.rm = T),
      first_obs_de = min(first_obs, na.rm = T), 
      last_obs_de = max(last_obs, na.rm = T), 
      gen_time = as.integer(max(last_obs_de) - min(first_obs_de))
    )
  return(sp_de_summary)
}
  
p_napi_de_summary <- get_sp_de_summary(p_napi_y_sum_o)
p_rapae_de_summary <- get_sp_de_summary(p_rapae_y_sum_o)
anth_c_de_summary <- get_sp_de_summary(anth_c_y_sum_o)
eu_chr_de_summary <- get_sp_de_summary(eu_chr_y_sum_o)
gone_rh_de_summary <- get_sp_de_summary(gone_rh_y_sum_o)

p_napi_de_summary <- get_sp_de_summary(p_napi_y_sum_i)
p_rapae_de_summary <- get_sp_de_summary(p_rapae_y_sum_i)
anth_c_de_summary <- get_sp_de_summary(anth_c_y_sum_i)
eu_chr_de_summary <- get_sp_de_summary(eu_chr_y_sum_i)
gone_rh_de_summary <- get_sp_de_summary(gone_rh_y_sum_i)

par(mfrow=c(2,2))
plot(mean_gen_time ~ year, data = p_napi_de_summary, pch = 20, main = "Pieris napi", ylab = "Mean generation length (days)")
plot(mean_gen_time ~ year, data = p_rapae_de_summary, pch = 20, main = "Pieris rapae", ylab = "Mean generation length (days)")
plot(mean_gen_time ~ year, data = anth_c_de_summary, pch = 20, main = "Anthocharis cardamines", ylab = "Mean generation length (days)")
plot(mean_gen_time ~ year, data = gone_rh_de_summary, pch = 20, main = "Genopteryx rhamni", ylab = "Mean generation length (days)")

#p napi and p rapi clear increase, others unstable 

get_state_summary <- function(sp_y_summary) {
  sp_state_summary <- sp_y_summary %>%
    group_by(year, stateProvince) %>%
    summarise(
      n = n(),
      mean_gen_time = as.integer(mean(as.numeric(last_obs - first_obs)), na.rm = T),
      first_obs_de = min(first_obs, na.rm = T), 
      last_obs_de = max(last_obs, na.rm = T), 
      gen_time = as.integer(max(last_obs_de) - min(first_obs_de)), 
      .groups = "drop"
    )
  return(sp_state_summary)
}

p_napi_state_summary <- get_state_summary(p_napi_y_sum_o)
p_rapae_state_summary <- get_state_summary(p_rapae_y_sum_o)
anth_c_state_summary <- get_state_summary(anth_c_y_sum_o)
eu_chr_state_summary <- get_state_summary(eu_chr_y_sum_o)
gone_rh_state_summary <- get_state_summary(gone_rh_y_sum_o)


#5.3. Flight peaks - country scale -----
get_de_peaks <- function(species_grid_5) {
  species_peaks <- species_grid_5 %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>%
    mutate(day = day(eventDate)) %>%
    group_by(year, month, day) %>%
    summarise(
      n = n(),  
      eventDate = eventDate,
      .groups = "drop") %>%
    arrange(month, day) 
  return(species_peaks)
}

p_napi_peaks_o <- get_de_peaks(p_napi_grid_o)
p_rapae_peaks_o <- get_de_peaks(p_rapae_grid_o)
anth_c_peaks_o <- get_de_peaks(anth_c_grid_o)
eu_chr_peaks_o <- get_de_peaks(eu_chr_grid_o)
gone_rh_peaks_o <- get_de_peaks(gone_rh_grid_o)


p_napi_peaks_i <- get_de_peaks(p_napi_grid_i)
p_rapae_peaks_i <- get_de_peaks(p_rapae_grid_i)
anth_c_peaks_i <- get_de_peaks(anth_c_grid_i)
eu_chr_peaks_i <- get_de_peaks(eu_chr_grid_i)
gone_rh_peaks_i <- get_de_peaks(gone_rh_grid_i)


#summary - only counts: 

get_counts_y <- function(sp_de_sum) { 
  spec_counts_y <- sp_de_sum %>% 
    filter(year %in% c(2016:2025)) %>%
    group_by(year) %>% 
    summarise(
     n = n(), 
     .groups = "drop") %>% 
       arrange(year)
     return(spec_counts_y)
}

p_napi_counts_y <- get_counts_y(p_napi_y_sum_o) 
p_rapae_counts_y <- get_counts_y(p_rapae_y_sum_o)
anth_c_counts_y <- get_counts_y(anth_c_y_sum_o)
eu_chr_counts_y <- get_counts_y(eu_chr_y_sum_o)
gone_rh_counts_y <- get_counts_y(gone_rh_y_sum_o)

p_napi_counts_y_i <- get_counts_y(p_napi_y_sum_i) 
p_rapae_counts_y_i <- get_counts_y(p_rapae_y_sum_i)
anth_c_counts_y_i <- get_counts_y(anth_c_y_sum_i)
eu_chr_counts_y_i <- get_counts_y(eu_chr_y_sum_i)
gone_rh_counts_y_i <- get_counts_y(gone_rh_y_sum_i)

par(mfrow=c(1,2))
plot(p_napi_counts_y, pch = 19, 
     ylab = "Total number of observations")
plot(p_napi_counts_y_i, pch = 19, 
     ylab = "Total number of observations")

plot(p_rapae_counts_y, main = "Pieris rapae - observation.org", pch = 19, 
     ylab = "Total number of observations")
plot(p_rapae_counts_y_i, main = "Pieris rapae - iNaturalist", pch = 19, 
     ylab = "Total number of observations")


plot(anth_c_counts_y, pch = 19, 
     ylab = "Total number of observations on observation.org (A. cardamines)")
plot(anth_c_counts_y_i, pch = 19, 
     ylab = "Total number of observations on iNaturalist (A. cardamines)")

plot(gone_rh_counts_y, pch = 19,
     ylab = "Total number of observations on observation.org (G. rhamni)")
plot(gone_rh_counts_y_i, pch = 19,
     ylab = "Total number of observations on iNaturalist (G. rhamni)")


#p napi shows some instability

#Monthly plots
get_de_peaks2 <- function(species_grid_5) {
  species_peaks2 <- species_grid_5 %>% 
    st_drop_geometry() %>%
    mutate(month = month(eventDate)) %>%
    mutate(year = year(eventDate)) %>% 
    filter(year %in% c(2016:2025)) %>% 
    group_by(month) %>%
    summarise(
      n = n(),
      .groups = "drop") %>%
    arrange(month) 
  return(species_peaks2)
}

p_napi_peaks2_o <- get_de_peaks2(p_napi_grid_o)
p_rapae_peaks2_o <- get_de_peaks2(p_rapae_grid_o)
anth_c_peaks2_o <- get_de_peaks2(anth_c_grid_o)
eu_chr_peaks2_o <- get_de_peaks2(eu_chr_grid_o)
gone_rh_peaks2_o <- get_de_peaks2(gone_rh_grid_o)

p_napi_peaks2_i <- get_de_peaks2(p_napi_grid_i)
p_rapae_peaks2_i <- get_de_peaks2(p_rapae_grid_i)
anth_c_peaks2_i <- get_de_peaks2(anth_c_grid_i)
eu_chr_peaks2_i <- get_de_peaks2(eu_chr_grid_i)
gone_rh_peaks2_i <- get_de_peaks2(gone_rh_grid_i)


plot(gone_rh_peaks_i$n ~ gone_rh_peaks_i$month)
plot(gone_rh_peaks_o$n ~ gone_rh_peaks_o$month)
plot(gone_rh_peaks2_i$n ~ gone_rh_peaks2_i$month)
plot(gone_rh_peaks2_i$n ~ gone_rh_peaks2_i$month)
#2-5 
#6-9

plot(p_napi_peaks_o$n ~ p_napi_peaks_o$eventDate)
plot(p_napi_peaks_o$n ~ p_napi_peaks_o$month) #4, 7
plot(p_napi_peaks2_o$n ~ p_napi_peaks2_o$month, pch = 19)
#1st peak: 3-5 
#2nd peak: 6-10

plot(p_rapae_peaks_o$n ~ p_rapae_peaks_o$eventDate)
plot(p_rapae_peaks_o$n ~ p_rapae_peaks_o$month) #small peak 3-5, bigger 6-11
plot(p_rapae_peaks2_o$n ~ p_rapae_peaks2_o$month, pch = 19)

plot(anth_c_peaks_o$n ~ anth_c_peaks_o$eventDate)
plot(anth_c_peaks_o$n ~ anth_c_peaks_o$month) #4
plot(anth_c_peaks2_o$n ~ anth_c_peaks2_o$month, pch = 19)

#Do the different apps capture the same phenology pattern? 
par(mfrow=c(1,3))
plot(p_napi_peaks2_o, pch = 19)
plot(p_napi_peaks2_i, pch = 19)

plot(p_rapae_peaks2_o, pch = 20)
plot(p_rapae_peaks2_i, pch = 20)

plot(anth_c_peaks2_o, pch = 20)
plot(anth_c_peaks2_i, pch = 20)

plot(gone_rh_peaks2_o, pch = 20)
plot(gone_rh_peaks2_i, pch = 20)


# 5.4. Peaks stateProvince daily --------
#not used in final analysis
get_state_peaks <- function(species_grid_5) {
  species_state_peaks <- species_grid_5 %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>%
    mutate(day = day(eventDate)) %>%
    group_by(stateProvince, year, month, day) %>%
    summarise(
      n = n(), 
      eventDate = eventDate,
      .groups = "drop") %>%
    arrange(month, day) 
  return(species_state_peaks)
}

p_napi_state_peaks <- get_state_peaks(p_napi_grid_o)
p_rapae_state_peaks <- get_state_peaks(p_rapae_grid_o)
anth_c_state_peaks <- get_state_peaks(anth_c_grid_o)
eu_chr_state_peaks <- get_state_peaks(eu_chr_grid_o)
gone_rh_state_peaks <- get_state_peaks(gone_rh_grid_o)


plot_function <- function(species_state_peaks) { 
  ggplot(species_state_peaks, aes(x = month, y = n)) + 
    geom_point() + 
    facet_wrap(~stateProvince) + 
    labs(x = "Month", y = "Observations")
}

p_napi_plot <- plot_function(p_napi_state_peaks) %>% plot()
p_rapae_plot <- plot_function(p_rapae_state_peaks) %>% plot()
anth_c_plot <- plot_function(anth_c_state_peaks) %>% plot()
eu_chr_plot <- plot_function(eu_chr_state_peaks) %>% plot()
gone_rh_plot <- plot_function(gone_rh_state_peaks) %>% plot()

# Yearly all DE plots

#good data: P napi: Bayern, Niedersachsen, Nordrhein West,
#Schleswig-Holstein, Rheinland-Pfalz, Baden Wüttemberg
#P rapae: Bayern, Schleswig Holstein, Nordr Wrst Niedersachsen

#Peaks state province monthly 
get_state_peaks_m <- function(species_grid_5) {
  species_state_peaks <- species_grid_5 %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>%
    group_by(stateProvince, year, month) %>%
    summarise(
      n = n(),  
      .groups = "drop") %>%
    arrange(year, month) 
  return(species_state_peaks)
}

p_napi_state_peaks_m <- get_state_peaks_m(p_napi_grid_o)
p_rapae_state_peaks_m <- get_state_peaks_m(p_rapae_grid_o)
anth_c_state_peaks_m <- get_state_peaks_m(anth_c_grid_o)
eu_chr_state_peaks_m <- get_state_peaks_m(eu_chr_grid_o)
gone_rh_state_peaks_m <- get_state_peaks_m(gone_rh_grid_o)


plot_function <- function(species_state_peaks) { 
  ggplot(species_state_peaks, aes(x = month, y = n)) + 
    geom_point() + 
    facet_wrap(~stateProvince) + 
    labs(x = "Month", y = "Observations")
}

p_napi_plot <- plot_function(p_napi_state_peaks) %>% plot()
p_rapae_plot <- plot_function(p_rapae_state_peaks) %>% plot()
anth_c_plot <- plot_function(anth_c_state_peaks) %>% plot()
eu_chr_plot <- plot_function(eu_chr_state_peaks) %>% plot()
gone_rh_plot <- plot_function(gone_rh_state_peaks_m) %>% plot()


# 5.5. splitting P.n., P.r., G.r., data: --------
get_peaks_pn_1 <- function(sp_grid) {
    sp_grid_new <- sp_grid %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>% 
    filter(month %in% c(3,4,5))
    return(sp_grid_new)
}

get_peaks_pn_2 <- function(sp_grid) {
  sp_grid_new <- sp_grid %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>% 
    filter(month %in% c(6,7,8,9,10))
    return(sp_grid_new)
}

p_napi_1_o <- get_peaks_pn_1(p_napi_grid_o)
p_napi_1_i <- get_peaks_pn_1(p_napi_grid_i)

p_napi_2_o <- get_peaks_pn_2(p_napi_grid_o)
p_napi_2_i <- get_peaks_pn_2(p_napi_grid_i)

p_rapae_1_o <- get_peaks_pn_1(p_rapae_grid_o)
p_rapae_1_i <- get_peaks_pn_1(p_rapae_grid_i)

p_rapae_2_o <- get_peaks_pn_2(p_rapae_grid_o)
p_rapae_2_i <- get_peaks_pn_2(p_rapae_grid_i)


#G.r. : different months

get_peaks_gr_1 <- function(sp_grid) {
  sp_grid_new <- sp_grid %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>% 
    filter(month %in% c(2,3,4,5))
  return(sp_grid_new)
}

get_peaks_gr_2 <- function(sp_grid) {
  sp_grid_new <- sp_grid %>% 
    st_drop_geometry() %>%
    mutate(year = year(eventDate)) %>%
    mutate(month = month(eventDate)) %>% 
    filter(month %in% c(6,7,8,9))
  return(sp_grid_new)
}

gone_rh_1_o <- get_peaks_gr_1(gone_rh_grid_o)
gone_rh_1_i <- get_peaks_gr_1(gone_rh_grid_i)

gone_rh_2_o <- get_peaks_gr_2(gone_rh_grid_o)
gone_rh_2_i <- get_peaks_gr_2(gone_rh_grid_i)

#6. Climate data -----
# 6.1. Climate data filtering ------
all_grids_o <- bind_rows(p_napi_y_sum_o %>% dplyr::select(grid_id, lon, lat),
                       p_rapae_y_sum_o %>% dplyr::select(grid_id, lon, lat), 
                       anth_c_y_sum_o %>% dplyr::select(grid_id, lon, lat), 
                       eu_chr_y_sum_o %>% dplyr::select(grid_id, lon, lat),
                       gone_rh_y_sum_o %>% dplyr::select(grid_id, lon, lat)) %>%
  distinct()
                       
all_grids_i <- bind_rows(p_napi_y_sum_i %>% dplyr::select(grid_id, lon, lat),
                         p_rapae_y_sum_i %>% dplyr::select(grid_id, lon, lat), 
                         anth_c_y_sum_i %>% dplyr::select(grid_id, lon, lat), 
                         eu_chr_y_sum_i %>% dplyr::select(grid_id, lon, lat),
                         gone_rh_y_sum_i %>% dplyr::select(grid_id, lon, lat)) %>%
  distinct()


all_grids <- bind_rows(all_grids_i, all_grids_o) %>% distinct()

all_prec_gr <- all_prec %>% 
  filter(grid_id %in% all_grids$grid_id) 
all_tmax_gr <- all_tmax %>% 
  filter(grid_id %in% all_grids$grid_id) 
all_tmin_gr <- all_tmin %>% 
  filter(grid_id %in% all_grids$grid_id) 
all_tmean_gr <- all_tmean %>% 
  filter(grid_id %in% all_grids$grid_id) 


# 6.2. Climate summaries ---------
prec_subset <- all_prec_gr %>% filter(month %in% c(1,2,3,5,6,7,11,12))

prec_subs_summary <- prec_subset %>%
group_by(year, grid_id, lon, lat) %>%
  summarise(
    total_jan_prec = sum(prec[month == 1], na.rm = TRUE), 
    total_feb_prec = sum(prec[month == 2], na.rm = TRUE), 
    total_mar_prec = sum(prec[month == 3], na.rm = TRUE),
    .groups = "drop" 
  )

prec_subs_summer_last <- prec_subset %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    total_may_prec = sum(prec[month == 5], na.rm = TRUE), 
    total_jun_prec = sum(prec[month == 6], na.rm = TRUE), 
    total_jl_prec = sum(prec[month == 7], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(year = year + 1)
#4 variables (x2)
    
prec_wint_last <- prec_subset %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    total_n_prec = sum(prec[month == 11], na.rm = TRUE), 
    total_d_prec = sum(prec[month == 12], na.rm = TRUE),
    .groups = "drop") %>% 
  mutate(year = year + 1)

prec_res <- inner_join(prec_subs_summer_last, prec_subs_summary,
                       by = c("year", "grid_id", "lon", "lat")) %>% 
  inner_join(prec_wint_last, by = c("year", "grid_id", "lon", "lat"))


tmin_subset <- all_tmin_gr %>% filter(month %in% c(1,2,3,11,12))
length(which(is.na(tmin_subset)))


tmin_subs_summary <- tmin_subset %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    av_jan_tmin = mean(tmin[month == 1], na.rm = TRUE),
    av_f_tmin = mean(tmin[month == 2], na.rm = TRUE),
    av_mar_tmin = mean(tmin[month == 3], na.rm = TRUE),
    abs_jan_tmin = min(tmin[month == 1]),
    abs_f_tmin = min(tmin[month == 2]),
    abs_mar_tmin = min(tmin[month == 3]),
    .groups = "drop"
  )


tmin_subs_summary_last <- tmin_subset %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    av_n_tmin = mean(tmin[month == 11], na.rm = TRUE), 
    av_d_tmin = mean(tmin[month == 12], na.rm = TRUE), 
    abs_n_tmin = min(tmin[month == 11]), 
    abs_d_tmin = min(tmin[month == 12]), 
    .groups = "drop"
  ) %>% 
  mutate(year = year + 1)
#14 variables (10 + 4)

tmin_res <- inner_join(tmin_subs_summary, tmin_subs_summary_last,
                       by = c("year", "grid_id", "lon", "lat"))

tmin_subset_last <- tmin_subset %>% filter(month %in% c(11,12)) %>%
  mutate(year = year + 1) %>% filter(year %in% c(2016:2024))

tmin_subset_late <- tmin_subset %>% filter(month %in% c(1,2,3)) %>%
  filter(year %in% c(2016:2024))


tmax_subs_summary <- all_tmax_gr %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    av_jan_tmax = mean(tmax[month == 1], na.rm = TRUE),
    av_f_tmax = mean(tmax[month == 2], na.rm = TRUE),
    av_mar_tmax = mean(tmax[month == 3], na.rm = TRUE),
    abs_jan_tmax = max(tmax[month == 1]),
    abs_f_tmax = max(tmax[month == 2]),
    abs_mar_tmax = max(tmax[month == 3]), 
    .groups = "drop"
  ) %>% 
  filter(year %in% 2016:2024)

tmax_subs_summary_last <- all_tmax_gr %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    av_may_tmax = mean(tmax[month == 5], na.rm = TRUE),
    av_jun_tmax = mean(tmax[month == 6], na.rm = TRUE),
    av_jl_tmax = mean(tmax[month == 7], na.rm = TRUE),
    av_aug_tmax = mean(tmax[month == 8], na.rm = TRUE), 
    av_s_tmax = mean(tmax[month == 9], na.rm = TRUE), 
    av_n_tmax = mean(tmax[month == 11], na.rm = TRUE),
    av_d_tmax = mean(tmax[month == 12], na.rm = TRUE),
    abs_may_tmax = max(tmax[month == 5]), 
    abs_jun_tmax = max(tmax[month == 6]),
    abs_jl_tmax = max(tmax[month == 7]),
    abs_aug_tmax = max(tmax[month == 8]), 
    abs_s_tmax = max(tmax[month == 9]),
    abs_n_tmax = max(tmax[month == 11]),
    abs_d_tmax = max(tmax[month == 12]),
    .groups = "drop"
  )  %>% 
  mutate(year = year + 1) %>% 
  filter(year %in% 2016:2024)

tmax_res <- inner_join(tmax_subs_summary, tmax_subs_summary_last,
                       by = c("year", "grid_id", "lon", "lat"))
  
#keep all in mean, select some for current and last year
tmean_subs_summary <- all_tmean_gr %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
    av_jan_tmean = mean(tmean[month == 1], na.rm = TRUE),
    av_f_tmean = mean(tmean[month == 2], na.rm = TRUE),
    av_mar_tmean = mean(tmean[month == 3], na.rm = TRUE),
    .groups = "drop"
  )

tmean_subs_summary_last <- all_tmean_gr %>%
  group_by(year, grid_id, lon, lat) %>%
  summarise(
av_may_tmean = mean(tmean[month == 5], na.rm = TRUE),
av_jun_tmean = mean(tmean[month == 6], na.rm = TRUE),
av_jl_tmean = mean(tmean[month == 7], na.rm = TRUE),
av_aug_tmean = mean(tmean[month == 8], na.rm = TRUE), 
av_s_tmean = mean(tmean[month == 9], na.rm = TRUE), 
av_oc_tmean = mean(tmean[month == 10], na.rm = TRUE), 
av_n_tmean = mean(tmean[month == 11], na.rm = TRUE),
av_d_tmean = mean(tmean[month == 12], na.rm = TRUE),
.groups = "drop" 
  ) %>% 
  mutate(year = year + 1)

#10 variables (8 + 5)
tmean_res <- inner_join(tmean_subs_summary, tmean_subs_summary_last,
                       by = c("year", "grid_id", "lon", "lat"))

all_clim_res <- prec_res %>% 
  inner_join(tmax_res, 
             by = c("year", "grid_id", "lon", "lat")) %>%
  inner_join(tmin_res, 
             by = c("year", "grid_id", "lon", "lat")) %>%
  inner_join(tmean_res, 
             by = c("year", "grid_id", "lon", "lat"))
#write.csv(all_clim_res)

#merged clim 1: merge last year: dec, nov, rename
#4 + 4 + 5 = 13
#mergeed clim 2: this year: jan-sep (use merging method that keeps all cols)
# 8 + 14 + 10 + 4 = 36
#merge 1 + 2 - all clim (3)
#49 var


# 6.3 Merging climate and species summaries ------

anth_c_obs_all <- inner_join(anth_c_y_sum_o, all_clim_res, by = c("year", "grid_id", "lon", "lat"))
eu_chr_obs_all <- inner_join(eu_chr_y_sum_o, all_clim_res, by = c("year", "grid_id", "lon", "lat"))
anth_c_inat_all <- inner_join(anth_c_y_sum_i, all_clim_res, by = c("year", "grid_id", "lon", "lat"))
eu_chr_inat_all <- inner_join(eu_chr_y_sum_i, all_clim_res, by = c("year", "grid_id", "lon", "lat"))
anth_c_ng_all <- inner_join(anth_c_y_sum_ng, all_clim_res, by = c("year", "grid_id", "lon", "lat"))
eu_chr_ng_all <- inner_join(eu_chr_y_sum_ng, all_clim_res, by = c("year", "grid_id", "lon", "lat"))

#save
write.csv(anth_c_obs_all, "Datasets_analysis/anth_c_obs")
write.csv(eu_chr_obs_all, "Datasets_analysis/eu_chr_obs")
write.csv(anth_c_inat_all, "Data_analysis_inat/anth_c_inat")
write.csv(eu_chr_inat_all, "Data_analysis_inat/eu_chr_inat")
write.csv(anth_c_ng_all, "Datasets_analysis_ng/anth_c_inat")
write.csv(eu_chr_ng_all, "Datasets_analysis_ng/eu_chr_inat")



#6.4. Flight peak calculation for P.n. P.r., G.r. ------
get_peak_part <- function(species_peak) {
  species_peak %>% 
  group_by(year, grid_id, lon, lat) %>% 
    mutate(
      eventDate = as.integer(yday(eventDate))) %>%
    summarise(
    first_obs = min(eventDate),
    last_obs = max(eventDate),
    peak = as.integer(mean(c(last_obs, first_obs))), 
    n = n(),
    stateProvince = stateProvince,
    gen_time = as.integer(max(last_obs) - min(first_obs)),
    .groups = "drop"
    ) %>% 
    arrange(year, grid_id, lon, lat) %>%
    distinct()
}
  
p_n_1_obs <- get_peak_part(p_napi_1_o)
p_n_2_obs <- get_peak_part(p_napi_2_o)

p_n_1_inat <- get_peak_part(p_napi_1_i)
p_n_2_inat <- get_peak_part(p_napi_2_i)

p_r_1_obs <- get_peak_part(p_rapae_1_o)
p_r_2_obs <- get_peak_part(p_rapae_2_o)

p_r_1_inat <- get_peak_part(p_rapae_1_i)
p_r_2_inat <- get_peak_part(p_rapae_2_i)

g_r_1_obs <- get_peak_part(gone_rh_1_o)
g_r_2_obs <- get_peak_part(gone_rh_2_o)

g_r_1_inat <- get_peak_part(gone_rh_1_i)
g_r_2_inat <- get_peak_part(gone_rh_2_i)



get_per_year_stats <- function(dataset) { 
  summary <- dataset %>% 
    filter(year %in% c(2016:2025)) %>%
    group_by(year) %>% 
    summarise(
      peak = mean(peak), 
      gen_length = mean(gen_time), 
      .groups = "drop"
    )
  return(summary)
}

p_n_1_sum <- get_per_year_stats(p_n_1_obs)
p_n_2_sum <- get_per_year_stats(p_n_2_obs)
p_r_1_sum <- get_per_year_stats(p_r_1_obs)
p_r_2_sum <- get_per_year_stats(p_r_2_obs)
g_r_1_sum <- get_per_year_stats(g_r_1_obs)
g_r_2_sum <- get_per_year_stats(g_r_2_obs)


p_n_1_sum_i <- get_per_year_stats(p_n_1_inat)
p_n_2_sum_i <- get_per_year_stats(p_n_2_inat)
p_r_1_sum_i <- get_per_year_stats(p_r_1_inat)
p_r_2_sum_i <- get_per_year_stats(p_r_2_inat)
g_r_1_sum_i <- get_per_year_stats(g_r_1_inat)
g_r_2_sum_i <- get_per_year_stats(g_r_2_inat)


### Peak 
#P napi
plot(peak ~ year, 
     data = p_n_1_sum, pch = 19)
plot(peak ~ year, 
     data = p_n_2_sum, pch = 19)
#clear decrease of first peak with some plasticity eg 2023 being late


#P rapae
plot(peak ~ year, 
     data = p_r_1_sum, pch = 19)
plot(peak ~ year, 
     data = p_r_2_sum, pch = 19)
#decrease with plasticity eg 2023 being late

plot(peak ~ year, 
     data = g_r_1_sum, pch = 19)
plot(peak ~ year, 
     data = g_r_2_sum, pch = 19)
#unclear pattern but decrease of first peak, 2023 late


p_n_1_obs_all <- inner_join(p_n_1_obs, all_clim_res, by = c("year", "lon", "lat"))
p_n_2_obs_all <- inner_join(p_n_2_obs, all_clim_res, by = c("year", "lon", "lat"))
p_n_1_inat_all <- inner_join(p_n_1_inat, all_clim_res, by = c("year", "lon", "lat"))
p_n_2_inat_all <- inner_join(p_n_2_inat, all_clim_res, by = c("year", "lon", "lat"))

p_r_1_obs_all <- inner_join(p_r_1_obs, all_clim_res, by = c("year", "lon", "lat"))
p_r_2_obs_all <- inner_join(p_r_2_obs, all_clim_res, by = c("year", "lon", "lat"))
p_r_1_inat_all <- inner_join(p_r_1_inat, all_clim_res, by = c("year", "lon", "lat"))
p_r_2_inat_all <- inner_join(p_r_2_inat, all_clim_res, by = c("year", "lon", "lat"))


g_r_1_obs_all <- inner_join(g_r_1_obs, all_clim_res, by = c("year", "lon", "lat"))
g_r_2_obs_all <- inner_join(g_r_2_obs, all_clim_res, by = c("year", "lon", "lat"))
g_r_1_inat_all <- inner_join(g_r_1_inat, all_clim_res, by = c("year", "lon", "lat"))
g_r_2_inat_all <- inner_join(g_r_2_inat, all_clim_res, by = c("year", "lon", "lat"))

#save 
write.csv(p_n_1_obs_all, "Datasets_analysis/p_n_1_obs_all")
write.csv(p_n_2_obs_all, "Datasets_analysis/p_n_2_obs_all")
write.csv(p_r_1_obs_all, "Datasets_analysis/p_r_1_obs_all")
write.csv(p_r_2_obs_all, "Datasets_analysis/p_r_2_obs_all")
write.csv(g_r_1_obs_all, "Datasets_analysis/g_r_1_obs_all")
write.csv(g_r_2_obs_all, "Datasets_analysis/g_r_2_obs_all")

write.csv(p_n_1_inat_all, "Data_analysis/p_n_1_inat_all")
write.csv(p_n_2_inat_all, "Data_analysis/p_n_2_inat_all")
write.csv(p_r_1_inat_all, "Data_analysis/p_r_1_inat_all")
write.csv(p_r_2_inat_all, "Data_analysis/p_r_2_inat_all")
write.csv(g_r_1_inat_all, "Data_analysis/g_r_1_inat_all")
write.csv(g_r_2_inat_all, "Data_analysis/g_r_2_inat_all")





