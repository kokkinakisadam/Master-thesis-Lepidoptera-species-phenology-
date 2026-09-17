#Extraction and processing of daily climate data of Germany

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

# 1 Germany grid ------
germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")
germany <- germany[6]
crs(germany) #4326
#plot(germany, col = NA)
#germany <- st_transform(germany, 4326)
extent <- ext(germany)
# transform to UTM Zone 32N to have a equal area grid 
germany_utm <- st_transform(germany, 32632)
#germany_proj <- st_transform(germany, 4326)

# 5 km grid for analysis 
grid_5_utm <- st_make_grid(germany_utm, cellsize = 5000, square = TRUE)
# adding unique grid IDs
grid_5_utm <- st_sf(grid_id = 1:length(grid_5_utm), geometry = grid_5_utm)
grid_5_utm <- st_intersection(grid_5_utm, germany_utm)
# extracting grid centroids
grid_5_utm$centroid <- st_centroid(grid_5_utm$geometry)

grid_5_proj <- st_transform(grid_5_utm, 4326)

grid_5_proj$centroid <- st_centroid(grid_5_proj$geometry)
grid_5_proj$lon <- st_coordinates(grid_5_proj$centroid)[, 1]
grid_5_proj$lat <- st_coordinates(grid_5_proj$centroid)[, 2]

# 2. Tmax  --------
start_date <- as.Date("2015-01-01") #from last year
end_date <- as.Date("2024-12-31") #maximum in the available NetCDF files


nc_file_tmax <- "Data_clim/tx_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc"
# output folder
output_folder <- file.path("Data_clim/extracted_results_tmax")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

nc <- nc_open(nc_file_tmax)
print(names(nc$var)) 
nc_close(nc)
tmax_rast <- rast(nc_file_tmax, subds = "tx")
tmax_rast <- project(tmax_rast, crs(grid_5_proj), method = "bilinear" )

nl <- nlyr(tmax_rast) #5114
time_vals_tmax <- time(tmax_rast) 
# loop over each NetCDF file (assumed to be annual)
#for (file in nc_files) {

time_mask_tmax <- (time_vals_tmax >= start_date & time_vals_tmax <= end_date)
time_ind_tmax <- which(time_mask_tmax)

tmax_subset <- tmax_rast[[time_ind_tmax]]
dates_subset_tmax <- time_vals_tmax[time_ind_tmax] #error

grid_coords <- grid_5_proj %>% 
  st_drop_geometry() %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

nlyr(tmax_subset)  #3288 (days) - 3653

chunk_size <- 100 
n_cells <- nrow(grid_coords)
n_chunks <- ceiling(n_cells/chunk_size)

tmax_results <- list ()
#message("Processing file...")


for(chunk in 1:n_chunks) { 
  
  message("Processing chunk", chunk, "of", n_chunks)
  
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_cells)
  
  chunk_coords <- grid_coords[start_idx:end_idx, ]
  chunk_ids <- grid_5_proj$grid_id[start_idx:end_idx]
  
  # extract time values from NetCDF (if available)
  #time_vals <- if (!is.null(getZ(rbrick))) getZ(rbrick) else paste0("Day_", 1:n_layers)
  
  # extract raster values at coordinates
  extracted_vals_tmax <- terra::extract(tmax_rast, chunk_coords)
  extracted_df_tmax <- as.data.frame(extracted_vals_tmax)
  
  names(extracted_df_tmax)[1] <- "point_id"
  extracted_df_tmax$grid_id <- chunk_ids
  extracted_df_tmax$lon <- chunk_coords[,1]
  extracted_df_tmax$lat <- chunk_coords[,2]
  
  # Convert to long format
  extracted_long_tmax <- extracted_df_tmax %>%
    pivot_longer(cols = -c(point_id, grid_id, lon, lat),
                 names_to = "layer_name",
                 values_to = "tmax")
  
  extracted_long_tmax <- extracted_long_tmax %>% mutate(
    layer_num = as.numeric(str_extract(layer_name, "\\d+")),
    date = dates_subset_tmax[layer_num],
    year = year(date),
    month = month(date),
    day = day(date)
  ) %>%
    dplyr::select(-point_id, -layer_name, -layer_num)
  
  
  tmax_results[[chunk]] <- extracted_long_tmax
  
  # Clean up
  rm(extracted_vals_tmax, extracted_df_tmax, extracted_long_tmax)
  gc()
}

all_tmax <- bind_rows(tmax_results)
saveRDS(all_tmax, "Data_clim/results_tmax")

all_tmax <- readRDS("Data_clim/results_tmax") 

head(all_tmax)
str(all_tmax) 

# 3. Tmin -----

nc_file_tmin <- "Data_clim/tn_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc"
# output folder
output_folder <- file.path("Data_clim/extracted_results_tmax")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

nc <- nc_open(nc_file_tmin)
print(names(nc$var)) 
nc_close(nc)
tmin_rast <- rast(nc_file_tmin, subds = "tn")

tmin_rast <- project(tmin_rast, crs(grid_5_proj), method = "bilinear" )

nl_tmin <- nlyr(tmin_rast) #5114
time_vals_tmin <- time(tmin_rast) 
# loop over each NetCDF file (assumed to be annual)
#for (file in nc_files) {

time_mask_tmin <- (time_vals_tmin >= start_date & time_vals_tmin <= end_date)
time_ind_tmin <- which(time_mask_tmin)

tmin_subset <- tmin_rast[[time_ind_tmin]]
dates_subset_tmin <- time_vals_tmin[time_ind_tmin] #error

grid_coords <- grid_5_proj %>% 
  st_drop_geometry() %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

nlyr(tmin_subset)  #3288 (days)

chunk_size <- 100 
n_cells <- nrow(grid_coords)
n_chunks <- ceiling(n_cells/chunk_size)

tmin_results <- list ()
#message("Processing file...")


for(chunk in 1:n_chunks) { 
  
  message("Processing chunk", chunk, "of", n_chunks)
  
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_cells)
  
  chunk_coords <- grid_coords[start_idx:end_idx, ]
  chunk_ids <- grid_5_proj$grid_id[start_idx:end_idx]
  
  # extract time values from NetCDF (if available)
  #time_vals <- if (!is.null(getZ(rbrick))) getZ(rbrick) else paste0("Day_", 1:n_layers)
  
  # extract raster values at coordinates
  extracted_vals_tmin <- terra::extract(tmin_rast, chunk_coords)
  extracted_df_tmin <- as.data.frame(extracted_vals_tmin)
  
  names(extracted_df_tmin)[1] <- "point_id"
  extracted_df_tmin$grid_id <- chunk_ids
  extracted_df_tmin$lon <- chunk_coords[,1]
  extracted_df_tmin$lat <- chunk_coords[,2]
  
  # Convert to long format
  extracted_long_tmin <- extracted_df_tmin %>%
    pivot_longer(cols = -c(point_id, grid_id, lon, lat),
                 names_to = "layer_name",
                 values_to = "tmin")
  
  extracted_long_tmin <- extracted_long_tmin %>% mutate(
    layer_num = as.numeric(str_extract(layer_name, "\\d+")),
    date = dates_subset_tmin[layer_num],
    year = year(date),
    month = month(date),
    day = day(date)
  ) %>%
    dplyr::select(-point_id, -layer_name, -layer_num)
  
  
  tmin_results[[chunk]] <- extracted_long_tmin
  
  # Clean up
  rm(extracted_vals_tmin, extracted_df_tmin, extracted_long_tmin)
  gc()
}

all_tmin <- bind_rows(tmin_results)
saveRDS(all_tmin, "Data_clim/results_tmin")

all_tmin <- readRDS("Data_clim/results_tmin") 

head(all_tmin)
str(all_tmin) 

# 4. Tmean ----------

nc_file_tmean <- "Data_clim/tg_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc"
# output folder
output_folder <- file.path("Data_clim/extracted_results_tg")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

nc <- nc_open(nc_file_tmean) #tg ok
print(names(nc$var)) 
nc_close(nc)
tmean_rast <- rast(nc_file_tmean, subds = "tg")

tmean_rast <- project(tmean_rast, crs(grid_5_proj), method = "bilinear" )

nl <- nlyr(tmean_rast) #5114
time_vals_tmean <- time(tmean_rast) 
# loop over each NetCDF file (assumed to be annual)
#for (file in nc_files) {

time_mask_tmean <- (time_vals_tmean >= start_date & time_vals_tmean <= end_date)
time_ind_tmean <- which(time_mask_tmean)

tmean_subset <- tmean_rast[[time_ind_tmean]]
dates_subset_tmean <- time_vals_tmean[time_ind_tmean] #error

grid_coords <- grid_5_proj %>% 
  st_drop_geometry() %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

nlyr(tmean_subset)  #3288 (days)

chunk_size <- 100 
n_cells <- nrow(grid_coords)
n_chunks <- ceiling(n_cells/chunk_size)

tmean_results <- list ()
#message("Processing file...")


for(chunk in 1:n_chunks) { 
  
  message("Processing chunk", chunk, "of", n_chunks)
  
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_cells)
  
  chunk_coords <- grid_coords[start_idx:end_idx, ]
  chunk_ids <- grid_5_proj$grid_id[start_idx:end_idx]
  
  # extract time values from NetCDF (if available)
  #time_vals <- if (!is.null(getZ(rbrick))) getZ(rbrick) else paste0("Day_", 1:n_layers)
  
  # extract raster values at coordinates
  extracted_vals_tmean <- terra::extract(tmean_rast, chunk_coords)
  extracted_df_tmean <- as.data.frame(extracted_vals_tmean)
  
  names(extracted_df_tmean)[1] <- "point_id"
  extracted_df_tmean$grid_id <- chunk_ids
  extracted_df_tmean$lon <- chunk_coords[,1]
  extracted_df_tmean$lat <- chunk_coords[,2]
  
  # Convert to long format
  extracted_long_tmean <- extracted_df_tmean %>%
    pivot_longer(cols = -c(point_id, grid_id, lon, lat),
                 names_to = "layer_name",
                 values_to = "tmean")
  
  extracted_long_tmean <- extracted_long_tmean %>% mutate(
    layer_num = as.numeric(str_extract(layer_name, "\\d+")),
    date = dates_subset_tmean[layer_num],
    year = year(date),
    month = month(date),
    day = day(date)
  ) %>%
    dplyr::select(-point_id, -layer_name, -layer_num)
  
  
  tmean_results[[chunk]] <- extracted_long_tmean
  
  # Clean up
  rm(extracted_vals_tmean, extracted_df_tmean, extracted_long_tmean)
  gc()
}

all_tmean <- bind_rows(tmean_results)
saveRDS(all_tmean, "Data_clim/results_tmean")

all_tmean <- readRDS("Data_clim/results_tmean") 

head(all_tmean)
str(all_tmean) 

# 5.  Precipitation ----------
start_date <- as.Date("2015-01-01") #from our analysis
end_date <- as.Date("2024-12-31") #maximum in the available NetCDF files


nc_file_prec <- "Data_clim/rr_ens_mean_0.1deg_reg_2011-2024_v31.0e.nc"
# output folder
output_folder <- file.path("Data_clim/extracted_results_prec")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

nc <- nc_open(nc_file_prec)
print(names(nc$var)) #rr
nc_close(nc)

prec_rast <- rast(nc_file_prec, subds = "rr")
prec_rast <- project(prec_rast, crs(grid_5_proj), method = "bilinear" )


nl <- nlyr(prec_rast) #5114
time_vals <- time(prec_rast) 
# loop over each NetCDF file (assumed to be annual)
#for (file in nc_files) {

time_mask <- (time_vals >= start_date & time_vals <= end_date)
time_ind <- which(time_mask)

prec_subset <- prec_rast[[time_ind]]
dates_subset <- time_vals[time_ind] #error

grid_coords <- grid_5_proj %>% 
  st_drop_geometry() %>%
  dplyr::select(lon, lat) %>%
  as.matrix()

nlyr(prec_subset)  #3288 (days)

chunk_size <- 100 
n_cells <- nrow(grid_coords)
n_chunks <- ceiling(n_cells/chunk_size)

prec_results <- list ()

for(chunk in 1:n_chunks) { 
  
  message("Processing chunk", chunk, "of", n_chunks)
  
  start_idx <- (chunk - 1) * chunk_size + 1
  end_idx <- min(chunk * chunk_size, n_cells)
  
  chunk_coords <- grid_coords[start_idx:end_idx, ]
  chunk_ids <- grid_5_proj$grid_id[start_idx:end_idx]
  
  # extract time values from NetCDF (if available)
  #time_vals <- if (!is.null(getZ(rbrick))) getZ(rbrick) else paste0("Day_", 1:n_layers)
  
  # extract raster values at coordinates
  extracted_vals <- terra::extract(prec_rast, chunk_coords)
  extracted_df <- as.data.frame(extracted_vals)
  
  names(extracted_df)[1] <- "point_id"
  extracted_df$grid_id <- chunk_ids
  extracted_df$lon <- chunk_coords[,1]
  extracted_df$lat <- chunk_coords[,2]
  
  # Convert to long format
  extracted_long <- extracted_df %>%
    pivot_longer(cols = -c(point_id, grid_id, lon, lat),
                 names_to = "layer_name",
                 values_to = "prec")
  
  extracted_long <- extracted_long %>% mutate(
    layer_num = as.numeric(str_extract(layer_name, "\\d+")),
    date = dates_subset[layer_num],
    year = year(date),
    month = month(date),
    day = day(date)
  ) %>%
    dplyr::select(-point_id, -layer_name, -layer_num)
  
  
  prec_results[[chunk]] <- extracted_long
  
  # Clean up
  rm(extracted_vals, extracted_df, extracted_long)
  gc()
}

all_prec <- bind_rows(prec_results)
saveRDS(all_prec, "Data_clim/results_prec")

all_prec <- readRDS("Data_clim/results_prec") 

head(all_prec)
str(all_prec) 


# 6. Analysis --------

all_tmax <- readRDS("Data_clim/results_tmax") 
all_tmin <- readRDS("Data_clim/results_tmin") 
all_tmean <- readRDS("Data_clim/results_tmean") 
all_prec <- readRDS("Data_clim/results_prec") 



