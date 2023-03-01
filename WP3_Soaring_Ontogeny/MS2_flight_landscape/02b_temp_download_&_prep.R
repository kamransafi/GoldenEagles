#script for downloading monthly temperature data to 1) append to step selection input data and 2) make predictions for the interaction plot and the Alpine region maps
#Jan. 20. 2022, Konstanz, DE
#Elham Nourani PhD

#follows from data_processing_&_annotation.R

library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(ncdf4)
library(reticulate)

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

# STEP 1: Prep specs for data download -----------------------------------------------------------------------------------------------

#open data. this has been already annotated with terrain features
load("alt_50_20_min_70_ind_static_time_ann.RData") #cmpl_ann

#open Alpine region to extract extent
Alps_extent <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(wgs)) %>% 
  st_bbox()

area <- round(as.numeric(Alps_extent[c(4,1,2,3)]),2)

#extract unique years
years <- unique(year(cmpl_ann$timestamp))

#unique months
mnths <- c(1:12) %>% str_pad(2,"left","0")

####################################### DOWNLOAD MONONTHLY DATA ###############################################


# STEP 2: Download data from CDS -----------------------------------------------------------------------------------------------

#prepare the connection
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

output_path <- "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/monthly_temp_cds/"

#try all months and years at once wont work

for(i in years){
  for(j in mnths){
    
    request <- r_to_py(list(
      product_type = "monthly_averaged_reanalysis",
      variable = "2m_temperature", 
      year = i,
      month = j,
      time = "00:00",
      area = area,  # North, West, South, East.
      format = "netcdf",
      dataset_short_name = "reanalysis-era5-land-monthly-means"
    ))
    
    server$retrieve("reanalysis-era5-land-monthly-means",
                    request,
                    target = paste0(output_path,"temp_data_", j, "_", i, ".nc")) 
    
  }
  
}

# STEP 3: process netcdf data: modeling -----------------------------------------------------------------------------------------------

ncdf_list <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/monthly_temp_cds", full.names = T)

#extract temperature values for each point of the ssf input dataset. do it by month-year 
mnth_yr_ls <- split(cmpl_ann, list(year(cmpl_ann$timestamp), month(cmpl_ann$timestamp)))

#remove empty elements
mnth_yr_ls <- mnth_yr_ls[lapply(mnth_yr_ls, nrow) > 0]


data_temp <- lapply(mnth_yr_ls, function(x){
  
  mnth_yr <- paste(unique(month(x$timestamp)) %>% str_pad(2,"left","0"),  unique(year(x$timestamp)) %>% as.character(), sep = "_")
  
  #open corresponding temp data as a raster
  temp_r <- ncdf_list[ncdf_list %>% 
    str_which(mnth_yr)] %>% 
    rast() %>% 
    as("SpatRaster")
  
  x_temp <- x %>% 
    st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
    mutate(month_temp = extract(temp_r,.)$t2m,
           location.long = st_coordinates(.)[,1],
           location.lat = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()

  x_temp
  
}) %>% 
  reduce(rbind)

saveRDS(data_temp, "alt_50_20_min_70_ind_static_time_ann_temp.rds") #data_temp


# STEP 4: download data: landscape maps -----------------------------------------------------------------------------------------------

#have a look at the distribution of data across months
cmpl_ann %>% filter(used == 1) %>% mutate(month = month(timestamp)) %>%  pull(month) %>%  hist()

#let's do June vs Oct for 10 years
years <- as.character(2013:2022)
mnths <- c("06", "10")

#area is the same as above
ext <- round(as.numeric(Alps_extent[c(4,1,2,3)]),2)

#prepare the connection
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

output_path <- "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/monthly_temp_cds_10yrs/"


for(i in years){
  for(j in mnths){
    
    request <- r_to_py(list(
      product_type = "monthly_averaged_reanalysis",
      variable = "2m_temperature", 
      year = i,
      month = j,
      time = "00:00",
      area = ext,  # North, West, South, East.
      format = "netcdf",
      dataset_short_name = "reanalysis-era5-land-monthly-means"
    ))
    
    server$retrieve("reanalysis-era5-land-monthly-means",
                    request,
                    target = paste0(output_path,"temp_data_", j, "_", i, ".nc")) 
    
  }
  
}

# STEP 5: process data: landscape maps -----------------------------------------------------------------------------------------------

ncdf_list <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/monthly_temp_cds_10yrs", full.names = T)

# I need one value per cell. So, open all 10 layer for each month, take the average per cell, then append to the alpine pred dataframe

#open corresponding temp data as a raster
Jun_r <- ncdf_list[ncdf_list %>% 
                      str_which("06")] %>% 
  rast() %>% 
  mean()

writeRaster(Jun_r, file = "avg_temp_10yr_Jun.tif")

Oct_r <- ncdf_list[ncdf_list %>% 
                     str_which("10")] %>% 
  rast() %>% 
  mean()

writeRaster(Oct_r, file = "avg_temp_10yr_Oct.tif")



####################################### DOWNLOAD DAILY DATA ###############################################

#STEP 2: download

#unique months
days <- c(1:31) %>% str_pad(2,"left","0")

#prepare the connection
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

output_path <- "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/daily_temp_cds/"

#try all months and years at once wont work

for(i in years){
  for(j in mnths){
    
    request <- r_to_py(list(
      product_type = "reanalysis",
      variable ="2m_temperature",
      year = i,
      month = j,
      time = "12:00",
      day = days,
      area = area,  # North, West, South, East.
      format = "netcdf",
      dataset_short_name = "reanalysis-era5-single-levels"
    ))
    
    server$retrieve("reanalysis-era5-single-levels",
                    request,
                    target = paste0(output_path,"temp_data_", j, "_", i, ".nc")) 
    
  }
  
}

#STEP 3: process

ncdf_list <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/daily_temp_cds", full.names = T)

#extract temperature values for each point of the ssf input dataset. do it by day-month-year 
mnth_yr_ls <- split(cmpl_ann, list(year(cmpl_ann$timestamp), month(cmpl_ann$timestamp)))

#remove empty elements
mnth_yr_ls <- mnth_yr_ls[lapply(mnth_yr_ls, nrow) > 0]


data_temp <- lapply(mnth_yr_ls, function(x){
  
  mnth_yr <- paste(unique(month(x$timestamp)) %>% str_pad(2,"left","0"),  unique(year(x$timestamp)) %>% as.character(), sep = "_")
  
  #open corresponding temp data as a raster
  temp_r <- ncdf_list[ncdf_list %>% 
                        str_which(mnth_yr)] %>% 
    rast() %>% 
    as("SpatRaster")
  
  x_temp <- x %>% 
    st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
   extract(x = temp_r, y = ., layer = paste0("t2m_", day(.$timestamp)), bind = T) %>% 
    as.data.frame(geom = "XY") %>% 
    rename(location.long = x,
           location.lat = y,
           t2m = value)
  
  x_temp
  
}) %>% 
  reduce(rbind)

saveRDS(data_temp, "alt_50_20_min_70_ind_static_time_ann_dailytemp.rds") #data_temp