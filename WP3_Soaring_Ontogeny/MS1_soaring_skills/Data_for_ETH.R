#send a couple of days of flight data to the ETH people, to annotate with COSMO data and see how promising it is.
#May 12. 2022. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(sf)
library(lubridate)
library(mapview)
library(rgdal)
library(move)


setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
source("/home/enourani/ownCloud/Work/Projects/wind_support_Kami.R")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# STEP 1: prepare sample data to send to eth ---------------------------------------
#open data
load("all_data_w_emig.RData") #data_w_info; from embc_segmentation.R

#subset two days
two_days <- data_w_info %>% 
  filter(individual.local.identifier == "Sinestra1 19 (eobs 7003)" &
           date(timestamp) %in% c(as.Date("2020-05-27"),as.Date("2021-04-17"))) %>% 
  mutate(day = day(timestamp)) %>% 
  drop_na("location.lat") %>% 
  dplyr::select(c("timestamp", "location.long", "location.lat", "height.above.ellipsoid", "id", "utm.easting", "utm.northing", "utm.zone"))

ind_sf <- two_days %>%  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)

mapview(ind_sf, zcol = date(two_days$timestamp))

write.csv(two_days, "sample_golden_eagles.csv")

data <- read.csv("sample_golden_eagles.csv") %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)

# STEP 2: investigate data from eth ---------------------------------------
#coordinates are rounded up, so the eth points end up overlapping. Match the order of points by time and append to original file
#if only col-binding based on the order of rows, no need to convert the time column to timestamp ;)

time_ref_row <- read.table("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/COSMO_wind/sample_golden_eagles-2020-05-27.uvw", nrow = 1)

time_ref <- as.POSIXct(strptime(paste(paste(str_sub(time_ref_row[,3], 1,4), str_sub(time_ref_row[,3], 5,6), str_sub(time_ref_row[,3], 7,8), sep = "-"), #the date
                  paste(str_sub(time_ref_row[,3], 10,11), str_sub(time_ref_row[,3], 12,13), "00", sep = ":"), sep = " "), #the time
                  format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")

wind <- read.table("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/COSMO_wind/sample_golden_eagles-2020-05-27.uvw", skip = 5)
colnames(wind) <- c("time", "lon", "lat", "z", "u", "v","w")

#wind <- wind %>% 
#  mutate(hr_min = paste(str_split(time, "\\.", simplify = T)[, 1] %>% str_pad(2, "left", "0"), #hour
#                        str_split(time, "\\.", simplify = T)[, 2] %>% str_pad(2, "right", "0"), sep = ":") %>%  hm()) %>% #minute
#  mutate(timestamp = time_ref + hr_min) 
  
wind_data <- read.csv("sample_golden_eagles.csv") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter(date(timestamp) == date(time_ref)) %>% #extract the day of interest :p
  bind_cols(wind %>% dplyr::select(5:7))



wind_sf <- wind_data %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)

mapview(wind_sf, zcol = "w") + mapview(data, color = "gray")

# STEP 3: estimate metrics using the data ---------------------------------------

#estimate wind speed, wind support and crosswind

#convert to move object to estimate heading
duplicated_times_to_remove <- unlist(lapply(getDuplicatedTimestamps(x = wind_data$id, timestamps = wind_data$timestamp), "[", 2))

wind_data <- wind_data[-duplicated_times_to_remove,]

mv <- move(x = wind_data$location.long, y = wind_data$location.lat, time = wind_data$timestamp, proj = wgs, animal = "id", data = wind_data)

mv$heading <- c(angle(mv), NA)
mv$gr_speed <- c(speed(mv), NA)   

wind_data <- mv %>%
  as("sf") %>% 
  mutate(ws = wind_support(u = u, v = v, heading = heading),
         cw = cross_wind(u = u, v = v, heading = heading),
         wind_speed = sqrt(u^2 + v^2)) %>% 
  mutate(airspeed = sqrt((gr_speed - ws)^2 + (cw)^2))

# STEP 4: compare with ERA 5 --------------------------------------

data <- wind_data %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
  dplyr::select(-1) %>% 
  as.data.frame()

#rename columns
colnames(data)[c(2,3)] <- c("location-long","location-lat")

write.csv(data, "ETH_for_ERA5_annotation.csv", row.names = FALSE)

era5 <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ETH_ERA5_annotation/ETH_forERA5_noheight.csv-6585178841830394962/ETH_forERA5_noheight.csv-6585178841830394962.csv") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  rename(era5_u_900 = ECMWF.ERA5.PL.U.Wind,
         era5_v_900 = ECMWF.ERA5.PL.V.Wind) %>% 
  full_join(wind_data)

#they are pretty different. but also I don't know if the 900 mbar is a good level to use
Us <- melt(era5[,c("u","era5_u_900")])
ggplot(Us, aes(x = value, fill = variable)) + geom_density(alpha = 0.25) + theme_bw()

Vs <-  melt(era5[,c("v","era5_v_900")])
ggplot(Vs, aes(x = value, fill = variable)) + geom_density(alpha = 0.25) + theme_bw()
