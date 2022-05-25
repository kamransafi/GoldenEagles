#send a couple of days of flight data to the ETH people, to annotate with COSMO data and see how promising it is.
#May 12. 2022. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(sf)
library(lubridate)
library(mapview)
library(rgdal)

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

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

#### open data from ETH and investigate
#coordinates are rounded up, so the eth points end up overlapping. Match the order of points by time and append to original file
#if only colbinding based on the order of rows, no need to convert the time column to timestamp ;)

time_ref_row <- read.table("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/COSMO_wind/sample_golden_eagles-2020-05-27.uvw", nrow = 1)

time_ref <- as.POSIXct(strptime(paste(paste(str_sub(time_ref_row[,3], 1,4), str_sub(time_ref_row[,3], 5,6), str_sub(time_ref_row[,3], 7,8), sep = "-"), #the date
                  paste(str_sub(time_ref_row[,3], 10,11), str_sub(time_ref_row[,3], 12,13), "00", sep = ":"), sep = " "), #the time
                  format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")

wind <- read.table("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/COSMO_wind/sample_golden_eagles-2020-05-27.uvw", skip = 5)
colnames(wind) <- c("time", "lon", "lat", "z", "u", "v","w")

wind <- wind %>% 
  mutate(hr_min = paste(str_split(time, "\\.", simplify = T)[, 1] %>% str_pad(2, "left", "0"), #hour
                        str_split(time, "\\.", simplify = T)[, 2] %>% str_pad(2, "right", "0"), sep = ":") %>%  hm()) %>% #minute
  mutate(timestamp = time_ref + hr_min) 
  
wind_data <- read.csv("sample_golden_eagles.csv") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter(date(timestamp) == date(time_ref)) %>% #extract the day of interest :p
  bind_cols(wind %>% dplyr::select(5:7))



wind_sf <- wind_data %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)

mapview(wind_sf, zcol = "w") + mapview(data, color = "gray")
