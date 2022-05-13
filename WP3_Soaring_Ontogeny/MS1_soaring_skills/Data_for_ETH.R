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
