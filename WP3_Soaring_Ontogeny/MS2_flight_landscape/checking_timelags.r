#code for looking up the data resolution of golden eagle data
#Jul. 24. 2023. Elham Nourani, PhD
#Canberra, AU

library(move2)
library(tidyverse)
library(lubridate)
library(mapview)
library(units)

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
GE_id <- 	282734839

info <- movebank_retrieve(entity_type = "individual", study_id =282734839)

movebank_download_study_info(282734839)

#downlaod a bit of data to check data collection frequency
gps <- movebank_retrieve(study_id = 	282734839, individual_local_identifier = c("Siat20 (eobs 7037)", "Mals2_20 (eobs 7579) "),
                         entity_type = "event", sensor_type_id= "gps", attributes = "all", 
                         timestamp_end = "20200915000000000",
                         timestamp_end = "20200925000000000") %>% 
  drop_na(location_long)

mv <- mt_as_move2(gps, coords = c("location_long", "location_lat"), time_column = "timestamp", track_id_column = "individual_local_identifier", track_attributes = "tag_id")

time_lags <- mt_time_lags(mv) #the median is 19.9 minutes!
