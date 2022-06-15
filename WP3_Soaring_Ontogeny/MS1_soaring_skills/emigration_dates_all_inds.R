## 06.13.2022 by Hester Bronnvik
## golden eagle emigration date estimation using the recurse package
## Elham's update to include all individuals: June 14. 2022

library(tidyverse)
library(move)
library(lubridate)
library(sf)
library(stringr)
require(recurse)
require(scales)
require(sp)

setwd("/home/enourani/Desktop/Hester_GE/recurse")


data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  #mutate(id = strsplit(individual.local.identifier, " ") %>% map_chr(., 1),
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))

#fledging dates from Hester
fledging_dates <- read.csv("dates.csv", header = F, stringsAsFactors = F, encoding = "latin1") %>% 
  rename(local.identifier = V1,
         fledging_dt = V2,
         FPT_radius = V3) %>% 
  drop_na(fledging_dt) %>% #remove NA values: 2 individuals
  mutate(fledging_dt = as.POSIXct(fledging_dt, tz = "UTC"))

#only keep individuals that we have fledging data for
data_70 <- data %>% 
  filter(individual.local.identifier %in% fledging_dates$local.identifier)

(b <- Sys.time())
emig_dates <- lapply(split(data_70, data_70$individual.local.identifier), function(x){
  
  #extract fledging dt
  fledging_dt <- fledging_dates %>% 
    filter(local.identifier == unique(x$individual.local.identifier))
  
  #keep only the first 300 days post fledging AND subsample to 15 min
  post_fledging_15 <- x %>% 
    filter(between(timestamp, fledging_dt$fledging_dt, fledging_dt$fledging_dt + days(300))) %>% 
    mutate(dt_15min = round_date(timestamp, "15 minutes")) %>% 
    group_by(dt_15min) %>% 
    slice(1) %>% 
    drop_na(location.long) %>% 
    as.data.frame()
  
  #create a move object
  ind_15_mv <- move(x = post_fledging_15$location.long, y = post_fledging_15$location.lat, time = post_fledging_15$timestamp, 
                    proj = "+proj=longlat +datum=WGS84 +no_defs", data = post_fledging_15)
  
  # project to UTM so that the units are in meters and the radius is easy to set
  # this UTM only works for Swiss birds
  ind_15_mv <- spTransform(ind_15_mv, CRSobj = "+init=EPSG:21781")#"+proj=utm +zone=32+datum=WGS84")
  
  
  #tryCatch({
    ind_visit <- getRecursions(ind_15_mv, radius =  7000, timeunits = "days", threshold = 14)
    #return the emigration time
    data.frame(individual.local.identifier = unique(x$individual.local.identifier),
                       emigration_dt = ind_visit$revisitStats$exitTime[1])

  
}) %>% 
  reduce(rbind)

Sys.time() - b #8 min

emig_dates



#bind fledging and emigration dates
emig_fledg_dates <- emig_dates %>% 
  full_join(fledging_dates, by = c("individual.local.identifier" = "local.identifier"))
  

save(emig_fledg_dates, file = "em_fl_dt_recurse_70ind.RData")
  
  