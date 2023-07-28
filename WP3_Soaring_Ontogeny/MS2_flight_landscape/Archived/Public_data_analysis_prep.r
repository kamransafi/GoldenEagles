#Script for step-selection analysis of golden eagle commuting flights as reported in Nourani et al. 2023
#Elham Nourani, PhD. 25.07.2023
#enourani@ab.mpg.de

library(tidyverse)
library(lubridate)

##### STEP 1: open tracking data#####

#open file with all data (March 15. 2023)
data <- readRDS("/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack_Golden_Eagle_Alps_1min.rds")

summary(data$location.long)


#20 min subset to get rid of super bursts and keep regular sampling
data <- data %>% 
  drop_na(location.long) %>% 
  mutate(dt_20min = round_date(timestamp, "20 minutes")) %>% 
  group_by(individual.local.identifier,dt_20min) %>% 
  slice(1) %>% 
  ungroup()

#only keep the 55 individuals
emig_dates <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/dispersal_dates.csv")
data <- data  %>% 
filter(individual.local.identifier %in% emig_dates$individual.local.identifier) 


#intermediate save

saveRDS(data, file = "/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack_Golden_Eagle_Alps_55_20min.rds")


data <- data %>% 
  left_join(emig_dates %>% dplyr::select(individual.local.identifier, dispersal_dt), by = "individual.local.identifier") %>% 
  mutate(weeks_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("weeks")) %>% as.numeric() %>% ceiling(),
         days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")) %>% as.numeric())

#keep the post dispersal only
data_pd <- data %>% 
  filter(between(weeks_since_emig, 1, 156)) %>% 
  select(timestamp, location.long, location.lat, eobs.horizontal.accuracy.estimate, height.above.ellipsoid, sensor.type, individual.taxon.canonical.name,
         tag.local.identifier, individual.local.identifier, dispersal_dt, days_since_emig, weeks_since_emig)

write.csv(data_pd, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/public_data/GPS_data.csv", row.names = F)

#the tracking data contains commuting flights for 55 juvenile golden eagles in the alps

emig_dates <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/fleding_emigration_timing_Mar2023.rds")
emig_dates[emig_dates$individual.local.identifier == "Matsch19 (eobs 7035)", "emigration_dt"] <- "2020-03-08 12:17:00" 
emig_dates[emig_dates$individual.local.identifier == "Güstizia18 (eobs 5942)", "emigration_dt"] <- "2019-03-24 11:18:31"

#prepare to save in public data folder
emig_dates <- emig_dates %>% 
  filter(individual.local.identifier %in% unique(data$individual.local.identifier)) %>% 
  select(individual.local.identifier, emigration_dt) %>% 
  mutate(recurse_radius_km = ifelse(individual.local.identifier %in% c("Matsch19 (eobs 7035)", "Güstizia18 (eobs 5942)"), 30, 7)) %>% 
  rename(dispersal_dt = emigration_dt) %>% 
  arrange(individual.local.identifier)

write.csv(emig_dates, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/dispersal_dates.csv")
#public data folder: "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/"

emig_dates <- read.csv("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/dispersal_dates.csv")

#open tracking data (after random step generation. 02_data_processing_&_annotation.R)
used_av_track <- readRDS("alt_50_60_min_55_ind2.rds") %>% 
  select(used, individual.local.identifier, timestamp, location.long, location.lat, height.above.ellipsoid, geoid, dem, height_msl, height_ground, 
         burst_id, step_id, stratum, turning_angle, step_length) %>% 
  data.frame() #this removes all attributes

write.csv(used_av_track, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/ssf_input_data.csv")



##### TRY2: Jul28#####

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#I redid all the steps for data prep with move2 and tidyverse, but the output file was different from the original one and model did not converge.
#Use the dataset used for the paper: 1) correct the week assignment, correct the dispersal date for two individuals

#open data. data was prepared in 03_energy_landscape_modeling_method1.R. 
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  mutate(emigration_dt = as.POSIXct(emigration_dt, tz = "UTC"),
         animal_ID = as.numeric(as.factor(ind1)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum))) %>% 
  select(c("timestamp" , "height.above.ellipsoid" , "individual.local.identifier",
           "emigration_dt" , "days_since_emig" ,"weeks_since_emig" , "dem" , "geoid",                      
           "height_msl", "height_ground" , "embc_clst", "embc_clst_smth" ,    "burst_id",                 
           "step_length" ,  "turning_angle" , "step_id" ,"used"  , "heading" ,                   
           "stratum" , "TRI_100", "ridge_100", "location.long" ,"location.lat",             
           "animal_ID" , "stratum_ID" ))

#remove extra data for the two individuals:
data_f <- data %>% 
  filter(case_when(individual.local.identifier == "Matsch19 (eobs 7035)" ~ timestamp >= as.POSIXct("2020-03-08 12:17:00", tz = "UTC"),
                   individual.local.identifier == "Güstizia18 (eobs 5942)" ~ timestamp >= as.POSIXct("2019-03-24 11:18:31", tz = "UTC"),
                   TRUE ~ timestamp >= as.POSIXct("2010-01-01 00:00:00", tz = "UTC"))) #case_when needs an else statement, which is specified by true. keep everything else

data_f[data_f$individual.local.identifier == "Matsch19 (eobs 7035)", "emigration_dt"] <- "2020-03-08 12:17:00" 
data_f[data_f$individual.local.identifier == "Güstizia18 (eobs 5942)", "emigration_dt"] <- "2019-03-24 11:18:31"

## fix the date columns
data <- data_f %>% 
  rename(dispersal_dt = emigration_dt) %>% 
  mutate(weeks_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("weeks")) %>% as.numeric() %>% ceiling() + 1,
         days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")) %>% as.numeric() + 1) %>% 
  filter(weeks_since_emig < 157)

write.csv(data, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/public_data/SSF_input_data.csv", row.names = F)
