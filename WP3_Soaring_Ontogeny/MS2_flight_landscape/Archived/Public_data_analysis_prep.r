#Script for step-selection analysis of golden eagle commuting flights as reported in Nourani et al. 2023
#Elham Nourani, PhD. 25.07.2023
#enourani@ab.mpg.de




##### STEP 1: open tracking data#####

#open file with all data (March 15. 2023)
data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S",tz = "UTC"))




#20 min subset to get rid of super bursts and keep regular sampling
data <- data %>% 
  mutate(dt_20min = round_date(timestamp, "20 minutes")) %>% 
  group_by(individual.local.identifier,dt_20min) %>% 
  slice(1) 

#only keep the 55 individuals
emig_dates <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/dispersal_dates.csv")
data <- data  %>% 
filter(individual.local.identifier %in% emig_dates$individual.local.identifier) 


#intermediate save

saveRDS(data, file = "/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack_Golden_Eagle_Alps_55_20min.rds")


data <- data %>% 
  left_join(dispersal_dt %>% dplyr::select(individual.local.identifier, dispersal_dt), by = "individual.local.identifier") %>% 
  mutate(weeks_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("weeks")) %>% as.numeric(),
         days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")) %>% as.numeric())

#add week since emigration
dd <- dd %>% 
  left_join(dispersal_dt %>% dplyr::select(individual.local.identifier, dispersal_dt), by = "individual.local.identifier") %>% 
  mutate(weeks_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("weeks")) %>% as.numeric(),
         days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")) %>% as.numeric())

#keep the post dispersal only
data_pd <- data %>% 
  filter(weeks_since_emig > -1)

saveRDS(data_pd, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/data_wks_20min.rds" )


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

