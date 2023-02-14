#script for estimating trends in flight behavior of golden eagles to quantify improvement in soaring flight.
#Elham Nourani, PhD. Konstanz, DE.
#Feb. 7. 2023

library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(mapview)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

#data prepared by Hester
flight <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/post-dispersal_classified_GPS_34IDs_230207.rds") #one element per ind

#emigration date info
load("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/em_fl_dt_recurse_70ind.RData") #emig_fledg_dates

# STEP 1: calculate ratios for each burst, (the previous version did this per day) -------------------------------------------------------------

flight_summary <- lapply(flight, function (x){
  
  #calculate week and day since emigration
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), funs(ceiling(as.numeric(.))))
  
  #calculate the ratio of each mode of flight to total flight time (as number of rows within a burst)
  ratios <- x %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID, thermalClust) %>% 
    summarize(n_rows = length(thermalClust)) %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    complete(thermalClust, fill = list(n_rows = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    mutate(total_flight_in_burst = sum(n_rows)) %>% 
    mutate(flight_type_to_all = n_rows/total_flight_in_burst,
           circling_to_linear = n_rows[thermalClust == "circular"]/ n_rows[thermalClust == "linear"]) #calculate the ratio of thermal soaring to linear soaring. 19 bursts have neither ciclring nor linear soaring 
 
  ratios
  
}) %>% 
  reduce(rbind)
  

saveRDS(flight_summary, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_ratios_n34.rds")


# STEP 2: exploration -------------------------------------------------------------

#all flight modes
ggplot(data = flight_summary %>%  filter(weeks_since_emig <= 135), aes(x = weeks_since_emig, y = flight_type_to_all)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~ thermalClust) +
  theme_classic()

#circling vs linear soaring
ggplot(data = flight_summary %>% filter(weeks_since_emig <= 135 & !is.na(circling_to_linear) & !is.infinite(circling_to_linear)) %>% group_by(local_identifier, days_since_emig,burstID) %>%  slice(1), 
       aes(x = weeks_since_emig, y = circling_to_linear)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_classic()

