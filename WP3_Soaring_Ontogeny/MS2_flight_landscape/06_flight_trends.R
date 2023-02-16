#script for estimating trends in flight behavior of golden eagles to quantify improvement in soaring flight.
#Elham Nourani, PhD. Konstanz, DE.
#Feb. 7. 2023
#to be submitted to the cluster Raven


library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(mapview)
library(parallel)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

#data prepared by Hester
flight <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/Golden_eagle_wind_behav_full_Nov_2022.rds") #one element per ind. 1.4 GB

#emigration date info
load("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/em_fl_dt_recurse_70ind.RData") #emig_fledg_dates

# STEP 1: Summarise flight metrics per burst -------------------------------------------------------------

#calculate total n-rows of each behavior category for each burst 

flight_summary <- lapply(flight, function (x){
  
  #calculate week and day since emigration and fledging
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), ~ ceiling(as.numeric(.)))
  
  #summary for behavior (flapping and not flapping). n of rows per category. there is probably a more elegant to do the following in one chunk of code. but I couln't figure out how to make count run on multiple columns 
  flapping <- x %>% 
    mutate(behavior = factor(behavior, levels = c("Flapping", "NotFlapping", "Unclassified"))) %>% #convert to a factor with three levels
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    count(behavior) %>% 
    complete(behavior, fill = list(n = 0)) %>%  #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = behavior, values_from = n, names_prefix = "n_") %>%  #widen the dataframe. have one column for each category of behavior. nrow is not the n of unique bursts
    arrange(local_identifier, weeks_since_emig, days_since_emig, burstID) #just an insurance policy. the rows are already ordered by these.
  
  thermClust <-  x %>% 
    mutate(thermalClust = factor(thermalClust, levels = c("circular", "linear", "other"))) %>% #convert to a factor with three levels
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    count(thermalClust) %>% 
    complete(thermalClust, fill = list(n = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = thermalClust, values_from = n, names_prefix = "n_") %>% 
    arrange(local_identifier, weeks_since_emig, days_since_emig, burstID) 
  
  soarClust <- x %>% 
    mutate(soarClust = factor(soarClust, levels = c("soar", "glide"))) %>% #convert to a factor with three levels
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    count(soarClust) %>% 
    complete(soarClust, fill = list(n = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = soarClust, values_from = n, names_prefix = "n_") %>% 
    arrange(local_identifier, weeks_since_emig, days_since_emig, burstID) 
  
  #append all above. this dataframe will have one row per burst.
  flight_summary <- x %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    summarise(n_rows_burst = n()) %>%  #nrow of the burst
    ungroup() %>% 
    bind_cols(flapping[,-c(1:4)], thermClust[,-c(1:4)], soarClust[,-c(1:4)])  #bind columns. the order of rows is the same. remove the overlapping columns.
  
  flight_summary
  
})
  

#saveRDS(flight_summary, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_ratios_n34.rds")


# STEP 2: exploration -------------------------------------------------------------

flight_summary <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_ratios_n34.rds")

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
