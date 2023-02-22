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

# STEP 1: Summarize flight metrics per burst -------------------------------------------------------------

#calculate total n-rows of each behavior category for each burst 

(b <- Sys.time())
flight_summary_ls <- lapply(flight, function (x){
  
  #b <- Sys.time()
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
    group_by(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) %>% 
    count(behavior) %>% 
    complete(behavior, fill = list(n = 0)) %>%  #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = behavior, values_from = n, names_prefix = "n_") %>%  #widen the dataframe. have one column for each category of behavior. nrow is not the n of unique bursts
    arrange(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) #just an insurance policy. the rows are already ordered by these.
  
  thermClust <-  x %>% 
    mutate(thermalClust = factor(thermalClust, levels = c("circular", "linear", "other"))) %>% #convert to a factor with three levels
    group_by(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) %>% 
    count(thermalClust) %>% 
    complete(thermalClust, fill = list(n = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = thermalClust, values_from = n, names_prefix = "n_") %>% 
    arrange(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) 
  
  soarClust <- x %>% 
    mutate(soarClust = factor(soarClust, levels = c("soar", "glide"))) %>% #convert to a factor with three levels
    group_by(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) %>% 
    count(soarClust) %>% 
    complete(soarClust, fill = list(n = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    ungroup() %>% 
    pivot_wider(names_from = soarClust, values_from = n, names_prefix = "n_") %>% 
    arrange(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) 
  
  #append all above. this dataframe will have one row per burst. Also calculate min, max & median of vertical speed, odba, and wind u and v estimated from the thermals.
  flight_summary <- x %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, burstID) %>% 
    summarise(med_odba = median(odbaMedian, na.rm = T),
              max_odba = max(odbaMedian, na.rm = T),
              min_odba = min(odbaMedian, na.rm = T),
              med_vertspd = median(vertSpeed_smooth, na.rm = T),
              max_vertspd = max(vertSpeed_smooth, na.rm = T),
              min_vertspd = min(vertSpeed_smooth, na.rm = T),
              med_windx = median(windX, na.rm = T),
              max_windx = max(windX, na.rm = T),
              min_windx = min(windX, na.rm = T),
              med_windy = median(windY, na.rm = T),
              max_windy = max(windY, na.rm = T),
              min_windy = min(windY, na.rm = T),
              n_rows_burst = n()) %>%  #nrow of the burst
    ungroup() %>% 
    bind_cols(flapping[,-c(1:5)], thermClust[,-c(1:5)], soarClust[,-c(1:5)]) %>%   #bind columns. the order of rows is the same. remove the overlapping columns.
    rename(n_flap_Unclassified = n_Unclassified, #rename vague columns for clarity
           n_soaring_other = n_other)
  
  #Sys.time() -b #59 secs per individual
  flight_summary

})

Sys.time() - b #23 min

saveRDS(flight_summary_ls, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_summary_per_burst_full.rds")


# STEP 2: calculate flight metrics for post-dispersal -------------------------------------------------------------


#transform to a dataframe for easier plotting
flight_summary_df <- flight_summary_ls %>% 
  reduce(rbind) %>% 
  filter(between(weeks_since_emig, 1, 135)) %>% 
  mutate(max_wind_speed = sqrt(max_windx^2 + max_windy ^2),
         flapping_ratio = n_Flapping/n_NotFlapping,
         flapping_to_all = n_Flapping/n_rows_burst,
         circling_to_all = n_circular/n_rows_burst,
         circling_to_linear = n_circular/n_linear) #calculate wind speed



ratios <- lapply(flight_summary_ls, function(x){
  
  #make plots for the whole period?
  
  #filter for post-dispersal and calculate metrics
  x_post_d <- x %>% 
    filter(days_since_emig > 0) %>% 
    
  
  
})



# STEP 2: exploration -------------------------------------------------------------

flight_summary_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_summary_per_burst_full.rds") %>% 
  reduce(rbind) %>% 
  filter(weeks_since_emig <= 135) #filter for weeks since dispersal to match the landscape scale

dd <- flight_summary_df %>%  filter(circling_to_all > 0)

ggplot(data = flight_summary_df %>%  filter(circling_to_all > 0), 
       aes(x = 1/weeks_since_emig, y = circling_to_all, group = local_identifier, color = local_identifier)) +
  geom_jitter() +
  geom_smooth(method = "lm") 


ggplot(data = flight_summary_df %>% drop_na(circling_to_all) %>% filter(is.finite(circling_to_all) ), 
       aes(x = weeks_since_emig, y = circling_to_all, group = local_identifier, color = local_identifier)) +
  geom_jitter() +
  geom_smooth(method = "lm") 

ggplot(data = flight_summary_df %>% drop_na(flapping_to_all) %>% filter(is.finite(flapping_to_all)), 
       aes(x = days_since_fled, y = flapping_to_all, group = local_identifier, color = local_identifier)) +
  geom_jitter() +
  geom_smooth(method = "lm") 


ggplot(data = flight_summary_df %>%  , aes(x = weeks_since_emig, y = circling_to_linear)) +
  geom_jitter() +
  geom_smooth(method = "gam") 

+
  facet_wrap(~ thermalClust) +
  theme_classic()

#circling vs linear soaring
ggplot(data = flight_summary %>% filter(weeks_since_emig <= 135 & !is.na(circling_to_linear) & !is.infinite(circling_to_linear)) %>% group_by(local_identifier, days_since_emig,burstID) %>%  slice(1), 
       aes(x = weeks_since_emig, y = circling_to_linear)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_classic()
