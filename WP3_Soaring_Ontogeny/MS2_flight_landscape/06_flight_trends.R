#script for estimating trends in flight behavior to quantify improvement in soaring flight.
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

#calculate ratios for each burst, instead of each day

ratios <- lapply(flight, function (x){
  
  #calculate week and day since emigration
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), funs(ceiling(as.numeric(.))))
  
  #estimate n of unique days and n of bursts per day
  burst_IDs <- x %>% 
    group_by(days_since_emig) %>% 
    summarise(n_bursts = n_distinct(burstID)) %>% 
    filter(n_bursts >= thr)
  
  #filter the data
  x_subset <- x %>% 
    filter(days_since_emig %in% burst_IDs$days_since_emig)
  
  x_subset
  
}) %>% 
  reduce(rbind)




########################## daily stuff. no clear pattern. very slight increase in the linear/circling trend ###################
# STEP 1: decide on the threshold of the number of bursts per day -------------------------------------------------------------

#If below the threshold, don't include the day.
n_ls <- lapply(flight, function (x){
  
  #calculate week and day since emigration
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), funs(ceiling(as.numeric(.))))
  
  #calculate the ratio of each mode of flight to total flight time (as number of rows)
  ratios <- x %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    mutate(total_flight_in_burst = n()) %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig,burstID, thermalClust) %>% 
    summarize(amount_in_burst = n(),
              total_flight_in_burst = head(total_flight_in_burst,1)) %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    complete(thermalClust, fill = list(amount_in_burst = 0)) %>% #add zero for instances where a factor level was not observed in a burst
    mutate(ratio = amount_in_burst/total_flight_in_burst)
  
  #plot----
  ggplot(data = ratios, aes(x = weeks_since_emig, y = ratio)) +
    geom_jitter() +
    geom_smooth(method = "lm") +
    facet_wrap(~ thermalClust) +
    theme_classic()
  #--------
 
  #calculate the ratio of thermal soaring to linear soaring (as number of rows)
  soaring_ratio <- x %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID, thermalClust) %>% #group by all these variables to make sure they are retained in the final data frame
    summarize(n_rows = length(thermalClust)) %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, burstID) %>% 
    complete(thermalClust, fill = list(n_rows = 0)) %>% #from library(tidyr). wherever one factor level had no rows, insert a zero
    summarize(soaring_ratio = n_rows[thermalClust == "circular"]/ n_rows[thermalClust == "linear"])
  
  #plot----
  ggplot(data = soaring_ratio, aes(x = weeks_since_emig, y = soaring_ratio)) +
    geom_jitter() +
    geom_smooth(method = "lm")
  #--------
  
  
})
  
  
  
  
  
  #------------------------------old stuff
  #calculate total flight and amount of flight in each category as the number of rows 
  flight_summary <- data %>% 
    group_by(local_identifier, days_since_emig) %>% 
    mutate(total_flight_in_day = n()) %>% 
    group_by(local_identifier, days_since_emig, flight_type) %>% 
    summarize(amount_in_day = n(),
              total_flight_in_day = head(total_flight_in_day,1),
              weeks_since_emig = head(weeks_since_emig,1))
  
  ratios <- flight_summary %>% 
    mutate(general_behaviors = ifelse(flight_type %in% c("soar_circular", "circular_no_ID"), "soaring", flight_type)) %>% 
    group_by(local_identifier, days_since_emig) %>%
    mutate(ratio = amount_in_day/total_flight_in_day)
  
  #just look at thermal clust: circular and linear (omit other. which is probably gliding)
  ratio_soaring <-  data %>% 
    group_by(local_identifier, days_since_emig) %>% 
    mutate(total_flight_in_day = n()) %>% 
    filter(thermalClust %in% c("circular", "linear")) %>% #get rid of the OTHER category
    group_by(local_identifier, days_since_emig, thermalClust) %>% 
    summarize(amount_in_day = n(),
              total_flight_in_day = head(total_flight_in_day,1),
              weeks_since_emig = head(weeks_since_emig,1)) %>% 
    arrange(thermalClust, .by_group = T) %>% 
    summarize(ratio = tail(amount_in_day,1)/head(amount_in_day,1),
              weeks_since_emig = head(weeks_since_emig,1)) #i.e. linear/circular
  
  
  
  
  
}) %>% 
  reduce(rbind)

#plot
hist(n_ls$n_bursts)

#take the third quartile as the threshold
thr <- summary(n_ls$n_bursts)[5] %>%  as.numeric() 

# STEP 2: subset the data and keep days that have more than 19 bursts -------------------------------------------------------------

data <- lapply(flight, function (x){
  
  #calculate week and day since emigration
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), funs(ceiling(as.numeric(.))))
  
  #estimate n of unique days and n of bursts per day
  burst_IDs <- x %>% 
    group_by(days_since_emig) %>% 
    summarise(n_bursts = n_distinct(burstID)) %>% 
    filter(n_bursts >= thr)
  
  #filter the data
  x_subset <- x %>% 
    filter(days_since_emig %in% burst_IDs$days_since_emig)
  
  x_subset
  
}) %>% 
  reduce(rbind)

saveRDS(data, "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_data_20bursts.rds")

# STEP 3: exploration -------------------------------------------------------------

#calculate total flight and amount of flight in each category as the number of rows 
flight_summary <- data %>% 
  group_by(local_identifier, days_since_emig) %>% 
  mutate(total_flight_in_day = n()) %>% 
  group_by(local_identifier, days_since_emig, flight_type) %>% 
  summarize(amount_in_day = n(),
            total_flight_in_day = head(total_flight_in_day,1),
            weeks_since_emig = head(weeks_since_emig,1))

ratios <- flight_summary %>% 
  mutate(general_behaviors = ifelse(flight_type %in% c("soar_circular", "circular_no_ID"), "soaring", flight_type)) %>% 
  group_by(local_identifier, days_since_emig) %>%
  mutate(ratio = amount_in_day/total_flight_in_day)

#just look at thermal clust: circular and linear (omit other. which is probably gliding)
ratio_soaring <-  data %>% 
  group_by(local_identifier, days_since_emig) %>% 
  mutate(total_flight_in_day = n()) %>% 
  filter(thermalClust %in% c("circular", "linear")) %>% #get rid of the OTHER category
  group_by(local_identifier, days_since_emig, thermalClust) %>% 
  summarize(amount_in_day = n(),
            total_flight_in_day = head(total_flight_in_day,1),
            weeks_since_emig = head(weeks_since_emig,1)) %>% 
  arrange(thermalClust, .by_group = T) %>% 
  summarize(ratio = tail(amount_in_day,1)/head(amount_in_day,1),
            weeks_since_emig = head(weeks_since_emig,1)) #i.e. linear/circular



#plots

ggplot(data = ratios, aes(x = days_since_emig, y = ratio)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  facet_wrap(~ general_behaviors) +
  theme_classic()

ggplot(data = ratio_soaring, aes(x = weeks_since_emig, y = ratio)) +
  #ylim(0,20) +
  geom_jitter() +
  geom_smooth(method = "lm")

# STEP 2b: if only looking at ratio of linear/circling, no need to filter for days with certain n of bursts. OR try to clac ratios per burst instead of per day.....