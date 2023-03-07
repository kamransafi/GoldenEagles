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

cov <- function(data){(sd(data, na.rm = T) / mean(data, na.rm = T)) * 100}

#emigration date info
load("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/em_fl_dt_recurse_70ind.RData") #emig_fledg_dates

# STEP 1: Summarize flight metrics per burst -------------------------------------------------------------

#calculate total n-rows of each behavior category for each burst 

(b <- Sys.time())
th_summary <- lapply(flight, function (x){
  
  #calculate week and day since emigration and fledging
  ind_dates <- emig_fledg_dates %>% 
    filter(individual.local.identifier == x$local_identifier[1])
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("days")),
           weeks_since_emig = difftime(timestamp,ind_dates$emigration_dt, units = c("weeks")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("days")),
           weeks_since_fled = difftime(timestamp,ind_dates$fledging_dt, units = c("weeks"))) %>% 
    mutate_at(vars(matches("_since_")), ~ ceiling(as.numeric(.)))
  
  #filter for only thermalling attempts
  th_summary <- x %>% 
    filter(!is.na(thermalID_burstID)) %>% 
    group_by(local_identifier, weeks_since_emig, days_since_emig, days_since_fled, thermalID_burstID) %>% 
    summarize(vert_sp_cov = cov(vert.speed),
              s_vert_sp_cov = cov(vertSpeed_smooth),
              ta_cov = cov(turn.angle),
              s_ta_cov = cov(turnAngle_smooth),
              max_wind_y = max(windY, na.rm = T),
              max_wind_x = max(windX, na.rm = T))

  th_summary

}) %>% 
  reduce(rbind)

Sys.time() - b #7 sec

saveRDS(th_summary, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/thermal_summary.rds")


# STEP 2: calculate flight metrics for post-dispersal -------------------------------------------------------------


#transform to a dataframe for easier plotting
th_to_plot <- th_summary %>% 
  #filter(between(weeks_since_emig, 1, 135)) %>% 
  mutate(max_wind_speed = sqrt(max_wind_x^2 + max_wind_y ^2),
         med_wind_speed = sqrt(max_wind_x^2 + max_wind_y^2))

# STEP 2:plots -------------------------------------------------------------

ggplot(data = th_to_plot, aes(x = weeks_since_emig, y = med_wind_speed, group = local_identifier, color = local_identifier)) +
  geom_point(show.legend = F) +
  geom_smooth(method = "gam", show.legend = F)

#overall trend

ggplot(data = th_to_plot, aes(x = weeks_since_emig, y = max_wind_speed)) +
  geom_point(show.legend = F) +
  ylim(0,1) +
  geom_smooth(method = "lm", show.legend = F)

# STEP xx: some old exploration -------------------------------------------------------------

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
