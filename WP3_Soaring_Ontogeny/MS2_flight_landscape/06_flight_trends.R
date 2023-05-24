#the OG script was for estimating trends in flight behavior of golden eagles to quantify improvement in soaring flight. 
#Now it looks at 1) geomorphological characteristics of the flight types, and 
# 2) the characteristics of the commuting flights.
#Elham Nourani, PhD. Konstanz, DE.
#Feb. 7. 2023
#to be submitted to the cluster Raven
#update April 18th 2023: using the new behavioral segmentation done by Martina in March 2023

library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(mapview)
library(parallel)
library(ggridges)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

# STEP 1: Open data and subset for inds in the landscape analysis -------------------------------------------------------------
#files <- list.files("/home/enourani/Desktop/GE-classification_Mar23", pattern = ".rds", full.names = T)

#behav_data <- lapply(files, readRDS)

#saveRDS(behav_data, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/Golden_eagle_wind_behav_full_Mar_2023.rds")

#individuals in the landscape analysis
inds_to_keep <- readRDS("all_inds_annotated_static_apr23.rds") %>% 
  distinct(individual.local.identifier)%>% 
  pull(individual.local.identifier)

#read in data
behav_data <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/Golden_eagle_wind_behav_full_Mar_2023.rds")
names(behav_data) <- sapply(behav_data, function(x) unique(x$individual.local.identifier))

behav_ls <- behav_data[names(behav_data) %in% inds_to_keep] #only 42 of the 55 individuals are in here

behav_df <- lapply(behav_ls, function(x){
  if("eobs.status" %in% names(x)){ x <- x %>%  dplyr::select(-"eobs.status")} #some elements of the list have one extra column
  x
}) %>% 
  reduce(rbind)

saveRDS(behav_df, file = "behav_segmented_42inds.rds")

#emigration date info
emig_dates <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/fleding_emigration_timing_Mar2023.rds") #this file only has the juvies.

# STEP 2: summarize geo characteristics of each flight type ----------------------------------------------------------------

# 2.1: extract flight categories of interest
soaring_only <- behav_df %>% 
  mutate(flightClust_smooth3 = as.character(flightClust_smooth3)) %>% 
  mutate(flight_mode = ifelse(behavior == "Flapping", "flapping",
                              flightClust_smooth3)) %>% 
  filter(flight_mode %in% c("linear soaring", "circular soaring"))


# 2.2: extract terrain values at gps points (25 m resolution terrain. prepared by Louise) 
ridge25 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/map_distance_ridge.tif")
TRI_25 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/TRI/TRI.sdat")
dem_25 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/region-alpes-dem.tif")
slope_TPI_25 <- terrain(terrain(dem_25, "slope"), "TPI", filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/slope_TPI_25.tif")
  
#create a stack using raster paths
topo <- c(ridge25, dem_25, TRI_25, slope_TPI_25)
names(topo) <- c("ridge25", "dem_25", "TRI_25", "slope_TPI_25")

#reproject tracking data to match topo, extract values from topo, convert back to wgs and save as a dataframe
topo_ann_df <- soaring_only %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_transform(crs = crs(topo)) %>% 
  extract(x = topo, y = ., method = "simple", bind = T) %>%
  project(wgs) %>% 
  data.frame(., geom(.)) %>% 
  dplyr::select(-c("geom", "part", "hole")) %>% 
  rename(location.long = x,
         location.lat = y)

saveRDS(topo_ann_df, file = "behav_segmented_terrain_42inds.rds")

# 2.3: summarize terrain characteristics at soaring locations

topo_ann_df <- readRDS("behav_segmented_terrain_42inds.rds") %>% 
  mutate(month = month(timestamp)) %>% 
  filter(month %in% c(5:8)) #only keep summer

##plots
to_plot <- topo_ann_df %>% 
  pivot_longer(cols = ends_with("25"), names_to = "variable", values_to = "value")


ggplot(to_model, aes(x = ridge25, y = as.factor(flight_mode))) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
  )


ggplot(data = to_plot, aes(x = variable, y = value, fill = flight_mode)) +
  geom_boxplot()

#models
to_model <- topo_ann_df %>% 
  mutate(circling = ifelse(flight_mode == "circular soaring", 1, 0)) %>%  #create new binomial response with 0 and 1
  drop_na(c("ridge25", "dem_25", "TRI_25", "slope_TPI_25"))

m <- glm(circling ~ scale(ridge25) + scale(dem_25) + scale(TRI_25) + scale(slope_TPI_25),
         data = to_model, family = binomial)

plot_summs(m)
#all vars: AIC: 7089395

m2 <- glm(circling ~ scale(ridge25) + scale(dem_25),
         data = to_model, family = binomial)

plot_summs(m2)
#AIC: 7089942


#effect plots:
effect_plot(m2, pred = ridge25, interval = TRUE, y.label = "prob of circling", x.label = "distance to ridge", ylim = c(0,1))
effect_plot(m2, pred = dem_25, interval = TRUE, y.label = "prob of circling", x.label = "elevation")
#effect_plot(m, pred = TRI_25, interval = TRUE, y.label = "prob of circling", x.label = "TRI")
#effect_plot(m, pred = slope_TPI_25, interval = TRUE, y.label = "prob of circling", x.label = "slope TPI")


# STEP 2: trends in commuting flight characteristics ----------------------------------------------------------------

data <- readRDS("all_inds_annotated_static_apr23.rds") 

summs <- data %>% 
  filter(used == 1) %>% 
  group_by(individual.local.identifier, date(timestamp), burst_id) %>% #each bout of continuous commuting flight has a unique ID. These are bursts with 1hr difference and no more.
  arrange(desc(timestamp)) %>% 
  summarize(duration = as.numeric(difftime(head(timestamp, 1), tail(timestamp, 1), units = "hours")),
            weeks_since_emig = head(weeks_since_emig, 1))
  
plot(x = summs$weeks_since_emig, y = summs$duration)
  
  
  
  
  
  
### old stuff


# STEP 1: Summarize flight metrics per burst -------------------------------------------------------------

#calculate total n-rows of each behavior category for each burst 

(b <- Sys.time())
flight_summary_ls <- lapply(behav_data, function (x){
  
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
flight_summary_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/flight_summary_per_burst_full.rds") %>% 
  reduce(rbind) %>% 
  filter(between(weeks_since_emig, 1, 135)) %>% 
  mutate(max_wind_speed = sqrt(max_windx^2 + max_windy ^2),
         med_wind_speed = sqrt(med_windx^2 + med_windy^2),
         flapping_ratio = n_Flapping/n_NotFlapping,
         flapping_to_all = n_Flapping/n_rows_burst,
         circling_to_all = n_circular/n_rows_burst,
         circling_to_linear = n_circular/n_linear)


#summaries per day OR week
wk_summaries <- flight_summary_df %>% 
  group_by(local_identifier, weeks_since_emig) %>%
  summarize_at(c("max_wind_speed","med_wind_speed","flapping_ratio", "flapping_to_all", "circling_to_all", "circling_to_linear"), mean, na.rm = T)

# STEP 2:plots -------------------------------------------------------------

ggplot(data = wk_summaries, aes(x = weeks_since_emig, y = circling_to_all, group = local_identifier, color = local_identifier)) +
  geom_point(show.legend = F) +
  geom_smooth(method = "lm", show.legend = F)

#overall trend

ggplot(data = wk_summaries, aes(x = weeks_since_emig, y = circling_to_all)) +
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
