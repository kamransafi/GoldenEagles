### estimate wind vectors for golden eagle data
### Hester Bronnvik

library(move)
library(moveWindSpeed)
library(geosphere)
# library(data.table)
library(sf)
library(tidyverse)

# get the functions for the MoveWindSpeed estimations (https://github.com/anflack/Updated-moveWindSpeed-)
source("/home/hbronnvik/Documents/chapter2/getWindEstimates_update.R")
# updated function using avg. eagle values
source("/home/hbronnvik/Documents/chapter4/thermallingFeaturesFunction.R")
source("/home/hbronnvik/Documents/chapter2/getTrackSegments_updated.R")
source("/home/hbronnvik/Documents/chapter2/getWindEstimate_update.R")
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 15), 
                                text = element_text(size = 17),
                                panel.background = element_rect(fill = "#F5F5F5"),
                                strip.background = element_rect(fill = "#F5F5F5", color = "#F5F5F5"),
                                strip.text = element_text(face = "bold")))
colfunc <- colorRampPalette(c("#6F101E", "#8F1D28", "#BF413E", "#FACEB6", "white", 
                              "#A4CCE3", "#3D84B4", "#1D5497", "#0A3162"))
wgs <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")# map projection

# functions:
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}
wind_direction <- function(u, v){(90-(atan2(v, u)*(180/pi)))%%360 }

# classified flight 
fls <- list.files("/home/hbronnvik/Documents/GE_data/seg24/", full.names = T)

## pick the thermals that last long enough for wind vector estimation
thermal_results <- lapply(fls, function(f){
  tr <- readRDS(f)
  # identify separate thermaling events
  tr <- tr %>% 
    filter(thermalClust == "circular" & soarClust == "soar") %>% 
    mutate(date = date(timestamp), 
           time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           time_diff = ifelse(is.na(time_diff), 1, time_diff),
           min_split = time_diff > 60,
           thermal_event = paste(individual.id, burstID, cumsum(min_split)+1, sep = "_")) %>% 
    dplyr::select(-time_diff, -min_split) %>% 
    group_by(thermal_event) %>% 
    mutate(thermal_duration = n(),
           vspeed_thermal = (height.above.ellipsoid[n()]-height.above.ellipsoid[1])/thermal_duration,
           turn_var_thermal = var(turn.angle, na.rm = T)/thermal_duration) %>% 
    ungroup()
  # select thermals that can be used for wind estimates
  tr <- tr %>% 
    # thermals that lasted at least 30 seconds and in which the bird climbed
    filter(thermal_duration >= 30 & vspeed_thermal > 0)
  return(tr)
})
# 52 individuals left
thermal_results <- thermal_results[!sapply(thermal_results, function(x) is.null(x))]

# estimate wind vectors according to Weinzierl et al. 2016,  https://doi.org/10.1002/ece3.2585
library(parallel)
wind_results <- lapply(1:length(thermal_results), function(n){
  print(n)
  tr <- thermal_results[[n]] %>%
    rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
    group_by(thermal_event) %>% 
    group_split()
  ind <- lapply(tr, function(burst){
    mv_burst <- move(x=burst$location.long, y=burst$location.lat, 
                     time=burst$timestamp,
                     proj="+proj=longlat +datum=WGS84 +no_defs",
                     animal=burst$individual.id,
                     data=burst)
    burst_summary <- burst %>%
      mutate(gap = c(1,timeLag(mv_burst)),
             consistent = gap == 1,
             event = cumsum(consistent==0)) %>%
      group_by(event) %>%
      summarize(obs = n()) %>% 
      mutate(sufficient = obs>29)
    if(T %in% unique(burst_summary$sufficient)){
      class_burst <- getWindEstimates(mv_burst)
      class_burst <- as.data.frame(class_burst)
      return(class_burst)
    }
  })
  ind <- ind[!sapply(ind, function(x) is.null(x))]
  ind <- data.table::rbindlist(ind, use.names = T)
  saveRDS(ind, file = paste0("/home/hbronnvik/Documents/GE_data/wind_estimates_25/", unique(ind$individual.id), "_seg_wind_",Sys.Date(),".rds"))
  return(ind)
})

summary(wind_speed(wind_results[[1]]$windX, wind_results[[1]]$windY))
