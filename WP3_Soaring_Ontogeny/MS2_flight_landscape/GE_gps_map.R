#make a map of the golden eagle data distribution
#May 18, 2023. Konstanz, De.
#Elham Nourani

library(tidyverse)
library(leaflet)
library(mapview)
library(sp)
library(oce)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")


#open ssf data to filter for individuals included in the study
ssf_inds <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  distinct(individual.local.identifier) %>% 
  pull(individual.local.identifier)

#open data
data <- readRDS("/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack_Golden_Eagle_Alps_1min.rds") %>% 
  filter(individual.local.identifier %in% ssf_inds) %>% 
  drop_na(location.lat) %>% 
  mutate(dt_1hr = round_date(timestamp, "1 hour")) %>% 
  group_by(individual.local.identifier,dt_1hr) %>% 
  slice(1) 

saveRDS(data, file = "GE_gps_1hr.rds")

data <- readRDS("GE_gps_1hr.rds")

data_sp <- data
coordinates(data_sp) <-~ location.long + location.lat

#ridge_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif")

# Note: if you do not already installed it, install it with:
# install.packages("leaflet")

clr <- oce::oceColorsPalette(100)[2]
clr <- oce::oceColorsPalette(100)[80]

# Background 1: NASA
# m <- leaflet() %>% 
#   addTiles() %>% 
#   setView( lng = 11, lat = 46, zoom = 5 ) %>% 
#   addProviderTiles("NASAGIBS.ViirsEarthAtNight2012")
# m
# 
# # Background 2: World Imagery
# m <- leaflet() %>% 
#   addTiles() %>% 
#   addScaleBar() %>% 
#   setView(lng = 11, lat = 46, zoom = 6.5) %>% 
#   addCircles(data = data_sp, color = clr, radius = 0.2, fillOpacity = 0.6, fillColor = clr) %>% 
#   addProviderTiles("OpenStreetMap")
# 
# m

m <- leaflet() %>% 
  addTiles() %>% 
  addScaleBar() %>% 
  setView(lng = 10.5, lat = 46.5, zoom = 7.5) %>% 
  addCircles(data = data_sp, color = "", radius = 0.05, fillOpacity = 0.8, fillColor = clr) %>% 
  addProviderTiles("Esri.WorldImagery")

m

# "#af46b7"
mapview::mapshot(m, file = "/GE_distr_map.png")

mapshot(m, file = "~/GE_distr_map.png")

