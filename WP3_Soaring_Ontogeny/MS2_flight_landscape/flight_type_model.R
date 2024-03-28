#modeling flight type as a function of the terrain. to get a better idea of what variables to include when building the energy landscape. Potentially can also include 
#day/week since fledging as an explanatory variable
#Elham Nourani, PhD. May 5, 2022. Konstanz, Germany.

#update: Feb, 12. 2024: use TRI and distance to ridgeline... make boxplots

library(tidyverse)
library(sf)
library(corrr)
library(mapview)
library(terra)

# ----------- STEP 1: open data with soaring flight assignments -----------------

#open terrain layers: terrain ruggedness index and distance to ridgeline

TRI_100 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/TRI_100_LF.tif")
ridge_100 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif")


#open bird flight data
flight <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/segmented_eagle_tracks/", pattern = ".rds", full.names = T)

#explore one individual first

load("/home/enourani/Desktop/Golden_Eagle_data/gps_acc_age/Almen19 (eobs 7001)_gps_acc_age.RData") #gaa_df 

one_ind <- flight[[10]] %>% 
  readRDS() %>% 
  #filter(flightClust_smooth3 %in% c("circular soaring", "linear soaring")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = "epsg:4326") %>% 
  as("SpatVector") %>% 
  terra::project(crs(ridge_100)) #it is faster to reproject this to the terrain and then back to wgs, than to project the terrain to wgs... 
  
#extract values for tracking points
terrain_ann <- one_ind %>% 
  extract(x = TRI_100, y = ., method = "simple", bind = T) %>% 
  extract(x = ridge_100, y = ., method = "simple", bind = T) 

terrain_df <- terrain_ann %>% 
  data.frame(., geom(.))

# ----------- STEP 2: box plots -----------------

ggplot(aes(y = TRI, x = flightClust_smooth3), data = terrain_df) + 
  geom_boxplot() +
 # geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) +
  theme_minimal()

ggplot(aes(y = distance_to_ridge_line_mask, x = flightClust_smooth3), data = terrain_df) + 
  geom_boxplot() +
  # geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.8) +
  theme_minimal()



#%%% old
# ----------- STEP 2: cluster/PCA to classify flight type based on the terrain -----------------

#make density plots comparing the two categories: all look very similar

ggplot(tr_df, aes(x = TRI)) +
  geom_density(aes(color = thermalClust)) +
  scale_color_manual(values = c("gray48", "gold2", "forestgreen"))+
  facet_wrap(~stage)
  
ggplot(tr_df, aes(x = dem)) +
  geom_density(aes(color = thermalClust)) +
  scale_color_manual(values = c("gray48", "gold2", "forestgreen"))+
  facet_wrap(~stage)

ggplot(tr_df, aes(x = slope_TPI)) +
  geom_density(aes(color = thermalClust)) +
  scale_color_manual(values = c("gray48", "gold2", "forestgreen"))+
  facet_wrap(~stage)

ggplot(tr_df, aes(x = aspect_TPI)) +
  geom_density(aes(color = thermalClust)) +
  scale_color_manual(values = c("gray48", "gold2", "forestgreen"))+
  facet_wrap(~stage)

ggplot(tr_df, aes(x = slope)) +
  geom_density(aes(color = thermalClust)) +
  scale_color_manual(values = c("gray48", "gold2", "forestgreen")) +
  facet_wrap(~stage)


#check for correlation
tr_df %>% 
  dplyr::select(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100")) %>% 
  correlate() #slope and TRI are correlated (.98); slope_tpi is correlated with slope and TRI

lm(as.factor(soarClust) ~ TRI_100 + TPI_100 + aspect_TPI_100, data = tr_df)
