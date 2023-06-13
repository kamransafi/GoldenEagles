#make a map of the golden eagle data distribution
#May 18, 2023. Konstanz, De.
#Elham Nourani

library(tidyverse)
library(leaflet)
library(mapview)
library(sf)
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
clr2 <- oce::oceColorsPalette(100)[80]

# Background 1: NASA


m <- leaflet() %>% 
  addTiles() %>% 
  addScaleBar() %>% 
  setView(lng = 10.5, lat = 46.5, zoom = 7.5) %>% 
  addCircles(data = data_sp, color = "", radius = 0.05, fillOpacity = 0.8, fillColor = clr) %>% 
  addProviderTiles("Esri.WorldImagery")

m


m2 <- leaflet() %>% 
  addTiles() %>% 
  addScaleBar() %>% 
  setView(lng = 10.5, lat = 46.5, zoom = 5.3) %>% 
  addProviderTiles("Esri.WorldImagery")

m2
# "#af46b7"
#mapview::mapshot(m, file = "/GE_distr_map.png")

#mapshot(m, file = "~/GE_distr_map.png")


#create a larger area map

#library(rworldmap)

world <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
 # st_crop(xmin = -10, xmax = 40, ymin = 30, ymax = 60) %>%
  st_union()

X11()

png("/home/enourani/ownCloud/Work/conferences/GRC_GRS_23/Europe_map.png", width = 12, height = 12, units = "in", res = 400)
plot(world, fill = "white", col = "black")
polygon(y = c(45,46.65, 46.65, 45.2), x = c(7.2, 7.2, 13.89, 13.89), border = "white", add = T)
dev.off()


library(marmap)
library(ggspatial)

png("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/GRC_GRS_23/Europe_map2.png", width = 6, height = 3.7, units = "in", res = 400)
ggplot() +
  geom_sf(data = ld, fill = "grey40", col = NA) +
  xlim(c(-6, 28)) + ylim(c(37, 52)) +
  geom_polygon(data = data.frame(y = c(45.2,46.65, 46.65, 45.2), x = c(7.2, 7.2, 13.89, 13.89)), aes(x = x, y = y), 
               color = clr2, fill = NA, linewidth = 1.8) +
  ggspatial::annotation_scale(
    location = "bl",
    bar_cols = c("black", "white"),
    text_family = "ArcherPro Book",
    text_col = "black") +
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"), style = north_arrow_fancy_orienteering(line_col = "black", line_width = 1.3)) +
  theme_void()

dev.off()


###########prep topographic maps for GRC poster
TRI_100 <- raster("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/TRI_100_LF.tif")
ridge_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif")

TRI_df <- as.data.frame(TRI_100, xy = T)

ridge_df <- as.data.frame(ridge_100, xy = T)

library(rasterVis)


levelplot(TRI_100, margin=FALSE) 

levelplot(ridge_100, margin=FALSE) 



png("/home/enourani/ownCloud/Work/conferences/GRC_GRS_23/TRI_map.png", width = 4.5, height = 3, res = 300, unit = "in")
ggplot(TRI_df) +
  geom_tile(aes(x = x, y = y, fill = TRI)) +
  scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,450), name = "") +
  theme_void() +
  theme(legend.position = "bottom")
dev.off()

png("/home/enourani/ownCloud/Work/conferences/GRC_GRS_23/ridge_map2.png", width = 4.5, height = 3, res = 300, unit = "in")
ggplot(ridge_df) +
  geom_tile(aes(x = x, y = y, fill = distance_to_ridge_line_mask)) +
  #scale_fill_gradientn2(colours = oce::oceColorsPalette(200), limits = c(0,9000), midpoint = 1500, name = "") +
  scale_fill_gradient2(low = clr, high = clr2, midpoint = 3000, name = "") +
  theme_void() +
  theme(legend.position = "bottom")
dev.off()


#exponential growth and plateau for the GRS talk
y <- c(1,1.5,2,2.5, 3,4,5,6,6,6,6)
x <- c(1:11)


