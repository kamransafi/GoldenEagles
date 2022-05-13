#create topographic layers for the Alps
#May13, 20222. Konstanz, DE
#Elham Nourani, PhD

library(tidyverse)
library(terra)
library(sf)

#download dem tiles E30N20 and E40N20 from EU-DEM: https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1?tab=mapview 

#estimate terrain values for French Alps

#open tiles
#dem1 <- rast("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF") %>% 
#  crop(ext(4023480,4834755,2188227,2849983)) #I did all of this procedure for this tile in data_

dem2 <- rast("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E30N20/eu_dem_v11_E30N20.TIF") %>% 
  crop(ext(3862443, 4017499, 2149463, 2614631))



slope <- terrain(dem2, v = "slope", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/slope.tif")
aspect <- terrain(dem2, v = "aspect", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/aspect.tif")
TRI <- terrain(dem2, v = "TRI", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/TRI.tif")
TPI <- terrain(dem2, v = "TPI", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/TPI.tif")

#estimate slope and aspect unevenness
sl_uneven <- terrain(slope, v = "TPI", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/slope_TPI.tif")
as_uneven <- terrain(aspect, v = "TPI", unit = "degrees", filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps.aspect_TPI.tif")

#aggregate raster cells to 100 m resolution
dem_100 <- raster::aggregate(x = dem2, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/dem_100.tif")
slope_100 <- raster::aggregate(x = slope, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/slope_100.tif")
aspect_100 <- raster::aggregate(x = aspect, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/aspect_100.tif")
sl_uneven_100 <- raster::aggregate(x = sl_uneven, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/slope_TPI_100.tif")
as_uneven_100 <- raster::aggregate(x = as_uneven, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/aspect_TPI_100.tif")
TRI_100 <- raster::aggregate(x = TRI, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/TRI_100.tif")
TPI_100 <- raster::aggregate(x = TPI, fact = 4, filename = "/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/TPI_100.tif")

#mask for alpine region only
#open dem
dem <- rast("/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/dem_100.tif")
tri <- rast("/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/TRI_100.tif")

stck <- c(dem,tri) 


#open topograptri#open topography layers
load("/home/enourani/Desktop/golden_eagle_static_layers/all_topo_100m_wgs.RData") #topo_wgs

#open Apline permitere layer
Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(stck)) %>% 
  as("SpatVector")

#mask topo layers with border of the Alps. it cuts the western corner a bit.... consider downloading the other tile as well....
topo_Alps <- stck %>% 
  mask(Alps) %>% 
  project("+proj=longlat +datum=WGS84 +no_defs")

saveRDS(topo_Alps, file = "Alps_dem_tri_wgs.rds")