#modeling flight type as a function of the terrain. to get a better idea of what variables to include when building the energy landscape. Potentially can also include 
#day/week since fleding as an explanatory variable
#Elham Nourani, PhD. May 5, 2022. Konstanz, Germany.

library(tidyverse)
library(sf)
library(corrr)
library(mapview)


# ----------- STEP 1: open data with soaring flight assignments -----------------

#build a binomial logistic regression, with thermal soaring as 0 and slope soaring as 1
#or, even better, make a multinomial model#open terrain layers and put onto one stack. there was no difference between flight types regarding terrain at 100 meters. try the higher res
ls <- list.files("/home/enourani/Desktop/golden_eagle_static_layers", pattern = ".tif", full.names = T)
High_res <- ls[grep(ls, pattern = "100", invert = T)]

dem <- rast("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/eu_dem_vdem_10011_E40N20.TIF")

stck <- lapply(High_res, rast) %>% 
  rast()


#open matched data with wind speed and flight type: one individual first
load("/home/enourani/Desktop/Golden_Eagle_data/gps_acc_age/Almen19 (eobs 7001)_gps_acc_age.RData") #gaa_df 
one_ind <- gaa_df %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
  as("SpatVector") %>% 
  project(crs(stck)) #it is faster to reproject this to the terrain and then back to wgs, than to project the terrain to wgs... also when I save
                     #the topo_wgs, then I can't open it anymore. some permission errors with the new POP OS version
  
#extract values for tracking points
tr_ann <- cbind(one_ind, extract(x = stck, y = one_ind, method = "bilinear"))

saveRDS(tr_ann, file = "one_ind_terrain_annotation.rds") #raster file (100m res)

tr_df <- tr_ann %>%
  project("+proj=longlat +datum=WGS84") %>% #project back to wgs 
  as.data.frame(geom = "XY") #keep coordinates in the df

saveRDS(tr_df, file = "one_ind_terrain_annotation_df.rds") #df file (100m res)
saveRDS(tr_df, file = "one_ind_terrain_annotation_df_hres.rds") #df file (original res)

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
