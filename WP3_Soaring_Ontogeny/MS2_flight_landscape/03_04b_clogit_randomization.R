#script for constructing the energy landscape for golden eagles. 
#build one model per ind per week. Use randomization to see if there is a trend in the coefficients over time.
#alternative to 03_04_clogit_workflow.R
#04.04.2023. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(magrittr)
library(corrr)
library(jtools) #to plot coeffs of clogit
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot)
library(terra)
library(gstat) #for interpolations
#library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(patchwork) #patching up interaction plots

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")
source("/home/enourani/ownCloud/Work/Projects/functions.R")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")


#STEP 1: prep annotated data ----------------------------------------------------------
#open annotated data and join together
PL_files <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/annotations/Mar_31_23/PL",
                        pattern = ".csv", full.names = T)

PL_df <- lapply(PL_files, read.csv) %>% 
  reduce(rbind)

surface_files <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/annotations/Mar_31_23/surface/",
                       pattern = ".csv", full.names = T)

ann_df <- lapply(surface_files, read.csv) %>% 
  reduce(rbind) %>% 
  bind_cols(PL_df %>% dplyr::select(c("ECMWF.ERA5.PL.U.Wind", "ECMWF.ERA5.PL.V.Wind"))) #better practice is to include row_id with the annotatin files and use that to join files. 

data <- readRDS("alt_50_60_min_55_ind_static_100.rds") %>% 
  bind_cols(ann_df %>% dplyr::select("ECMWF.ERA5.SL.Boundary.Layer.Height", "ECMWF.ERA5.SL.Mean.Surface.Sensible.Heat.Flux", 
                                     "ECMWF.ERA5.SL.Temperature..2.m.above.Ground.", "ECMWF.ERA5.SL.Instantaneous.Moisture.Flux",
                                     "ECMWF.ERA5.PL.U.Wind", "ECMWF.ERA5.PL.V.Wind")) %>% 
  rename(u_wind = ECMWF.ERA5.PL.U.Wind,
         v_wind = ECMWF.ERA5.PL.V.Wind,
         blh = ECMWF.ERA5.SL.Boundary.Layer.Height,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.) %>% 
  mutate(wind_support= wind_support(u = u_wind, v = v_wind, heading = heading),
         cross_wind= cross_wind(u = u_wind, v = v_wind, heading = heading),
         wind_speed = sqrt(u_wind^2 + v_wind^2),
         abs_cross_wind = abs(cross_wind(u = u_wind, v = v_wind, heading = heading))) %>% 
  mutate_at(c("step_length", "dem_100", "TRI_100", "slope_TPI_100", "t2m", "wind_support", "wind_speed", "cross_wind", "abs_cross_wind", "weeks_since_emig"), list(z = ~(scale(.))))  #calculate the z scores

#save data
saveRDS(data, file = "all_inds_annotated_apr23.rds")

#correlations
data %>% 
  dplyr::select(c("step_length", "dem_100_z", "TRI_100_z", "slope_TPI_100_z", "t2m_z", "wind_support_z", "wind_speed_z", "cross_wind_z", "abs_cross_wind_z", "weeks_since_emig_z")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.6) #wind speed and absolute cross wind are correlated (r = 0.665)! nothing else!!

# STEP 1: ssf modeling per ind per week ----------------------------------------------------------------

inds_ls <- split(data, data$individual.local.identifier) #there are 55 individuals. no filter was run on the n of weeks available per ind. it doesnt' matter :p

#for one model per week, no need to include temperature because there is no seasonality at that scale.

form1 <- used ~ dem_100_z + TRI_100_z + wind_support_z + abs_cross_wind_z + 
  strata(stratum)

form1 <- used ~ dem_100_z  + TRI_100_z + slope_TPI_100_z +
  strata(stratum)


coeffs_ls <- lapply(inds_ls, function(ind){
  
  all_wks <- lapply(split(ind, ind$weeks_since_emig), function(wk){
    
    ssf <- clogit(form1, data = wk)
    
    summary(ssf)$coefficients %>%  
      as.data.frame() %>%
      rownames_to_column( "variable") %>% 
      mutate(individual.local.identifier = unique(wk$individual.local.identifier),
             weeks_since_emig = unique(wk$weeks_since_emig))
    
  }) %>% 
    purrr::reduce(rbind)
  
  #run a linear model for the changes in each coefficient over time. (per ind)
  #for dem
  dem_coeffs <- all_wks %>% 
    filter(variable == "dem_100_z")
  
  dem_lm <- lm(coef ~ weeks_since_emig, data = all_wks %>% filter(variable == "dem_100_z"))
  TRI_lm <- lm(coef ~ weeks_since_emig, data = all_wks %>% filter(variable == "TRI_100_z"))
  slope_TPI_lm <- lm(coef ~ weeks_since_emig, data = all_wks %>% filter(variable == "slope_TPI_100_z"))
  
  
  
})



# formula3 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + dem_100_z * cos(turning_angle) * weeks_since_emig_z +
#   TRI_100_z * step_length_z * weeks_since_emig_z + TRI_100_z * cos(turning_angle) * weeks_since_emig_z +
#   slope_TPI_100_z * step_length_z * weeks_since_emig_z + slope_TPI_100_z * cos(turning_angle) * weeks_since_emig_z +
#   strata(stratum)
# 
# ssf3 <- clogit(formula3, data = data)
# summary(ssf3)
# plot_summs(ssf3)
# modelsummary(ssf3)

formula4 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100_z * step_length_z * weeks_since_emig_z + 
  slope_TPI_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf4 <- clogit(formula4, data = data)
summary(ssf4)
plot_summs(ssf4)



#---------------

form1a <- used ~ dem_100_z * weeks_since_emig_z * t2m_z + 
  TRI_100_z * weeks_since_emig_z + 
  strata(stratum)
ssf <- clogit(form1a, data = data)
summary(ssf)
plot_summs(ssf)

form2 <- used ~ dem_100_z + t2m_z + slope_TPI_100_z +
  TRI_100_z + weeks_since_emig_z + wind_support_z + cross_wind_z + abs_cross_wind_z +
  strata(stratum)

ssf2 <- clogit(form2, data = data)
summary(ssf2)
plot_summs(ssf2)

form3 <- used ~ dem_100_z + slope_TPI_100_z + TRI_100_z + weeks_since_emig_z * wind_support_z + abs_cross_wind_z +
  strata(stratum)

ssf3 <- clogit(form3, data = data)
summary(ssf3)
plot_summs(ssf3)


form4 <- used ~ dem_100_z * weeks_since_emig_z + slope_TPI_100_z * weeks_since_emig_z + TRI_100_z * weeks_since_emig_z + weeks_since_emig_z * abs_cross_wind_z +
  strata(stratum)

ssf4 <- clogit(form4, data = data)
summary(ssf4)
plot_summs(ssf4)

