#code for using the SSF built in energy_landscape_modeling_method1.R to make flight suitability maps for regular peirods of time for the entire Alpine region
#May 16. 2022. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(terra)
library(sf)
library(mapview)

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")


#I need to make predictions with the INLA model for the entire Alpine region. i.e. I need to append enough rows with NA values to the dataset 
# 1) only make predictions for unique combos of dem and tri, then associate the prediction to the rest of the cells
#with the same values. 3) run on the cluster, maybe in batches
#AND, let's not forget that I need to do this for different timestamps (weeks)!! :p

# STEP 0: prep topo layers. do 200 m instead of 100

dem <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/dem_100.tif")
TRI <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/TRI_100.tif")

dem_200 <- terra::aggregate(dem, fact = 2, fun = "mean", filename = "/home/enourani/Desktop/golden_eagle_static_layers/whole_region/dem_200.tif")
tri_200 <- terra::aggregate(TRI, fact = 2, fun = "mean", filename = "/home/enourani/Desktop/golden_eagle_static_layers/whole_region/tri_200.tif")

#open Apline perimeter layer
Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(tri_200)) %>% 
  as("SpatVector")

stck <- c(dem_200,tri_200) %>% 
  mask(Alps)

names(stck) <- c("dem_200", "tri_200")

# STEP 1: prep new data ----------------------------------------------------------------

#unique dem-tri combos: 4,755,351 lol
alps_topo_df <- as.data.frame(stck, xy = T) 

saveRDS(alps_topo_df, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_df.rds")

#alps_unique <- alps_topo_df %>% 
#  group_by(dem_200,tri_200) %>% 
#  slice(n = 1)

alps_unique <- alps_topo_df %>% 
  distinct(dem_200, tri_200)

#how many unique combinations do we have
n_unique <- nrow(alps_unique)

#open eagle data
all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_inlaready_wks.rds") %>% 
  select(c("location.long", "location.lat", "track", "stratum", "step_length", "turning_angle", "used", "dem_100", "TRI_100", "ind1", "ind2", "weeks_since_emig_n", "weeks_since_emig_n_z"))

# STEP 2: create new datasets: one for each interval ----------------------------------------------------------------

#do every week for the first month, then every 6 months and 2 yrs

set.seed(500)
n <- n_unique #unique combo of tri and dem values 

#prep data for week one, then (on the cluster) loop over weeks and just change the week number and run the model 
#although the new terrain info for the alps is in 200 m resolution, don't change the column name in the input data, because
#it will mess up the model and predictions!!
#(/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/alps_preds/order_of_business_alps)

alps_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = NA, 
         dem_100 = alps_unique$dem_200, #add these right here, so I won't need to back-transform later
         TRI_100 = alps_unique$tri_200,
         dem_100_z = (alps_unique$dem_200 - mean(all_data$dem_100))/sd(all_data$dem_100), #convert these to z-scores based on the mean and variance of the tracking data.
         TRI_100_z = (alps_unique$tri_200 - mean(all_data$TRI_100))/sd(all_data$TRI_100)) %>% 
  full_join(all_data)

saveRDS(alps_data, file = "inla_preds_for_cluster/generic_alt_50_20_min_48_ind_wmissing.rds")


#also prepare a vector with the sd and mean of the weeks since fledging. to be used ltr for the modeling

weeks_since_z_info <- data.frame(mean_wks = mean(all_data$weeks_since_emig_n), #center
                                 sd_wks = sd(all_data$weeks_since_emig_n)) #scale

saveRDS(weeks_since_z_info, file = "inla_preds_for_cluster/weeks_since_z_info.rds")

#run the weekly predictions on the cluster

# STEP 3: make sure model coefficients match the original model's ----------------------------------------------------------------

#list coefficient files to make sure they match the original model. check! :D
graph_files <- list.files("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results_alps/GE_ALPS/", pattern = "graph", full.names = T)

#(write a ftn?) plot the coefficients
lapply(graph_files, function(wk){
  
  graph <- readRDS(wk)
  graph <- graph[graph$Factor != "weeks_since_emig_n_z",]
  if(class(graph$Factor) != "character"){droplevels(graph$Factor)}
  VarOrder <- rev(unique(graph$Factor))
  VarNames <- VarOrder
  
  graph$Factor <- factor(graph$Factor, levels = VarOrder)
  levels(graph$Factor) <- VarNames
  
  #min <- min(graph$Lower,na.rm = T)
  #max <- max(graph$Upper,na.rm = T)
  
  graph$Factor_n <- as.numeric(graph$Factor)
  
  #plot in ggplot2
  #X11(width = 4.7, height = 2.7)
  
  coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray", size = 0.5) +
    geom_point(color = "cornflowerblue", size = 2)  +
    xlim(-0.1,0.5) +
    scale_y_discrete(name = "",
                     labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM")) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", size = 1) +
    theme_classic()
  
  
  ggsave(plot = coefs, 
         filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/weekly_alps_preds_jul26/inla_coeffs_48ind_", str_sub(wk,-7,-5),".png"), 
         width = 4.7, height = 2.7, dpi = 300) #they are all the same. so should be OK
  
})


# STEP 4: create the prediction maps: one for each week ----------------------------------------------------------------

#was thinking to also include a map of sd for each week, but I didnt save the sd info in the prediction results. if planning to
#repeat for more weeks, save the sd in the results.
#figure out what to do with values of TRI and TPI that are outside the range of the training data
#option 1: create a new raster layer with these regions. overlay this with the prediction raster with a gray fill


#open original data
all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_inlaready_wks.rds")

#list prediction files
pred_files <- list.files("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results_alps/GE_ALPS/", 
                         pattern = "preds", full.names = T)

#open the alpine region dataframe.
alps_topo_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_df.rds")

tri_200 <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/tri_200.tif")

#open Apline perimeter layer
Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(tri_200)) %>% 
  as("SpatVector")


#add one column per week, and cbind the prediction values with the original alps_topo_df based on TRI and TPI
lapply(pred_files, function(wk){
  
  preds <- readRDS(wk) %>% 
    mutate(wk_n = str_sub(wk,-7,-5),
           #backtransform the z-scores
           tri_200 = (TRI_100_z * sd(all_data$TRI_100)) + mean(all_data$TRI_100),   #round the values for easier matching with the original alps data
           dem_200 = (dem_100_z * sd(all_data$dem_100)) + mean(all_data$dem_100),
           weeks_since_emig_n = (weeks_since_emig_n_z * sd(all_data$weeks_since_emig_n)) + mean(all_data$weeks_since_emig_n)) #%>% #this will have one value. because each lapply round includes infor for one week
    #match the prediction values to the alps_df
    #full_join(alps_topo_df, by = c("dem_200", "tri_200"))

  alps_preds <- alps_topo_df %>% 
  left_join(preds, by = c("dem_200", "tri_200"))
  
  
})


#for plotting, create a raster stack, with one layer per week.



