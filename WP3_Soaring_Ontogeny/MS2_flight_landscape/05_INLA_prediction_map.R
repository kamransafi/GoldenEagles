#code for using the SSF built in 03_energy_landscape_modeling_method1.R to make flight suitability maps for regular periods of time for the entire Alpine region
#May 16. 2022. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(terra)
library(sf)
library(mapview)

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

#I need to make predictions with the INLA model for the entire Alpine region. i.e. I need to append enough rows with NA values to the dataset 
#previously, I 1) only made the predictions for unique combos of dem and tri, then associate the prediction to the rest of the cells, and
# 2) used 200 m instead of 100 m DEM. But I don't think I need to do these, because I'm using the cluster anyway

# STEP 0A: prep topo layers ----------------------------------------------------------------

dem_100 <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/dem_100.tif")
TRI_100 <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/TRI_100.tif")

#open Apline perimeter layer
Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(TRI_100)) %>% 
  as("SpatVector")

Alps_buffer <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(TRI_100)) %>% 
  st_buffer(1000) %>% 
  as("SpatVector") #use buffered alps to make sure all temp values along the borders get transferred to the topo stack

stck_topo <- c(dem_100,TRI_100) %>% 
  mask(Alps_buffer)

names(stck_topo) <- c("dem_100", "TRI_100")

# STEP 0B: prep temp layers ----------------------------------------------------------------
#decade-long temp data based on 02b_temo_download_prep.R: use monthly temp for June and October
Jun_temp <- rast("avg_temp_10yr_Jun.tif")
Oct_temp <- rast("avg_temp_10yr_Oct.tif")

#temp and topo have different extents and resolutions. make separate stacks
stck_temp <- c(Jun_temp, Oct_temp) %>%
  project(Alps) %>% 
  resample(stck_topo, method="bilinear", threads = T) %>%  #match the resolution to stck_topo
  mask(Alps_buffer)

names(stck_temp) <- c("Jun_temp", "Oct_temp")

stck_all <- c(stck_temp, stck_topo) %>% 
  mask(Alps) #after putting everything in a stack, mask with the non-buffered Alps layer


# STEP 1: prep new data ----------------------------------------------------------------

alps_df <- as.data.frame(stck_all, xy = T) %>% 
  rename(location.long = x,
         location.lat = y)
#there are 33,000 NAs. just omit them. they are at the border of the alps
#nas_sf <- alps_df %>%  filter(is.na(Jun_temp)) %>% st_as_sf(coords = c("location.long","location.lat"), crs = wgs)

alps_df_no_na <- alps_df %>% 
  drop_na(Jun_temp) %>% 
  as.data.frame()
  
saveRDS(alps_df_no_na, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")



# STEP 2: create new dataset for predictions ----------------------------------------------------------------

#open alps topo and temp data
alps_df_no_na <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")

#open eagle data
all_data <- readRDS("alt_50_20_min_48_ind_static_temp_inlaready_wks.rds") %>% 
  select(c("location.long", "location.lat", "track", "stratum", "step_length", "turning_angle", "used", "ind1", "ind2", "weeks_since_emig_n", "weeks_since_emig_n_z", 
           "dem_100", "dem_100_z", "TRI_100", "TRI_100_z", "month_temp", "month_temp_z"))


set.seed(500)
n <- nrow(alps_df_no_na) 

#unique combo of tri and dem values 
n_unique <- alps_df_no_na %>% 
  distinct(dem_100, TRI_100) %>% 
  nrow()
  #group_by(dem_100, TRI_100) %>% 
  #slice(1) %>% 
  #ungroup() %>% 
  #nrow()
  
#prep data for week one, then (on the cluster) loop over weeks and just change the week number and run the model. one file per month

for(i in c("Jun_temp", "Oct_temp")){
  
  alps_data <- all_data %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = T) %>% 
    mutate(used = NA,
           weeks_since_emig_n = NA, 
           dem_100 = alps_df_no_na %>%  pull(dem_100), #add these right here, so I won't need to back-transform later
           TRI_100 = alps_df_no_na %>% pull(TRI_100),
           month_temp = alps_df_no_na %>%  pull(i), #include temperature for one of the two months
           dem_100_z = (alps_df_no_na$dem_100 - mean(all_data$dem_100))/sd(all_data$dem_100), #convert these to z-scores based on the mean and variance of the tracking data.
           TRI_100_z = (alps_df_no_na$TRI_100 - mean(all_data$TRI_100))/sd(all_data$TRI_100),
           month_temp_z = (alps_df_no_na %>%  pull(i) - mean(all_data$month_temp))/sd(all_data$month_temp)) %>% 
    bind_rows(all_data)
  
  saveRDS(alps_data, file = paste0("inla_preds_for_cluster/alps_alt_50_20_min_48_ind_wmissing_", i, ".rds")) 
  
}

#also prepare a vector with the sd and mean of the weeks since fledging. to be used ltr for the modeling
weeks_since_z_info <- data.frame(mean_wks = mean(all_data$weeks_since_emig_n), #center
                                 sd_wks = sd(all_data$weeks_since_emig_n)) #scale

saveRDS(weeks_since_z_info, file = "inla_preds_for_cluster/weeks_since_z_info.rds")

#run the weekly predictions on the cluster: cluster_prep/order_of_business_alps.txt

# STEP 3: make sure model coefficients match the original model's ----------------------------------------------------------------

#list coefficient files to make sure they match the original model. check! :D
graph_files <- list.files("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results_alps/alps_preds_Jan23/GE_ALPS/", pattern = "graph", full.names = T)

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
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "gray", size = 0.5) +
    geom_point(color = "cornflowerblue", size = 2)  +
    #xlim(-0.2,0.6) +
    scale_y_discrete(name = "",
                     labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "Montly temperature", "TRI", "DEM")) +
    geom_linerange(aes(xmin = Lower, xmax = Upper), color = "cornflowerblue", size = 1) +
    theme_classic()
  
  
  ggsave(plot = coefs, 
         filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/weekly_alps_preds_jan23/inla_coeffs_48ind_", str_sub(wk,-7,-5),".png"), 
         width = 4.7, height = 2.7, dpi = 300) #they are all the same. so should be OK
  
})


# STEP 4: create the prediction maps: one for each week ----------------------------------------------------------------

#was thinking to also include a map of sd for each week, but I didnt save the sd info in the prediction results. if planning to
#repeat for more weeks, save the sd in the results.
#figure out what to do with values of TRI and TPI that are outside the range of the training data
#option 1: create a new raster layer with these regions. overlay this with the prediction raster with a gray fill


#open original data
all_data <- readRDS("alt_50_20_min_48_ind_static_temp_inlaready_wks.rds") 
alps_df_no_na <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")
#list prediction files
 pred_files <- list.files("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results_alps/alps_preds_Jan23/GE_ALPS/", 
                          pattern = "preds", full.names = T)

tri_200 <- rast("/home/enourani/Desktop/golden_eagle_static_layers/whole_region/tri_200.tif")
 
#open Apline perimeter layer
 Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
   st_transform(crs(tri_200)) %>% 
   as("SpatVector")
# 
# #explore the variation in the predictions for each week
# lapply(pred_files, function(wk){
#   
#   preds <- readRDS(wk) %>% 
#     mutate(wk_n = str_sub(wk,-7,-5),
#            #backtransform the z-scores
#            tri_200 = (TRI_100_z * sd(all_data$TRI_100)) + mean(all_data$TRI_100),   #round the values for easier matching with the original alps data
#            dem_200 = (dem_100_z * sd(all_data$dem_100)) + mean(all_data$dem_100),
#            weeks_since_emig_n = (weeks_since_emig_n_z * sd(all_data$weeks_since_emig_n)) + mean(all_data$weeks_since_emig_n)) #%>% #this will have one value. because each lapply round includes infor for one week
#   
#   #append the predictions to the original alps topo data. 
#   alps_preds <- alps_topo_df %>% 
#     left_join(preds, by = c("dem_200", "tri_200"))
# 
#   png(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/weekly_alps_preds_dec2/histograms/wk_", str_sub(wk,-7,-5),".png"),
#       width = 4.7, height = 2.7, units = "in", res = 300)  
#   hist(alps_preds$prob_pres)
#   dev.off()
#   
#     })

 
#open new data
new_data <- readRDS("inla_preds_for_cluster/alps_alt_50_20_min_48_ind_wmissing_Jun_temp.rds")

used_na <- which(is.na(new_data$used)) #the NA rows are at the beginning of new data

#create one map per week
lapply(pred_files, function(wk){

  preds <- readRDS(wk) %>% 
    mutate(wk_n = str_sub(wk,-7,-5),
           #backtransform the z-scores
           TRI_100 = (TRI_100_z * sd(all_data$TRI_100)) + mean(all_data$TRI_100),  
           dem_100 = (dem_100_z * sd(all_data$dem_100)) + mean(all_data$dem_100),
           month_temp = (month_temp_z * sd(all_data$month_temp)) + mean(all_data$month_temp),
           weeks_since_emig_n = (weeks_since_emig_n_z * sd(all_data$weeks_since_emig_n)) + mean(all_data$weeks_since_emig_n)) #%>% #this will have one value. because each lapply round includes infor for one week
    #match the prediction values to the alps_df
    #full_join(alps_topo_df, by = c("dem_200", "tri_200"))

  #append the predictions to the original alps topo data. 
  alps_preds <- alps_df_no_na %>% 
  left_join(preds, by = c("dem_100", "TRI_100"))
  
  r <- rast(alps_preds[,c("x", "y", "prob_pres")], crs = crs(Alps)) #use the ggplot code at the end to plot it. reorder the columns
  
})


#for plotting, create a raster stack, with one layer per week.



