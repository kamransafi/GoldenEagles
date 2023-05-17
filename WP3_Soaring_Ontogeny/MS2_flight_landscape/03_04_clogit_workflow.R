#script for constructing the energy landscape for golden eagles
#trial with clogit. INLA predictions are all close to 0
#follows from 03_energy_landscape_modeling_method1.R and 04_INLA_output_plots.R
#02.03.2023. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(magrittr)
library(corrr)
library(jtools) #to plot coeffs of clogit
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
library(modelsummary) #to get AIC for the clogit models
library(TwoStepCLogit)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")


setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open data. data was prepared in 03_energy_landscape_modeling_method1.R
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds")

# STEP 1: ssf modeling ----------------------------------------------------------------

f <-  used ~ TRI_100_z * step_length_z * weeks_since_emig_z + 
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf <- clogit(f, data = data)
summary(ssf)
plot_summs(ssf)
modelsummary(ssf) #AIC = 347,895.0

# STEP 2: predict using the ssf ----------------------------------------------------------------

#data created in 03_energy_landscape_modeling_method1.r
new_data <- readRDS("new_data_only_ssf_preds.rds")

#predict using the model
preds <- predict(ssf, newdata = new_data, type = "risk")
preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise %>% 
  mutate(probs = preds/(preds+1))

#prepare for plotting
y_axis_var <- c("TRI_100", "ridge_100", "step_length")
x_axis_var <- "weeks_since_emig"
labels <- data.frame(vars = y_axis_var, 
                     label = c("Terrain Ruggedness Index", "Distance to ridge (km)", "Step length (m)"))
#older versions of this script contained backtransforming the z-scores in a for loop.

for (i in y_axis_var){

  label <- labels %>% 
    filter(vars == i) %>% 
    pull(label)
  
  #interaction to be plotted
  interaction_term <- paste0("wk_", i)
  pred_r <- preds_pr %>% 
    filter(interaction == interaction_term) %>%  #only keep the rows that contain data for this interaction term
    dplyr::select(c(which(names(.) %in% c(x_axis_var, i)), "probs")) %>% 
    terra::rast(type = "xyz") %>%
    #focal(w = 3, fun = mean, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) #%>% 
    #rename(probs = focal_mean)
  
  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = probs)) +
   scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                        na.value = "white", name = "Intensity of use")+
    labs(x = "", y = label) +
    theme_classic()
  
  
  #save the plot
  ggsave(plot = pred_p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_", interaction_term,".png"), 
         width = 11, height = 3, dpi = 400)
}


# STEP 3: Alpine predictions ----------------------------------------------------------------

#open data 
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") #this is limited to 3 yrs post dispersal

#create a dataframe fro the Alps. layers made by Louise from 02_data_processing_&_annotation.R
#the dimensions and extents are different. so crop based on ridge to match the extents
ridge_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif")
TRI_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/TRI_100_LF.tif") %>% 
  crop(ridge_100)
dem_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/dem_100_LF.tif") %>% 
  crop(ridge_100)

ext(ridge_100) <- ext(dem_100)

topo <- c(dem_100, TRI_100, ridge_100)
names(topo) <- c("dem_100", "TRI_100", "ridge_100")

topo_df <- as.data.frame(topo, xy = T) %>% 
  rename(location.long = x,
         location.lat = y) %>% 
  mutate(dem_100_z = (dem_100 - attr(data[,colnames(data) == "dem_100_z"],'scaled:center'))/attr(data[,colnames(data) == "dem_100_z"],'scaled:scale'), #calculate z scores using the mean and sd of the modeling dataset
         TRI_100_z = (TRI_100 - attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'))/attr(data[,colnames(data) == "TRI_100_z"],'scaled:scale'),
         ridge_100_z = (ridge_100 - attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'))/attr(data[,colnames(data) == "ridge_100_z"],'scaled:scale'))

rm(ridge_100, dem_100, TRI_100, topo)
gc(gc())

#use the same stratum for all predictions
set.seed(7)
stratum_for_pred <- sample(data$stratum, 1)

#median step length for each week
step_length_wks <- data %>% 
  group_by(weeks_since_emig) %>% 
  reframe(step_length = median(step_length, na.rm = T)) #calculate median step length for each week

lapply(c(1:nrow(step_length_wks)), function(one_week){ #one pred for each week in the dataset (until year 3)
  
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  #extract the mean of step lengths within the nth week (for all inds)
  step_length <- step_length_wks %>% 
    filter(weeks_since_emig == one_week) %>% 
    pull(step_length)
  
  new_data <- topo_df %>% 
    mutate(step_length = step_length,
           step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'),
           weeks_since_emig = one_week,
           weeks_since_emig_z =  (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
           stratum = stratum_for_pred)
  
  #predict using the model
  preds <- predict(ssf, newdata = new_data, type = "risk")
  preds_prob <- preds/(preds+1)
  
  preds_pr <- new_data %>% 
    mutate(preds = preds,
           probs = preds_prob)
  
  saveRDS(preds_pr, paste0("/media/enourani/Ellham's HDD/Elham_GE/clogit_alpine_preds/alpine_preds_wk_", week_i, ".rds"))
  
  #progress report
  print(paste0("Week ", week_i, " done!"))
})



pred_ls <- list.files("/media/enourani/Ellham's HDD/Elham_GE/clogit_alpine_preds", 
                      pattern = ".rds", full.names = T)

#make the maps
lapply(pred_ls, function(x){
  
  week_i <- str_sub(x, 76, -5) %>% 
    str_pad(3,"left","0")
  
  png(paste0("/media/enourani/Ellham's HDD/Elham_GE/clogit_alpine_maps/alps_wk_", week_i, ".png"),
      width = 13, height = 9, units = "in", res = 400)
  print(readRDS(x) %>% 
          ggplot() +
          geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
          scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                               na.value = "white", name = "Intensity of use") +
          labs(x = "", y = "", title = paste0("Week ", week_i, " since dispersal")) +
          theme_void()
  )
  
  dev.off()
  
})

#calculate suitable areas
areas_.7 <- data.frame()
#list the prediction files
pred_ls <- list.files("/media/enourani/Ellham's HDD/Elham_GE/clogit_alpine_preds/", 
                      pattern = ".rds", full.names = T)

#areas_.7 <- lapply(pred_ls, function(x){
  
for(x in pred_ls){
  
  week_i <- str_sub(x, 76, -5) %>% 
    str_pad(3,"left","0")
  
  #create a raster layer
  area_.7 <- readRDS(x) %>%
    filter(probs >= 0.7) %>% 
    summarize(pixels = n()) %>% #count the 
    mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3),
           week_since_dispersal = week_i)
  
  areas_.7 <- rbind(areas_.7, area_.7)
}

saveRDS(areas_.7, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/suitable_areas_wkly.rds")

#}) %>% 
#  reduce(rbind)



#plot the trends in the suitable areas over time

plot(as.numeric(areas_.7$week_since_dispersal), areas_.7$area_km2)
