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


#try twostep clogit
f <-  used ~ TRI_100_z * step_length_z * weeks_since_emig_z + 
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum) + cluster(ind1)

ssf <- Ts.estim(f, data = data, random = ~ TRI_100_z + ridge_100_z, D = "UN") #throws an error because week is constant within strata.

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
y_axis_var <- c("dem_100", "TRI_100", "ridge_100", "step_length")
x_axis_var <- "weeks_since_emig"

#older versions of this script contained backtransforming the z-scores in a for loop.

for (i in y_axis_var){

  #interaction to be plotted
  interaction_term <- paste0("wk_", i)
  
  #create a raster and run a moving window to make the pattern easier to see.
  pred_r <- preds_pr %>% 
    filter(interaction == interaction_term) %>%  #only keep the rows that contain data for this interaction term
    dplyr::select(c(which(names(.) %in% c(x_axis_var, i)), "probs")) %>% 
    terra::rast(type = "xyz") %>%
    focal(w = 3, fun = median, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) %>% 
    rename(probs = focal_median)
  
  #plot
  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = probs)) +
    scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                         na.value = "white", name = "Intensity of use") +
    labs(x = "", y = i) +
    theme_classic()
  
  #save the plot
  ggsave(plot = pred_p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_", interaction_term,".png"), 
         width = 9, height = 4, dpi = 500)
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

#use the same stratum for all predictions
set.seed(7)
stratum_for_pred <- sample(data$stratum, 1)

#median step length for each week
step_length_wks <- data %>% 
  group_by(weeks_since_emig) %>% 
  reframe(step_length = median(step_length, na.rm = T)) #calculate median step length for each week

wkly_preds <- lapply(c(1:nrow(step_length_wks)), function(one_week){ #one pred for each week in the dataset (until year 3)
  
  week_i <- one_week %>% 
    str_pad(2,"left","0")
  
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
    
  })


#have a look at the distribution of probabilities
lapply(wkly_preds, function(x) summary(x$probs))

#list of pred files
#pred_ls <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/clogit_alpine_preds/", 
#                      pattern = ".rds", full.names = T)


pred_ls <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/clogit_alpine_preds/", 
                      pattern = ".rds", full.names = T)

#make the maps
lapply(pred_ls, function(x){
  
  week_i <- str_sub(x, 145, -5) %>% 
    str_pad(2,"left","0")
  
  png(paste0("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/weekly_alp_preds_clogit_May11/alps_wk_", week_i, ".png"),
      width = 13, height = 9, units = "in", res = 400)
  print(readRDS(x) %>% 
          ggplot() +
          geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
          scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                               na.value = "white", name = "Intensity of use") +
          labs(x = "", y = "", title = paste0("Week ", str_sub(x, 145, -5), " since dispersal")) +
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
    str_pad(2,"left","0")
  
  #create a raster layer
  area_.7 <- readRDS(x) %>%
    filter(probs >= 0.7) %>% 
    summarize(pixels = n()) %>% #count the 
    mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3),
           week_since_dispersal = week_i)
  
  areas_.7 <- rbind(areas_.7, area_.7)
}

saveRDS(areas_.7, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/suitable_areas_wks_1_71.rds")

#}) %>% 
#  reduce(rbind)



######### old code below:
################################## predictions for the alps ###########
data <- readRDS("alt_50_20_min_48_ind_static_100_daytemp_inlaready_wks.rds") %>% 
  dplyr::select(c("location.long", "location.lat", "track", "stratum", "step_length", "turning_angle", "used","weeks_since_emig_n", "weeks_since_emig_n_z", 
                  "dem_100", "dem_100_z", "TRI_100", "TRI_100_z", "t2m_z", "t2m"))

alps_df_no_na <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")



n <- nrow(alps_df_no_na)

#first to 4th week, after 6 months, after  a year, after 1.5 years
lapply(c(1:24,52,76), function(i){
  
  wk <- i %>% str_pad(2,"left","0")
  
  new_data <- data %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = T) %>% 
    dplyr::select(-c("dem_100", "TRI_100","location.long", "location.lat")) %>%  #remove the rows that will be replaced by the alpine data. location.lat and location.long will only be added for the alpine rows, and not the original
    bind_cols(alps_df_no_na) %>% 
    mutate(used = NA,
           weeks_since_emig_n = i,
           weeks_since_emig_n_z = (i - mean(data$weeks_since_emig_n))/sd(data$weeks_since_emig_n),
           dem_100_z = (dem_100 - mean(data$dem_100))/sd(data$dem_100), #convert these to z-scores based on the mean and variance of the tracking data.
           TRI_100_z = (TRI_100 - mean(data$TRI_100))/sd(data$TRI_100))
  
  #predict using the model
  preds <- predict(ssf, newdata = new_data, type = "risk")
  preds_prob <- preds/(preds+1)
  
  preds_pr <- new_data %>% 
    mutate(preds = preds,
           probs = preds_prob,
           week = wk)
  
  

  saveRDS(preds_pr, file = paste0("clogit_alp_pred_week_", wk, ".rds"))
  
  #r <- preds_pr %>% 
  #  dplyr::select(c("location.long", "location.lat", "probs")) %>% 
  #  rast(crs = crs(TRI_100))
  
})

file_ls <- list.files(pattern = "clogit_alp_pred", full.names = T) 

#save plots
lapply(file_ls, function(x){

  p <- readRDS(x) %>% 
  ggplot() +
  geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
   labs(x = "", y = "", title = paste0("week ", str_sub(x,24,-5), " since dispersal")) +
  theme_classic()
  
  ggsave(plot = p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_alps_wk_", str_sub(x, 24, -5), ".png"), 
         width = 9, height = 6, dpi = 300)

})

################################## trends in suitable areas ###########

#empty raster stack to store reclassified rasters
p_stack <- NULL

#empty list to store frequency of each suitability clss
freq_ls <- list()
  
#list the weekly prediction files
file_ls <- list.files(pattern = "clogit_alp_pred_week", full.names = T) 

#create matrix for reclassification. The first two columns are "from" "to" of the input values, and the third column "becomes" 
classes <- matrix(c(c(0, .3, .5, .7, .85), 
                  c(.3, .5, .7, .85, 1),
                  c(1, 2, 3, 4, 5)), ncol = 3, nrow = 5, byrow = F)


for (i in file_ls){ 
 
  r_p <- readRDS(x) %>%
    dplyr::select(c("location.long", "location.lat", "probs")) %>% 
    rast(type = "xyz", crs = crs(TRI_100)) %>% 
    classify(rcl = classes, include.lowest = T, right = T)  #classify the raster into suitability categories
  
  names(r_p) <- "suitability"
  
  #concatenate to the raster stack
  p_stack <- c(p_stack, r_p)
  
 #estimate frequency of each class
   freq_r <- r_p %>% 
     freq() #in the next update, add week number to this data.frame
 
   freq_ls[[length(freq_ls) + 1]] <- freq_r 
}

  
wk <- c(1:24,52,76) %>% str_pad(2,"left","0")
  
names(freq_ls) <- wk 
names(p_stack) <- wk

#prep data for plotting
data <- freq_ls %>% 
  bind_rows(.id = "week") %>% 
  rename(suitability = value)

#make a stacked bar plot of different suitability classes

ggplot(data, aes(fill = suitability, y = count, x = week)) + 
  geom_bar(position = "stack", stat = "identity")





