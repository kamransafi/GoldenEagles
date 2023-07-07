#script for making Alpine prediction maps for golden eagles
#21.06.2023. Canberra, AU.
#Elham Nourani, PhD.

library(dplyr)
library(ggplot2)
library(glmmTMB)
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(ggnewscale)


#STEP 1: open data. data was prepared in 03_energy_landscape_modeling_method1.R.
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  mutate(animal_ID = as.numeric(as.factor(ind1)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum)))

#read in the ssf model
TMB_M <- readRDS( "TMB_model.rds")

#Alpine df prepared in 03_04_clogit_workflow.R
topo_df <- readRDS("topo_df_100_LF.rds")

#prepare dist to ridge to be used as the base layer for plotting
ridge_100 <- rast("ridge_100_LF.tif") %>% 
  as.data.frame(xy = T) %>% 
  drop_na(distance_to_ridge_line_mask)
ridge_100[ridge_100$distance_to_ridge_line_mask > 5000, "distance_to_ridge_line_mask"] <- NA #do this for the plot to look nicer... NA values will be white


#STEP 2: make predictions for each week and save the output

wks_ls <- split(data, data$weeks_since_emig)

#I will run this four times for each batch of weeks
for(x in wks_ls[c(1:39)]){
  #for(x in wks_ls[c(40:79)]){
  #for(x in wks_ls[c(80:119)]){
  #for(x in wks_ls[c(120:156)]){
  
  #extract week number
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
  #create week ID to be used for naming files and plots
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  #generate a new dataset
  topo_df <- topo_df %>%
    mutate(step_length = mean(x$step_length),
           step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'),
           weeks_since_emig = one_week,
           weeks_since_emig_z = (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
           stratum_ID = NA,
           animal_ID = sample(x$animal_ID, nrow(topo_df), replace = T)) #set the grouping variables to NA as per glmm.TMB recommendation
  
  #predict using the model
  topo_df$preds <- predict(TMB_M, newdata = topo_df, type = "link")
  
  topo_df$probs <- gtools::inv.logit(topo_df$preds)
  
  #calculate suitable areas
  area_.6 <- topo_df %>% 
    filter(probs >= .6) %>% 
    summarize(pixels = n()) %>% #count the 
    mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3),
           week_since_dispersal = week_i)
  
  #save area as a file
  save(area_.6, file = paste0("area_alps_wk_", week_i, ".r"))
  
  
  #plot the raw map
  t <- topo_df %>% 
    ggplot() +
    geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
    scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                         na.value = "white", name = "Intensity of use") +
    labs(x = "", y = "", title = paste0("Week ", one_week, " since dispersal")) +
    theme_void()
  
  ggsave(plot = t, filename = paste0("alps_wk_", week_i, "_raw.tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400)
  
  rm(t)
  
  
  #plot the density map
  p <- ggplot() +  
    geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = topo_df %>% filter(probs >= .6) %>% dplyr::select("location.lat", "location.long", "probs"), 
                    aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2)) +
    labs(title = paste0("Week ", one_week, " since dispersal"), x = "", y = "") +
    theme_void()
  
  ggsave(plot = p, filename = paste0("alps_wk_", week_i, "_density.tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400) #3min
  
  rm(p)
  
}



