#script for constructing the energy landscape for golden eagles
#trial with clogit. INLA predictions are all close to 0
#follows from 03_energy_landscape_modeling_method1.R and 04_INLA_output_plots.R
#02.03.2023. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)
library(magrittr)
library(corrr)
library(jtools) #to plot coeffs of clogit
library(fields)
library(survival)
library(terra)
library(gstat) #for interpolations
#library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(modelsummary) #to get AIC for the clogit models
library(hrbrthemes)
library(ggnewscale)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open data. data was prepared in 03_energy_landscape_modeling_method1.R
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds")

#define colors
clr <- oce::oceColorsPalette(100)[9] #was [2] before
clr_light <- oce::oceColorsPalette(100)[10]

# STEP 1: ssf modeling ----------------------------------------------------------------

f <- used ~ TRI_100_z * step_length_z * weeks_since_emig_z + 
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf <- clogit(f, data = data)

#save to use on Raven
saveRDS(ssf, file = "ssf_model_for_raven.rds")

summary(ssf)
plot_summs(ssf)
modelsummary(ssf) #AIC = 347,895.0

#proper plot
X11(width = 7, height = 2) 
p_coeffs <- plot_summs(ssf, colors = clr, point.size = 1.5, point.alpha = .5, point.shape = 19) +
  scale_y_discrete(labels = rev(c("TRI", "Step length", "Distance to ridge", "TRI: Step length", "TRI: Week",
                                  "Step length: Week", "Step length: Distance to ridge", "Distance to ridge: Week",
                                  "TRI: Step length: Week", "Distance to ridge: Step length: Week"))) +
  labs(x = "Estimate", y = "") +
  xlim(-.71, .25) +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

#save the plot
pdf("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/figs/clogit_coeffs.pdf", 
    width = 7, height = 2)
p_coeffs
dev.off()


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
  
  #X11(width = 6.9, height = 3.5)
  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = probs)) +
   scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                        na.value = "white", name = "Intensity of use")+
    #guides(fill = guide_colourbar(title.position = "top")) +
    guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
    labs(x = "Week since dispersal", y = label) + #add label for GRC plot
    theme_minimal() +
    theme(plot.margin = margin(0, 15, 0, 0, "pt"),
          legend.direction="horizontal",
          legend.position = "bottom",
          legend.key.width=unit(.7,"cm"),
          legend.key.height=unit(.25,"cm"),
          text = element_text(size = 8), #font size should be between 6-8
          axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
          axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  
  assign(paste0(i, "_p"), pred_p)
}

#plot all interaction plots together
X11(width = 7, height = 3)
combined <- TRI_100_p + ridge_100_p + step_length_p & theme(legend.position = "bottom")
combined + plot_layout( guides = "collect")

pdf("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/figs/clogit_interactions.pdf", 
     width = 7, height = 3)
combined + plot_layout( guides = "collect")
dev.off()

# STEP 3: Alpine predictions ----------------------------------------------------------------

#open data 
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") #this is limited to 3 yrs post dispersal
ssf <- readRDS("ssf_model_for_raven.rds")

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

saveRDS(topo_df, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/topo_df_100_LF.rds")

rm(ridge_100, dem_100, TRI_100, topo)
gc();gc()


#In the previous version, I used the same stratum for all predictions. In this round, take a random startum per map.
# set.seed(7)
# stratum_for_pred <- sample(data$stratum, 1)

#in the previous version, I used the median step length for each week to make the predictions. Now, take 10 random step length values per week. 

wks_ls <- split(data, data$weeks_since_emig)
wks_ls2 <- wks_ls[42:156]

(start <- Sys.time())
areas <- lapply(wks_ls2, function(x){
  
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
     week_i <- one_week %>% 
      str_pad(3,"left","0")
    
     #generate a new dataset
    new_data <- topo_df %>%
      mutate(step_length = mean(x$step_length), #try again with the average step length
             #step_length = sample(x$step_length, nrow(topo_df), replace = T),
             #step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'),
             step_length_z = 0,
             weeks_since_emig = one_week,
             weeks_since_emig_z =  (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
             stratum = sample(x$stratum, nrow(topo_df), replace = T))
    
    #predict using the model
    preds <- predict(ssf, newdata = new_data, type = "risk")
    preds_prob <- preds/(preds+1)
    
    preds_pr <- new_data %>%
      mutate(preds = preds,
             probs = preds_prob)
    
    #skip saving this file. it's 1 GB. just get out of it what I need
    saveRDS(preds_pr, paste0("/home/enourani/Documents/alpine_preds/alpine_preds_wk_", week_i, ".rds"))
  
      
      #calculate suitable areas
      area_.7 <- preds_pr %>%
        filter(probs >= 0.7) %>% 
        summarize(pixels = n()) %>% #count the 
        mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
               area_km2 = round(area_m2/1e6,3),
               week_since_dispersal = week_i)
      
       area_.7
       
})

Sys.time() - start # 3.366299 hours for 42-156


## add areas of weeks 1 - 41 to the dataframe
files <- list.files("/home/enourani/Documents/alpine_preds", full.names = T, pattern = ".rds")[1:41]

areas_1_41 <- data.frame()

for (i in files){
  preds_pr <- readRDS(i)
  
  week_i <- str_sub(i, -7, -5)
  
  area_.7 <- preds_pr %>%
    filter(probs >= 0.7) %>% 
    summarize(pixels = n()) %>% #count the 
    mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3),
           week_since_dispersal = week_i)
  
  areas_1_41 <- rbind(areas_1_41, area_.7)
  
}


areas_df <- areas %>% 
  reduce(rbind) %>% 
  bind_rows(areas_1_41,.) %>% 
  mutate(week_since_dispersal = as.numeric(week_since_dispersal))

saveRDS(areas_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/suitable_areas_meansl_wkly.rds")



#plot the trends in the suitable areas over time

areas_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/suitable_areas_meansl_wkly.rds")

X11(width = 3.4, height = 3.4) 
p <- ggplot(areas_df, aes(x = week_since_dispersal, y = area_km2)) +
  geom_smooth(method = "lm", alpha = .1, level = .95, color = clr, fill = clr_light, lwd = 1) + #95% standard error
  geom_point(size = 1.5,  alpha = .5, color = clr, fill = clr) +
  labs(x = "Weeks since dispersal",
       y = bquote("Flyable area " (km^2))) +
  theme_minimal() +
  theme(plot.margin = margin(0, 5, 0, 0, "pt"),
        text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

pdf("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/figs/area_over_wks_meansl.pdf", 
    width = 3.4, height = 3.4)
p
dev.off()

cor.test(areas_df$week_since_dispersal, areas_df$area_km2)

#calculate the percentage increase in suitable areas
#calc area of the Alpine region

#open alpine df
alps <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/topo_df_100_LF.rds")

alps_area <- alps %>%
  dplyr::select(ridge_100) %>% 
  drop_na(ridge_100) %>% 
  summarize(pixels = n()) %>% #count the 
  mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
         area_km2 = round(area_m2/1e6,3))

#190,544.7 km2

#compare with the flyable ares
areas_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/suitable_areas_meansl_wkly.rds")

#add a ratio column 
areas <- areas_df %>% 
  mutate(alp_ratio = area_km2/alps_area$area_km2)

#percentage change over time
year_1 <- (areas %>% filter(week_since_dispersal == 52) %>%  pull(area_km2)) / (alps_area %>% slice(1) %>% pull(area_km2))
year_3 <- (areas %>% filter(week_since_dispersal == 156) %>%  pull(area_km2)) / (alps_area %>% slice(1) %>% pull(area_km2))

year_1_3 <- (areas %>% slice(1) %>% pull(area_km2)) / (areas %>% filter(week_since_dispersal == 156) %>%  pull(area_km2)) 



#68% increase


# STEP 4: Alpine maps ----------------------------------------------------------------

pred_files <- list.files("/home/enourani/Documents/alpine_preds", full.names = T, pattern = ".rds")
  
#put dist to ridge as the base layer
ridge_100 <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif") %>% 
  as.data.frame(xy = T)
ridge_100[ridge_100$distance_to_ridge_line_mask > 5000, "distance_to_ridge_line_mask"] <- NA #do this for the plot to look nicer... NA values will be white

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/figs/alpine_preds/")

(b <- Sys.time())

for(wk in pred_files[79:156]){
#lapply(pred_files[4:6], function(wk){
  
  week_i <- str_sub(wk, -7, -5)
  
  r_.7 <- readRDS(wk) %>% 
    filter(probs >= .7) %>% 
    dplyr::select(c("location.long", "location.lat", "probs"))
  
  p <- ggplot() +
    geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = r_.7, aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2), guide = "none") +
    labs(title = paste0("Week ", as.numeric(week_i), " since dispersal"), x = "", y = "") +
    theme_void()
  
  
  ggsave(plot = p, filename = paste0(week_i, "_alpine_pred.png"), width = 7, height = 5, dpi = 300)


  print(paste0("week ", week_i, " done"))
#})
}

Sys.time() - b

#save a raster
#preds_r <- preds_pr %>%  #only keep the rows that contain data for this interaction term
#  dplyr::select(c("location.long", "location.lat", "probs")) %>% 
#  terra::rast(type = "xyz")

#make the maps
X11(width = 13, height = 9)

r_.7 %>%
  ggplot() +
  geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
  scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                       na.value = "white", name = "Intensity of use") +
  labs(x = "", y = "") +
  theme_void()


X11(width = 13, height = 9); ggplot(r_.7, aes(x = location.long, y = location.lat)) + geom_density2d()

ggplot(r_.7, aes(x = location.long, y = location.lat)) +
  stat_density_2d(aes(color = stat(level)))+
  scale_color_gradientn(colours = oce::oceColorsPalette(100))

ggplot(r_.7, aes(x = location.long, y = location.lat)) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon")+
  scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[80:100], alpha = .3)) #set transparency

ggplot(r_.7, aes(x = location.long, y = location.lat)) +
  stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE) +
  scale_fill_gradientn(colours = oce::oceColorsPalette(100), na.value = "white") +
  labs(x = "", y = "") +
  theme_void()

ggplot(r_.7, aes(x = location.long, y = location.lat)) + geom_bin2d(bins = 100)



#this one is good!
X11(width = 13, height = 9); ggplot() +
geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
       scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
  new_scale_fill() +
  geom_density2d(data = r_.7, aes(x = location.long, y = location.lat, fill = stat(density)))
  #stat_density_2d(r_.7, aes(x = location.long, y = location.lat, fill = stat(density)), geom = "polygon", contour = FALSE) +
  #scale_fill_gradientn(colours = oce::oceColorsPalette(100), na.value = "white") +
  labs(x = "", y = "") +
  theme_void()
  
  
  X11(width = 13, height = 9); ggplot() +
    geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = r_.7, aes(x = location.long, y = location.lat, fill = stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2), guide = "none") +
    theme_void()

#create animation of the maps. run the following code in the terminal
#ffmpeg -framerate 25 -pattern_type glob -i "*.png" output.mp4