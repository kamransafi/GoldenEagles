#Script for step-selection analysis of golden eagle commuting flights as reported in Nourani et al. 2023
#this script contains the code for reproducing the analyses and the plots in the paper
#Elham Nourani, PhD. 25.07.2023
#enourani@ab.mpg.de

library(tidyverse)
library(lubridate)
library(corrr)
library(glmmTMB)
library(performance)
library(oce)
library(terra)
library(ggnewscale)
library(sf)
library(move2)
library(stars)

#setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/public_data/")


##### STEP 1: Open annotated data #####

data <- read.csv("SSF_input_data.csv") %>%  #this is the post-emigration data, including all the available and alternative steps
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

#define colors and variables to use later on
clr <- oce::oceColorsPalette(100)[9] 
clr2 <- oce::oceColorsPalette(100)[80]
crs_wgs84 <- st_crs(4326) 

##### STEP 2: Check for autocorrelation #####

data %>% 
  dplyr::select(c("step_length", "TRI_100", "ridge_100", "weeks_since_emig")) %>% 
  correlate() 

##### STEP 3: Z-transform the predictor variables #####

data <- data %>% 
  mutate_at(c("step_length", "TRI_100", "ridge_100", "weeks_since_emig"), list(z = ~(scale(.))))

##### STEP 4: Run the step selection model #####
#set up the model structure based on Muff et al: 
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

TMB_struc <- glmmTMB(used ~ -1 + TRI_100_z * step_length_z * weeks_since_emig_z + 
                       ridge_100_z * step_length_z * weeks_since_emig_z + (1|stratum_ID) + 
                       (0 + ridge_100_z | animal_ID) + 
                       (0 + TRI_100_z | animal_ID), 
                     family = poisson, data = data, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA, 1:2))), #2 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3), 0, 0))) #add a 0 for each random slope. in this case, 2


TMB_M <- glmmTMB:::fitTMB(TMB_struc) #the model is provided in the public data repository
summary(TMB_M)

#saveRDS(TMB_M, file = "TMB_model.rds")

##### STEP 5: Model validation #####
performance_rmse(TMB_M)

##### STEP 6: PLOT Fig.1 - coefficient estimates #####

graph <- confint(TMB_M) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!(Factor %in% c("weeks_since_emig_z", "Std.Dev.ridge_100_z|animal_ID", "Std.Dev.TRI_100_z|animal_ID.1")))

colnames(graph)[c(2,3)] <- c("Lower", "Upper") 

labels <- rev(c("TRI", "Step length", "Distance to ridge", "TRI: Step length", "TRI: Week",
                "Step length: Week", "Step length: Distance to ridge", "Distance to ridge: Week",
                "TRI: Step length: Week", "Distance to ridge: Step length: Week"))

VarOrder <- rev(unique(graph$Factor))
graph$Factor <- factor(graph$Factor, levels = VarOrder)

#plot
X11(width = 7, height = 2) 

ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(color = clr, size = 1.7)  +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = clr, linewidth = .7) +
  labs(x = "Estimate", y = "") +
  scale_y_discrete(labels = labels) +
  xlim(-.71, .25) +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t = 6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin = margin(r = 6)))

#ggsave(coefs, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/coeffs.pdf", 
#       width = 7, height = 2, dpi = 300)

##### STEP 7: PLOT Fig.S1 - individual-specific coefficients #####

rnd <- ranef(TMB_M) %>% 
  as.data.frame() %>% 
  filter(grpvar == "animal_ID") %>% 
  mutate(ID = rep(levels(as.factor(data$individual.local.identifier)),2),
         variable = c(rep("Distance_to_ridge", 55), rep("TRI", 55)),
         Lower = condval - condsd,
         Upper = condval + condsd)  #assign individuals' names

#order the observations based on distance to ridge line
order <- rnd %>% 
  filter(variable == "Distance_to_ridge") %>% 
  arrange(desc(condval)) %>% 
  pull(ID)

rnd$ID <- factor(rnd$ID, levels = order)

cols <- c(TRI = clr2,
          Distance_to_ridge = clr)

X11(width = 7, height = 9)
ggplot(rnd, aes(x = condval, y = ID, color = variable)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(size = 2, position = position_dodge(width = .7))  +
  geom_linerange(aes(xmin = Lower, xmax = Upper), size = 0.8, position = position_dodge(width = .7)) +
  scale_color_manual(values = cols) + 
  scale_y_discrete(labels = order) +
  labs(x = "Difference from fixed effect estimate", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

#ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/ind_coefs.pdf", 
#       width = 7, height = 9, dpi = 300)

##### STEP 8: PLOT Fig. 2 - interaction terms #####

## 8.1: create new data: make predictions for values that fall on a regular grid for optimal visualization -------------------------------

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. The max of the variables might be outliers, so use the 90% quantile instead
grd_tri <- expand.grid(x = (1:156),
                       y = seq(from = min(data$TRI_100, na.rm = T), to = quantile(data$TRI_100, .9, na.rm = T),  length.out = 15)) %>% 
  rename(weeks_since_emig = x,
         TRI_100 = y)  %>% 
  mutate(dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'), #set other variables to their mean
         ridge_100 = attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'),
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         interaction = "wk_TRI_100")

grd_dist <- expand.grid(x = (1:156),
                        y = seq(from = min(data$ridge_100, na.rm = T), to = quantile(data$ridge_100, .9, na.rm = T), length.out = 15)) %>%
  rename(weeks_since_emig = x,
         ridge_100 = y) %>% 
  mutate(TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'), #set other variables to their mean
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'),
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         interaction = "wk_ridge_100")

grd_sl <- expand.grid(x = (1:156),
                      y = seq(from = min(data$step_length, na.rm = T), to = quantile(data$step_length, .9, na.rm = T), length.out = 15)) %>%  # 226
  rename(weeks_since_emig = x,
         step_length = y) %>% 
  mutate(TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'), #set other variables to their mean
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'),
         ridge_100 = attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'),
         interaction = "wk_step_length")

#merge all together and assigned NA as used
grd_all <- bind_rows(grd_tri, grd_dist, grd_sl) %>% 
  mutate(used = NA)

set.seed(777)
n <- nrow(grd_all)

#complete the new_data dataframe with the other variables used for modeling. The resulting dataframe is provided in the public repository
new_data_only <- data %>%
  group_by(stratum_ID) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  #only keep the columns that I need
  dplyr::select(c("stratum_ID", "animal_ID")) %>% 
  bind_cols(grd_all) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation for consistency
  mutate(TRI_100_z = (TRI_100 - attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'))/attr(data[,colnames(data) == "TRI_100_z"],'scaled:scale'),
         ridge_100_z = (ridge_100 - attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'))/attr(data[,colnames(data) == "ridge_100_z"],'scaled:scale'),
         step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'), 
         weeks_since_emig_z = (weeks_since_emig - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'))

#write.csv(new_data_only, file = "new_data_for_pred.csv", row.names = F)

## 8.2: make predictions with the model -------------------------------

#read in the model file and the new data
TMB_M <- readRDS("TMB_model.rds")
new_data <- read.csv("new_data_for_pred.csv") %>% 
  mutate(stratum_ID = NA) ##accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
#(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."

preds <- predict(TMB_M, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #convert the log odds to a probability. read more: https://rpubs.com/crossxwill/logistic-poisson-prob

## 8.3: plot the predictions for interaction terms -------------------------------

#prepare for plotting
y_axis_var <- c("TRI_100", "ridge_100")
x_axis_var <- "weeks_since_emig"
labels <- data.frame(vars = y_axis_var, 
                     label = c("Topographic Ruggedness Index", "Distance to ridge line (km)"))

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
  if(i == "ridge_100"){
    pred_p <- pred_r %>% 
      ggplot() +
      geom_tile(aes(x = x, y = y, fill = probs)) +
      scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                           na.value = "white", name = "Energy Availability Index")+
      guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
      labs(x = "Week since emigration", y = label) + #add label for GRC plot
      theme_minimal() +
      theme(plot.margin = margin(0, 15, 0, 0, "pt"),
            legend.direction="horizontal",
            legend.position = "bottom",
            legend.key.width=unit(.7,"cm"),
            legend.key.height=unit(.25,"cm"),
            text = element_text(size = 8), #font size should be between 6-8
            axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
            axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  } else{
    pred_p <- pred_r %>% 
      ggplot() +
      geom_tile(aes(x = x, y = y, fill = probs)) +
      scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                           na.value = "white", name = "Energy Availability Index")+
      guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
      labs(x = "", y = label) + #add label for GRC plot
      theme_minimal() +
      theme(plot.margin = margin(0, 15, 0, 0, "pt"),
            legend.direction="horizontal",
            legend.position = "bottom",
            legend.key.width=unit(.7,"cm"),
            legend.key.height=unit(.25,"cm"),
            text = element_text(size = 8), #font size should be between 6-8
            axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
            axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  }
  assign(paste0(i, "_p"), pred_p)
  
}

#plot all interaction plots together
X11(width = 4.5, height = 4.8)
combined <- TRI_100_p + ridge_100_p & theme(legend.position = "bottom")
p <- combined + plot_layout(guides = "collect", nrow = 2)

##### STEP 9: Predictions for the Alps #####

# Open the Alpine dataframe. This file contains all combinations of TRI and distance to ridgeline for each lat-lon grid in the Alpine region 
topo_df <- readRDS("topo_df_100_LF.rds")

# Prepare dist to ridge to be used as the base layer for plotting 
ridge_100 <- rast("ridge_100_LF.tif") %>% 
  as.data.frame(xy = T) %>% 
  drop_na(distance_to_ridge_line_mask)
ridge_100[ridge_100$distance_to_ridge_line_mask > 5000, "distance_to_ridge_line_mask"] <- NA #do this for the plot to look nicer... NA values will be white

#here we will make one prediction for the Alpine region per week. For every week, we will save the raw predictions for supplementary video AND will save the area with a predicted value > 0.7
#For four weeks c(1,4,24,52), we will also save a map of the hotspots of energy availability (Fig. 3).
#this whole process is quite memory hungry. 
#to save the plots, create an output directory and call the path "output_path"

#output_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/alpine_preds_7/raw_preds/"

#read in model rds
TMB_M <- readRDS("TMB_model.rds")

wks_ls <- split(data, data$weeks_since_emig)[43:156] #make a list with one element per week.

#regularly clean up the memory
gc()

for(x in wks_ls){
  
  #extract week number
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
  #create week ID to be used for naming files and plots
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  #generate a new dataset for predictions
  topo_df <- topo_df %>%
    mutate(step_length = mean(x$step_length),
           step_length_z = 0,
           weeks_since_emig = one_week,
           weeks_since_emig_z = (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
           stratum_ID = NA, #set the grouping variables to NA as per glmm.TMB recommendation
           animal_ID = sample(x$animal_ID, nrow(topo_df), replace = T)) 
  
  #predict using the model
  topo_df$preds <- predict(TMB_M, newdata = topo_df, type = "link")
  
  topo_df$probs <- gtools::inv.logit(topo_df$preds) #convert the log odds to a probability
  
  #clean up the RAM
  gc()
  
  saveRDS(topo_df, file = paste0("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/weekly_preds_alps_sep23/",
                                 "alps_wk_", week_i, ".rds"))
  
  #calculate suitable areas
  # area_.7 <- topo_df %>%
  #   filter(probs >= .7) %>%
  #   summarize(pixels = n()) %>% #count the
  #   mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
  #          area_km2 = round(area_m2/1e6,3),
  #          week_since_emig = week_i)
  # 
  # #save area as a file
  # save(area_.7, file = paste0(output_path, "area_alps_wk_", week_i, ".r"))
  # 
  
  #plot the raw map
  # t <- topo_df %>%
  #   ggplot() +
  #   geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
  #   scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
  #                        na.value = "white", name = "Perceived Energy\nAvailability Index") +
  #   labs(x = "", y = "", title = paste0("Week ", one_week, " since emigration"))  +
  #   guides(colour=guide_legend(title.vjust = 3)) +
  #   theme(legend.text = element_text(size = 10),
  #         legend.title = element_text(size = 12)) +
  #   theme_bw()
  # 
  # ggsave(plot = t, filename = paste0(output_path, "alps_wk_", week_i, ".tiff"), device = "tiff", width = 8, height = 4.8, dpi = 400)
  # 
  # rm(t)
  gc()
  
  #density plots for four weeks
  # if(one_week %in% c(1, 4, 24, 52)){
  #   p <- ggplot() +
  #     geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
  #     scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
  #     new_scale_fill() +
  #     stat_density_2d(data = topo_df %>% filter(probs >= .7) %>% dplyr::select("location.lat", "location.long", "probs"),
  #                     aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
  #     scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2)) +
  #     labs(title = paste0("Week ", one_week, " since emigration"), x = "", y = "") +
  #     theme_void()
  # 
  #   #dev.off()
  #   ggsave(plot = p, filename = paste0(week_i, "_alpine_pred.tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400) #3min
  #   rm(p)
  #   gc()
  # }
  
  print(paste0("week ", week_i, " of 156 done!"))
  
}

##### STEP 10: Video S1 energy landscape animation #####
#to create an mp4 using the tiff files, use the following command in the terminal (make sure your working directory contains your tiff files):
#sudo apt install ffmpeg
#ffmpeg -framerate 15 -pattern_type glob -i "*.tiff" output.mp4

##### STEP 11: PLOT Fig. 4  Flyable areas #####

#list the files containing the areas with predicted value >=7 (created in step 9)
files <- list.files(output_path,
                    pattern = "area_alps_wk_", full.names = T)

#load all the files into one data frame
areas_df <- files %>%
  map_df(~ get(load(file = .x)))  %>%
  mutate(week_since_emig = as.numeric(week_since_emig))

X11(width = 3.4, height = 2) 
p <- ggplot(areas_df, aes(x = week_since_emig, y = area_km2)) +
  geom_point(size = 1.5,  alpha = .4, color = clr, fill = clr) +
  labs(x = "Weeks since emigration",
       y = bquote("Flyable area " (km^2))) +
  theme_minimal() +
  theme(plot.margin = margin(0, 5, 0, 0, "pt"),
        text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t = 6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin = margin(r = 6)))

#ggsave(p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/area_over_wks_meansl.pdf", 
#       width = 3.4, height = 2, dpi = 400)

#compare with the flyable areas wit the Alpine region
#calculate the area of the Alps: calculate the number of cells within the ridge line raster that has been masked using the perimeter of the Alpine
alps <- rast("ridge_100_LF.tif")

alps_area <- alps %>%
  as.data.frame() %>% 
  drop_na(distance_to_ridge_line_mask) %>% 
  summarize(pixels = n()) %>% #count the number of cells
  mutate(area_m2 = pixels * 100 * 100, #the resolution of the raster is 100 * 100 m
         area_km2 = round(area_m2/1e6,3)) #the area of the Alpine region is 190,544.7 km2

#add a column for the ratio of the Alpine region being classified as flyable 
areas <- areas_df %>% 
  mutate(alp_ratio = area_km2/alps_area$area_km2)

#calculate total area used by the birds per week
data_used <- data %>% 
  filter(used == 1) %>% 
  mutate(ind_step_id = paste(individual.local.identifier, burst_id, sep = "_")) #assign an ID to each burst of uninterrupted flight. to make sure that non-flight trajectories aren't added to the data

#overlap the tracking data with the alpine raster and count the overlapping cells.
#but first, convert the tracking points to trajectories
mv_w <- mt_as_move2(data_used %>% 
                      st_as_sf(coords = c("location.long", "location.lat"), crs = crs_wgs84),
                    time_column = "timestamp", track_id_column = "ind_step_id")

#match the projection to the Alpine raster
mv_w <- mv_w %>% sf::st_transform(crs = st_crs(alps))

#estimate used area per week. This should be cumulative area
mv_ls <- split(mv_w, mv_w$weeks_since_emig)

used_area_wk <- lapply(1:length(mv_ls), function(x){
  
  #weeks prior to this one + this
  w_lines <- mv_ls[1:x] %>% 
    reduce(rbind) %>% 
    mt_track_lines() #create a LINESTRING object
  
  #create raster object, convert to a dataframe and count the cells.
  r_lines <- w_lines %>% 
    st_rasterize(dx = 100, dy = 100) %>% 
    as.data.frame() %>% 
    drop_na(ID) %>% 
    summarize(pixels = n()) %>% #count n of rows
    mutate(weeks_since_emig = x,
           area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3)) # 9228.63 area_km2
  
  r_lines
  
}) %>% 
  reduce(rbind)



##### Session info #####

# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.2 LTS
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] stars_0.6-2        abind_1.4-5        move2_0.2.2        sf_1.0-14         
# [5] ggnewscale_0.4.9   terra_1.7-39       oce_1.8-1          gsw_1.1-1         
# [9] performance_0.10.4 corrr_0.4.4        lubridate_1.9.2    forcats_1.0.0     
# [13] stringr_1.5.0      dplyr_1.1.2        purrr_1.0.1        readr_2.1.4       
# [17] tidyr_1.3.0        tibble_3.2.1       ggplot2_3.4.2      tidyverse_2.0.0   
# 
# loaded via a namespace (and not attached):
# [1] rstudioapi_0.15.0  magrittr_2.0.3     vctrs_0.6.3        KernSmooth_2.23-22
# [5] TMB_1.9.5          lifecycle_1.0.3    pkgconfig_2.0.3    Matrix_1.6-0      
# [9] R6_2.5.1           colorspace_2.1-0   lwgeom_0.2-8       fansi_1.0.4       
# [13] timechange_0.2.0   compiler_4.3.1     proxy_0.4-26       bit64_4.0.5       
# [17] withr_2.5.0        DBI_1.1.3          MASS_7.3-60        classInt_0.4-9    
# [21] gtools_3.9.4       tools_4.3.1        units_0.8-0        glue_1.6.2        
# [25] nlme_3.1-162       grid_4.3.1         generics_0.1.3     gtable_0.3.3      
# [29] tzdb_0.4.0         class_7.3-22       hms_1.1.3          utf8_1.2.3        
# [33] pillar_1.9.0       vroom_1.6.3        splines_4.3.1      lattice_0.21-8    
# [37] bit_4.0.5          tidyselect_1.2.0   stringi_1.7.12     boot_1.3-28       
# [41] codetools_0.2-19   cli_3.6.1          munsell_0.5.0      Rcpp_1.0.11       
# [45] parallel_4.3.1     assertthat_0.2.1   scales_1.2.1       e1071_1.7-9       
# [49] insight_0.19.3     crayon_1.5.2       rlang_1.1.1 
