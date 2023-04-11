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



formula3 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + dem_100_z * cos(turning_angle) * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z + TRI_100_z * cos(turning_angle) * weeks_since_emig_z +
  strata(stratum)

ssf3 <- clogit(formula3, data = data)
summary(ssf3)
plot_summs(ssf3)

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


X11(width = 14, height = 12)
coefs <- plot_summs(ssf, point.size = 7, point.shape = 16) +
  labs(y = NULL) +
  labs(x = NULL) +
  scale_y_discrete(labels = c("Elevation:Weeks since dispersal:Temperature", "Weeks since dispersal:Ruggedness", "Week since dispersal: Temperature","Elevation:Temperature",
                               "Elevation:Weeks since dispersal","Ruggedness", "Temperature", "Elevation")) +
  theme_classic() +
  theme(axis.text.y = element_text(face="bold", size = 20),
                axis.text.x = element_text(face="bold", size = 20))
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_coeffs.png", 
       width = 14, height = 12, dpi = 500)


# STEP 2: predict using the ssf: terrain and week since emig interactions ----------------------------------------------------------------

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations
grd_dem <- expand.grid(x = seq(from = min(all_data$weeks_since_emig_n_z), to = max(all_data$weeks_since_emig_n_z), by = 0.1),
                   y = seq(from = min(all_data$dem_100_z), to = max(all_data$dem_100_z), by = 0.1))

grd_tri <- expand.grid(x = seq(from = min(all_data$weeks_since_emig_n_z), to = max(all_data$weeks_since_emig_n_z), by = 0.1),
                       y = seq(from = min(all_data$TRI_100_z), to = max(all_data$TRI_100_z), by = 0.1))

grd_temp <- expand.grid(x = seq(from = min(all_data$t2m_z), to = max(all_data$t2m_z), by = 0.1),
                       y = seq(from = min(all_data$dem_100_z), to = max(all_data$dem_100_z), by = 0.1))

n <- nrow(grd_dem) #the dem grid has more number of rows, so use it

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = grd_dem$x, 
         dem_100_z = grd_dem$y,
         TRI_100_z = sample(grd_tri$y, n, replace = T),
         t2m_z = 0) #for the terrain interaction plot, keep temp at mean level


#predict using the model
preds <- predict(ssf, newdata = new_data, type = "risk")
preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise %>% 
  mutate(probs = preds/(preds+1))

#prepare for plotting
y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

for (i in y_axis_var){
  
  i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
  i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')
  
  #summarize values, so each (x,y) combo has one probability value
  avg_pred <- preds_pr %>% 
    group_by_at(c(which(names(preds_pr) == i),which(names(preds_pr) == x_axis_var))) %>%  #group by weeks since emigration and i
    summarise(avg_pres = mean(probs)) %>% 
    ungroup() %>% 
    mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
           dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
    rename(x = which(names(.) == x_axis_var), #weeks since emig
           y = which(names(.) == i)) %>% #y axis variables
    dplyr::select(c("x","y","avg_pres")) %>%  #reorder the columns for making xyz raster
    as.data.frame()
  
  #saveRDS(avg_pred, file = paste0("inla_pred_clogit_not_smoothed_", i,".rds"))
  
  #comment out below to save the unsmoothed values
  r <- avg_pred %>%
    rast(type = "xyz") %>%
    focal(w = 7, fun = mean, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) %>%
    rename(avg_pres = focal_mean)

  saveRDS(r, file = paste0("inla_pred_clogit_", i,".rds"))

}

#create plots
p_dem <- readRDS("inla_pred_clogit_dem_100_z.rds") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "", y = "Elevation \n (m)") +
  theme_classic()

p_rugg <- readRDS("inla_pred_clogit_TRI_100_z.rds") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "Weeks since dispersal", y = "Terrain Ruggedness \n Index") +
  theme_classic()

#put both plots in one device
X11(width = 9, height = 4)
combined <- p_dem + p_rugg & theme(legend.position = "right")
(p_2  <- combined + plot_layout(guides = "collect", nrow = 2))

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_tr_interactions.png", 
       width = 9, height = 4, dpi = 500)

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_tr_interactions_raw.png", 
       width = 9, height = 4, dpi = 500)

# STEP 3: predict using the ssf: dem and temperature interaction ----------------------------------------------------------------

grd_temp <- expand.grid(x = seq(from = min(all_data$t2m_z), to = max(all_data$t2m_z), by = 0.1),
                        y = seq(from = min(all_data$dem_100_z), to = max(all_data$dem_100_z), by = 0.1))

n <- nrow(grd_temp) #the dem grid has more number of rows, so use it

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = 0, 
         dem_100_z = grd_temp$y,
         TRI_100_z = 0,
         t2m_z = grd_temp$x) 


#predict using the model
preds <- predict(ssf, newdata = new_data, type = "risk")
preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise %>% 
  mutate(probs = preds/(preds+1))

#prepare for plotting
y_axis_var <- "dem_100_z"
x_axis_var <- "t2m_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

i <- y_axis_var
  
  i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
  i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')
  
  #summarize values, so each (x,y) combo has one probability value
  avg_pred <- preds_pr %>% 
    group_by_at(c(which(names(preds_pr) == i),which(names(preds_pr) == x_axis_var))) %>%  #group by weeks since emigration and i
    summarise(avg_pres = mean(probs)) %>% 
    ungroup() %>% 
    mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
           dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
    rename(x = which(names(.) == x_axis_var), #weeks since emig
           y = which(names(.) == i)) %>% #y axis variables
    dplyr::select(c("x","y","avg_pres")) %>%  #reorder the columns for making xyz raster
    as.data.frame()
  
  #saveRDS(avg_pred, file = paste0("inla_pred_clogit_not_smoothed_", i,".rds"))
  
  #comment out below to save the unsmoothed values
  r <- avg_pred %>%
    rast(type = "xyz") %>%
    focal(w = 7, fun = mean, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) %>%
    rename(avg_pres = focal_mean)
  

#create plots
X11(width = 5, height = 5)

p <- r %>% 
  ggplot() +
  geom_tile(aes(x = x-273.15, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "Temperature (°C)", y = "Elevation (m)") +
  theme_classic()

ggsave(plot = p, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_temp_interaction.png", 
       width = 5, height = 5, dpi = 500)

################################## predictions for the alps ###########
all_data <- readRDS("alt_50_20_min_48_ind_static_100_daytemp_inlaready_wks.rds") %>% 
  dplyr::select(c("location.long", "location.lat", "track", "stratum", "step_length", "turning_angle", "used","weeks_since_emig_n", "weeks_since_emig_n_z", 
                  "dem_100", "dem_100_z", "TRI_100", "TRI_100_z", "t2m_z", "t2m"))

alps_df_no_na <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")



n <- nrow(alps_df_no_na)

#first to 4th week, after 6 months, after  a year, after 1.5 years
lapply(c(1:24,52,76), function(i){
  
  wk <- i %>% str_pad(2,"left","0")
  
  new_data <- all_data %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = T) %>% 
    dplyr::select(-c("dem_100", "TRI_100","location.long", "location.lat")) %>%  #remove the rows that will be replaced by the alpine data. location.lat and location.long will only be added for the alpine rows, and not the original
    bind_cols(alps_df_no_na) %>% 
    mutate(used = NA,
           weeks_since_emig_n = i,
           weeks_since_emig_n_z = (i - mean(all_data$weeks_since_emig_n))/sd(all_data$weeks_since_emig_n),
           dem_100_z = (dem_100 - mean(all_data$dem_100))/sd(all_data$dem_100), #convert these to z-scores based on the mean and variance of the tracking data.
           TRI_100_z = (TRI_100 - mean(all_data$TRI_100))/sd(all_data$TRI_100))
  
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





