#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#This version tries the method of fitting one model per week using the iSSA framework (selection-free movement kernel + habitat selection kernel)
#update Apr 20: there are upto 80 weeks since emigration (mean 40). not all inds have data for all weeks. So, making weekly models is not a good idea. So, repeat method2,
#but with weeks since instead of days since. ALSO, try weeks since fledging.... ALSO, include month as a proxy for weather conditions
#update May 5: add landform categories to the models... or try....
#Apr 20. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot) #to plot coeffs of INLA
library(jtools) #to plot coeffs of clogit
library(terra)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open annotated data (static variables and days & weeks since fledging and emigration)
load("alt_50_20_min_25_ind_static_time_ann_weeks.RData") #cmpl_ann

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig)),#round up
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_")) #forgot to include stratum id in the previous code ) %>%  

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# STEP 1: data exploration ----------------------------------------------------------------

weeks_per_ind <- cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(n_weeks = n_distinct(weeks_since_emig_n)) #the median is 40, the 3rd quart is 60

#how many points and strata are there per week per indiiividual?
n_per_week <- cmpl_ann %>% 
  group_by(individual.local.identifier,weeks_since_emig_n) %>% 
  summarize(n_pts = n(),
            n_str = n_distinct(stratum))

#for example, how many unique strata are there in the second week since emig? There is a small number and   0 for some individuals....how about including weeks since emig in method 2 instead of days since?
n_per_week %>% 
  filter(weeks_since_emig_n == 3) %>% 
  summarise(n_str = sum(n_str))

#terrain ~ days since emigration ..... no patterns in the plots

plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "dem_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "slope_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "aspect_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "slope_TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "TRI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "aspect_TPI_100")])


ggplot(cmpl_ann[cmpl_ann$used == 1,], aes(as.numeric(weeks_since_emig_n), dem_100)) +
  geom_point() +
  stat_smooth(aes(group = individual.local.identifier), method = "lm") +
  theme_minimal() +
  theme(legend.position = "none")


# STEP 2: summary plots ----------------------------------------------------------------

#one set of boxplots for a few weeks

data_int <- cmpl_ann %>%
  filter(weeks_since_emig_n %in% c(1, 5,10, 25, 30)) %>% 
  mutate(weeks_f = as.factor(weeks_since_emig_n))


variables <- c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100" )

X11(width = 6, height = 9)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/box_plots_weeks.png", width = 6, height = 9, units = "in", res = 300)

par(mfrow= c(4,2), 
    oma = c(0,0,3,0), 
    las = 1,
    mgp = c(0,1,0))
for(i in 1:length(variables)){
  boxplot(data_int[,variables[i]] ~ data_int[,"weeks_f"], data = data_int, boxfill = NA, border = NA, main = variables[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(data_int[data_int$used == 1, variables[i]] ~ data_int[data_int$used == 1,"weeks_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9),  lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1, "weeks_f"])) - 0.15)
  boxplot(data_int[data_int$used == 0, variables[i]] ~ data_int[data_int$used == 0, "weeks_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1 , "weeks_f"])) + 0.15)
  
}
dev.off()


#also consider making the plots for periods of time and not only one day..indeed :p

# STEP 3: check for collinearity ----------------------------------------------------------------

cmpl_ann %>% 
  dplyr::select(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "days_since_emig_n")) %>% 
  correlate() #slope and TRI are correlated (.98)

# STEP 4: standardize variables ----------------------------------------------------------------

all_data <- cmpl_ann %>% 
  mutate_at(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "weeks_since_emig_n", "eobs.temperature"),
            list(z = ~(scale(.)))) %>% 
  mutate(month = month(timestamp))

# STEP 5: ssf modeling ----------------------------------------------------------------

#Apr13 note: prediction rasters don't need to have the same origin and resolution. That can be fixed using terra::resample ftn.
#model formula for all individuals will be the same
#model with only habitat selection 
formula <- used ~ dem_100_z * weeks_since_emig_n_z + 
  slope_TPI_100_z * weeks_since_emig_n_z + 
  aspect_TPI_100_z * weeks_since_emig_n_z  + 
  TRI_100_z * weeks_since_emig_n_z + 
  TPI_100_z * weeks_since_emig_n_z + 
  strata(stratum)

#add month or temperature to control for env conditions. also remove aspect tpi
formula2 <- used ~ dem_100_z * weeks_since_emig_n_z + log(step_length) * eobs.temperature_z + 
  slope_TPI_100_z * weeks_since_emig_n_z + 
  TRI_100_z * weeks_since_emig_n_z + 
  TPI_100_z * weeks_since_emig_n_z + 
  strata(stratum)

ssf2 <- clogit(formula2, data = all_data)
plot_summs(ssf2)

n <- 200 #define n for new data generation

#name the variables
y_axis_var <- c("dem_100_z", "slope_TPI_100_z", "aspect_TPI_100_z", "TRI_100_z", "TPI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the lapply call
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

path <- paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs",
               strsplit(as.character(Sys.Date()), "22")[[1]][2])
dir.create(path)

#for all individuals
lapply(split(all_data, all_data$individual.local.identifier), function(ind){
  #the model
  ssf <- clogit(formula, data = ind)
  #save model output
  saveRDS(ssf, file = paste0(path, "/", ind$individual.local.identifier[1], "_model.rds"))
  
  #save coeff plot
  png(paste0(path, "/", ind$individual.local.identifier[1], "_coeffs.png"), width = 7, height = 5, units = "in", res = 300)
  plot_summs(ssf)
  dev.off()
  
  #generate new data
  new_data <- ind %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = F) %>% 
    mutate(used = NA,  #regular intervals for all variables, so we can make a raster later on
           weeks_since_emig_n = sample(seq(min(ind$weeks_since_emig_n),max(ind$weeks_since_emig_n), length.out = 20), n, replace = T), 
           weeks_since_emig_n_z = sample(seq(min(ind$weeks_since_emig_n_z),max(ind$weeks_since_emig_n_z), length.out = 20), n, replace = T), #more levels for time than the other variables
           dem_100_z = sample(seq(min(ind$dem_100_z),max(ind$dem_100_z), length.out = 20), n, replace = T),
           slope_TPI_100_z = sample(seq(min(ind$slope_TPI_100_z),max(ind$slope_TPI_100_z), length.out = 20), n, replace = T),
           aspect_TPI_100_z = sample(seq(min(ind$aspect_TPI_100_z),max(ind$aspect_TPI_100_z), length.out = 20), n, replace = T),
           TRI_100_z = sample(seq(min(ind$TRI_100_z),max(ind$TRI_100_z), length.out = 20), n, replace = T),
           TPI_100_z = sample(seq(min(ind$TPI_100_z),max(ind$TPI_100_z), length.out = 20), n, replace = T)) 
  
  #predict using the model
  preds <- predict(ssf, newdata = new_data, type = "risk")
  preds_pr <- new_data %>% 
    mutate(preds = preds) %>% 
    rowwise %>% 
    mutate(probs = preds/(preds+1))
  
  lapply(y_axis_var, function(y_var){
    
    #extract scale and center values needed to back transform the scaled values. do this using the whole dataset.
    var_scale <- attr(all_data[,colnames(all_data) == y_var],'scaled:scale')
    var_center <- attr(all_data[,colnames(all_data) == y_var],'scaled:center')
    
    #summarize values, so each (x,y) combo has one probability value
    avg_pred <- preds_pr %>% 
      group_by_at(c(y_var,x_axis_var)) %>%  #group by days since emigration and y_var
      summarise(avg_pres = mean(probs)) %>% 
      ungroup() %>% 
      mutate(dplyr::select(.,all_of(y_var)) * var_scale + var_center, #back-transform the values for plotting
             dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
      #rename(x_backtr = 1, #days since fledging
      #       i_backtr = 2) %>%  #y axis variables
      as.data.frame()
    
    #create a raster
    coordinates(avg_pred) <- c(x_axis_var, y_var)
    gridded(avg_pred) <- TRUE
    r <- raster(avg_pred)
    
    #save rater
    saveRDS(r, file = paste0(path, "/", ind$individual.local.identifier[1], "_", y_var, "_pred.rds"))
    
    #save as image
    proj4string( r) <- wgs
    
    #save raster plot
    png(paste0(path, "/", 
               y_var, "_pred.png"), width = 7, height = 5, units = "in", res = 300)
    plot(r, main = y_var)
    dev.off()
    
  })
  
})

# STEP 6: merge interaction plots ----------------------------------------------------------------

#all rasters have the same name, so use this function to be able to assign new name to them

y_vars <- c("dem_100", "slope_TPI_100", "aspect_TPI_100", "TRI_100", "TPI_100")

lapply(y_vars, function(x){
  rasters <- list.files(path, pattern = paste0(x ,".+rds"), full.names = TRUE) %>% #list only the rds files
    map(readRDS)
  
  #resample to wgs, just so they look nice and have a crs
  ls_r <- lapply(rasters, function(r){
    crs(r) <-  "+proj=longlat +datum=WGS84 +no_defs"
    return(as(r,"SpatRaster"))
  })
  
  #set the second raster as the base
  base <- ls_r[[2]]
  
  #resample rasters to the same origin and resolution
  ls_r_r <- lapply(ls_r, function(r){
    rs <- r %>% 
      resample(base, method = "near")
  })
  
  #average all into one raster layer
  avg_r <- ls_r_r %>% 
    reduce(c) %>% 
    app(fun = "mean", na.rm = T)
  

  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(avg_r, xy = T)[, c(1,2)], col = 2), as.data.frame(avg_r, xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = ext(avg_r)[1],to = ext(avg_r)[2],by = 2),
                     y = seq(from = ext(avg_r)[3],to = ext(avg_r)[4],by = 2))
  
  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords, surf.1.pred)
  
  colnames(interpdf) <- c("x_backtr","i_backtr","prob_pres")
  
  coordinates(interpdf) <- ~ x_backtr + i_backtr
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  
  saveRDS(interpr, file = paste0(path, "/avg_raster_preds/", x, "_avg_pred.rds"))
  
  proj4string(interpr) <- wgs
  
  #save raster plot
  png(paste0(path, "/avg_raster_preds/", x, "_avg_pred.png"), width = 7, height = 5, units = "in", res = 300)
  plot(interpr, xlim = c(1,60), main = x)
  dev.off()
  
})




#TO DO:
#why do I have negative values for days in the prediction rasters???..... the actual intensity of use values are NAs



# STEP 7: merge interaction plots ----------------------------------------------------------------

load("alt_50_20_min_25_ind_static_time_ann_weeks.RData") #cmpl_ann

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig)),#round up
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_")) #forgot to include stratum id in the previous code ) %>%  

all_data <- cmpl_ann %>% 
  mutate_at(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "weeks_since_emig_n", "eobs.temperature"),
            list(z = ~(scale(.)))) %>% 
  mutate(month = month(timestamp))


#open the landform layer
ll <- rast("/home/enourani/Desktop/LandformClassification/meybeck05/w001001.adf") %>% 
  as.factor()
lf <- rast("/home/enourani/Desktop/LandformClassification/iwahashi2/w001001.adf") %>% 
  as.factor()

#open landform legend
meta_d <- read_excel("/home/enourani/Desktop/LandformClassification/legend/LEGEND.xlsx", sheet = "Meybeck", range = "A1:B16")
meta_i <- read_excel("/home/enourani/Desktop/LandformClassification/legend/LEGEND.xlsx", sheet = "Iwahashi", range = "A1:B17")

#extract landform values at tracking points
data_lf <- all_data %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(extract(x = ll, y = st_coordinates(.))) %>% #extract landforms: Meybeck
  full_join(meta_d, by = c("COUNT_" = "ID")) %>% 
  rename(landform_class_M = COUNT_,
         landform_legend_M = Legend) %>% 
  mutate(extract(x = lf, y = st_coordinates(.))) %>% #extract landforms: Iwahashi
  full_join(meta_i, by = c("COUNT_" = "ID")) %>% 
  rename(landform_class_I = COUNT_,
         landform_legend_I = Legend) 


#ssf analysis. use landform Iwahashi, there is a more even distribution of values across character levels
formula1 <- used ~ as.factor(landform_class_I) * weeks_since_emig_n_z +
  strata(stratum)

ssf1 <- clogit(formula1, data = data_lf)
plot_summs(ssf1)
