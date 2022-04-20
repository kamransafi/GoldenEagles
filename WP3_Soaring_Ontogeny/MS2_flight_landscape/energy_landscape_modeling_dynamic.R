#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R and weather_annotation.R
#the first attempt only included static variables: energy_landscape_modeling_method2.R.
#Apr 19. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot) #to plot coeffs of INLA
library(jtools) #to plot coeffs of clogit
library(terra)

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open annotated data (static variables and time since fledging and emigration)
dyn_ann <- readRDS("alt_50_60_min_25_ind_dyn_time_ann.RData")

dyn_ann <- dyn_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)), #round up
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_")) #forgot to include stratum id in the previous code ) %>%  


Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# STEP 1: data exploration ----------------------------------------------------------------

#number of days available for each individual
dyn_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% #round up days since emigration. we don't want zeros
  ggplot(aes(x = individual.local.identifier, y = max_day)) +
  geom_col()

dyn_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% 
  summarize(mode = Mode(max_day), # 633
            mean = mean(max_day)) #431


# no patterns in the plots

plot(dyn_ann[dyn_ann$used == 1, c("days_since_emig", "wind_speed")])
plot(dyn_ann[dyn_ann$used == 1, c("days_since_emig", "wind_support")])
plot(dyn_ann[dyn_ann$used == 1, c("days_since_emig", "abs_cross_wind")])


ggplot(dyn_ann[dyn_ann$used == 1,], aes(as.numeric(days_since_emig), wind_speed)) +
  geom_point() +
  stat_smooth(aes(group = individual.local.identifier), method = "lm") +
  theme_minimal() +
  theme(legend.position = "none")


# STEP 2: summary plots ----------------------------------------------------------------

#one set of boxplots for select days: 10, 50 , 100, 500

data_int <- dyn_ann %>%
  filter(days_since_emig_n %in% c(1, 10, 30 , 100, 300, 500)) %>% 
  mutate(days_f = as.factor(days_since_emig_n))


variables <- c("wind_support", "wind_speed", "abs_cross_wind", "cross_wind")

X11(width = 6, height = 7)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/dyn_box_plots.png", width = 6, height = 7, units = "in", res = 300)

par(mfrow= c(2,2), 
    oma = c(0,0,3,0), 
    las = 1,
    mgp = c(0,1,0))
for(i in 1:length(variables)){
  boxplot(data_int[,variables[i]] ~ data_int[,"days_f"], data = data_int, boxfill = NA, border = NA, main = variables[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(data_int[data_int$used == 1, variables[i]] ~ data_int[data_int$used == 1,"days_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9),  lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1, "days_f"])) - 0.15)
  boxplot(data_int[data_int$used == 0, variables[i]] ~ data_int[data_int$used == 0, "days_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1 , "days_f"])) + 0.15)
  
}
dev.off()


#also consider making the plots for periods of time and not only one day

# STEP 3: check for collinearity ----------------------------------------------------------------

dyn_ann %>% 
  dplyr::select(c("wind_support", "wind_speed", "abs_cross_wind", "cross_wind", "days_since_emig_n", "step_length")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #abs crosswind and wind speed are correlated (.67)

# STEP 4: standardize variables ----------------------------------------------------------------

all_data <- dyn_ann %>% 
  mutate_at(c("wind_support", "wind_speed", "abs_cross_wind", "cross_wind", "days_since_emig_n"),
            list(z = ~(scale(.)))) 

# STEP 5: ssf modeling ----------------------------------------------------------------

#run one model for all data to see general patterns
formula <- used ~ wind_support_z * days_since_emig_n_z + 
  wind_speed_z * days_since_emig_n_z +
  strata(stratum)

ssf_all <- clogit(formula, data = all_data)

formula2 <- used ~ wind_support_z * days_since_emig_n_z + 
  abs_cross_wind_z * days_since_emig_n_z +
  strata(stratum) 

ssf_all2 <- clogit(formula2, data = all_data) #both interaction terms are significant :P

#include the movement kernel
formula3 <- used ~ log(step_length) * days_since_emig_n_z + #the lowest AIC of all the other combos... scale the step-length?
  #cos(turning_angle) * days_since_emig_n_z +
  wind_support_z * days_since_emig_n_z + 
  abs_cross_wind_z * days_since_emig_n_z +
  strata(stratum) 

ssf_all3 <- clogit(formula3, data = all_data)

plot_summs(list(ssf_all2, ssf_all3))


#include the movement kernel
formula4 <- used ~ log(step_length) * days_since_emig_n_z * wind_support_z + #the lowest AIC of all the other combos... scale the step-length?
  strata(stratum) 

ssf_all4 <- clogit(formula4, data = all_data)

plot_summs(list(ssf_all2, ssf_all3))




n <- 200 #define n for new data generation

#name the variables
y_axis_var <- c("wind_support_z", "abs_cross_wind_z")
x_axis_var <- "days_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the lapply call
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

#for all individuals
lapply(split(all_data, all_data$individual.local.identifier), function(ind){
  #the model
  ssf <- clogit(formula2, data = ind)
  #save model output
  saveRDS(ssf, file = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_dyn_model_outputs_Apr19/", 
                          ind$individual.local.identifier[1], "_dynmodel.rds"))
  #save coeff plot
  png(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_dyn_model_outputs_Apr19/", 
             ind$individual.local.identifier[1], "_coeffs.png"), width = 7, height = 5, units = "in", res = 300)
  plot_summs(ssf)
  dev.off()
  
  #generate new data
  new_data <- ind %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = T) %>% #i need replacement to be true, because this dataset doesn't have > 200 unique strata 
    mutate(used = NA,  #regular intervals for all variables, so we can make a raster later on
           days_since_emig_n = sample(seq(min(ind$days_since_emig_n),max(ind$days_since_emig_n), length.out = 20), n, replace = T), 
           days_since_emig_n_z = sample(seq(min(ind$days_since_emig_n_z),max(ind$days_since_emig_n_z), length.out = 20), n, replace = T), #more levels for time than the other variables
           wind_support_z = sample(seq(min(ind$wind_support_z),max(ind$wind_support_z), length.out = 20), n, replace = T),
           abs_cross_wind = sample(seq(min(ind$abs_cross_wind),max(abs_cross_wind), length.out = 20), n, replace = T)) 
  
  #predict using the model
  preds <- predict(ssf, newdata = new_data, type = "risk")
  preds_pr <- new_data %>% 
    mutate(preds = complete.cases(preds)) %>% 
    rowwise %>% 
    mutate(probs = preds/(preds+1))
  
  #lapply(y_axis_var, function(y_var){
    
    #extract scale and center values needed to back transform the scaled values. do this using the whole dataset.
    var_scale <- attr(all_data[,colnames(all_data) == y_var],'scaled:scale')
    var_center <- attr(all_data[,colnames(all_data) == y_var],'scaled:center')
    
    #summarize values, so each (x,y) combo has one probability value
    avg_pred <- preds_pr %>% 
      group_by_at(c(y_var,x_axis_var)) %>%  #group by days since emigration and y_var
      summarise(avg_pres = mean(probs, na.rm = T)) %>% 
      ungroup() %>% 
      mutate(dplyr::select(.,all_of(y_var)) * var_scale + var_center, #back-transform the values for plotting
             dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
      #rename(x_backtr = 1, #days since fledging
      #       i_backtr = 2) %>%  #y axis variables
      as.data.frame()
    
    #create a raster
    coordinates(avg_pred) <- c(x_axis_var, y_var)
    gridded(avg_pred) <- TRUE
    r <- raster(avg_pred) #many NA values....
    
    #save rater
    saveRDS(r, file = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_dyn_model_outputs_Apr19/", 
                          ind$individual.local.identifier[1], "_", y_var, "_pred.rds"))
    
    
  #})
  
})

# STEP 6: merge interaction plots ----------------------------------------------------------------

#all rasters have the same name, so use this function to be able to assign new name to them

y_vars <- c("dem_100", "slope_TPI_100", "aspect_TPI_100", "TRI_100", "TPI_100")

lapply(y_vars, function(x){
  rasters <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/RDS_files/", pattern = x, full.names = TRUE) %>% 
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
  
  saveRDS(interpr, file = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/prediction_raster_avg_Apr19/", 
                             x, "_avg_pred.rds"))
  
  proj4string(interpr) <- wgs
  
  #save raster plot
  png(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/prediction_raster_avg_Apr19/", 
             x, "_avg_pred.png"), width = 7, height = 5, units = "in", res = 300)
  plot(interpr, main = x)
  dev.off()
  
})




#TO DO:
#why do I have negative values for days in the prediction rasters???..... the actual intensity of use values are NAs


