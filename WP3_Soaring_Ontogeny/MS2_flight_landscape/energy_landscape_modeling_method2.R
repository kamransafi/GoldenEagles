#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#This version tries the method of fitting one model per individual. using the iSSA framework (selection-free movement kernel + habitat selection kernel)
#Apr 11. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot) #to plot coeffs of INLA
library(jtools) #to plot coeffs of clogit

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open annotated data (static variables and time since fledging and emigration)
load("alt_50_20_min_25_ind_static_time_ann.RData") #cmpl_ann

cmpl_ann <- cmpl_ann %>% 
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
cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% #round up days since emigration. we don't want zeros
  ggplot(aes(x = individual.local.identifier, y = max_day)) +
  geom_col()

cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% 
  summarize(mode = Mode(max_day), # 537
            mean = mean(max_day)) #444


#terrain ~ days since emigration ..... no patterns in the plots

plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "dem_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "slope_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "aspect_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "slope_TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "TRI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "aspect_TPI_100")])


ggplot(cmpl_ann[cmpl_ann$used == 1,], aes(as.numeric(days_since_emig), aspect_TPI_100)) +
  geom_point() +
  stat_smooth(aes(group = individual.local.identifier), method = "lm") +
  theme_minimal() +
  theme(legend.position = "none")


# STEP 2: summary plots ----------------------------------------------------------------

#one set of boxplots for select days: 10, 50 , 100, 500

data_int <- cmpl_ann %>%
  filter(days_since_emig_n %in% c(1, 10, 30 , 100, 300)) %>% 
  mutate(days_f = as.factor(days_since_emig_n))


variables <- c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100" )

X11(width = 6, height = 9)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/box_plots.png", width = 6, height = 9, units = "in", res = 300)

par(mfrow= c(4,2), 
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

cmpl_ann %>% 
  dplyr::select(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "days_since_emig_n")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #slope and TRI are correlated (.98)

# STEP 4: standardize variables ----------------------------------------------------------------

all_data <- cmpl_ann %>% 
  mutate_at(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "days_since_emig_n"),
            list(z = ~(scale(.)))) 

# STEP 5: ssf modeling ----------------------------------------------------------------

#Apr13 note: prediction rasters don't need to have the same origin and resolution. That can be fixed using terra::resample ftn.
#model formula for all individuals will be the same
#model with only habitat selection 
formula <- used ~ dem_100_z * days_since_emig_n_z + 
  slope_TPI_100_z * days_since_emig_n_z + 
  aspect_TPI_100_z * days_since_emig_n_z  + 
  TRI_100_z * days_since_emig_n_z + 
  TPI_100_z * days_since_emig_n_z + 
  strata(stratum)

n <- 200 #define n for new data generation

#name the variables
y_axis_var <- c("dem_100_z", "slope_TPI_100_z", "aspect_TPI_100_z", "TRI_100_z", "TPI_100_z")
x_axis_var <- "days_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the lapply call
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

#for all individuals
lapply(split(all_data, all_data$individual.local.identifier), function(ind){
  #the model
  ssf <- clogit(formula, data = ind)
  #save model output
  saveRDS(ssf, file = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/RDS_files/", 
                          ind$individual.local.identifier[1], "_model.rds"))
  #save coeff plot
  png(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/RDS_files/", 
             ind$individual.local.identifier[1], "_coeffs.png"), width = 7, height = 5, units = "in", res = 300)
  plot_summs(ssf)
  dev.off()
  
  #generate new data
  new_data <- ind %>%
    group_by(stratum) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = F) %>% 
    mutate(used = NA,  #regular intervals for all variables, so we can make a raster later on
           days_since_emig_n = sample(seq(min(ind$days_since_emig_n),max(ind$days_since_emig_n), length.out = 20), n, replace = T), 
           days_since_emig_n_z = sample(seq(min(ind$days_since_emig_n_z),max(ind$days_since_emig_n_z), length.out = 20), n, replace = T), #more levels for time than the other variables
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
  
  lapply(y_axis_var, function(var){
    
    #extract scale and center values needed to back transform the scaled values. do this using the whole dataset.
    var_scale <- attr(all_data[,colnames(all_data) == var],'scaled:scale')
    var_center <- attr(all_data[,colnames(all_data) == var],'scaled:center')
    
    #summarize values, so each (x,y) combo has one probability value
    avg_pred <- preds_pr %>% 
      group_by_at(c(var,x_axis_var)) %>%  #group by days since emigration and var
      summarise(avg_pres = mean(probs)) %>% 
      ungroup() %>% 
      mutate(dplyr::select(.,all_of(var)) * var_scale + var_center, #back-transform the values for plotting
             dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
      #rename(x_backtr = 1, #days since fledging
      #       i_backtr = 2) %>%  #y axis variables
      as.data.frame()
    
    #create a raster
    coordinates(avg_pred) <- c(x_axis_var, var)
    gridded(avg_pred) <- TRUE
    r <- raster(avg_pred)
    
    #save rater
    saveRDS(r, file = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/RDS_files/", 
                          ind$individual.local.identifier[1], "_", var, "_pred.rds"))
    
    
  })
  
})

# STEP 6: merge interaction plots ----------------------------------------------------------------

#all rasters have the same name, so use this function to be able to assign new name to them

lapply(y_axis_var, function(x){
  rasters <- list.files("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ind_model_outputs_Apr13/RDS_files/", pattern = x, full.names = TRUE) %>% 
    map(readRDS)
  
  #decide on which raster to resample all other ones to

  
  #extract resolutions
  cell_size <- lapply(rasters, res) %>% 
    reduce(rbind) %>% 
    as.data.frame() %>% 
    rename(days = 1, var = 2) %>% 
    summarize(avg = mean(var),
              md = Mode(var)) %>% 
    pull(md)
  
  #extract min value for y
  days_origin <- lapply(rasters, function(x) extent(x)[1]) %>% 
    reduce(rbind) %>% 
    as.data.frame() %>% 
    summarize(min(V1)) %>% 
    pull()
  
  #resample rasters to the same origin and resolution
  lapply(rasters, function(x){
    r <- as(x,"SpatRaster") %>% 
      resample()
    })
  
  
  
  
  
  
})


#TO DO:
#why do I have negative values for days in the prediction rasters???.....


