#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R (and temp_download_&prep.R)
#data preparation was done in 03_energy_landscape_modeling.R. This scripts builds one ssf per week using static variables, instead of adding week_since_emig as an interaction with the terrain variables
#Feb 22. 2022. Elham Nourani. Konstanz, DE

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
library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
#setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# STEP 0: open data ----------------------------------------------------------------
#open prepared data. use 100m res throughout 
all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_temp_inlaready_wmissing_wks_n1500.rds") %>% 
  filter(TRI_100 < quantile(TRI_100,.99)) %>% #remove TRI outliers 
  mutate(TRI_100_z = scale(TRI_100)) #recalculate the z score
  
#open alpine tri_dem data. From 05_INLA_prediction_map.R
#alps_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alps_topo_100m_temp_df.rds")


# STEP 1: ssf modeling: one per week ----------------------------------------------------------------

#investigate the amount of data per week
wk_ids <- all_data %>% 
  group_by(weeks_since_emig_n) %>% 
  summarize(n = n()) %>% 
  filter(n > 1000) %>% 
  pull(weeks_since_emig_n) #121 weeks remain

#define variables and formula
formulaM <- used ~ -1 +
  dem_100_z + TRI_100_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, TRI_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

mean.beta <- 0
prec.beta <- 1e-4 

lapply(wk_ids, function(wk){
  
  #----------------------------------------------------- STEP 1: filter for week number wk
  data <- all_data %>% 
    filter(weeks_since_emig_n == wk)
  
  #----------------------------------------------------- STEP 2: model! 
  (b <- Sys.time())
  
  M_main <- inla(formulaM, family = "Poisson",
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = data,
                 num.threads = 20, #if making alpine preds, set to 1 to have enough memory allocated to the thread
                 control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                 control.compute = list(openmp.strategy =  "pardiso", config = TRUE, cpo = T), #consider deactivating cpo to save computing power
                 control.inla(strategy = "adaptive", int.strategy = "eb"),
                 inla.mode="experimental", verbose = F)
  
  Sys.time() - b # 3 seconds
  
  #quick coeffs plot: Efxplot(M_main)
  
  #----------------------------------------------------- STEP 3: extract model evaluation
  
  
  #----------------------------------------------------- STEP 4: extract model coefficients
  graph <- as.data.frame(summary(M_main)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph$Factor <- rownames(graph)
  
  saveRDS(graph, file = paste0("graph_M_main100_wk_", wk, ".rds"))

  
  #----------------------------------------------------- STEP 5: extract individual variation
  
  summ_rndm <- M_main$summary.random
  saveRDS(summ_rndm, file = paste0("rnd_coeff_M_main100_wk_", wk, ".rds"))
  
  
  #----------------------------------------------------- STEP 6: clean up
  
  rm(M_main, data, graph, summ_rndm)
  
  gc(gc())
  
})


### script started for running the models and alpine preds simultaneously. but then decided to only run the model on training dataset as above.
# lapply(wk_ids, function(wk){
#   
#   #----------------------------------------------------- STEP 1: filter for week number wk
#   data <- all_data %>% 
#     filter(weeks_since_emig_n == wk)
#   
#   #----------------------------------------------------- STEP 2: add new data for alpine predictions
#   set.seed(400)
#   
#   #n needs to be large enough to cover the whole range of 
#   n <- nrow(alps_df) 
#   
#   new_data <- data %>%
#     group_by(stratum) %>% 
#     slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
#     ungroup() %>% 
#     slice_sample(n = n, replace = T) %>% 
#     dplyr::select(-c("dem_100", "TRI_100", "location.long", "location.lat")) %>%  #remove the columns that will be replaced by the alpine data. location.lat and location.long will only be added for the alpine rows, and not the original
#     bind_cols(alps_df) %>%  
#     mutate(used = NA,
#            dem_100_z = (dem_100 - mean(all_data$dem_100))/sd(all_data$dem_100), #convert these to z-scores based on the mean and variance of the tracking data.
#            TRI_100_z = (TRI_100 - mean(all_data$TRI_100))/sd(all_data$TRI_100)) %>% 
#     bind_rows(data %>%  mutate(location.lat = NA, location.long = NA)) %>% #location.lat and location.long will only be added for the alpine rows, and not the original.
#     as.data.frame(data)
#   
#   #----------------------------------------------------- STEP 3: model! 
#   (b <- Sys.time())
#   
#   M_pred <- inla(formulaM, family = "Poisson",
#                  control.fixed = list(
#                    mean = mean.beta,
#                    prec = list(default = prec.beta)),
#                  data = new_data,
#                  num.threads = 1, #to have enough memory allocated to the thread, set to 1
#                  control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
#                  control.compute = list(openmp.strategy =  "pardiso", config = TRUE, cpo = T), #consider deactivating cpo to save computing power
#                  control.inla(strategy = "adaptive", int.strategy = "eb"),
#                  inla.mode="experimental", verbose = F)
#   
#   Sys.time() - b # 3 seconds
#   
#   #quick coeffs plot: Efxplot(M_main)
#   
#   #----------------------------------------------------- STEP 4: extract model coefficients
#   graph <- as.data.frame(summary(M_pred)$fixed)
#   colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
#   colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
#   colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
#   
#   graph <- graph %>%
#     mutate(Factor = rownames(graph),
#            week = wk)
#   
#   #----------------------------------------------------- STEP 5: extract predictions
#   used_na <- which(is.na(new_data$used))
#   
#   #extract information for rows that had NAs as response variables. append to the original new_data (with NAs as used)
#   preds <- new_data %>% 
#     filter(is.na(used)) %>% 
#     dplyr::select(c("location.long", "location.lat", "dem_100", "TRI_100", "dem_100_z", "TRI_100_z", "weeks_since_emig_n")) %>% 
#     mutate(preds = M_pred$summary.fitted.values[used_na,"mean"],
#            preds_sd = M_pred$summary.fitted.values[used_na,"sd"]) %>% 
#     mutate(prob_pres = exp(preds)/(1+exp(preds)))
#   
#   #----------------------------------------------------- STEP 5: extract individual variation
#   
#   summ_rndm <- M_pred$summary.random
#   
#   
# })