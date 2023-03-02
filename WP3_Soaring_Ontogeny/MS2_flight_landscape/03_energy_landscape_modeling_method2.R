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
all_data <- readRDS("alt_50_20_min_70_ind_static_time_ann_dailytemp.rds") %>% 
  filter(TRI_100 < quantile(TRI_100,.99)) %>% #remove TRI outliers 
  mutate_at(c("dem_100", "slope_100", "aspect_100", "slope_TPI_100", "aspect_TPI_100", "TRI_100", "t2m"),
            list(z = ~(scale(.)))) %>%  #calculate the z scores
  mutate(ind1 = factor(individual.local.identifier),
         ind2 = factor(individual.local.identifier),
         ind3 = factor(individual.local.identifier),
         days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig)))
  
  
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
  dem_100_z * t2m_z + TRI_100_z * t2m_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, TRI_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

mean.beta <- 0
prec.beta <- 1e-4 

graph_ls <- lapply(wk_ids, function(wk){
  
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
                 num.threads = 8, #if making alpine preds, set to 1 to have enough memory allocated to the thread
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
  graph$week <- wk
  graph
  #saveRDS(graph, file = paste0("graph_M_main100_wk_", wk, ".rds"))

  
  #----------------------------------------------------- STEP 5: extract individual variation
  
  #summ_rndm <- M_main$summary.random
  #saveRDS(summ_rndm, file = paste0("rnd_coeff_M_main100_wk_", wk, ".rds"))
  
  
  #----------------------------------------------------- STEP 6: clean up
  
  #rm(M_main, data, graph, summ_rndm)

  
})



