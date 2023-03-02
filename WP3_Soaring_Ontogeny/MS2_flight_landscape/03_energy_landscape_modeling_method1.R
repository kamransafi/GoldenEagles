#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R (and temp_download_&prep.R)
#only includes static variables. 
#try adding daily temp to see. seems like dem and temp are interacting!.. from 03_method2
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
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#open annotated data
all_data <- readRDS("alt_50_20_min_70_ind_static_time_ann_dailytemp.rds") %>% 
  filter(TRI_100 < quantile(TRI_100,.99)) %>% #remove TRI outliers 
  mutate(ind1 = factor(individual.local.identifier), # INLA formula using interaction terms for time and predictor variables. repeat for random effects
         ind2 = factor(individual.local.identifier),
         ind3 = factor(individual.local.identifier),
         days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig))) %>% 
  mutate_at(c("dem_100", "TRI_100", "t2m", "weeks_since_emig_n"), list(z = ~(scale(.))))  #calculate the z scores
  
saveRDS(all_data, file = "alt_50_20_min_48_ind_static_100_daytemp_inlaready_wks.rds") #this has the limit on TRI range

# STEP 2: ssf modeling ----------------------------------------------------------------

## mid_step: investigate using clogit
# control for monthly temperature....

form1a <- used ~ dem_100_z * weeks_since_emig_n_z *t2m_z + 
  TRI_100_z * weeks_since_emig_n_z + 
  strata(stratum)

ssf <- clogit(form1a, data = all_data)
summary(ssf)
plot_summs(ssf)

form1a <- used ~ dem_200_z * TRI_200_z + 
  strata(stratum)

form1a <- used ~ dem_200_z * TRI_200_z * weeks_since_emig_n + 
  strata(stratum)

form1a <- used ~ weeks_since_emig_n + 
  strata(stratum)

#I go more into detail of using this model in 03_04_clogit_workflow.R

# STEP 3: create new data for prediction ----------------------------------------------------------------

#add one new row to unique strata instead of entire empty copies of strata. assign week since emigration and terrain values on a regular grid, so we can make a raster later on
set.seed(7777)

#n needs to be large enough to cover the whole range of 
n <- 1000

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 10), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 10), n, replace = T),
         t2m_z = sample(seq(min(all_data$t2m_z),max(all_data$t2m_z), length.out = 10), n, replace = T)) %>% #set to the mean of temperature. it is zero, because we are working with a z-score 
  full_join(all_data)

saveRDS(new_data,"alt_50_20_min_48_ind_static_daytemp_100_inlaready_wmissing_wks_n1000.rds")


#the model will be run on the cluster. see cluster_prep/order_of_business_main_model.txt
