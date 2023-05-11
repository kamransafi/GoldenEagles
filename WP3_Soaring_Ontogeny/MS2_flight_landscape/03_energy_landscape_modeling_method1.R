#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R (and temp_download_&prep.R)
#only includes static variables. 
#try adding daily temp to see. seems like dem and temp are interacting!.. from 03_method2
#Feb 22. 2022. Elham Nourani. Konstanz, DE
#revisiting this code on April 11. 2023:
# - use hourly steps instead of 20 minutes
# - include step length in the analysis
# - consider using the 25 resolution for the terrain


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
#library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
#library(oce) #color palette for interaction plots
library(modelsummary)

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

#open annotated data. prepared in 03_04b_clogit_randomization.R
data <-  readRDS("alt_50_60_min_55_ind_static_r_100.rds") %>% 
  filter(weeks_since_emig <= 156) %>%  #cap the week since at 156 weeks (3 yrs).
  mutate_at(c("step_length", "dem_100", "TRI_100", "ridge_100", "weeks_since_emig"), list(z = ~(scale(.)))) %>%   #calculate the z scores
  mutate(ind1 = individual.local.identifier,
         ind2 = individual.local.identifier,
         ind3 = individual.local.identifier,
         ind4 = individual.local.identifier,
         month = month(timestamp)) 

saveRDS(data, file = "all_inds_annotated_static_3yrs_apr23.rds") #for a version with dynamic variables, see all_inds_annotated_apr23.rds

### look at the correlations
data %>% 
  dplyr::select(c("step_length", "dem_100", "TRI_100", "ridge_100", "weeks_since_emig")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.6) #no autocorrelation

# STEP 2: CLOGIT- ssf modeling exploration  ----------------------------------------------------------------

data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds")

## mid_step: investigate using clogit
f <-  used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100 * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf1 <- clogit(f, data = data)
summary(ssf1)
plot_summs(ssf1)
modelsummary(ssf1) #AIC = 347,895.0

f2 <-  used ~ TRI_100 * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf2 <- clogit(f2, data = data)
summary(ssf2)
plot_summs(ssf2)
modelsummary(ssf2) #AIC = 351,242.1

#I go more into detail of using this model in 03_04_clogit_workflow.R

# STEP 2: INLA code for Raven ----------------------------------------------------------------
##full model with seasonality (to run on raven in oder_of_business_main_model.txt)
F_full <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z + #make sure all vars are scaled
  TRI_100_z * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", cyclic = T, scale.model = TRUE)) +
  f(ind2, TRI_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", cyclic = T, scale.model = TRUE)) +
  f(ind3, ridge_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", cyclic = T, scale.model = TRUE))

## model without seasonality to compare
F_OG <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, ridge_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

## model with seasonality only for dem. this model crashed on Raven. I've tried to run the model with only tri and ridge, which dont have seasonality
F_full <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))), #seasonality only for dem
    group = month, control.group = list(model = "ar1", scale.model = TRUE, cyclic = T)) +
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, ridge_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


# STEP 3: create new data for predictions to plot interaction terms----------------------------------------------------------------

data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") #this is limited to 3 yrs post dispersal

#the model will be run on the cluster.

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. one per interaction term. merge later on
grd_dem <- expand.grid(x = (1:104), #make the preds up to week 50, which is the mean of weeks_since and is almost one year
                       y = seq(from = min(data$dem_100, na.rm = T), to = quantile(data$dem_100, .9, na.rm = T), by = 200)) %>% # n = 1050
  rename(weeks_since_emig = x,
         dem_100 = y) %>% 
  mutate(TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'), #set other variables to their mean
         ridge_100 = attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'),
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         interaction = "wk_dem_100")

grd_tri <- expand.grid(x = (1:104),
                       y = seq(from = min(data$TRI_100, na.rm = T), to = quantile(data$TRI_100, .9, na.rm = T), by = 20)) %>%  #n = 675
  rename(weeks_since_emig = x,
         TRI_100 = y)  %>% 
  mutate(dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'), #set other variables to their mean
         ridge_100 = attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'),
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         interaction = "wk_TRI_100")

grd_dist <- expand.grid(x = (1:104),
                        y = seq(from = min(data$ridge_100, na.rm = T), to = quantile(data$ridge_100, .9, na.rm = T), by = 25)) %>%  # 525
  rename(weeks_since_emig = x,
         ridge_100 = y) %>% 
  mutate(TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'), #set other variables to their mean
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'),
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         interaction = "wk_ridge_100")

grd_sl <- expand.grid(x = (1:104),
                      y = seq(from = min(data$step_length, na.rm = T), to = quantile(data$step_length, .9, na.rm = T), by = 200)) %>%  # 226
  rename(weeks_since_emig = x,
         step_length = y) %>% 
  mutate(TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'), #set other variables to their mean
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'),
         ridge_100 = attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'),
         interaction = "wk_step_length")


grd_all <- bind_rows(grd_dem, grd_tri, grd_dist, grd_sl) %>% 
  mutate(used = NA,
         month = 6) #make all preds for month 6


set.seed(777)
n <- nrow(grd_all)

new_data_only <- data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  #only keep the columns that I need
  dplyr::select(c("stratum", "ind1", "ind2", "ind3")) %>% 
  bind_cols(grd_all) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate( dem_100_z = (dem_100 - attr(data[,colnames(data) == "dem_100_z"],'scaled:center'))/attr(data[,colnames(data) == "dem_100_z"],'scaled:scale'),
          TRI_100_z = (TRI_100 - attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'))/attr(data[,colnames(data) == "TRI_100_z"],'scaled:scale'),
          ridge_100_z = (ridge_100 - attr(data[,colnames(data) == "ridge_100_z"],'scaled:center'))/attr(data[,colnames(data) == "ridge_100_z"],'scaled:scale'),
          step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'), 
          weeks_since_emig_z = (weeks_since_emig - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'))

saveRDS(new_data_only, file = "new_data_only_ssf_preds.rds")

new_data <- data %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only)

saveRDS(new_data, file = "new_data_ssf_inla_preds.rds")

