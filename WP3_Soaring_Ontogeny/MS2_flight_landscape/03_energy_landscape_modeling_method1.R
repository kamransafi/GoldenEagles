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


#open annotated data. prepared in 03_04b_clogit_randomization.R
data <-  readRDS("alt_50_60_min_55_ind_static_100.rds") %>% 
  mutate_at(c("step_length", "dem_100", "TRI_100", "slope_TPI_100", "weeks_since_emig"), list(z = ~(scale(.)))) %>%   #calculate the z scores
  mutate(ind1 = individual.local.identifier,
         ind2 = individual.local.identifier,
         ind3 = individual.local.identifier,
         month = month(timestamp))

saveRDS(data, file = "all_inds_annotated_static_apr23.rds") #for a version with dynamic variables, see all_inds_annotated_apr23.rds

# STEP 2: CLOGIT- ssf modeling exploration  ----------------------------------------------------------------

## mid_step: investigate using clogit

form1 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100_z * step_length_z * weeks_since_emig_z 
  slope_TPI_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf1 <- clogit(form1, data = data)
summary(ssf1)
plot_summs(ssf1)
modelsummary(ssf1) #AIC = 196770.1

form2 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100_z * step_length_z * weeks_since_emig_z + 
  slope_TPI_100_z * step_length_z * weeks_since_emig_z + 
  t2m +
  strata(stratum)

ssf2 <- clogit(form2, data = data)
summary(ssf2)
plot_summs(ssf2)
modelsummary(ssf2) #AIC = 189219.6

#I go more into detail of using this model in 03_04_clogit_workflow.R

# STEP 3: INLA parameterization ----------------------------------------------------------------
#the aim is to take seasonality into account

#do a test with a small sample of the data
sample <- data %>% 
  filter(stratum %in% sample(unique(data$stratum), 300, replace = F))


mean.beta <- 0
prec.beta <- 1e-4

#define variables and formula. use one variable for now to test the model 

#model without seasonality
F1 <- used ~ -1 + dem_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(3,0.05))))

M2 <- inla(F1, family = "Poisson", 
          control.fixed = list(
            mean = mean.beta,
            prec = list(default = prec.beta)),
          data = sample, 
          num.threads = 7,
          control.predictor = list(compute = TRUE, link = 1), 
          control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))

Efxplot(M1)

##model with seasonality
F2 <- used ~ -1 +
  dem_100_z * weeks_since_emig_z * step_length_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE))

M2 <- inla(F2, family = "Poisson",
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = sample,
           num.threads = 5)

Efxplot(list(M1,M2))

# STEP 3: INLA code for Raven ----------------------------------------------------------------
##full model with seasonality (to run on raven in oder_of_business_main_model.txt)
F_full <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE)) +
  f(ind2, TRI_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE)) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE))

## model without seasonality to compare
F_OG <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

# STEP 3: create new data for prediction: step length * dem ----------------------------------------------------------------

data <- readRDS("all_inds_annotated_static_apr23.rds")

#the model will be run on the cluster.
#create one separate dataset for each interaction plot, and keep other variables at their mean

#add one new row to unique strata instead of entire empty copies of strata. assign values on a regular grid, to make a raster later on
#n needs to be large enough to cover the whole range of vars

set.seed(7)
n <- 1000

new_data <- data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  mutate(used = NA,
         dem_100 = sample(seq(min(data$dem_100, na.rm = T), quantile(data$dem_100, 0.9, na.rm = T), length.out = 10), n, replace = T), #use the 90% quantile instead of max. get rid of outliers
         step_length = sample(seq(min(data$step_length, na.rm = T), quantile(data$step_length, 0.9, na.rm = T), length.out = 10), n, replace = T),
         TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'),
         slope_TPI_100 = attr(data[,colnames(data) == "slope_TPI_100_z"],'scaled:center'),
         weeks_since_emig = attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'), 
         month = 6)  %>%  #pick June for all preds.
  mutate(dem_100_z = (dem_100 - attr(data[,colnames(data) == "dem_100_z"],'scaled:center'))/attr(data[,colnames(data) == "dem_100_z"],'scaled:scale'),
         step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'), #get the mean(center) and sd(scale) from the previous z transformation. for consistency
         TRI_100_z = 0,
         slope_TPI_100_z = 0,
         weeks_since_emig_z = 0) %>% 
  full_join(data) %>% 
  as.data.frame()
  
saveRDS(new_data,"inla_preds_input_sl_dem.rds")

# STEP 4: create new data for prediction: dem * week since dispersal ----------------------------------------------------------------

set.seed(77)
n <- 1000

new_data <- data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  mutate(used = NA,
         dem_100 = sample(seq(min(data$dem_100, na.rm = T), quantile(data$dem_100, 0.9, na.rm = T), length.out = 10), n, replace = T), #use the 90% quantile instead of max. get rid of outliers
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'),
         slope_TPI_100 = attr(data[,colnames(data) == "slope_TPI_100_z"],'scaled:center'),
         weeks_since_emig = sample(seq(min(data$weeks_since_emig, na.rm = T), quantile(data$weeks_since_emig, 0.9, na.rm = T), length.out = 10), n, replace = T), 
         month = 6)  %>%  #pick June for all preds.
  mutate(dem_100_z = (dem_100 - attr(data[,colnames(data) == "dem_100_z"],'scaled:center'))/attr(data[,colnames(data) == "dem_100_z"],'scaled:scale'),
         step_length_z = 0, 
         TRI_100_z = 0,
         slope_TPI_100_z = 0,
         weeks_since_emig_z = (weeks_since_emig - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale')) %>% 
  full_join(data) %>% 
  as.data.frame()

saveRDS(new_data,"inla_preds_input_dem_wk.rds")

# STEP 5: create new data for prediction: tri * week since dispersal ----------------------------------------------------------------

set.seed(777)
n <- 1000

new_data <- data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  mutate(used = NA,
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'), 
         step_length = attr(data[,colnames(data) == "step_length_z"],'scaled:center'),
         TRI_100 = sample(seq(min(data$TRI_100, na.rm = T), quantile(data$TRI_100, 0.9, na.rm = T), length.out = 10), n, replace = T),
         slope_TPI_100 = attr(data[,colnames(data) == "slope_TPI_100_z"],'scaled:center'),
         weeks_since_emig = sample(seq(min(data$weeks_since_emig, na.rm = T), quantile(data$weeks_since_emig, 0.9, na.rm = T), length.out = 10), n, replace = T), 
         month = 6)  %>%  #pick June for all preds.
  mutate(dem_100_z = 0,
         step_length_z = 0, 
         TRI_100_z = (TRI_100 - attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'))/attr(data[,colnames(data) == "TRI_100_z"],'scaled:scale'),
         slope_TPI_100_z = 0,
         weeks_since_emig_z = (weeks_since_emig - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale')) %>% 
  full_join(data) %>% 
  as.data.frame()

saveRDS(new_data,"inla_preds_input_tri_wk.rds")

# STEP 6: create new data for prediction: sl * week since dispersal ----------------------------------------------------------------

set.seed(7777)
n <- 1000

new_data <- data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  mutate(used = NA,
         dem_100 = attr(data[,colnames(data) == "dem_100_z"],'scaled:center'), 
         step_length = sample(seq(min(data$step_length, na.rm = T), quantile(data$step_length, 0.9, na.rm = T), length.out = 10), n, replace = T),
         TRI_100 = attr(data[,colnames(data) == "TRI_100_z"],'scaled:center'),
         slope_TPI_100 = attr(data[,colnames(data) == "slope_TPI_100_z"],'scaled:center'),
         weeks_since_emig = sample(seq(min(data$weeks_since_emig, na.rm = T), quantile(data$weeks_since_emig, 0.9, na.rm = T), length.out = 10), n, replace = T), 
         month = 6)  %>%  #pick June for all preds.
  mutate(dem_100_z = 0,
         step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'), 
         TRI_100_z = 0,
         slope_TPI_100_z = 0,
         weeks_since_emig_z = (weeks_since_emig - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale')) %>% 
  full_join(data) %>% 
  as.data.frame()

saveRDS(new_data,"inla_preds_input_sl_wk.rds")


