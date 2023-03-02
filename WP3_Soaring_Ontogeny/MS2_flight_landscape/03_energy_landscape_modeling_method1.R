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

# STEP 5: ssf modeling ----------------------------------------------------------------

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

#use the clogit to make predictions. 
n <- 1000

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 20), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 20), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 20), n, replace = T),
         t2m_z = sample(seq(min(all_data$t2m_z),max(all_data$t2m_z), length.out = 20), n, replace = T))


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
  
  #saveRDS(avg_pred, file = paste0("inla_pred_clogit_", i,".rds"))
  
  r <- avg_pred %>% 
    rast(type = "xyz") %>% 
    focal(w = 7, fun = max, na.policy = "only", na.rm = T) %>% 
    as.data.frame(xy = T) %>%
    rename(avg_pres = focal_mean)
  
  saveRDS(r, file = paste0("inla_pred_clogit_", i,".rds"))
  
  #coordinates(avg_pred) <-~ x + y
  #gridded(avg_pred) <- TRUE
  #r <- raster(avg_pred) #many NA values....
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






#######
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
