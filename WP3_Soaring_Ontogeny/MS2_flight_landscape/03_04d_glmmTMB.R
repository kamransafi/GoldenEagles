#script for constructing the energy landscape for golden eagles
#using glmmTMB, which allows for random effects to be included
#follows from 03_energy_landscape_modeling_method1.R, 04_INLA_output_plots.R
#procedure similar to 03_04_clogit_workflow.R
#21.06.2023. Canberra, AU.
#Elham Nourani, PhD.

library(tidyverse)
library(glmmTMB)

#code is based on Muff et al:
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open data. data was prepared in 03_energy_landscape_modeling_method1.R. 
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  mutate(animal_ID = as.numeric(as.factor(ind1)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum)))

TMB_struc <- glmmTMB(used ~ -1 + TRI_100_z * step_length_z * weeks_since_emig_z + 
                     ridge_100_z * step_length_z * weeks_since_emig_z + (1|stratum_ID) + 
                     (0 + ridge_100_z | animal_ID) + 
                     (0 + TRI_100_z | animal_ID), 
                   family = poisson, data = data, doFit = FALSE,
                   #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                   map = list(theta = factor(c(NA,1:2))), #2 is the n of random slopes
                   #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                   start = list(theta = c(log(1e3),0,0))) #add a 0 for each random slope. in this case, 2

#TMBStruc$parameters$theta[1] = log(1e3) 
#TMBStruc$mapArg = list(theta=factor(c(NA,1:2)))


TMB_m <- glmmTMB:::fitTMB(TMB_struc)
summary(TMB_m)

#extract coefficient estimates and confidence intervals
confint(TMB_m)

#extract individual-specific random effects:
ranef(TMB_m)[[1]]$animal_ID



