#Script for step-selection analysis of golden eagle commuting flights as reported in Nourani et al. 2023
#this script contains the code for reproducing all the plots in the paper
#Elham Nourani, PhD. 25.07.2023
#enourani@ab.mpg.de

library(tidyverse)
library(lubridate)
library(corrr)
library(glmmTMB)

##### STEP 0: Open annotated data #####

data <- read.csv("GPS_data_annotated.csv") %>%  #this is the post-dispersal data. weeks_since_emig is the column representing weeks since dispersal
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"),
         animal_ID = as.numeric(as.factor(individual.local.identifier)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum)))

##### STEP 2: Check for autocorrelation #####

data %>% 
  dplyr::select(c("step_length", "TRI_100", "ridge_100", "weeks_since_emig")) %>% 
  correlate() 

##### STEP 3: Scale the predictor variables #####

data <- data %>% 
  mutate_at(c("step_length", "TRI_100", "ridge_100", "weeks_since_emig"), list(z = ~(scale(.))))

##### STEP 4: Run the ssf model #####
#set up the model structure based on Muff et al: 
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

TMB_struc <- glmmTMB(used ~ -1 + TRI_100_z * step_length_z * weeks_since_emig_z + 
                       ridge_100_z * step_length_z * weeks_since_emig_z + (1|stratum_ID) + 
                       (0 + ridge_100_z | animal_ID) + 
                       (0 + TRI_100_z | animal_ID), 
                     family = poisson, data = data, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA,1:2))), #2 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3),0,0))) #add a 0 for each random slope. in this case, 2


TMB_M <- glmmTMB:::fitTMB(TMB_struc) #this throws an error. redoing the data prep didnt produce the same dataset as previously done. 
summary(TMB_M)

saveRDS(TMB_M, file = "TMB_model.rds")

##### STEP 5: Model validation #####

##### STEP 6: PLOT Fig.1 - coefficient estimates #####


##### STEP 7: PLOT Fig.S1 - individual-specific coefficients #####

##### STEP 8: PLOT Fig. 2 - interaction terms #####

##### STEP 9: Predictions for the Alps #####

##### STEP 10: PLOT Fig.3 flyability hotspot maps #####

##### STEP 11: PLOT video S1 energy landscape maps #####

##### System info #####

#the tracking data contains commuting flights for 55 juvenile golden eagles in the alps

tracking_data <- read.csv("ssf_input_data.csv") %>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

dispersal_dt <- read.csv("dispersal_dates.csv") %>% 
  mutate(dispersal_dt = as.POSIXct(dispersal_dt, tz = "UTC"))

##### STEP 2: annotate with week since dispersal#####

tracking_data <- tracking_data %>% 
  left_join(dispersal_dt %>% dplyr::select(individual.local.identifier, dispersal_dt), by = individual.local.identifier) %>% 
  mutate(days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")) %>% as.numeric() %>%  ceiling(),
         weeks_since_emig = difftime(timestamp,emigration_dt, units = c("weeks")) %>% as.numeric() %>%  ceiling())

dd <- dd %>% 
  left_join(dispersal_dt %>% dplyr::select(individual.local.identifier, dispersal_dt), by = "individual.local.identifier") %>% 
  mutate(days_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("days")),
         weeks_since_emig = difftime(date(timestamp), date(dispersal_dt), units = c("weeks")))

