#Script for step-selection analysis of golden eagle commuting flights as reported in Nourani et al. 2023
#this script contains the code for reproducing all the plots in the paper
#Elham Nourani, PhD. 25.07.2023
#enourani@ab.mpg.de

library(tidyverse)
library(lubridate)

##### STEP 0: Open ssf input data #####


##### STEP 2: Check for autocorrelation #####

##### STEP 3: Scale the predictor variables #####

##### STEP 4: Run the ssf model #####

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

