## 06.20.2022
## golden eagle excursion date and length estimation
## Hester Brønnvik (mpiab)

# required libraries
library(move)
library(lubridate)
library(sf)
library(stringr)
require(recurse)
require(scales)
require(sp)
library(dplyr)
# required files:
# the login data for Movebank
load("loginStored.RData")
# the unique ID of the LifeTrack Golden Eagle Alps study on Movebank
eagleStudyId <- 282734839
# basic information on the birds in the study
all_inds <- getMovebank("individual", login=loginStored, study_id=eagleStudyId)
# the previously estimated emigration dates
load("em_fl_dt_recurse_70ind.RData")
# the names of the individuals in those data
inds_to_do <- emig_fledg_dates$individual.local.identifier

start_time <- Sys.time()
# a loop to go through each individual and estimate emigration time
excursions <- lapply(inds_to_do, function(x){
  print(x)
  # select the first time
  ts <- all_inds$timestamp_start[which(all_inds$local_identifier == x)]
  # select the previously estimated dispersal date for the given individual
  et <- paste0(emig_fledg_dates$emigration_dt[which(emig_fledg_dates$individual.local.identifier == x)], "000")
  # download the data between those days from Movebank
  ind <- getMovebankData(study = "LifeTrack Golden Eagle Alps", 
                         animalName =  x, removeDuplicatedTimestamps=T,
                         login = loginStored, timestamp_start = gsub("\\.|\\-|\\:| ","", ts), timestamp_end = gsub("\\.|\\-|\\:| ","", et))
  print(paste0("Downloaded data for ", x, "."))
  # project to UTM so that the units are in meters and the radius is easy to set
  ind <- spTransform(ind, CRSobj = "+init=EPSG:32632")
  # select only the first location as a nest proxy
  loc1 <- data.frame(x = ind$location_long[1], y = ind$location_lat[1])
  # find all the times that an animal crossed the circle with center at loc1 and radius 7km
  ind_visit <- getRecursionsAtLocations(ind, locations = loc1, radius =  7000, timeunits = "days")

  # record excursion times (times when the bird left the radius, however briefly, before dispersal)
  recursions <- ind_visit$revisitStats
  
  ifelse( # if there are not excursions:
                  nrow(recursions) == 1,
                   excursions <- data.frame(local_identifier = x,
                                            excursion_start = NA,
                                            excursion_end = NA,
                                            excursion_length_in_days = NA), 
                   # if there are excursions:
                   excursions <- lapply(1:(nrow(recursions)-1), function(y){
                     data.frame(local_identifier = x, 
                                excursion_start = recursions$exitTime[y], 
                                excursion_end = recursions$entranceTime[y+1],
                                excursion_length_in_days = recursions$timeSinceLastVisit[y+1])})
                   %>% reduce(rbind)
  )
 excursions
}) %>% reduce(rbind)
# run time:
Sys.time() - start_time
