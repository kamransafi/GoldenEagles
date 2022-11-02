### Juvenile golden eagle soaring performance code
### Hester Brønnvik
### 25.03.2022
### hbronnvik@ab.mpg.de


## This file contains code for three tasks: 1) calculating time spent soaring, 2) calculating 
## maximum wind speed in thermals per day, and 3) calculating ODBA per time in each thermal.

## Packages and functions:
library(move)
library(lubridate)
library(data.table)
library(tidyverse)
library(ggpubr)
#library(plyr)
library(doParallel)
is.error <- function(x) inherits(x, "try-error")
source("C:/Users/hbronnvik/Documents/golden_eagle_tools/associate_ACCinfoToGPS.R") #For ACCtoGPS_parallel function
load("C:/Users/hbronnvik/Documents/storkSSFs/Laptop/loginStored.RData") # The MoveBank login
times <- read.csv("C:/Users/hbronnvik/Documents/golden_eagle_tools/fledging_emig_times.csv") %>% 
  select(-"X")
eagleStudyId <- 282734839 # The golden eagle study on MoveBank

## Task 1: ratio of soaring modes
##############################################################################################
## 21.02.22

## 1. Calculate the total time spent flying in these data (gliding, soaring, flapping)
## 2. Identify all soaring events
## 3. Measure how much time is spent in each
## 4. Sum these over each day
## 5. Calculate the ratio of soaring to flight
## 6. Plot the ratio over days since fledging

### Classify behaviors
# birds to remove:
files <- list.files("D:/goldenEagles_fromMartina/accGPS_behavClass", full.names = T)[-c(3,4,5,23,24)]
# list the data with wind speed classifications removing the adults 
wind_files <- list.files("D:/goldenEagles_wind", full.names = T)
# the adult birds that should be removed
#ind_id <- c("Aosta1_20 (eobs 7511)", "Aosta2_20 (eobs 7558)", "Aosta21 (eobs7590)", "Mellau21 (eobs 6988)", "Memmingen20 (eobs 7507)")
# the estimated dates


file_names <- str_sub(files, 47, -24)
files <- files[-which(!file_names %in% times$individual.local.identifier)]

invisible(lapply(files, function(x){
  # read in the data for this ID containing flapping classifications 
  ind <- readRDS(x)
  # extract ID
  ID <- unique(ind$local_identifier)
  # read in the data for this ID containing wind estimates 
  load(wind_files[grepl(ID, wind_files, fixed = T)]) # burstsWindDF
  
  date <- times %>% filter(individual.local.identifier == ID)
  
  pfdp_wind <- burstsWindDF%>% 
    filter(between(timestamp, date$fledging_dt, date$emigration_dt)) %>% 
    rename(location_long = location.long,
           location_lat = location.lat)
  
  pfdp_behav <- ind %>% 
    filter(between(timestamp, date$fledging_dt, date$emigration_dt))
  
  pfdp <- pfdp_behav %>% 
    left_join(pfdp_wind, by = intersect(colnames(pfdp_behav), colnames(pfdp_wind))) %>% 
    #filter(between(timestamp, date$fledging_dt, date$emigration_dt)) %>% 
    mutate(flight_type = ifelse(behavior == "Flapping",
                                "flap",
                                ifelse(soarClust == "soar" & thermalClust == "circular",
                                       "soar_circular", 
                                       ifelse(thermalClust == "linear",
                                              "soar_linear",
                                              ifelse(soarClust == "glide",
                                                     "glide",
                                                     NA))))) %>% 
    mutate(flight_type = ifelse(is.na(thermalID) & flight_type == "soar_circular",
                                "circular_no_ID",
                                flight_type))
  
  
  # wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  # pfdpsf <- st_as_sf(pfdp, coords = c("location_long", "location_lat"), crs = wgs)
  # check <- pfdpsf[1:500,]
  # mapview(check, zcol = "flight_type")
  
  #saveRDS(pfdp, file = paste0("D:/Golden_eagle_wind_behav_June_2022/", ID, ".rds"))
  
}))

# 15 individuals have 0 rows of data due to missing wind estimates, possibly because
# they fledged in the fall and data transmission was poor during the pfdp

 
### Ratio of soaring to flight

pfdp_files <- list.files("C:/Users/hbronnvik/Documents/golden_eagle_tools/Golden_eagle_wind_behav_June_2022", full.names = T)
# remove individuals without data
pfdp_files <- pfdp_files[sapply(pfdp_files, file.size) > 2000]

# determine the total time spent in each behavior for each individual
flight_durations <- lapply(pfdp_files, function(x){
  ind <- readRDS(x)  
  ID <- unique(ind$local_identifier)
  
  # determine the center point of each day to allow environmental annotation
  ind <- ind %>% 
    mutate(date = date(timestamp)) %>% 
    group_by(date) %>% 
    mutate(poi = centroid(as.matrix(cbind(location_long, location_lat)))) %>%
    ungroup()

  # find the time spent in each behavior per day
  behavior_time <- ind %>% 
    drop_na(flight_type) %>% 
    group_by(date, flight_type) %>% 
    summarize(duration = difftime(tail(timestamp, 1), head(timestamp,1), units = "sec")) %>% 
    summarize(total_time = sum(duration),
              circle_time = sum(duration[flight_type %in% c("soar_circular", "circular_no_ID")]),
              linear_time = sum(duration[flight_type == "soar_linear"]),
              glide_time = sum(duration[flight_type == "glide"]),
              flap_time = sum(duration[flight_type == "flap"])) %>% 
    mutate(local_identifier = ID)
  
  # append the center points
  behavior_time <- full_join(ind[,c("date", "poi")], behavior_time) %>% 
    group_by(date) %>% 
    slice(1) %>% 
    # reduce the matrix
    mutate(center_long = poi[,1],
           center_lat = poi[,2]) %>% 
    ungroup() %>%
    # clean the unneeded column
    select(-poi) %>% 
    # remove dates for which behavior classification did not happen
    drop_na(local_identifier)
  
  return(behavior_time)
}) %>% reduce(rbind)

# calculate the ratios of each flight type to total flight time
flight_ratios <- flight_durations %>% 
  mutate(circles = as.numeric(circle_time)/as.numeric(total_time),
         lines = as.numeric(linear_time)/as.numeric(total_time),
         glides = as.numeric(glide_time)/as.numeric(total_time),
         flaps = as.numeric(flap_time)/as.numeric(total_time),
         soar_to_flight = (as.numeric(circle_time) + as.numeric(linear_time))/as.numeric(total_time)) %>% 
  rename(date = 1) %>% 
  left_join(times[c("individual.local.identifier", "emigration_dt", "fledging_dt")], by = c("local_identifier" = "individual.local.identifier")) %>% 
  # days since fledging for each individual
  mutate(days_since_fledging = difftime(date, fledging_dt, units = "days")) %>% 
  filter(total_time > 0)

# prep for movebank annotation
flight_ratios <- flight_ratios %>% 
  mutate(timestamp = paste0(date, " 12:00:00.000")) %>% 
  rename("location-long" = center_long,
         "location-lat" = center_lat)

# write.csv(flight_ratios, file = paste0("C:/Users/hbronnvik/Documents/golden_eagle_tools/behav_class_GE_", Sys.Date(), ".csv"))


library(reshape2)
plot_ratios <- flight_ratios %>% 
  rename("Thermal soaring" = circles,
         "Linear soaring" = lines,
         "Flapping" = flaps,
         "Gliding" = glides) %>% 
  melt(id.vars = c("local_identifier", "date", "days_since_fledging")) %>% 
  filter(variable %in% c("Thermal soaring", "Linear soaring", "Flapping", "Gliding")) %>% 
  mutate(days_since_fledging = as.numeric(days_since_fledging),
         value = as.numeric(value))

ggplot(plot_ratios, aes(days_since_fledging, value))+
  geom_smooth(method = "lm", se = F, color = "black")+
  geom_point(aes(color = variable))+
  labs(x = "Days since fledging", y = "Ratio of time spent per flight time")+
  theme_classic()+
  facet_wrap(~variable, scales = "free")+
  theme(legend.position="none", axis.text = element_text(color = "black"))



### Ratio of gliding to flapping between thermals

# find consecutive, unique thermals
# measure time spent flapping 
# measure time spent gliding
# calculate the ratio

### Wind speed exposure during thermal soaring events
pfdp_files <- list.files("D:/Golden_eagle_wind_behav_June_2022", full.names = T)
# remove individuals without data
pfdp_files <- pfdp_files[sapply(pfdp_files, file.size) > 2000] 

wind_estimates <- lapply(pfdp_files, function(x){
  ind <- readRDS(x) %>% 
    drop_na(location_lat)
  # extract wind speeds per thermal
  wind_exposure <- ind %>% 
    drop_na(thermalID) %>% 
    group_by(local_identifier, thermalID) %>% 
    arrange(timestamp) %>% 
    mutate(wind_speed = sqrt((windX)^2 + (windY)^2)) %>% 
    summarize(max_wind = max(wind_speed, na.rm = T),
              min_wind = min(wind_speed, na.rm = T),
              mean_wind = mean(wind_speed, na.rm = T), 
              timestamp = head(timestamp, 1),
              location_lat = head(location_lat, 1),
              location_long = head(location_long, 1))
  return(wind_exposure)
}) %>% reduce(rbind)
  
wind_times <- wind_estimates %>% 
  left_join(times[c(-4)], by = c("local_identifier" = "individual.local.identifier")) %>% 
  mutate(hours_since_fledging = difftime(timestamp, fledging_dt, units = "hours"))

ggplot(wind_times, aes(hours_since_fledging, max_wind))+
  geom_point()+
  geom_smooth(method = "lm", se = F, aes(color = local_identifier))+
  theme_classic()+
  theme(legend.position="none", axis.text = element_text(color = "black"))


wind_daily <- wind_times %>% 
  mutate(date = date(timestamp)) %>% 
  group_by(local_identifier, date) %>% 
  summarise(wind_speed = mean(mean_wind))

flight_ratios <- full_join(flight_ratios, wind_daily, by = c("local_identifier", "date"))

# check the correlations of the variables
pairs(~ circles + lines + glides + flaps + wind_speed, data = flight_ratios, upper.panel = NULL)

library(lme4)
mod_data <- flight_ratios %>% 
  select(circles, lines, glides, flaps, local_identifier, days_since_fledging) %>% 
  mutate(days_since_fledging = as.numeric(days_since_fledging), 
         circles = scale(circles), 
         lines = scale(lines),
         glides = scale(glides),
         flaps = scale(flaps))

mod <- lmer(days_since_fledging ~ circles + lines + glides + flaps + (1|local_identifier), data = mod_data)
qqnorm(resid(mod))

##############################################################################################

## Task 2: wind speed exposure within thermals
##############################################################################################
## 10.03.2022


## 1. Identify the maximum wind speed experienced each day
## 2. Plot against days since fledging


# calculate the wind speed at each location
# mod(180 + 180/pi * atan2(u,v), 360)
# sqrt(u^2 + v^2)
cfls$wind_speed <- sqrt((cfls$windX)^2 + (cfls$windY)^2)


# extract the maximum wind speed on each day
pd_ws <- data.frame()

for (i in unique(cfls$id_date)) {
  # for each unique individual on each day
  df <- cfls[which(cfls$id_date == i),]
  # identify the maximum wind speed encountered
  df2 <- df[which.max(df$wind_speed),]
  # save the date and maximum for plotting
  pd_ws <- rbind(pd_ws, df2)
  # signal
  print(paste0("Separated wind speed for ", i, "."), quote = F)
  rm(df);rm(df2);rm(i)
}

## 2. Plot
ggplot(pd_ws, aes(x = dsf, y = wind_speed, color = individual_local_identifier)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Maximum wind speed") +
  scale_x_continuous(breaks = seq.int(0, max(pd_ws$dsf), by = 40)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
##############################################################################################

## Task3: DBA during thermaling
##############################################################################################
## 14.03.2022


## 1. Associate accelerometry data and calculate ODBA
## 2. Identify unique thermals
## 3. Sum ODBA within each thermal
## 4. Find the ratio of ODBA to time spent 
## 5. Plot the ratio against days since fledging

## 1. Associate accelerometry data and calculate ODBA (code from Martina)

## Download all of the acc data
# List names of the individuals included in the study and filter only those having ACC data
allInds <- getMovebank("individual", login=loginStored, study_id=eagleStudyId)
cfls <- data.frame(individual_local_identifier = str_sub(pfdp_files, 82, -5))
allInds <- allInds[which(allInds$local_identifier %in% unique(cfls$individual_local_identifier)),]
# Keep only individuals that have ACC information and were not "tests"
allInds_acc <- allInds[grep("gps|acceleration", allInds$sensor_type_ids),]
allInds_acc <- allInds_acc[which(!allInds_acc$local_identifier %in% grep("test|Test", allInds_acc$local_identifier, value=T)),]
# For the remaining individuals, download the ACC data one by one
cl <- makeCluster(detectCores()-4, export = list("loginStored","eagleStudyId", "cfls"), 
                  lib = list("move", "tidyverse", "lubridate"))
clusterExport(cl, list=c("cfls", "eagleStudyId", "loginStored", "allInds_acc", "ACCtoGPS_parallel", "is.error")) 

library(plyr)
library(dplyr)
start_time <- Sys.time()
if(nrow(allInds_acc) > 0){
  indInfo <- data.frame(indName=as.character(allInds_acc$local_identifier[!allInds_acc$local_identifier%in%c(NA,"")]),
                        indId=as.numeric(allInds_acc$id[!allInds_acc$local_identifier%in%c(NA,"")]))
  accLs <- llply(indInfo$indId, function(ind)try({  #with llply
    # library(move)
    # library(lubridate)
    # library(tidyverse)
    # load("C:/Users/hbronnvik/Documents/storkSSFs/Laptop/loginStored.RData")
    # is.error <- function(x) inherits(x, "try-error")
    # source("C:/Users/hbronnvik/Documents/golden_eagle_tools/associate_ACCinfoToGPS.R") #For ACCtoGPS_parallel function
    # load("C:/Users/hbronnvik/Documents/storkSSFs/Laptop/loginStored.RData") # The MoveBank login
    # eagleStudyId <- 282734839 # The golden eagle study on MoveBank
    # allInds <- getMovebank("individual", login=loginStored, study_id=eagleStudyId)
    # pfdp_files <- list.files("C:/Users/hbronnvik/Documents/golden_eagle_tools/Golden_eagle_wind_behav_June_2022", full.names = T)
    # # remove individuals without data
    # pfdp_files <- pfdp_files[sapply(pfdp_files, file.size) > 2000]
    # cfls <- data.frame(individual_local_identifier = str_sub(pfdp_files, 82, -5))
    # allInds <- allInds[which(allInds$local_identifier %in% unique(cfls$individual_local_identifier)),]
    # # Keep only individuals that have ACC information and were not "tests"
    # allInds_acc <- allInds[grep("gps|acceleration", allInds$sensor_type_ids),]
    # allInds_acc <- allInds_acc[which(!allInds_acc$local_identifier %in% grep("test|Test", allInds_acc$local_identifier, value=T)),]
    # Download only the timestamps of the acc data
    Trange <- getMovebank(entity="event", study_id=eagleStudyId, individual_id=ind, sensor_type_id=2365683, login=loginStored,
                          attributes=c("timestamp"))
    Trange <- as.vector(Trange$timestamp)
    print(paste0("Time range created."))
    if(length(Trange)>0){
      # Split the timestamps in chunks of 100k lines per download (otherwise the download breaks)
      dlsq <- cbind(seq(0, length(Trange), by=100000)+1, c(seq(0, length(Trange), by=100000), length(Trange))[-1])
      # Create start and end timestamps for download
      DLDtimeSeq <- data.frame(startTime=Trange[dlsq[,1]], endTime=Trange[dlsq[,2]])
      print(paste0("Start and end times created."))
      # Download the chunks one by one
      ACCdwld <- do.call(rbind, lapply(1:nrow(DLDtimeSeq), function(j){
        chunkACC <- getMovebankNonLocationData(study=eagleStudyId, animalName=ind, sensorID="Acceleration", login=loginStored, 
                                               timestamp_start=gsub(".", "", gsub("-|:| ", "", DLDtimeSeq$startTime[j]), fixed=T), 
                                               timestamp_end=gsub(".", "", gsub("-|:| ", "", DLDtimeSeq$endTime[j]), fixed=T))
        chunkACC$timestamp <- as.POSIXct(strptime(chunkACC$timestamp, "%Y-%m-%d %H:%M:%OS", tz="UTC"))
        print(paste0("Chunk created."))
        return(chunkACC)
      }))
      return(ACCdwld)
    }
    print(paste0("Downloaded ACC data for bird #", which(allInds_acc$id == ind), ", ", allInds_acc$local_identifier[which(allInds_acc$id == ind)], "."), quote = F)
  }), .parallel=F)
  # Remove potential individuals that returned errors during download
  accLs <- accLs[!vapply(accLs, is.error, logical(1))]
  # Exclude empty elements from the list, bind and save
  accLs <- accLs[which(!sapply(accLs, is.null))]
  if(length(accLs)>0){
    accDf <- do.call(rbind, accLs)
    save(accDf, file="C:/Users/hbronnvik/Documents/golden_eagle_tools/goldenEagles_movebankDownload_onlyAcc.rdata")
  }
}
#stopCluster(cl)
Sys.time() - start_time

#load("goldenEagles_movebankDownload_onlyAcc.rdata")

# split the acc data into files for each individual
#dir.create("accData")
accLs <- split(accDf, accDf$individual_local_identifier)
lapply(1:length(accLs), function(acc){
  save(acc, file=paste0("C:/Users/hbronnvik/Documents/golden_eagle_tools/accData/", unique(acc$individual_local_identifier),"_onlyAcc.RData"))
})

# for (i in 1:length(accLs)) {
#   acc <- accLs[[i]]
#   save(acc, file=paste0("C:/Users/hbronnvik/Documents/golden_eagle_tools/accData/", unique(accLs[[i]]$individual_local_identifier),"_onlyAcc.RData"))
# }

#list segmented GPS files
cfls_ls <- list.files("C:/Users/hbronnvik/Documents/golden_eagle_tools/Golden_eagle_wind_behav_June_2022", full.names = T)
# remove individuals without data
cfls_ls <- pfdp_files[sapply(pfdp_files, file.size) > 2000]
#list acc files
acc_ls <- list.files("C:/Users/hbronnvik/Documents/golden_eagle_tools/accData", full.name=T) 

# extract IDs from each
gpsInds <- gsub(".*/|\\(.*","",cfls_ls)
accInds <- gsub(".*/|\\(.*","",acc_ls)

start_time <- Sys.time()
gpsAcc_ls <- llply(gpsInds, function(ind)try({
  print(ind)
  
  HRdf <- readRDS(grep(ind, cfls_ls, value=T, fixed=T)) #load gps file (named HRdf)
  load(grep(ind, acc_ls, value=T, fixed=T)) #load acc file (named acc)
  # subset only gps data corresponding to thermals
  # for only the thermal data, this takes ~5 hours, but the data here do not have IDs appended
  # only the thermal soaring data
  gps <- HRdf[which(HRdf$thermalClust == "circular" & HRdf$soarClust == "soar"),]
  # Order both datasets by timestamp
  acc <- acc[order(acc$timestamp),]
  gps <- gps[order(gps$timestamp),]
  # Subset only acc locations within the time range of the gps for the interpolation
  acc <- acc[acc$timestamp > min(gps$timestamp) & acc$timestamp < max(gps$timestamp),]
  if(nrow(acc) > 0 & length(grep("accelerations_raw", names(acc)))==1){  
    # Extract column names (sometimes they differ depending on tags)
    axesCol <- grep("acceleration_axes", names(acc), value=T)
    accRawCol <- grep("accelerations_raw", names(acc), value=T)
    sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
    # Exclude ACC data that have only X or Y axis
    acc <- acc[which(!as.character(acc[, axesCol]) %in% c("X","Y")),]
    # Exclude rows with missing acc data
    acc <- acc[which(complete.cases(acc[, accRawCol])),]
    # Associate to each gps location the ACC EVENT ID of the acc values closest in time (within 5 min)
    names(acc)[names(acc) == "event_id"] <- "acc_event_id"
    acc$acc_event_id <- as.character(acc$acc_event_id)
    # nrow(acc)==length(unique(acc$acc_event_id))
    ColsToAssociate <- c("acc_event_id")
    timeTolerance <- 5*60 # 5 mins in seconds
    gps$timestamp <- as.POSIXct(as.character(gps$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    #create the empty columns that I want to fill in during the loop
    #start_time <- Sys.time()
    gpsAcc <- ACCtoGPS_parallel(ACCdata = acc, GPSdata = gps, ColsToAssociate, timeTolerance = timeTolerance)
    #Sys.time() - start_time
    
    # Subset the acc dataset only to the acc_event_id associated to the gps information
    accSub <- acc[which(acc$acc_event_id %in% gpsAcc$acc_event_id),]
    # Calculate acc info (vedba etc) for this subset of acc data (split acc strings and calculate mean and cumulative vedba)
    if(nrow(accSub)>0 & length(unique(accSub[, axesCol]))==1){ #Continue only if number of acc axes doesn't vary within the same individual
      accDf_vedba <- data.frame(acc_event_id=accSub$acc_event_id, n_samples_per_axis=NA, acc_burst_duration_s=NA, 
                                meanVedba=NA, cumVedba=NA, meanODBA=NA, cODBA=NA)
      for(j in 1:nrow(accSub)){
        Naxes <- nchar(as.character(accSub[j, axesCol]))
        accMx <- matrix(as.integer(unlist(strsplit(as.character(accSub[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
        n_samples_per_axis <- nrow(accMx)
        acc_burst_duration_s <- n_samples_per_axis/accSub[j, sampFreqCol]
        if(nchar(accSub[j, axesCol])<3){stop("The ACC data have fewer than 3 axes.")}
        vedba <- sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2)
        ODBA <- (accMx[,1]-mean(accMx[,1])) + (accMx[,2]-mean(accMx[,2])) + (accMx[,3]-mean(accMx[,3]))
        
        accDf_vedba[j, c("n_samples_per_axis", "acc_burst_duration_s", "meanVedba", "cumVedba", "meanODBA", "cODBA")] <- c(n_samples_per_axis, acc_burst_duration_s, 
                                                                                                                                   mean(vedba, na.rm=T), sum(vedba, nrm=T),
                                                                                                                                   mean(ODBA, na.rm=T), sum(ODBA, na.rm=T))
      }
      # Merge the resulting columns (vedba etc) to the gps data based on acc_event_id
      if(nrow(accDf_vedba)>0){
        gpsAcc <- merge(gpsAcc, accDf_vedba, by="acc_event_id", all.x=T)
        return(gpsAcc)
      }
    }
  }
}), .parallel=F)
is.error <- function(x) inherits(x, "try-error")
# Remove potential individuals that returned errors during download
gpsAcc_ls <- gpsAcc_ls[!vapply(gpsAcc_ls, is.error, logical(1))]
# Exclude empty elements from the list, bind and save
gpsAcc_ls <- gpsAcc_ls[which(!sapply(gpsAcc_ls, is.null))]
if(length(gpsAcc_ls)>0){
  gpsAccDf <- rbindlist(gpsAcc_ls, use.names = T)
  save(gpsAccDf, file="C:/Users/hbronnvik/Documents/golden_eagle_tools/goldenEagles_movebankDownload_gps&acc.rdata")
}
Sys.time() - start_time


# sticky end: find mean VeDBA for each thermal for each bird, add PFDP duration, plot

gpsAcc_ls2 <- lapply(gpsAcc_ls, function(x){
  x$ring_id <- as.character(x$ring_id)
  return(x)
})

check <- bind_rows(gpsAcc_ls2)

thermDBA <- gpsAccDf %>% 
  group_by(local_identifier, thermalID) %>% 
  summarize(avg_thermal_VeDBA = mean(meanVedba)) %>% 
  drop_na(thermalID, avg_thermal_VeDBA) %>% 
  summarize(avg_VeDBA = mean(avg_thermal_VeDBA),
            max_VeDBA = max(avg_thermal_VeDBA))

dates <- read.csv("C:/Users/hbronnvik/Documents/golden_eagle_tools/fledging_emig_times.csv")
dates <- dates %>% 
  mutate(pfdp_duration = difftime(emigration_dt, fledging_dt, units = "days")) %>% 
  filter(individual.local.identifier %in% thermDBA$local_identifier)

thermDBA <- thermDBA %>% 
  mutate(pfdp_duration = dates$pfdp_duration[which(unique(dates$individual.local.identifier) == unique(local_identifier))])

ggplot(thermDBA, aes(avg_VeDBA, pfdp_duration)) +
  geom_smooth(method = "lm", color = "black", size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Mean VeDBA in thermals", y = "PFDP duration (days)") +
  theme_classic()

ggplot(thermDBA, aes(max_VeDBA, pfdp_duration)) +
  geom_smooth(method = "lm", color = "black", size = 1.2) +
  geom_point(size = 3) +
  labs(x = "Mean VeDBA in thermals", y = "PFDP duration (days)") +
  theme_classic()


thermDBA <- thermDBA %>% gather(avg_VeDBA, max_VeDBA, key = "measure", value = "VeDBA")

ggplot(thermDBA, aes(VeDBA, pfdp_duration)) +
  geom_smooth(method = "lm", color = "black", size = 1.2) +
  geom_point(size = 3) +
  ylim(0,300) +
  labs(x = "VeDBA in thermals", y = "PFDP duration (days)") +
  theme_classic() +
  facet_wrap(~measure, scales = "free")

## 2. Identify unique thermals

# The method I used relies on T/F, so using the data pared down to thermals (above) failed
# I used the assignments from task 1, pared those down to thermals as well, and appended them
# based on shared ID & timestamp

# load("gps_acc/goldenEagles_movebankDownload_gps&acc_ids.RData")


# measure the time spent in each soaring "event"
# using only the data for thermal soaring
gpsAccDf <- gpsAccDf[which(gpsAccDf$thermalClust == "circular"),]

tempt <- data.frame()

# for each day
for (i in unique(gpsAccDf$id_date)) {
  df <- gpsAccDf[which(gpsAccDf$id_date == i),]
  # and or each thermal soaring event
  circle <- unique(df$thermalID)
  for (k in circle){
    df2 <- df[which(df$thermalID == k),]
    # the amount of time in that event
    df2$circle_dt <- as.numeric(difftime(df2$timestamp[nrow(df2)], df2$timestamp[1], units = "sec"))
    tempt <- rbind(tempt, df2)
    rm(df2)
  }
  print(paste0("Calculated time spent in each thermal soaring event for ", i, "."), quote = F)
  rm(df)
}

gpsAccDf <- tempt; rm(tempt)
#save(gpsAccDf, file = "df_with_soaring_times.RData")

hist(gpsAccDf$circle_dt, breaks = 400)

## 3. Sum ODBA within each thermal

# create a column of unique event IDs
gpsAccDf$thermal_event <- paste0(gpsAccDf$id_date, "_", gpsAccDf$thermalID)

dbaDF <- gpsAccDf %>% 
  # for each ID, day, and thermal with a unique acc_event_id (unique nearest timestamp, unique cumulative ODBA)
  group_by(thermal_event, acc_event_id) %>% 
  # take only the first value (just one mean ODBA per acc event)
  slice(1) %>%
  ungroup() %>% 
  # then within each thermal
  group_by(thermal_event) %>%
  # average the DBA 
  mutate(mean_ODBA_thermal = mean(meanODBA),
         mean_VeDBA_thermal = mean(meanVedba)) %>% 
  # finally, take just the one mean ODBA per thermal
  slice(1) %>% 
  ungroup()


## 4. Find the ratio of average ODBA to time spent 
dbaDF$ratio_odba_sec <- dbaDF$mean_ODBA_thermal/dbaDF$circle_dt
dbaDF$ratio_vedba_sec <- dbaDF$mean_VeDBA_thermal/dbaDF$circle_dt



##############################################################################################










