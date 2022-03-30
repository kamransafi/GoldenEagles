### Juvenile golden eagle soaring performance code
### Hester Brønnvik
### 25.03.2022
###


## This file contains code for three tasks: 1) calculating time spent soaring, 2) calculating 
## maximum wind speed in thermals per day, and 3) calculating ODBA per time in each thermal.

## Packages and functions:
library(move)
library(lubridate)
library(data.table)
library(tidyverse)
library(ggpubr)
library(plyr)
library(doParallel)
cl <- makeCluster(detectCores()-3, type='PSOCK')
registerDoParallel(cl)
is.error <- function(x) inherits(x, "try-error")
source("associate_ACCinfoToGPS.R") #For ACCtoGPS function
load("loginStored.RData") # The MoveBank login
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

### Assign life stages to the birds for which we have fledging dates

# the file with the names of the birds in different formats
# load("eagle_names.RData")
# 
# # the birds that have ACC data associated
# objs <- list.files("thermals")#, full.names = T)
# objs <- sub("\\) ", "\\)", objs) # correct Stürfis20
# objs <- str_sub(objs, 1,-48)
# objs <- sub(" ", "\\.", objs)
# objs <- sub("Art.San ", "", objs) # correct Art San Romerio18
# 
# # the file with the fledging and emigration dates (from Svea)
# stages <- read.csv("Goldeneagles10_2021.csv", stringsAsFactors = F)
# stages <- stages[order(stages$id),]
# stages <- stages[which(stages$id %in% objs),]
# # get the fledging date column re-arranged
# stages$date_fledging <- paste0(stages$date_fledging, " ", stages$time_fledging)
# stages$date_fledging <- as.POSIXct(gsub("\\.", "\\-", stages$date_fledging), format = "%d-%m-%Y %H:%M:%OS")
# # get the emigration date column re-arranged
# stages$date_emigration <- paste0(stages$date_emigration, " ", stages$time_emigration)
# stages$date_emigration <- as.POSIXct(gsub("\\.", "\\-", stages$date_emigration), format = "%d-%m-%Y %H:%M:%OS")
# # set all the missing emigration dates to today so that it is impossible for a bird to have timestamps > emigration date
# # to assign life stage, we look at timestamps < emigration date, and therefore cannot use NA
# stages$date_emigration[which(is.na(stages$date_emigration))] <- Sys.time()
# 
# # open the segmented files
# # give each a life_stage column
# # either add data or leave NA
# gps_fls <- list.files("thermals", full.names = T)
# 
# gps_ls <- list()
# 
# for (i in gps_fls) {
#   load(i)
#   if(unique(HRdf$local_identifier) %in% eagle_names$local_identifier){
#     gps_ls[[length(gps_ls) + 1]] <- HRdf
#     print(paste0("Successfully loaded ", unique(HRdf$local_identifier), "."), quote = F)
#   } else{rm(i);rm(HRdf)}
# }
# 
# 
# # loop through each individual and label its data with life stage
# for (i in 1:length(gps_ls)) {
#   # get the segmented GPS data from a single id
#   gps_df <- gps_ls[[i]]
#   # add the name for life stage data to the data frame with gps
#   gps_df$sname <- NA
#   gps_df$sname <- eagle_names$age_name[which(eagle_names$local_identifier == unique(gps_df$local_identifier))]
#   # for each id, match stages to move,
#   gps_df$stage <- NA
#   # then say before fledge, after fledge, and after emigration
#   if (unique(gps_df$local_identifier) %in% eagle_names$local_identifier){ #our_Inds$eobs
#     # timestamps less than date of fledging are classified as 1 (pre-fledging)
#     gps_df$stage[which(gps_df$timestamp < stages$date_fledging[which(stages$id == unique(gps_df$sname))])] <- 1
#     # timestamps greater than date of fledging & less than date of emigration are classfied as 2 (fledgling)
#     gps_df$stage[which(gps_df$timestamp > stages$date_fledging[which(stages$id == unique(gps_df$sname))] &
#                          gps_df$timestamp < stages$date_emigration[which(stages$id == unique(gps_df$sname))])] <- 2
#     # timesamps greater than date of emigration are classfied as 3 (emigrant)
#     gps_df$stage[which(gps_df$timestamp > stages$date_emigration[which(stages$id == unique(gps_df$sname))])] <- 3
#   }
#   # save the file again
#   save(gps_df, file = paste0("gps_age/", unique(gps_df$local_identifier), gsub("-",".", Sys.Date()), "_gps_age.RData"))
#   # signal
#   print(paste0("Added life stage information to ", unique(gps_df$local_identifier), "."), quote = F)
# }
# 
# # finally, check that there are life stage data associated with each individual
# for (i in 1:length(gps_ls)) {
#   df <- gps_ls[[i]]
#   if(NA %in% unique(df$stage)){
#     print(paste0(unique(gps_ls[[i]]$local_identifier), " has no stages."))
#   }
# }

# list the files that have the segmented gps data
cfls_ls <- list.files("thermals", full.names = T)

cfls <- data.frame()

# read in the data
for (i in cfls_ls) {
  load(i)
  gps_df <- gps_df[which(gps_df$stage == 2),] # select only the ones post-fledging
  cfls <- rbind(cfls,gps_df)
  rm(gps_df)
  print(paste0("Successfully loaded ", str_sub(sub("_.RData", "", sub("\\/", "", gsub("gps_age", "", i))), 1, -11), "."), quote = F)
}

# reclassify the thermal soaring events so that there is no circular gliding
cfls$thermalClust <- "other"
cfls$thermalClust[which(cfls$gr.speed >= 2 & cfls$soarClust == "soar" & cfls$turnAngle_smooth >= 300)] <- "circular"
cfls$thermalClust[which(cfls$gr.speed >= 2 & cfls$soarClust=="soar" & cfls$thermalClust != "circular")] <- "linear"
cfls$thermalClust <- factor(cfls$thermalClust, levels=c("circular","linear","other"))

## 1. Calculate the total time spent flying in these data (gliding, soaring, flapping)

# select only flight data
cfls <- cfls[which(cfls$gr.speed >= 2),]
# create a year-month-day column
cfls$ymd <- format(cfls$timestamp, format = "%Y-%m-%d")
# create a unique label for each day for each individual
cfls$id_date <- paste0(cfls$local_identifier," ", cfls$ymd)

cfls_t <- data.frame()

# calculate the time spent flying on each day in seconds
for (i in unique(cfls$id_date)) {
  df <- cfls[which(cfls$id_date == i),]
  df <- df[order(df$timestamp),]
  df$flight_time <- as.numeric(difftime(df$timestamp[nrow(df)], df$timestamp[1], units = "sec"))
  cfls_t <- rbind(cfls_t, df)
  rm(df)
  print(paste0("Calculated flight time per day for ", i, "."), quote = F)
}

cfls <- cfls_t; rm(cfls_t)


## 2. give IDs to each unique soaring event (Martina probably has faster code for this)

cfls2 <- data.frame()

# each orographic event
for (i in id_date) {
  df <- cfls[which(cfls$id_date == i), ]
  df <- df[order(df$timestamp),]
  df$linear <- NA
  # where the thermal cluster is linear soaring, True
  df$linear[which(df$thermalClust == "linear")] <- T
  # where it is circular or "other", False
  df$linear[which(df$thermalClust != "linear")] <- F
  # assign a new number to each slope soaring event 
  df$linearID <- inverse.rle(within.list(rle(df$linear), 
                                         values[values] <- seq_along(values[values])))
  df$linearID[which(df$linearID == 0)] <- NA
  cfls2 <- rbind(cfls2, df)
  rm(df)
  print(paste0("Assigned orographic event IDs for ", i, "."), quote = F)
}
cfls <- cfls2; rm(cfls2)

# a data frame to hold the IDs
cfls3 <- data.frame()


# each thermal event
for (i in unique(cfls$id_date)) {
  df <- cfls[which(cfls$id_date == i), ]
  df <- df[order(df$timestamp),]
  df$thermal <- NA
  # where the thermal cluster is thermal soaring, True
  df$thermal[which(df$thermalClust == "circular")] <- T
  # where it is linear or "other", False
  df$thermal[which(df$thermalClust != "circular")] <- F
  # assign a new number to each slope soaring event 
  df$thermalID <- inverse.rle(within.list(rle(df$thermal), 
                                          values[values] <- seq_along(values[values])))
  df$thermalID[which(df$thermalID == 0)] <- NA
  cfls3 <- rbind(cfls3, df)
  rm(df)
  print(paste0("Assigned thermal event IDs for ", i, "."), quote = F)
}
cfls <- cfls3; rm(cfls3)

## 3. measure the time spent in each soaring "event"

# pare down the data to relieve computing stress
# using only the data for linear or thermal soaring
cfls <- cfls[which(cfls$thermalClust != "other"),]

templ <- data.frame()
tempt <- data.frame()

# for each day
for (i in id_date) {
  df <- cfls[which(cfls$id_date == i),]
  # and each slope soaring event (as identified above)
  line <- unique(df$linearID)
  for (j in line) {
    df2 <- df[which(df$linearID == j),]
    # the last timestamp minus the first
    df2$line_dt <- as.numeric(difftime(df2$timestamp[nrow(df2)], df2$timestamp[1], units = "sec"))
    # bound to a new data frame
    templ <- rbind(templ, df2)
    rm(df2)
  }
  print(paste0("Calculated time spent in each linear soaring event for ", i, "."), quote = F)
  # and the same for each thermal soaring event
  circle <- unique(df$thermalID)
  for (k in circle){
    df2 <- df[which(df$thermalID == k),]
    df2$circle_dt <- as.numeric(difftime(df2$timestamp[nrow(df2)], df2$timestamp[1], units = "sec"))
    tempt <- rbind(tempt, df2)
    rm(df2)
  }
  print(paste0("Calculated time spent in each thermal soaring event for ", i, "."), quote = F)
  rm(df)
}

## 4. Sum these over each day

thermal_times <- data.frame()

# the thermal events in a day for an individual
for (i in unique(tempt$id_date)) {
  # single individual and day
  df <- tempt[which(tempt$id_date == i),]
  # its thermal events
  circle <- unique(df$thermalID)
  all_obs <- data.frame()
  for (j in circle) {
    # single thermal
    event <- df[which(df$thermalID == j),]
    # the first observation of it (all observations hold the same amount of time)
    obs <- event[1,]
    # added to the first observation of all other thermals on that day
    all_obs <- rbind(all_obs, obs)
  }
  # sum the time spent on thermals per day and add it to the id_date df
  df$tpd <- sum(all_obs$circle_dt)
  thermal_times <- rbind(thermal_times, df)
  print(paste0("Calculated time spent in thermals for ", i, "."), quote = F)
}

linear_times <- data.frame()

for (i in unique(templ$id_date)) {
  # single individual and day
  df <- templ[which(templ$id_date == i),]
  # its thermal events
  line <- unique(df$linearID)
  all_obs <- data.frame()
  for (j in line) {
    # single thermal
    event <- df[which(df$linearID == j),]
    # the first observation of it (all observations hold the same amount of time)
    obs <- event[1,]
    # added to the first observation of all other thermals on that day
    all_obs <- rbind(all_obs, obs)
  }
  # sum the time spent on thermals per day and add it to the id_date df
  df$tpd <- sum(all_obs$line_dt)
  linear_times <- rbind(linear_times, df)
  print(paste0("Calculated time spent linear soaring for ", i, "."), quote = F)
}

rm(templ);rm(tempt)

thermal_times <- thermal_times %>% 
  group_by(id_date) %>% 
  slice(1)
linear_times <- linear_times %>% 
  group_by(id_date) %>% 
  slice(1)

## 5. Calculate the ratio of soaring time to flight time

# calculate the ratio of soaring to total time spent flying per day
thermal_times$ttf <- thermal_times$tpd/thermal_times$flight_time
linear_times$ltf <- linear_times$tpd/linear_times$flight_time


# if fledging dates are appended in a column called "fledging_date",
# calculate the number of days between each observation and fledging
thermal_times$dsf <- as.numeric(as.Date(thermal_times$ymd) - as.Date(paste(year(thermal_times$fledging_date), month(thermal_times$fledging_date), day(thermal_times$fledging_date), sep = "-")))
linear_times$dsf <- as.numeric(as.Date(linear_times$ymd) - as.Date(paste(year(linear_times$fledging_date), month(linear_times$fledging_date), day(linear_times$fledging_date), sep = "-")))

## 6. Plot the ratio of soaring:flight over days since fledging

thermal_plot <- ggplot(thermal_times, aes(x = dsf, y = ttf)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Relative time spent thermal soaring (s)") +
  scale_x_continuous(breaks = seq.int(0, 500, by = 50)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
linear_plot <- ggplot(linear_times, aes(x = dsf, y = ltf)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Relative time spent slope soaring (s)") + 
  scale_x_continuous(breaks = seq.int(0, 500, by = 50)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
ggarrange(thermal_plot,linear_plot, common.legend = F)

###########################################################################################

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
allInds <- allInds[which(allInds$local_identifier %in% unique(cfls$individual_local_identifier)),]
# Keep only individuals that have ACC information and were not "tests"
allInds_acc <- allInds[grep("GPS.*Acceleration", allInds$sensor_type_ids),]
allInds_acc <- allInds_acc[which(!allInds_acc$local_identifier %in% grep("test|Test", allInds_acc$local_identifier, value=T)),]
# For the remaining individuals, download the ACC data one by one
clusterExport(cl, varlist=c("cfls", "eagleStudyId", "loginStored", "allInds_acc", "ACCtoGPS", "is.error")) 

start_time <- Sys.time()
if(nrow(allInds_acc) > 0){
  indInfo <- data.frame(indName=as.character(allInds_acc$local_identifier[!allInds_acc$local_identifier%in%c(NA,"")]),
                        indId=as.numeric(allInds_acc$id[!allInds_acc$local_identifier%in%c(NA,"")]))
  accLs <- llply(indInfo$indId, function(ind)try({  #with llply
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
    save(accDf, file="goldenEagles_movebankDownload_onlyAcc.rdata")
  }
}
Sys.time() - start_time

#load("goldenEagles_movebankDownload_onlyAcc.rdata")

# split the acc data into files for each individual
#dir.create("accData")
accLs <- split(accDf, accDf$individual_local_identifier)
lapply(1:length(accLs), function(acc){
  save(acc, file=paste0("accData/",unique(acc$individual_local_identifier),"_onlyAcc.RData"))
})

#list segmented GPS files
cfls_ls <- list.files("thermals", full.name=T) 
#list acc files
acc_ls <- list.files("accData", full.name=T) 

# extract IDs from each
gpsInds <- gsub(".*/|\\(.*","",cfls_ls)
accInds <- gsub(".*/|\\(.*","",acc_ls)

start_time <- Sys.time()
gpsAcc_ls <- llply(gpsInds, function(ind)try({
  print(ind)
  
  load(grep(ind, cfls_ls, value=T, fixed=T)) #load gps file (named HRdf)
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
        
        accDf_vedba[j, c("n_samples_per_axis", "acc_burst_duration_s", "meanVedba", "cumVedba", "ODBA", "meanODBA", "cODBA")] <- c(n_samples_per_axis, acc_burst_duration_s, 
                                                                                                                                   mean(vedba, na.rm=T), sum(vedba, nrm=T), ODBA,
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
  gpsAccDf <- rbindlist(gpsAcc_ls)
  save(gpsAccDf, file="gps_acc/goldenEagles_movebankDownload_gps&acc.rdata")
}
Sys.time() - start_time


## 2. Identify unique thermals

# The method I used relies on T/F, so using the data pared down to thermals (above) failed
# I used the assignments from task 1, pared those down to thermals as well, and appended them
# based on shared ID & timestamp

# load("gps_acc/goldenEagles_movebankDownload_gps&acc_ids.RData")


# measure the time spent in each soaring "event"
# using only the data for thermal soaring
cfls <- cfls[which(cfls$thermalClust == "circular"),]

tempt <- data.frame()

# for each day
for (i in id_date) {
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

## 3. Sum ODBA within each thermal

# create a column of unique event IDs
gpsAccDf$thermal_event <- paste0(gpsAccDf$id_date, "_", gpsAccDf$thermalID)

odbaDF <- gpsAccDf %>% 
  # for each ID, day, and thermal with a unique acc_event_id (unique nearest timestamp, unique cumulative ODBA)
  group_by(thermal_event, acc_event_id) %>% 
  # take only the first value (just one cumulative ODBA)
  slice(1) %>%
  ungroup() %>% 
  # then within each thermal
  group_by(thermal_event) %>%
  # sum the DBA 
  mutate(sum_ODBA = sum(cODBA)) %>% 
  # finally, take just the one sum ODBA per thermal
  slice(1)
ungroup()


## 4. Find the ratio of ODBA to time spent 
odbaDF$ratio_dba_sec <- odbaDF$sumODBA/odbaDF$circle_dt


## 5. Plot the ratio against days since fledging

# # while working with the subset of birds that Svea aged:
# 
# # append the fledging dates
# # the file with the fledging and emigration dates
# stages <- read.csv("Goldeneagles10_2021.csv", stringsAsFactors = F)
# # correct encoding
# stages$id <- gsub("\xfc", "ü", stages$id)
# stages$id <- gsub("\xf6", "ö", stages$id)
# stages <- stages[order(stages$id),]
# # get the fledging date column re-arranged
# stages$date_fledging <- paste0(stages$date_fledging, " ", stages$time_fledging)
# stages$date_fledging <- as.POSIXct(gsub("\\.", "\\-", stages$date_fledging), format = "%d-%m-%Y %H:%M:%OS")
# # get the emigration date column re-arranged
# stages$date_emigration <- paste0(stages$date_emigration, " ", stages$time_emigration)
# stages$date_emigration <- as.POSIXct(gsub("\\.", "\\-", stages$date_emigration), format = "%d-%m-%Y %H:%M:%OS")
# # set all the missing emigration dates to today so that it is impossible for a bird to have timestamps > emigration date
# # to assign life stage, we look at timestamps < emigration date, and therefore cannot use NA
# stages$date_emigration[which(is.na(stages$date_emigration))] <- Sys.time()
# 
# # the file containing the names in different formats
# load("eagle_names.RData")
# 
# # append the fledging dates
# odbaDF$fledging_date <- NA
# sub_df <- odbaDF[which(odbaDF$local_identifier %in% eagle_names$local_identifier),]
# tcfls <- sub_df
# sub_df <- data.frame()
# 
# for (i in unique(tcfls$local_identifier)) {
#   df <- tcfls[which(tcfls$local_identifier == i),]
#   df$fledging_date <- stages$date_fledging[which(stages$id == 
#                                                    eagle_names$age_name[which(eagle_names$local_identifier == i)])]
#   sub_df <- rbind(sub_df, df)
#   
# }
# 
# rm(tcfls)

ggplot(sub_df, aes(x = dsf, y = ratio_dba_sec))+
  geom_point() +
  geom_smooth(method = "lm", se = F, aes(color = sub_df$local_identifier)) +
  labs(x = "Days since fledging", y = "ODBA per time") +
  scale_x_continuous(breaks = seq.int(0, 1100, by = 200)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
