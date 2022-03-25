## Task3: DBA during thermaling
## Hester Brønnvik
## hbronnvik@ab.mpg.de
## 14.03.2022

## 1. Identify unique thermals
## 2. Calculate ODBA
## 3. Sum ODBA within each thermal
## 4. Find the ratio of ODBA to time spent 
## 5. Plot the ratio against days since fledging

### Get and clean the data
load("loginStored.RData")
eagleStudyId <- 282734839
library(plyr)
library(doParallel)
cl <- makeCluster(detectCores()-3, type='PSOCK')
registerDoParallel(cl)
source("associate_ACCinfoToGPS.R") #For ACCtoGPS function
library(data.table)
library(tidyverse)
is.error <- function(x) inherits(x, "try-error")


cfls_ls <- list.files("thermals/") 
cfls_ls <- paste0("thermals/", cfls_ls)

cfls <- list()

for (i in cfls_ls) {
  load(i)
  # only the thermal soaring data
  HRdf <- HRdf[which(HRdf$thermalClust == "circular" & HRdf$soarClust == "soar"),]
  cfls[[which(cfls_ls == i)]] <- HRdf
  print(paste0("Successfully loaded ",  unique(HRdf$name), "."), quote = F)
  rm(HRdf)
}



# create a unique id for each day
cfls$id_date <- paste(str_sub(cfls$individual_local_identifier, 1, -13), paste(year(cfls$timestamp), month(cfls$timestamp), day(cfls$timestamp), sep = "-"), sep = "_")


### 1. 

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

#dir.create("accData")
accLs <- split(accDf, accDf$individual_local_identifier)
lapply(1:length(accLs), function(acc){
  save(acc, file=paste0("accData/",unique(acc$individual_local_identifier),"_onlyAcc.RData"))
})

#list GPS files
cfls_ls <- list.files("thermals", full.name=T) 
#list acc files
acc_ls <- list.files("accData", full.name=T) 

gpsInds <- gsub(".*/|\\(.*","",cfls_ls)
accInds <- gsub(".*/|\\(.*","",acc_ls)

table(gpsInds %in% accInds)
gpsInds_sub <- gpsInds[gpsInds %in% accInds]

start_time <- Sys.time()
gpsAcc_ls <- llply(gpsInds_sub, function(ind)try({
  print(ind)
  
  load(grep(ind, cfls_ls, value=T, fixed=T)) #load gps file (named HRdf)
  load(grep(ind, acc_ls, value=T, fixed=T)) #load acc file (named acc)
  # subset only gps data corresponding to thermals
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
# Time difference of 5.306388 hours
registerDoSEQ()

### 3. Sum ODBA within each thermal
gpsAccDf$id_date <- paste0(gpsAccDf$local_identifier, "_", month(gpsAccDf$timestamp), "-", day(gpsAccDf$timestamp), "-", year(gpsAccDf$timestamp))
tidDF <- data.frame()

# assign unique IDs to each thermal event
for (i in unique(gpsAccDf$id_date)) {
  df <- gpsAccDf[which(gpsAccDf$id_date == i), ]
  df <- df[order(df$timestamp),]
  df$thermal <- NA
  # where the thermal cluster is thermal soaring, True
  df$thermal[which(df$thermalClust == "circular")] <- T
  # where it is linear or "other", False
  df$thermal[which(df$thermalClust != "circular")] <- F
  # assign a new number to each thermal soaring event 
  df$thermalID <- inverse.rle(within.list(rle(df$thermal), 
                                          values[values] <- seq_along(values[values])))
  df$thermalID[which(df$thermalID == 0)] <- NA
  tidDF <- rbind(tidDF, df)
  rm(df)
  print(paste0("Assigned thermal event IDs for ", i, "."), quote = F)
}

gpsAccDf <- tidDF; rm(tidDF)

# load("gps_acc/goldenEagles_movebankDownload_gps&acc_ids.RData")

# create a column of unique event IDS
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

# while working with the subset of birds that Svea aged:
  
# append the fledging dates
# the file with the fledging and emigration dates
stages <- read.csv("Goldeneagles10_2021.csv", stringsAsFactors = F)
# correct encoding
stages$id <- gsub("\xfc", "ü", stages$id)
stages$id <- gsub("\xf6", "ö", stages$id)
stages <- stages[order(stages$id),]
# get the fledging date column re-arranged
stages$date_fledging <- paste0(stages$date_fledging, " ", stages$time_fledging)
stages$date_fledging <- as.POSIXct(gsub("\\.", "\\-", stages$date_fledging), format = "%d-%m-%Y %H:%M:%OS")
# get the emigration date column re-arranged
stages$date_emigration <- paste0(stages$date_emigration, " ", stages$time_emigration)
stages$date_emigration <- as.POSIXct(gsub("\\.", "\\-", stages$date_emigration), format = "%d-%m-%Y %H:%M:%OS")
# set all the missing emigration dates to today so that it is impossible for a bird to have timestamps > emigration date
# to assign life stage, we look at timestamps < emigration date, and therefore cannot use NA
stages$date_emigration[which(is.na(stages$date_emigration))] <- Sys.time()

# the file containing the names in different formats
load("eagle_names.RData")

# append the fledging dates
odbaDF$fledging_date <- NA
sub_df <- odbaDF[which(odbaDF$local_identifier %in% eagle_names$local_identifier),]
tcfls <- sub_df
sub_df <- data.frame()

for (i in unique(tcfls$local_identifier)) {
  df <- tcfls[which(tcfls$local_identifier == i),]
  df$fledging_date <- stages$date_fledging[which(stages$id == 
                      eagle_names$age_name[which(eagle_names$local_identifier == i)])]
  sub_df <- rbind(sub_df, df)
  
}

rm(tcfls)

# calculate the time difference between fledging and the data
sub_df$dsf <- as.numeric(sub_df$timestamp - sub_df$fledging_date)

# plot
ggplot(sub_df, aes(x = dsf, y = sum_ODBA))+
  geom_point() +
  geom_smooth(method = "lm", se = F, aes(color = sub_df$local_identifier)) +
  labs(x = "Days since fledging", y = "Sum of ODBA in a thermal") +
  scale_x_continuous(breaks = seq.int(0, 1100, by = 200)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))

