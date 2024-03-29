## Download large amounts of ACC by getting the data for each individual at a time
## Check for your system whether doMC works. Otherwise can be replaced eventually with foreach
## extract ACC function extracts acc data from the file and accounts for different sampling frequencies and 
## number of axis. This is a crude function and should be replaced by functions from the moveACC
## package. 

library(move)
library(data.table)
library(foreach)
library(doParallel)
register(cores=6)
options(digits.secs=4)


extractACC <- function(timestamp, eobs_acceleration_axes, eobs_acceleration_sampling_frequency_per_axis, eobs_accelerations_raw){
  acc <- matrix(as.integer(unlist(strsplit(eobs_accelerations_raw, " "))), ncol=nchar(eobs_acceleration_axes), byrow = T)
  acc <- data.table(acc)
  if(eobs_acceleration_axes=="X"){
    names(acc) <- "accX"
  }
  if(eobs_acceleration_axes=="Y"){
    names(acc) <- "accY"
  }
  if(eobs_acceleration_axes=="Z"){
    names(acc) <- "accZ"
  }
  if(eobs_acceleration_axes=="XY"){
    names(acc) <- c("accX", "accY")
  }
  if(eobs_acceleration_axes=="XZ"){
    names(acc) <- c("accX", "accZ")
  }
  if(eobs_acceleration_axes=="YZ"){
    names(acc) <- c("accY", "accZ")
  }
  if(eobs_acceleration_axes=="XYZ"){
    names(acc) <- c("accX", "accY", "accZ")
  }
  timestamp <- timestamp + seq(0, nrow(acc)-1)/eobs_acceleration_sampling_frequency_per_axis
  Sp10thS <- as.integer(10*eobs_acceleration_sampling_frequency_per_axis)
  return(cbind(timestamp=timestamp, acc, Sp10thS))
}

ACCpartialDLD <- function(animalName){
  # get the time range for ACC data
  Trange <- getMovebank(entity="event", study_id=studyID, sensor_type_id=2365683 , animalName=animalName, login=creds, 
                        attributes=c("timestamp"))
  Trange <- as.vector(Trange$timestamp)
  # split the download in chuncks of 10k lines per download
  dlsq <- cbind(seq(0, length(Trange), by=10000)+1, c(seq(0, length(Trange), by=10000), length(Trange))[-1])
  # create start and end timestamps for download
  DLDtimeSeq <- data.frame(startTime=Trange[dlsq[,1]], endTime=Trange[dlsq[,2]])
  #download the chunks in parallel #nrow(DLDtimeSeq)
  ACCdld <- foreach(j=1:3, .combine="rbind") %do% {
    IndACC <- getMovebank(entity="event", study_id=studyID, sensor_type_id=2365683 , animalName=Inds$IndID[1], login=creds, 
                          attributes=c("timestamp", "eobs_acceleration_axes", "eobs_acceleration_sampling_frequency_per_axis", "eobs_accelerations_raw" ),
                          timestamp_start=gsub(".", "", sub(" ", "", gsub(":", "", gsub("-", "", DLDtimeSeq$startTime[j]))), fixed=T), 
                          timestamp_end=gsub(".", "", sub(" ", "", gsub(":", "", gsub("-", "", DLDtimeSeq$endTime[j]))), fixed=T))
    IndACC$timestamp <- as.POSIXct(strptime(IndACC$timestamp, "%Y-%m-%d %H:%M:%OS", tz="UTC"))
    foreach(x=1:nrow(IndACC), .combine="rbind") %dopar% (extractACC(IndACC[x,1], IndACC[x,2], IndACC[x,3], IndACC[x,4]))
  }
  return(ACCdld)
}


study <- "LifeTrack Golden Eagle Alps"
creds <- movebankLogin()
studyID <- getMovebankID(study, creds)
Animals <- getMovebankAnimals(study, creds)
Inds <- data.frame(IndName=Animals$local_identifier, IndID=Animals$individual_id)
Inds <- Inds[!duplicated(Inds),]
data <- NULL
for(i in 1:nrow(Inds)){
  SensorIDs <- Animals$sensor_type_id[Animals$individual_id==Inds$IndID[i]]
  IndGPS <- getMovebankData(study, animalName=Inds$IndName[i], creds, removeDuplicatedTimestamps=T)
  attribs <- getMovebankSensorsAttributes(studyID, creds)
  if(2365683 %in% SensorIDs)
  {
    IndACC <- ACCpartialDLD(animalName=Inds$IndID[i])
  }
  data[[Inds$IndName[i]]] <- list(IndGPS, IndACC)
}

