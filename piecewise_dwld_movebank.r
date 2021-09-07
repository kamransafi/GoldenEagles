library(move)
library(data.table)
library(plyr)
library(doMC)
registerDoMC(cores=6)
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

ACCpartialDLD <- function(startTime, endTime){
  IndACC <- getMovebank(entity="event", study_id=studyID, sensor_type_id=2365683 , animalName=Inds$IndID[1], login=creds, 
                        attributes=c("timestamp", "eobs_acceleration_axes", "eobs_acceleration_sampling_frequency_per_axis", "eobs_accelerations_raw" ),
                        timestamp_start=gsub(".", "", sub(" ", "", gsub(":", "", gsub("-", "", startTime))), fixed=T), 
                        timestamp_end=gsub(".", "", sub(" ", "", gsub(":", "", gsub("-", "", endTime))), fixed=T))
  IndACC$timestamp <- as.POSIXct(strptime(IndACC$timestamp, "%Y-%m-%d %H:%M:%OS", tz="UTC"))
  return(rbindlist(llply(1:nrow(IndACC), function(x) extractACC(IndACC[x,1], IndACC[x,2], IndACC[x,3], IndACC[x,4]), .parallel=T)))
}
  

study <- "LifeTrack Golden Eagle Alps"
creds <- movebankLogin()
studyID <- getMovebankID(study, creds)
Animals <- getMovebankAnimals(study, creds)
Inds <- data.frame(IndName=Animals$local_identifier, IndID=Animals$individual_id)
Inds <- Inds[!duplicated(Inds),]

for(i in 1:nrow(Inds)){
  SensorIDs <- Animals$sensor_type_id[Animals$individual_id==Inds$IndID[i]]
#  IndGPS <- getMovebankData(study, animalName=Inds$IndName[i], creds, removeDuplicatedTimestamps=T)
  attribs <- getMovebankSensorsAttributes(studyID, creds)
  if(2365683 %in% SensorIDs)
  {
    Trange <- getMovebank(entity="event", study_id=studyID, sensor_type_id=2365683 , animalName=Inds$IndID[i], login=creds, 
                          attributes=c("timestamp"))
    Trange <- as.vector(Trange$timestamp)
    dlsq <- cbind(seq(0, length(Trange), by=10000)+1, c(seq(0, length(Trange), by=10000), length(Trange))[-1])
    DLDtimeSeq <- data.frame(startTime=Trange[dlsq[,1]], endTime=Trange[dlsq[,2]])
    t0 <- Sys.time()
    acc_table <- rbindlist(lapply(1:3, function(x) ACCpartialDLD(DLDtimeSeq[x,1], DLDtimeSeq[x,2]))) #, .parallel = T nrow(DLDtimeSeq)4
    Sys.time()-t0
    }
  
}

