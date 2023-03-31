
## Script written by Martina in March 2023
## ACC data downloaded on March 23rd

#_________________________
## To run the same script on the cluster do the following:

# # After accessing draco, create an "eagles" folder in the cluster:
# mkdir /draco/ptmp/mscacco/eagles/
# # And copy the two necessary files to run the script (individual names and ACCToGps function (type the following in terminal where I am not logged in to draco):
# scp /home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/ContributionToElhamsPaper/individualNames_accIndsElham.RData mscacco@draco.mpcdf.mpg.de:/draco/ptmp/mscacco/eagles/
# scp /home/mscacco/ownCloud/Martina/PHD/R_functions/functionAccToGps_March2022.R mscacco@draco.mpcdf.mpg.de:/draco/ptmp/mscacco/eagles/
# 
# # Once enetering R on cluster, type this to install the moveACC package
# library('devtools')
# install_git('https://gitlab.com/anneks/moveACC.git', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes=T)

  
#_________________________________
## Download ACC data from Movebank ----

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
pathToHD <- "/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/"

library(move)
creds <- movebankLogin("TeamWikelski")

studyId <- getMovebankID("LifeTrack Golden Eagle Alps", creds)

# List names of the individuals included in the study and filter only those having ACC data
allInds <- getMovebank("individual", login=creds, study_id=studyId)
allInds <- allInds[grep("acceleration|gps", allInds$sensor_type_ids),]
# Subset only individuals with ACC that also have high resolution GPS (look at classified data)
gpsFls <- list.files("classifiedData", pattern="classifiedBursts_df")
table(gsub("_classifiedBursts_df.rds", "", gpsFls) %in% allInds$local_identifier)
table(allInds$local_identifier %in% gsub("_classifiedBursts_df.rds", "", gpsFls))
indsToDownload <- allInds$local_identifier[allInds$local_identifier %in% gsub("_classifiedBursts_df.rds", "", gpsFls)]

# # Import names of the individual included in Elham project
# load("individualNames_accIndsElham.RData") #object localIdentifier_accInds
# all(localIdentifier_accInds %in% allInds$local_identifier)
# elhamInds <- readRDS("burst_sample_sizes_1122.rds")
# all(elhamInds$local_identifier %in% allInds$local_identifier)
# indsToDownload <- elhamInds$local_identifier[!elhamInds$local_identifier %in% localIdentifier_accInds]

# # Import fledging and dispersal dates
# dispDate <- readRDS("fledging_emigration_timing.rds")
# table(dispDate$dispersal_day > yday(dispDate$fledging_dt))

## Download ACC on external hard drive
resLs <- lapply(indsToDownload, function(ind)try({
  print(ind)
  acc <- getMovebankNonLocationData(study=studyId, sensorID="Acceleration", login=creds, 
                                    animalName=ind)
  if(nrow(acc)>0){
    saveRDS(acc, file=paste0(pathToHD, ind, "_onlyAcc.rds"))
  }
}))
# for some the connection to movebank timed out, we find them and download them again
indsToDownload <- indsToDownload[grep("Error", resLs)]

#_____________________________________________________________
## Import and standardise ACC data to same burst duration ----

library(data.table)
library(plyr)
library(doParallel)
doParallel::registerDoParallel(6) 

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
pathToHD <- "/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/"

fls <- list.files(pathToHD, pattern="onlyAcc.rds", full.names = T)
# done <- list.files(pathToHD, pattern="equalDurationBursts", full.names = T)
# table(fls %in% gsub("_equalDurationBursts","",done))
# fls <- fls[!fls %in% gsub("_equalDurationBursts","",done)]

# Set decimals seconds (important for creating new acc timestamps after splitting the long bursts)
options(digits.secs=2)

lapply(fls, function(f){
  print(f)
  acc <- readRDS(f)
  animalName <- unique(acc[,grep("individual_local_identifier|individual.local.identifier", names(acc), value=T)])
  
  axesCol <- grep("acceleration_axes|acceleration.axes", names(acc), value=T)
  rawAccCol <- grep("accelerations_raw|accelerations.raw", names(acc), value=T)
  freqCol <- grep("acceleration_sampling_frequency_per_axis|acceleration.sampling.frequency.per.axis", names(acc), value=T)
  eventCol <- grep("event_id|event.id", names(acc), value=T)
  
  # transform factor columns into characters
  idx <- sapply(acc, is.factor)
  acc[idx] <- lapply(acc[idx], as.character)
  # Make sure timestamp column is in posixct
  acc$timestamp <- as.POSIXct(acc$timestamp, format="%Y-%m-%d %H:%M:%OS")
  # Calculate infos about sampling schedule
  acc$numberSamplesPerAxis <- sapply(strsplit(acc[,rawAccCol], " "), length)/nchar(acc[,axesCol])
  acc$burstDurationSecs <- acc$numberSamplesPerAxis/acc$eobs_acceleration_sampling_frequency_per_axis
  # There are two settings, 1.2 sec at 20 Hz and 7.9 sec at 33.33 Hz
  shortBursts <- acc[acc$burstDurationSecs < 2,]
  longBursts <- acc[acc$burstDurationSecs > 2,]
  #Split the long bursts in 1.2 secs chunks to make them more comparable
  newDurationSec <- 1.2
  #newWishedFreq <- 20
  if(nrow(longBursts) > 0){
    longBursts_cut <- as.data.frame(rbindlist(llply(1:nrow(longBursts), function(i){
      nAxes <- nchar(longBursts[i,axesCol])
      samplFreq <- longBursts[i,freqCol]
      rawSamples <- as.numeric(unlist(strsplit(longBursts[i, rawAccCol], " ")))
      tot_nSamplesInBurst <- length(rawSamples)
      burstDurationSec <- tot_nSamplesInBurst/samplFreq/nAxes
      new_nSamplesPerBurst <- round(tot_nSamplesInBurst/burstDurationSec) * newDurationSec #number of samples in 1.2 s bursts (same as the short bursts)
      # we could think of vary newDurationSec and check which duration, multiplied by the above, returns a number that is multiple of 3 (for which number%%3==0) in a separate function "findPossibleDuration"
      if(new_nSamplesPerBurst%%nAxes != 0){stop("The number of samples in the burst is not a multiple of the number of axes. A burst cannot have a variable number of samples per axis.")}
      nNewShortBursts <- floor(tot_nSamplesInBurst/new_nSamplesPerBurst) #number of new shorter bursts we are going to create from the original long burst (we take floor because the extra observations that are not enough for a new burst gets excluded)
      rawSamples_splitNewBursts <- split(rawSamples[1:(new_nSamplesPerBurst*nNewShortBursts)], rep(1:nNewShortBursts, each=new_nSamplesPerBurst))
      # Potentially here reduce frequency, function "findPossibleFrequency"
      # reduceBy <- round(samplFreq/newWishedFreq) #first it has to be split per axis, and then for each axis take every second value
      # lapply(rawSamples_splitNewBursts, function(b){
      # 
      #})
      rawSamples_splitNewBursts_format <- unlist(lapply(rawSamples_splitNewBursts, paste, collapse=" ")) #collapse them into the typical eobs format
      names(rawSamples_splitNewBursts_format) <- NULL
      # Replicate the columns of this burst as many times as the new splitted bursts
      newDf_splitBurst <- do.call(rbind,lapply(1:nNewShortBursts, function(x) return(longBursts[i,]))) #return the exact same line 6 times
      newDf_splitBurst[,rawAccCol] <- NA
      # replace the old duration and samples per axis the new values
      newDf_splitBurst$numberSamplesPerAxis <- new_nSamplesPerBurst/nAxes
      newDf_splitBurst$burstDurationSecs <- newDurationSec
      # assign new timestamp based on new burst duration
      for(j in 2:nNewShortBursts){(newDf_splitBurst$timestamp[j]=newDf_splitBurst$timestamp[j-1]+newDurationSec)}
      # assign a new acc event id to the split bursts
      for(k in 1:nNewShortBursts){newDf_splitBurst[k,eventCol] <- as.numeric(paste0(newDf_splitBurst[k,eventCol],k))}
      # replace the old raw data (all one burst) with the new split bursts
      newDf_splitBurst[,rawAccCol] <- as.vector(rawSamples_splitNewBursts_format)
      return(newDf_splitBurst)
    }, .parallel=T)))
    
    acc_format <- rbind(shortBursts, longBursts_cut)
    acc_format <- acc_format[order(acc_format$timestamp),]
  }else if(nrow(longBursts)==0){
    acc_format <- shortBursts
  }
  # table(acc_format$burstDurationSecs)
  # table(acc_format[,freqCol])
  #nrow(acc_format)==nrow(shortBursts)+(nrow(longBursts)*6)
  saveRDS(acc_format, file=paste0(pathToHD, animalName, "_onlyAcc_equalDurationBursts.rds"))
})


#____________________________________________
## Classify flapping on standardised ACC ----

# library('devtools')
# install_git('https://gitlab.com/anneks/moveACC.git', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes=T)
library(moveACC)
library(ggplot2)
library(plyr)
library(doParallel)
doParallel::registerDoParallel(6) 

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
pathToHD <- "/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/"

fls <- list.files(pathToHD, pattern="equalDurationBursts", full.names = T)
# done <- list.files(pathToHD, pattern="flappingClassification_DfPerBurst", full.names = T)
# table(!gsub("_onlyAcc_equalDurationBursts","",fls) %in% gsub("_flappingClassification_DfPerBurst","",done))
# fls <- fls[!gsub("_onlyAcc_equalDurationBursts","",fls) %in% gsub("_flappingClassification_DfPerBurst","",done)]

llply(fls, function(f){
  print(f)
  acc <- readRDS(f)
  # Speed test
  # t=Sys.time()
  # waveDf <- ACCwave(acc[1:50000,],transformedData=F)
  # Sys.time()-t
  
  # Calculate amplitude and wing beat frequency using FFT
  waveDf <- ACCwave(acc,transformedData=F)
  
  # clusterPlot(waveDf, forclustering= c("amplitude","odbaAvg"), cluster=T)
  # ggsave(filename=paste0("accData_flapping/",unique(waveDf$individualID),"_flappingClustMoveAcc.png"), device="png", dpi="print", height=12, width=15, units = "cm")
  # head(waveDf)
  # table(waveDf$burstDurationSecs)
  # clusterPlot(waveDf, cluster=F)
  # hist(waveDf$beatsSec, breaks="FD")
  # hist(waveDf$amplitude, breaks="FD")
  # wingBeatsPlot(dfw=waveDf, forclustering= c("amplitude","odbaAvg"))
  # wingBeatsHist(dfw=waveDf, forclustering= c("amplitude","odbaAvg"))
  
  # Classify flapping based on wind beat frequency and additionally amplitude and odba, and save dataset
  flapDf <- WingBeatsSelection(waveDf, forclustering= c("amplitude","odbaAvg"), minbeat=0, maxbeat=max(waveDf$beatsSec))
  saveRDS(flapDf, file=paste0(pathToHD,"classifiedData/",unique(waveDf$individualID),"_flappingClassification_DfPerBurst.rds"))

}, .parallel=T)    

#copy file to Martina's computer
pathToM <- "/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/" 
fls <- list.files(pathToHD, pattern="flappingClassification", full.name=F)
file.copy(from=paste0(pathToHD,fls), to=paste0(pathToM,"classifiedData/",fls), overwrite = T)

#________________
## Some plotting:

library(ggplot2)
library(gridExtra)

dir.create(paste0(pathToM,"ACCsegmentationPlots"))
dir.create(paste0(pathToHD,"ACCsegmentationPlots"))

accFls <- list.files(pathToHD, pattern="equalDurationBursts", full.names = T)
flapFls <- list.files(paste0(pathToHD,"classifiedData"), pattern="flappingClassification_DfPerBurst", full.names = T)

set.seed(1512)
for(i in sample(length(accFls), 15)){
  print(i)
  
  acc <- readRDS(accFls[i])
  flapDf <- readRDS(flapFls[i])
  
  # clusterPlot(flapDf, forclustering= c("amplitude","odbaAvg"), cluster=T)
  # ggsave(filename=paste0("accData_flapping/",unique(flapDf$individualID),"_flappingClustMoveAcc.png"), device="png", dpi="print", height=12, width=15, units = "cm")
  
  # table(flapDf$behavior)
  flaps <- flapDf[which(flapDf$behavior=="Flapping"),]
  noFlaps <- flapDf[which(flapDf$behavior!="Flapping"),]
  # Plot a sample of flapping events
  #pdf(paste0("accData_flapping/",unique(flapDf$individualID),"_flappingEvents.pdf"), height=12, width=20)
  if(nrow(flaps)>80){
    rowix <- sample(1:nrow(flaps), 45)
    rowix <- split(rowix, ceiling(seq_along(rowix)/15))
    pFlap <- lapply(rowix, function(x) PlotAccData(acc, bursts=flaps[x,"burstID"], interactivePlot = F))
  }
  #dev.off()
  ggsave(filename = paste0(pathToHD,"ACCsegmentationPlots/",unique(flapDf$individualID),"_flappingEvents.pdf"), 
    plot = marrangeGrob(pFlap, nrow=1, ncol=1), 
    width = 20, height = 12)
  # Plot a sample of non-flapping events
  rowix <- sample(1:nrow(noFlaps), 75)
  rowix <- split(rowix, ceiling(seq_along(rowix)/15))
  pNof <- lapply(rowix, function(x) PlotAccData(acc, bursts=noFlaps[x,"burstID"], interactivePlot = F))
  ggsave(filename = paste0(pathToHD,"ACCsegmentationPlots/",unique(flapDf$individualID),"_non-flappingEvents.pdf"), 
         plot = marrangeGrob(pNof, nrow=1, ncol=1), 
         width = 20, height = 12)
}

  # Plot some non flapping events
  # q <- quantile(noFlaps$amplitude, seq(0,1, 0.001))
  # plot(q)
  # e <- ecdf(noFlaps$amplitude)
  # plot(e)
  # e(100)
  # plot(noFlaps$amplitude, noFlaps$odbaMedian)

#________________________________________________
## Associate classified ACC data to GPS data ----

## Note that the GPS classification was run on Hester's computer.
## So the files with pattern "classifiedBursts_df" were sent by Hester

# Import the function
library(data.table)
library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(5) 

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
#setwd("/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/")
source("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/ContributionToElhamsPaper/GoldenEagles_ElHeMa_scripts/WP3_Soaring_Ontogeny/MS1_soaring_skills/functionAccToGps_March2022.R")

accFls <- list.files("classifiedData", pattern="flappingClassification_DfPerBurst", full.names = T)
gpsFls <- list.files("classifiedData", pattern="classifiedBursts_df", full.names = T)

accNames <- gsub(".*/|_flapping.*.","", accFls)
gpsNames <- gsub(".*/|_classif.*.","", gpsFls)

table(gpsNames %in% accNames)
table(accNames %in% gpsNames)
gpsNames[!gpsNames %in% accNames]
indDone <- gsub(".*/|_gps.*.","", list.files("classifiedData", pattern="gps&acc", full.names = T))
accNames_toDo <- accNames[!accNames %in% indDone]

llply(accNames_toDo, function(ind){
  print(ind)
  accDf <- readRDS(grep(ind, accFls, fixed=T, value=T))
  gpsDf <- readRDS(grep(ind, gpsFls, fixed=T, value=T))
  # Create new ACC event column (using the tag id and the burst id assigned during the acc classification)
  accDf$new_accEventId <- as.numeric(paste0(accDf$tagID, accDf$burstID))
  if(nrow(accDf)==length(unique(accDf$new_accEventId))){
    #format timestamps
    # gpsDf$timestamp <- as.POSIXct(as.character(gpsDf$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    # accDf$timestamp <- as.POSIXct(as.character(accDf$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    # Order both datasets by timestamp
    accDf <- accDf[order(accDf$timestamp),]
    gpsDf <- gpsDf[order(gpsDf$timestamp),]
    # Extract column names to associate to the GPS dataset
    accColsToAssociate <- c("beatsSec","amplitude","odbaAvg","odbaMedian","accAxes","numberSamplesPerAxis","burstDurationSecs","samplingFreqPerAxis","behavior")
    # Define time tolerance (e.g. 5 mins in seconds) to look for the closest ACC information and associate it to the gps data
    timeTolerance <- 5*60 
    accGps <- ACCtoGPS(GPSdata = gpsDf, ACCdata = accDf, accEventCol="new_accEventId",
                       timeTolerance = timeTolerance,
                       ColsToAssociate = accColsToAssociate) #columns of the acc data that you want to associate to the gps data
    # Save the resulting dataset (merging gps and acc) for each individual
    saveRDS(accGps, file=paste0("classifiedData/",ind,"_gps&acc_behavClass.rds"))
  }else{warning("ACC event id are not unique!")}
  #return(accGps)
}, .parallel=T)

# Bind all individuals and save final dataset
accGps_ls <- lapply(list.files("classifiedData", pattern="gps&acc"), readRDS)
accGps_class <- as.data.frame(rbindlist(accGps_ls))
saveRDS(accGps_class, file=paste0("classifiedData/goldenEagles_allIndividuals_gps&acc_behavClass.rds"))




## Optional (not done yet):
#__________________________________________________________________
## Filter out the time spent in the nest (take time after fledging) ----

# library(data.table)
# 
# setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/ContributionToElhamsPaper")
# # setwd("/draco/ptmp/mscacco/eagles/") #in cluster
# all_GpsAcc <- readRDS("accGps_behavClass/goldenEagles_allIndividuals_gps&acc_behavClass.rds")
# 
# # Import the file with the fledging date per individual
# load("life_stage_estimates.RData") #object life_stage
# life_stage <- as.data.frame(life_stage)
# 
# # Transform timestamps in POSIXct
# life_stage$FPT_fledge_date <- as.POSIXct(life_stage$FPT_fledge_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")
# all_GpsAcc$timestamp <- as.POSIXct(all_GpsAcc$timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
# # Split by individual id
# ind_ls <- split(all_GpsAcc, all_GpsAcc$name)
# # All individual names match
# table(names(ind_ls) %in% unique(life_stage$local_identifier))
# # Check time range before and after filtering by fledging date
# table(is.na(life_stage$FPT_fledge_date)) #There are NAs
# ind_ls_fledged <- lapply(ind_ls, function(ind){
#   indFledgeTime <- life_stage$FPT_fledge_date[life_stage$local_identifier == unique(ind$name)]
#   if(!is.na(indFledgeTime)){
#     ind <- ind[ind$timestamp > indFledgeTime,]
#   }
#   return(ind)
# })
# 
# all_GpsAcc_postFledge <- as.data.frame(rbindlist(ind_ls_fledged))
# 



