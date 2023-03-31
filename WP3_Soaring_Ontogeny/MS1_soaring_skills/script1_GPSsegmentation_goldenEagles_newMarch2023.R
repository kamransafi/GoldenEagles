
## Script written by Martina in March 2023
# input data is a csv with all individuals taken from Hester's and Elham's external hard drive (input data that are NOT on Martina's computer)
# data were downloaded mid-March 2023
# The output file that will be used for segmentation is animalName_gpsNoDup_moveObj (saved on both hard drive and Martina's computer)
# These data get classified into 4 categories (circular soaring, linear soaring, gliding and other).
# The names of the second output data files are animalID_classifiedBursts_df (saved on both hard drive and Martina's computer)
# Output plots are also saved on Martina's computer

#_____________________________________________________________________
### Import data, remove duplicates and calculate track variables ####
#_____________________________________________________________________

library(move)
library(data.table)
library(plyr)
library(doParallel)
doParallel::registerDoParallel(3) 

# Import golden eagle data from csv downloaded by Elham
df <- fread("/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/LifeTrack Golden Eagle Alps.csv")
names(df) <- gsub("-|:", ".", names(df))
df <- df[df$sensor.type=="gps",] 
# Check for missing infos in time and coords
# These usually correspond to eobs status B to D: 
#A = position and time within accuracy masks
#B = only time of week and weeknumber valid
#C = only weeknumber valid
#D = no valid data
anyNA(df$timestamp)
anyNA(df$individual.local.identifier)
anyNA(df$location.long); anyNA(df$location.lat)
table(df$eobs.status)
df <- df[df$eobs.status %in% c("","A"),]
df <- df[!is.na(df$location.long),]

# remove unnecessary columns
colsToKeep <- c("event.id", "timestamp", "location.long", "location.lat", "height.above.ellipsoid",
                "gps.satellite.count","gps.dop","eobs.status",
                "ground.speed","eobs.horizontal.accuracy.estimate",
                "study.name", "tag.local.identifier", "individual.local.identifier","individual.taxon.canonical.name")
df <- data.frame(df)[,colsToKeep]
df$event.id <- as.character(df$event.id) #important for solving the error: 'names' attribute [13] must be the same length as the vector [12]

## split into list of individuals and save each individual separately
length(unique(df$individual.local.identifier))
df_ls <- split(df, df$individual.local.identifier)
df_ls <- df_ls[sapply(df_ls, nrow)>=30] #keep only individuals with > 30 locations (needed for segmentation)

lapply(df_ls, function(ind){
  print(unique(ind$individual.local.identifier))
  saveRDS(ind, file = paste0("/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/",unique(ind$individual.local.identifier),"_gpsClean.rds"))
})

# Calculate track geometry infos for each individual
pathToHD <- "/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/"
fls <- list.files(pathToHD, pattern="gpsClean")
done <- list.files(pathToHD, pattern="moveObj")
fls <- fls[! gsub("_.*", "", fls) %in% gsub("_.*", "", done)]

llply(fls, function(f){
  
  print(f)
  ind <- readRDS(paste0(pathToHD,f))
  animalName <- unique(ind$individual.local.identifier)
  # Check for duplicates
  dup_ls <- getDuplicatedTimestamps(x=ind$individual.local.identifier,
                                    timestamps=ind$timestamp)
  #ind[dup_ls[[1]],] check duplicates, one of them is full of NA and no active status
  dupRows <- ind[unlist(dup_ls),] #how do they differ, remove those without Active status
  eventsToDrop <- dupRows$event.id[! dupRows$eobs.status=="A"]
  ind_nodup <- ind[! ind$event.id %in% eventsToDrop,]
  dup_ls <- getDuplicatedTimestamps(x=ind_nodup$individual.local.identifier,
                                    timestamps=ind_nodup$timestamp)
  if(length(dup_ls)>0){
    warning(paste0(animalName, " has still duplicated timestamps, only the first occurrence is kept."))
    eventsToDrop <- ind_nodup$event.id[sapply(dup_ls, "[", -1)] #remove all duplicates except the first one
    ind_nodup <- ind_nodup[! ind_nodup$event.id %in% eventsToDrop,]
  }
  
  # calculate track variables
  mv <- move(x=ind_nodup$location.long, y=ind_nodup$location.lat, 
             time=ind_nodup$timestamp,
             proj=crs("+proj=longlat +ellps=WGS84"),
             animal=ind_nodup$individual.local.identifier,
             data=ind_nodup)
  mv$timelag.sec <- c(NA,timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA,(mv$height.above.ellipsoid[-1] - mv$height.above.ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA,move::distance(mv))
  #mv$step.length <- distm(x=coordinates(mv), fun=distVincentyEllipsoid)
  #mv$step.length <- c(NA, raster::pointDistance(mv[-sum(n.locs(mv)),], mv[-1,], longlat=isLonLat(mv)))
  #mv$step.length2 <- distm(x=coordinates(mv), fun=pointDistance)
  mv$gr.speed <- c(NA, speed(mv))
  #mv$gr.speed <- c(NA, mv$step.length/mv$timelag.sec) #Error in `[[<-.data.frame`(`*tmp*`, name, value = c(NA, NA, 0.0119871263928527, : replacement has 1320202 rows, data has 1149
  
  # save on the hard drive
  save(mv, file = paste0(pathToHD,animalName,"_gpsNoDup_moveObj.rdata"))
  
}, .parallel=T)

# copy the files on Martina's computer
setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/data/") 
fls <- list.files(pathToHD, pattern="moveObj")
file.copy(from=paste0(pathToHD,fls), to=fls)


#___________________________________________________________________
### SEGMENTATION on bursts of continuous 1 sec resolution data ####
#__________________________________________________________________

library(move)
library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(6) 

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/")
pathToHD <- "/media/mscacco/Ellham's HDD/GE_all_gps_Mar23/"
dir.create("GPSsegmentationPlots")
dir.create("classifiedData")
dir.create(paste0(pathToHD,"GPSsegmentationPlots"))
dir.create(paste0(pathToHD,"classifiedData"))

# Import the files saved in the previous step in your working directory
fls <- list.files("data", pattern="gpsNoDup_moveObj", full.names = T)

# Define segmentation parameters
minResol <- 2 # 1 to max 2 sec timelag
minBurstDuration <- 30 # we want bursts of at least 30 secs

swV <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
swT <- 12 #smoothing window of 29 seconds for thermalling behaviour (according to Rolf/Bart's paper) (we changed it to smooth window of 21 sec)
circlDegrees <- 200 #degrees of rotation to be achieved in the time defined by swT*2
minBehavDuration <- 5 #minimum duration in seconds of a specific behaviour, when less than this and if in between two segments of a different behaviour it will be incorporated in the previous and follwoing segment 
minThermalDuration <- 20 #minimum duration for a circling event to be considered as thermalling
minGroundSpeed <- 3.5 #minimum ground speed expected during flight (in golden eagles good bursts of flight had a minimum gr speed = 5 m/s)
minN_classifiedObs <- 30 #minimum number of classified observations (non-other behaviour) in a burst, for keeping a burst
# golden eagles have bursts of 15 min = 900 secs, we could consider bursts where at least 30% of observations are classified

#llply(fls[1:5], function(f){ #for parallel
results <- lapply(fls, function(f)try({
  print(f)
  load(f) #object mv
  animalId <- mv@idData$individual.local.identifier
  # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
  mv$burstID <- c(0, cumsum(mv$timelag.sec[2:n.locs(mv)] > minResol))  #from row 2 (since the first is NA)
  # with table we can count all locations with the same ID (each one is a separate burst) and keep those high resolution bursts that have at least a certain number of locations (= minBurstDuration)
  burstDuration <- as.data.frame(table(mv$burstID))
  burstsToKeep <- burstDuration$Var1[which(burstDuration$Freq >= minBurstDuration)]
  # use those to subset the move obj and keep only the high resolution bursts
  HRmv <- mv[which(mv$burstID %in% burstsToKeep),]
  if(nrow(HRmv)>0){
    HRdf_bursts <- as.data.frame(HRmv)
    # Remove unnecessary columns
    HRdf_bursts <- HRdf_bursts[,-grep("mag|orientation|coords|timestamps|start.timestamp|optional|import|visible|algorithm|battery|decoding|accuracy|manually|activity|checksum|acceleration",
                                      colnames(HRdf_bursts))]
    # Split each individual dataframe by burst ID
    burst_ls_corr <- split(HRdf_bursts, HRdf_bursts$burstID)
    # Keep only bursts with minBurstDuration (30 of smoothing window will be NA) 
    burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= minBurstDuration)]
    # Compute smoothed turning angle separately for each burst
    HRdf <- as.data.frame(rbindlist(
      llply(burst_ls_corr_sub, function(b){
        b$vertSpeed_smooth <- NA
        b$turnAngle_smooth <- NA
        for(i in (swV+1):(nrow(b)-swV)){
          b$vertSpeed_smooth[i] <- mean(b$vert.speed[(i-swV):(i+swV)], na.rm=T)}
        for(i in (swT+1):(nrow(b)-swT)){
          b$turnAngle_smooth[i] <- max(abs(cumsum(b$turn.angle[(i-swT):(i+swT)])))}
        return(b) # return df with smoothed variables
      }, .parallel=T) # to make it run in parallel use llply (instead of lapply) with .parallel=T)
    ))
    # Classify soaring only based on vertical speed
    HRdf <- HRdf[complete.cases(HRdf$vertSpeed_smooth),]
    kmeanV <- kmeans(HRdf$vertSpeed_smooth, 2)   #Get the two clusters
    soarId <- which.max(aggregate(HRdf$vertSpeed_smooth~kmeanV$cluster, FUN=mean)[,2]) # which one is the soaring one?
    soarClust <- rep("glide", length(kmeanV$cluster))
    soarClust[which(kmeanV$cluster==soarId)] <- "soar"
    HRdf$soarClust <- factor(soarClust, levels=c("soar","glide"))  
    # Now classify thermalling only based on turning angle (cumulated to a 25 s time window in previous step)
    HRdf$flightClust <- "other"
    HRdf$flightClust[which(HRdf$gr.speed >= minGroundSpeed & HRdf$soarClust=="soar" & HRdf$turnAngle_smooth >= circlDegrees)] <- "circular soaring" #complete 150 degrees in 15 sec
    HRdf$flightClust[which(HRdf$gr.speed >= minGroundSpeed & HRdf$soarClust=="soar" & HRdf$flightClust != "circular soaring")] <- "linear soaring"
    HRdf$flightClust[which(HRdf$gr.speed >= minGroundSpeed & HRdf$soarClust=="glide")] <- "gliding"
    HRdf$flightClust <- factor(HRdf$flightClust, levels=c("circular soaring","linear soaring","gliding","other"))
    
    # Add some steps of smoothing based on duration of behaviours:
    burst_ls_class <- split(HRdf, HRdf$burstID)
    burst_ls_class_smooth <- llply(burst_ls_class, function(b){
      # We assign a unique ID to each consecutive flight segment based on a rule (the ID increases if the class of each obs is different from the previous one)
      print(unique(b$burstID))
      b <- b[order(b$timestamp),]
      b$flightNum <- c(0, cumsum(b$flightClust[-1] != b$flightClust[-nrow(b)]))
      
      # we calculate the duration and unique class of each behavioural segment
      behavDuration <- merge(aggregate(timelag.sec~flightNum, data=b, FUN=sum), 
                             aggregate(flightClust~flightNum, data=b, FUN=unique), by="flightNum")
      # we create a new category ID, which is the same as the original one, 
      # unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
      behavDuration$flightNum_smooth <- behavDuration$flightNum
      behavDuration$flightClust_smooth <- behavDuration$flightClust
      if(nrow(behavDuration)>2){
        for(i in 2:(nrow(behavDuration)-1)){
          if(behavDuration$timelag.sec[i] <= minBehavDuration & behavDuration$flightClust_smooth[i-1] == behavDuration$flightClust_smooth[i+1]){
            behavDuration$flightNum_smooth[c(i, i+1)] <- behavDuration$flightNum_smooth[i-1]
            behavDuration$flightClust_smooth[c(i, i+1)] <- behavDuration$flightClust_smooth[i-1]
          }else if(behavDuration$timelag.sec[i] <= minBehavDuration & behavDuration$flightClust_smooth[i-1] != behavDuration$flightClust_smooth[i+1]){
            behavDuration$flightNum_smooth[c(i)] <- behavDuration$flightNum_smooth[i-1]
            behavDuration$flightClust_smooth[c(i)] <- behavDuration$flightClust_smooth[i-1]
          }
        }
      if(nrow(behavDuration)==2){
        if(any(behavDuration$timelag.sec <= minBehavDuration)){
          longestBehav <- which.max(behavDuration$timelag.sec)
          behavDuration$flightNum_smooth <- behavDuration$flightNum_smooth[longestBehav]
          behavDuration$flightClust_smooth <- behavDuration$flightClust_smooth[longestBehav]
        }
      }}
      b <- merge(b, behavDuration[c("flightNum","flightNum_smooth","flightClust_smooth")], by="flightNum", all.x=T)
      
      # recalculate segment duration based on smoothed classification and reclassify as linear soaring all circling that lasts < 30 seconds
      behavDuration_smooth <- merge(aggregate(timelag.sec~flightNum_smooth, data=b, FUN=sum), 
                                    aggregate(flightClust_smooth~flightNum_smooth, data=b, FUN=unique), by="flightNum_smooth")
      behavDuration_smooth$flightNum_smooth2 <- behavDuration_smooth$flightNum_smooth
      behavDuration_smooth$flightClust_smooth2 <- behavDuration_smooth$flightClust_smooth
      if(nrow(behavDuration_smooth)>2){
        for(i in 2:(nrow(behavDuration_smooth)-1)){
          if(behavDuration_smooth$flightClust_smooth2[i] == "circular soaring"){
            if(behavDuration_smooth$timelag.sec[i] <= minThermalDuration & behavDuration_smooth$flightClust_smooth2[i-1] == "linear soaring" & behavDuration_smooth$flightClust_smooth2[i+1] == "linear soaring"){
              behavDuration_smooth$flightNum_smooth2[c(i, i+1)] <- behavDuration_smooth$flightNum_smooth2[i-1]
              behavDuration_smooth$flightClust_smooth2[c(i, i+1)] <- behavDuration_smooth$flightClust_smooth2[i-1]
            }}}
      }
      b <- merge(b, behavDuration_smooth[c("flightNum_smooth","flightNum_smooth2","flightClust_smooth2")], by="flightNum_smooth", all.x=T)
      # finally check the classification of the time window at the start and end of the track
      # these first and last points can only be classified as either gliding or linear, as their classification was only based on vertical speed but not turning angle
      # so if at the start or end there is linear soaring, but they are preceded or followed by circular soaring, they become circular
      behavDuration_smooth2 <- merge(aggregate(timelag.sec~flightNum_smooth2, data=b, FUN=sum), 
                                     aggregate(flightClust_smooth2~flightNum_smooth2, data=b, FUN=unique), by="flightNum_smooth2")
      behavDuration_smooth2$flightNum_smooth3 <- behavDuration_smooth2$flightNum_smooth2
      behavDuration_smooth2$flightClust_smooth3 <- behavDuration_smooth2$flightClust_smooth2
      if(nrow(behavDuration_smooth2) >= 2){
        if(behavDuration_smooth2$timelag.sec[1]==1){
          behavDuration_smooth2$flightNum_smooth3[1] <- behavDuration_smooth2$flightNum_smooth3[2]
        }else if(behavDuration_smooth2$timelag.sec[1] <= swT & 
                 behavDuration_smooth2$flightClust_smooth2[1] == "linear soaring" & behavDuration_smooth2$flightClust_smooth2[2] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[1] <- behavDuration_smooth2$flightNum_smooth3[2]
          behavDuration_smooth2$flightClust_smooth3[1] <- "circular soaring"
        }
        if(behavDuration_smooth2$timelag.sec[nrow(behavDuration_smooth2)]==1){
          behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)] <- behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)-1]
        }else if(behavDuration_smooth2$timelag.sec[nrow(behavDuration_smooth2)] <= swT & 
                 behavDuration_smooth2$flightClust_smooth2[nrow(behavDuration_smooth2)] == "linear soaring" & behavDuration_smooth2$flightClust_smooth2[nrow(behavDuration_smooth2)-1] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)] <- behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)-1]
          behavDuration_smooth2$flightClust_smooth3[nrow(behavDuration_smooth2)] <- "circular soaring"
        }
      }
      # merge with burst
      b <- merge(b, behavDuration_smooth2[c("flightNum_smooth2","flightNum_smooth3","flightClust_smooth3")], by="flightNum_smooth2", all.x=T)
      
      # Assign unique ID to the behavioural segment based on the final smoothest classification
      b$track_flight_id <- paste0(unique(b$individual.local.identifier),"_",unique(b$burstID),"_segm_",b$flightNum_smooth3) 
      
      return(b) #return each classified and smoothed burst to a list
    }, .parallel=T) # to make it run in parallel use llply (instead of lapply) with .parallel=T
    
    # Here we could decide to set a minimum number of classified observations (!= other) for a burst to be considered
    #summary(sapply(burst_ls_class_smooth,nrow)) #usual duration of a burst
    nClassObs <- sapply(burst_ls_class_smooth, function(b) length(which(b$flightClust_smooth3 != "other")))
    burst_ls_class_smooth <- burst_ls_class_smooth[nClassObs >= minN_classifiedObs]
    
    # Rbind all bursts and save classified and smoothed dataframe per individual
    HRdf_smooth <- as.data.frame(rbindlist(burst_ls_class_smooth))
    saveRDS(HRdf_smooth, file = paste0("classifiedData/",animalId,"_classifiedBursts_df.rds"))
  }
}))

#some individuals don't get process because they don't have high resolution GPS bursts. Check how many
results <- lapply(fls, function(f){
  load(f) #object mv
  # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
  mv$burstID <- c(0, cumsum(mv$timelag.sec[2:n.locs(mv)] > minResol))  #from row 2 (since the first is NA)
  # with table we can count all locations with the same ID (each one is a separate burst) and keep those high resolution bursts that have at least a certain number of locations (= minBurstDuration)
  burstDuration <- as.data.frame(table(mv$burstID))
  burstsToKeep <- burstDuration$Var1[which(burstDuration$Freq >= minBurstDuration)]
  return(length(burstsToKeep))
})
table(unlist(results)==0)

#copy file to HD
setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
fls <- list.files("classifiedData", pattern="classifiedBursts", full.name=T)
file.copy(from=fls, to=paste0(pathToHD,fls))


#______________________________________________________________
### PLOTTING sampled bursts to check segmentation results ####
#______________________________________________________________

library(rgl)
options(rgl.printRglwidget = TRUE)
library(webshot2)
library(scales)
library(ggplot2)
library(plotly)
library(data.table)

setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
fls <- list.files("classifiedData", pattern="classifiedBursts", full.name=T)

lapply(fls, function(f){
  
  print(f)
  HRdf_smooth <- readRDS(f) #obj HRdf_smooth
  animalID <- unique(HRdf_smooth$individual.local.identifier)
  burst_ls_class_smooth <- split(HRdf_smooth, HRdf_smooth$burstID)
  
  set.seed(1512)
  if(length(burst_ls_class_smooth)>20){n=20}else{n=length(burst_ls_class_smooth)}
  randomBursts <- sample(1:length(burst_ls_class_smooth), n)
  
  lapply(burst_ls_class_smooth[randomBursts], function(b){
    #print(unique(b$burstID))
    b <- b[order(b$timestamp),]
    # calculate aspect ratio for plot
    rangeLong <- max(b$location.long)-min(b$location.long)
    rangeLat <- max(b$location.lat)-min(b$location.lat)
    
    # Plot results
    png(paste0("GPSsegmentationPlots/",animalID,"_burst",unique(b$burstID),".png"))
    par(mfrow=c(1,2))
    plot(b$timestamp, b$height.above.ellipsoid, type="l", col="darkgrey", lwd=2,
         xlab="Timestamp", ylab="Height above ellipsoid (m)")
    points(b$timestamp, b$height.above.ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
    legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
    
    plot(b$location.long, b$location.lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2,
         xlab="Longitude", ylab="Latitude")
    points(b$location.long, b$location.lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
    mtext(paste0(unique(b$individual.local.identifier)," - burst ",unique(b$burstID), "\n raw classification"),                   # Add main title
          side = 3, line = - 2.5, outer = TRUE)
    dev.off()
    
    png(paste0("GPSsegmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth.png"))
    par(mfrow=c(1,2))
    plot(b$timestamp, b$height.above.ellipsoid, type="l", col="darkgrey", lwd=2,
         xlab="Timestamp", ylab="Height above ellipsoid (m)")
    points(b$timestamp, b$height.above.ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
    legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
    
    plot(b$location.long, b$location.lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2,
         xlab="Longitude", ylab="Latitude")
    points(b$location.long, b$location.lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
    mtext(paste0(unique(b$individual.local.identifier)," - burst ",unique(b$burstID), "\n smoothed classification"),                   # Add main title
          side = 3, line = - 2.5, outer = TRUE)
    dev.off()
    
    # this 3D plots are only useful when interactive, exporting them does not make much sense
    # plot3d(b[,c("location_long","location_lat","height_above_ellipsoid")], type="l", col="darkgrey")
    # points3d(b[,c("location_long","location_lat","height_above_ellipsoid")], col=c("red","darkgreen","blue","grey")[b$flightClust_smooth3], size=5)
    # aspect3d(x=rangeLong/rangeLat, y=1, z=1)
    # snapshot3d(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth3D.png"), width = 600, height = 600)
  })
})

#copy file to HD
setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/goldenEagles_march2023/") 
fls <- list.files("GPSsegmentationPlots", pattern="png", full.name=T)
file.copy(from=fls, to=paste0(pathToHD,fls), overwrite = T)

#__________________________________
# 3D plotting interactive
## Investigate specific bursts in 3D (these plots are not exported)
someBursts <- c("6","32", "2983", '3026', "3071", '6975', "7100", "7164", "3023", "3029", "3163")
# interactive 3D plots with plot_ly just for a few bursts
lapply(burst_ls_class_smooth[someBursts], function(b){
  animalID <- unique(b$individual.local.identifier)
  b <- b[order(b$timestamp),]
  # calculate aspect ratios for plot along 3 axes
  rangeLong <- max(b$location.long)-min(b$location.long)
  rangeLat <- max(b$location.lat)-min(b$location.lat)
  rangeLat_m <- (rangeLat*111.139)*1000 #transform range in metres (approximated) to compare with elevation
  rangeElev <- max(b$height.above.ellipsoid)-min(b$height.above.ellipsoid)
  ratioLat <- 1
  ratioLong <- rangeLong/rangeLat
  ratioElev <- rangeElev/rangeLat_m
  aspects <- c(x=ratioLong, y=ratioLat, z=ratioElev)
  aspects <- aspects/aspects[which.max(aspects)]
  # set colors
  pal <- c("circular soaring"="red",
           "linear soaring"="darkgreen",
           "gliding"="blue",
           "other"="grey")
  # interactive 3d plot
  p <- plot_ly(b, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, colors=pal,
               type="scatter3d", mode="lines", name=~burstID,
               line = list(color = 'darkgrey'),
               showlegend = F) %>%
    add_trace(data = b, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, 
              colors=pal, color=~flightClust_smooth3,
              type="scatter3d", mode="markers",
              text=~paste0("vert.speed: ", round(vert.speed,2),"\n",
                           "turn.angle: ",round(turnAngle_smooth,2),"\n",
                           "ground.speed: ",round(gr.speed,2)),
              marker = list(size = 5),
              showlegend = T, inherit = F) %>%
    layout(title = paste0("Animal ",unique(b$individual.local.identifier)," - burst ",unique(b$burstID)),
           scene=list(xaxis = list(title = "Longitude"), 
                      yaxis = list(title = "Latitude"), 
                      zaxis = list(title = "Height above ellipsoid (m)"),
                      aspectmode = "manual", aspectratio=list(x=aspects["x"], y=aspects["y"], z=aspects["z"]))) #preserve long/lat aspect ratio
  print(p)
  readline(prompt="Press [enter] to continue")
})
