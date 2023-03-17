
## Script written by Martina in March 2023
# data are taken from Hester's and Elham's external hard drive (input data that are NOT on Martina's computer)
# The output of the script are GPS data per individual classified into 4 categories (circular soaring, linear soaring, gliding and other).
# The names of the output data files are animalID_classifiedBursts_df and were saved on Martina's computer and copied also to the external hard drive.
# Output plots are also saved on Martina's computer

#____________________________________________________________
### Transform in movestack and calculate track variables ####
#____________________________________________________________

library(move)

# Set the directory (will also be used to save the output file) 
setwd("C:/Users/Tess Bronnvik/Desktop/Improvement_and_Golden_Eagles") 
# Import golden eagle movestack
load("eagle_ms.RData") # my object is called mvstk
# In my script the movestack object containing all eagles is calles mvstk, in the next line change it to the name of your object
mv_ls <- split(eagle_ms)
# The rest of the code is run individual by individual
mv_ls <- mv_ls[sapply(mv_ls, n.locs)>=30] #keep only move objects (individuals) with more than 30 locations (needed for segmentation)
lapply(names(mv_ls), function(animalName){
  mv <- mv_ls[[animalName]]
  mv$timelag.sec <- c(NA, timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA, (mv$height_above_ellipsoid[-1] - mv$height_above_ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA, distance(mv))
  #mv$step.length <- distm(x=coordinates(mv), fun=distVincentyEllipsoid)
  #mv$step.length <- c(NA, raster::pointDistance(mv[-sum(n.locs(mv)),], mv[-1,], longlat=isLonLat(mv)))
  #mv$step.length2 <- distm(x=coordinates(mv), fun=pointDistance)
  mv$gr.speed <- c(NA, speed(mv))
  #mv$gr.speed <- c(NA, mv$step.length/mv$timelag.sec) #Error in `[[<-.data.frame`(`*tmp*`, name, value = c(NA, NA, 0.0119871263928527, : replacement has 1320202 rows, data has 1149
  animalID <- mv@idData$local_identifier#strsplit(animalName, "\\.")[[1]]
  #animalID <- paste(animalID[c(1,length(animalID))], collapse="_")
  #if(substr(animalName,1,1)=="X"){animalName <- gsub("X","",animalName)} #When you split a movestack, if the id starts with a number R adds a X, which we don't want
  save(mv, file = paste0("tracks/", animalID,"_gpsNoDup_moveObj_", Sys.Date(), ".rdata"))
})

for (i in 1:length(eagle_ls)) {
  mv <- mv_ls[[i]]
  mv$timelag.sec <- c(NA, timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA, (mv$height_above_ellipsoid[-1] - mv$height_above_ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA, distance(mv))
  mv$gr.speed <- c(NA, speed(mv))
  animalID <- mv@idData$local_identifier
  save(mv, file = paste0("tracks/", animalID,"_gpsNoDup_moveObj_", Sys.Date(), ".rdata"))
  print(paste0("Calculated track metrics for animal ", i, ", ", animalID, "."), quote = F)
}


#___________________________________________________________________
### SEGMENTATION on bursts of continuous 1 sec resolution data ####
#__________________________________________________________________

library(move)
library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(detectCores()-1) 
#dir.create("prepped")

#setwd("/home/mscacco/ownCloud/Martina/ProgettiVari/GoldenEagles_WindMapsFromThermals/ContributionToElhamsPaper")
dir.create("newGPSsegmentation_March2023/segmentationPlots")
dir.create("newGPSsegmentation_March2023/classifiedData")

# Import the files saved in the previous step in your working directory
fls <- list.files("...hard drive/tracks/", pattern="gpsNoDup_moveObj", full.names = T)

# fls <- list.files("D:/Hester_GE/mar16/tracks/", pattern="gpsNoDup_moveObj", full.names = T)
# inds <- sapply(strsplit(fls, "/|_"), "[",6)
# flsDone <- list.files("D:/Hester_GE/mar16/thermals/", full.names = T)
# indsDone <- sapply(strsplit(flsDone, "/|_"), "[",6)
# indsToDO <- inds[!inds %in% indsDone]
# flsToDO <- sapply(indsToDO, function(ind){grep(ind, fls, value=T, fixed=T)})

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

#llply(fls, function(f){ #for parallel
lapply(fls, function(f){
  print(f)
  load(f) #object mv
  animalId <- paste0(mv@idData$individual_id,"_",mv@idData$local_identifier)
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
      b$track_flight_id <- paste0(unique(b$individual_id),"_",unique(b$burstID),"_segm_",b$flightNum_smooth3) 
      
      return(b) #return each classified and smoothed burst to a list
    }, .parallel=T) # to make it run in parallel use llply (instead of lapply) with .parallel=T
    
    # Here we could decide to set a minimum number of classified observations (!= other) for a burst to be considered
    #summary(sapply(burst_ls_class_smooth,nrow)) #usual duration of a burst
    nClassObs <- sapply(burst_ls_class_smooth, function(b) length(which(b$flightClust_smooth3 != "other")))
    burst_ls_class_smooth <- burst_ls_class_smooth[nClassObs >= minN_classifiedObs]
    
    # Rbind all bursts and save classified and smoothed dataframe per individual
    HRdf_smooth <- as.data.frame(rbindlist(burst_ls_class_smooth))
    save(HRdf_smooth, file = paste0("newGPSsegmentation_March2023/classifiedData/animal_",animalId,"_classifiedBursts_df.rdata"))
  }
})


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

burst_ls_class_smooth <- split(HRdf_smooth, HRdf_smooth$burstID)
set.seed(1512)
randomBursts <- sample(1:length(burst_ls_class_smooth), 100)

lapply(burst_ls_class_smooth[randomBursts], function(b){
  #b=burst_ls_class_smooth[["7006"]]
  
  print(unique(b$burstID))
  animalID <- unique(b$local_identifier)
  
  b <- b[order(b$timestamp),]
  # cbind(as.character(b$flightClust),as.character(b$flightClust_smooth),as.character(b$flightClust_smooth2), as.character(b$flightClust_smooth3))
  
  # calculate aspect ratio for plot
  rangeLong <- max(b$location_long)-min(b$location_long)
  rangeLat <- max(b$location_lat)-min(b$location_lat)
  
  # Plot results
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),".png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust], pch=19)
  dev.off()
  
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth.png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flightClust_smooth3], pch=19)
  dev.off()
  
  # this 3D plots are only useful when interactive, exporting them doesn not make much sense
  # plot3d(b[,c("location_long","location_lat","height_above_ellipsoid")], type="l", col="darkgrey")
  # points3d(b[,c("location_long","location_lat","height_above_ellipsoid")], col=c("red","darkgreen","blue","grey")[b$flightClust_smooth3], size=5)
  # aspect3d(x=rangeLong/rangeLat, y=1, z=1)
  # snapshot3d(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burstID),"_smooth3D.png"), width = 600, height = 600)
})

someBursts <- c("6","32", "2983", '3026', "3071", '6975', "7100", "7164", "3023", "3029", "3163")
# interactive 3D plots with plot_ly just for a few bursts
lapply(burst_ls_class_smooth[someBursts], function(b){
  animalID <- unique(b$local_identifier)
  b <- b[order(b$timestamp),]
  # calculate aspect ratios for plot along 3 axes
  rangeLong <- max(b$location_long)-min(b$location_long)
  rangeLat <- max(b$location_lat)-min(b$location_lat)
  rangeLat_m <- (rangeLat*111.139)*1000 #transform range in metres (approximated) to compare with elevation
  rangeElev <- max(b$height_above_ellipsoid)-min(b$height_above_ellipsoid)
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
  p <- plot_ly(b, x=~location_long, y=~location_lat, z=~height_above_ellipsoid, colors=pal,
               type="scatter3d", mode="lines", name=~burstID,
               line = list(color = 'darkgrey'),
               showlegend = F) %>%
    add_trace(data = b, x=~location_long, y=~location_lat, z=~height_above_ellipsoid, 
              colors=pal, color=~flightClust_smooth3,
              type="scatter3d", mode="markers",
              text=~paste0("vert.speed: ", round(vert.speed,2),"\n",
                           "turn.angle: ",round(turnAngle_smooth,2),"\n",
                           "ground.speed: ",round(gr.speed,2)),
              marker = list(size = 5),
              showlegend = T, inherit = F) %>%
    layout(title = paste0("Animal ",unique(b$local_identifier)," - burst ",unique(b$burstID)),
           scene=list(xaxis = list(title = "Longitude"), 
                      yaxis = list(title = "Latitude"), 
                      zaxis = list(title = "Height above ellipsoid (m)"),
                      aspectmode = "manual", aspectratio=list(x=aspects["x"], y=aspects["y"], z=aspects["z"]))) #preserve long/lat aspect ratio
  print(p)
  readline(prompt="Press [enter] to continue")
})
