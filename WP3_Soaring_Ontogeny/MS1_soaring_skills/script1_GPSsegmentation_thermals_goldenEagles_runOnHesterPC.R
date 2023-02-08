
## This is an example script wrote by Martin but run on Hester's computer using data that are not on Martina's computer
## The output of the script are GPS data per individual classified into soaring and gliding and circular and linear.
## The names of the output files are animalID_classifiedBursts_df and were sent by Hester to Martina.

## To run the script you only need to change the working directory (lines 11 and 41) and the name of your movestack object at line 15 :)

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
#_________________________________________________________
### Select bursts of continuous 1 sec resolution data ####
#_________________________________________________________
library(move)
library(plyr)
library(doParallel)
detectCores()
doParallel::registerDoParallel(detectCores()-1) 
#dir.create("prepped")

# Import the files saved in the previous step in your working directory
fls <- list.files("tracks/", pattern="gpsNoDup_moveObj", full.names = T)

# fls <- list.files("D:/Hester_GE/mar16/tracks/", pattern="gpsNoDup_moveObj", full.names = T)
# inds <- sapply(strsplit(fls, "/|_"), "[",6)
# flsDone <- list.files("D:/Hester_GE/mar16/thermals/", full.names = T)
# indsDone <- sapply(strsplit(flsDone, "/|_"), "[",6)
# indsToDO <- inds[!inds %in% indsDone]
# flsToDO <- sapply(indsToDO, function(ind){grep(ind, fls, value=T, fixed=T)})

minResol <- 1 # only 1 sec timelag
minBurstDuration <- 30 # we want bursts of at least 30 secs

swV <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
swT <- 14 #smoothing window of 29 seconds for thermalling behavior (according to Rolf/Bart's paper)

llply(flsToDO, function(f){
  print(f)
  load(f) #object mv
  animalID <- mv@idData$local_identifier
  # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
  mv$burstID <- c(0, cumsum(mv$timelag.sec[2:n.locs(mv)] != minResol))  #from row 2 (since the first is NA)
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
    # Keep only bursts with at least 40 s (30 of smoothing window will be NA) 
    burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= 40)]
    # Compute smoothed turning angle separately for each burst
    HRdf <- do.call(rbind, lapply(burst_ls_corr_sub, function(b){
      b$vertSpeed_smooth <- NA
      b$turnAngle_smooth <- NA
      for(i in (swV+1):(nrow(b)-swV)){
        b$vertSpeed_smooth[i] <- mean(b$vert.speed[(i-swV):(i+swV)], na.rm=T)}
      for(i in (swT+1):(nrow(b)-swT)){
        b$turnAngle_smooth[i] <- max(abs(cumsum(b$turn.angle[(i-swT):(i+swT)])))}
      return(b) # return df with smoothed variables
    }))
    # Classify soaring only based on vertical speed
    HRdf <- HRdf[complete.cases(HRdf$vertSpeed_smooth),]
    kmeanV <- kmeans(HRdf$vertSpeed_smooth, 2)   #Get the two clusters
    soarId <- which.max(aggregate(HRdf$vertSpeed_smooth~kmeanV$cluster, FUN=mean)[,2]) # which one is the soaring one?
    soarClust <- rep("glide", length(kmeanV$cluster))
    soarClust[which(kmeanV$cluster==soarId)] <- "soar"
    HRdf$soarClust <- factor(soarClust, levels=c("soar","glide"))  
    # Now classify thermaling only based on turning angle (cumulated to a 30 s time window in previous step)
    HRdf$thermalClust <- "other"
    HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$turnAngle_smooth >= 300)] <- "circular"
    HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$soarClust=="soar" & HRdf$thermalClust != "circular")] <- "linear"
    HRdf$thermalClust <- factor(HRdf$thermalClust, levels=c("circular","linear","other"))
    HRdf$name <- animalID
    # Save classified dataframe per individual
    save(HRdf, file = paste0("thermals/",animalID,"_classifiedBursts_df.RData"))
  }
}, .parallel=F)

load("prepped/Almen19_7001_classifiedBursts_df.RData")

# associated ACC and the ODBA should show whether a linear is high (flapping) or low (gliding) energy
# distance to thermal centroid and whether a second thermaling event is distant or near to the first
# total displacement within thermals should show the horizontal length and if it is associated with 
# wind speed it should be a blown thermal