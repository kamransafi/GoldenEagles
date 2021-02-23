setwd("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis")
birds <- list.files("./Data", pattern=".csv", full.names = T)

referenceDat <- read.csv(birds[grep("LifeTrack Golden Eagle Alps", birds)], header=T, sep=",", as.is=T)
birds <- birds[-grep("LifeTrack Golden Eagle Alps", birds)]

library(move)
library(lubridate)

biRDS <- list.files("./Data", pattern=".rds", full.names = T)
birds <- data.frame(file_orig=unlist(birds), file_rds= unlist(biRDS), newID= paste0("Eagle.", sprintf("%02d", 1:59)))

referenceDat$animal.id <- gsub("?", "", as.character(referenceDat$animal.id), fixed=T)
referenceDat$animal.sex[referenceDat$animal.sex==""] <- NA
referenceDat$animal.ring.id[referenceDat$animal.ring.id==""] <- NA
Names <- referenceDat$animal.id[referenceDat$animal.id!=""]
birds$animal.id <- NULL
birds$sex <- NULL
birds$ring.id <- NULL
for(i in 1:length(Names))
{
  if(any(grepl(Names[i], birds$file_orig, fixed=T)))
  {
    rn <- which(grepl(Names[i], birds$file_orig, fixed=T))
    birds[rn,"animal.id"] <- gsub(" ", "", Names[i], fixed=T) 
    birds[rn,"sex"] <- referenceDat$animal.sex[i]
    birds[rn,"ring.id"] <- gsub(" ", "", referenceDat$animal.ring.id[i], fixed=T) 
    print(paste("Matched ", Names[i], " in reference data to row ", rn, " in bird table.", sep=""))
  }else{
    print(paste("No match for ", Names[i], ".", sep=""))
  }
  rn <- NULL
}
