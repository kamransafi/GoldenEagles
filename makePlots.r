setwd("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/Data/")
# creds <- movebankLogin("Kami", "Wehile85")
# study <- "Life Track Golden Eagle Alps"
# animalInfo <- getMovebankAnimals(study, login=creds)
# AnimalNames <- unique(animalInfo$local_identifier)
# 
# simpleThin <- function(x){
#   eagle <- getMovebankData(study, creds, animal=x, removeDuplicatedTimestamps=T)
#   eagle1h <- eagle[!duplicated(cbind(year(timestamps(eagle)), month(timestamps(eagle)), day(timestamps(eagle)), hour(timestamps(eagle)))),]
#   message(paste("Done with ", x, ".", sep=""))
#   return(eagle1h)
# }
# allEagles <- lapply(AnimalNames, simpleThin)
# saveRDS(allEagles, "./allEagles.rds")
library(move)
library(lubridate)

Eagles <- readRDS("./allEagles.rds")
creds <- movebankLogin("Kami", "Wehile85")
study <- "Life Track Golden Eagle Alps"
animalInfo <- getMovebankAnimals(study, login=creds)
AnimalNames <- unique(animalInfo$local_identifier)
AnimalNames <- unique(animalInfo$local_identifier)
studyArea <- bbox(extent(c(6.147962, 15.512771, 45.388163, 48.016078)))

DEM <- raster("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/SRTM/DEM.grd")

# Project <- function(x)
# {
#   B <- move(x, removeDuplicatedTimestamps=TRUE)
#   Bproj <- spTransform(B, CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs"))
#   saveRDS(Bproj, paste(strsplit(x, ".csv"), ".rds", sep=""))
# }
# lapply(birds, Project)


lakes1 <- readOGR("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/GLWD-level1/glwd_1.shp", layer = "glwd_1")
lakes1 <- crop(lakes1, studyArea)

library(rnaturalearth)
library(scales)
#world countries
wrld <- ne_countries(scale=10)
wrld <- crop(wrld, studyArea)
slope <- terrain(DEM, opt='slope')
aspect <- terrain(DEM, opt='aspect')
hill <- hillShade(slope, aspect, 60, 270)
r <- raster(ncols=13, nrows=13, xmn=0)
gf <- focalWeight(r, 13, "Gauss")
hillBlur <- focal(hill, w=gf^2)


DEMlr <- aggregate(DEM, 3)
hillBlurLR <- aggregate(hillBlur, 3)

EaglesSim <- lapply(Eagles, function(x){
  ts <- timestamps(x)
  newT <- ts - years(min(year(ts)-2017))
  if(any(is.na(as.numeric(newT))))
  {
    x <- x[which(!is.na(as.numeric(newT))),]
    newT <- newT[!is.na(as.numeric(newT))]
  }
  timestamps(x) <- newT
  print(x@idData$local_identifier)
  return(x)
})

redEagle <- EaglesSim[c(-14,-37,-38,-49)]
OyEagles <- lapply(redEagle, function(x) {
  x[timestamps(x)>=as.POSIXct("2017-08-01", format="%Y-%m-%d", tz="UTC") & timestamps(x)<as.POSIXct("2018-08-01", format="%Y-%m-%d", tz="UTC"),]
})

frame <- seq(as.POSIXct("2017-08-01", format="%Y-%m-%d", tz="UTC"), as.POSIXct("2018-08-01", format="%Y-%m-%d", tz="UTC") , 3600)[-1]
pT <- lapply(1:length(frame), function(x) {
  ((day(frame[x])-1)*24+hour(frame[x]))/(24*max(day(frame)[month(frame)==month(frame[x])]))
})
pT <- unlist(pT)
frame <- cbind(frame, pT)

writeRaster(DEMlr, "/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/DEMlr.grd", overwrite=TRUE)
writeRaster(hillBlurLR, "/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/hillBlurLR.grd", overwrite=TRUE)
save("wrld", "lakes1", "EaglesSim", "redEagle", "OyEagles", "frame", "studyArea", file="/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/Items4Plots.Rdata")
