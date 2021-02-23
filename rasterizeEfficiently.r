library(move)
library(fields)
library(sf)
library(stars)
setwd("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/")
biRDS <- list.files("./Data/MBDL20201205", pattern=".rds", full.names = T)
#sizeDat <- sum(unlist(lapply(biRDS, function(x) {n.locs(readRDS(x))})))
t0 <- Sys.time()
studyArea <- st_bbox(c(xmin = 6.5, xmax = 15, ymax = 48.2, ymin = 45.4), crs=CRS("+proj=longlat +datum=WGS84"))
grd <- st_as_stars(studyArea, dx = 0.001, dy = 0.001, values = 0)
sf_extSoftVersion()["GDAL"]
Rastrack <- function(x, grid){
  ls <- st_sf(a = 1, st_sfc(st_linestring(coordinates(x))), crs=CRS("+proj=longlat +datum=WGS84"))
  tmp <- st_rasterize(ls, grd, options="ALL_TOUCHED=TRUE")
  return(tmp)
  message(paste("Done with ", x@idData$local_identifier, ".", sep=""))
}
test <- lapply(biRDS, function(x) {Rastrack(readRDS(x), grid = grd)})
sumRas <- test[[1]]
for(i in 2:length(test)){
  sumRas <- sumRas+test[[i]]
  print(i)
}
sumRas[sumRas==0] <- NA
writeRaster(sumRas, file="./trackDen.tif", format="GTiff", overwrite=T)
jpeg("./denplot.jpg", height=2600, width=8500)
plot(as(sumRas, "Raster"), col=tim.colors(24), maxpixels=15e9,  axes=F, border=NA, legend.cex=3)
dev.off()
print(paste("Needed ", round(difftime(Sys.time(), t0, units="secs"),1) , " seconds for loanding and rasterizing ", sizeDat, " locations at a resolution of 30 arc sec."))
