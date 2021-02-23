## SRTM download options for movement data
#srtm15p extracts terrestrial topography and marine bathymetry data from: Tozer, B, Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., & Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: SRTM15+. Earth and Space Science, 6, 1847– 1864. https://doi.org/10.1029/2019EA000658
#srtm90m extracts terrestrial data from: Jarvis A., H.I. Reuter, A. Nelson, E. Guevara, 2008, Hole-filled seamless SRTM data V4, International Centre for Tropical Agriculture (CIAT), available from http://srtm.csi.cgiar.org.
#citation: Reuter H.I, A. Nelson, A. Jarvis, 2007, An evaluation of void filling interpolation methods for SRTM data, International Journal of Geographic Information Science, 21:9, 983-1008.

library(RCurl)
library(raster)
srtm15p <- function(x){
  bb <- as.vector(extent(x)*1.1)
  txt <- postForm("https://topex.ucsd.edu/cgi-bin/get_srtm15.cgi",
           submitButton = "get data",
           north=bb[4],
           west=bb[1],
           east=bb[2],
           south=bb[3],
           style = "POST")
  xyz <- data.frame(matrix(as.numeric(unlist(strsplit(unlist(strsplit(txt, "\n")), "\t"))), ncol=3, byrow=T))
  names(xyz) <- c("long", "lat", "alt")
  coordinates(xyz) <- ~long+lat
  gridded(xyz) <- TRUE
  return(raster(xyz))
}

srtm90m <- function(x, method="wget", quiet=F){
  bb <- as.vector(extent(x)*1.1)
  xtilellc <- ceiling((bb[1]+180)/5)
  ytilellc <- 24 - floor((bb[3]+60)/5)
  xtileurc <- ceiling((bb[2]+180)/5)
  ytileurc <- 24 - floor((bb[4]+60)/5)
  combi <- expand.grid(xtilellc:xtileurc, ytilellc:ytileurc)
  tile <- apply(combi, 1, function(x) paste("srtm_", sprintf("%02d",x[1]), "_" , sprintf("%02d", x[2]), ".zip", sep=""))
  tileURL <- lapply(tile, function(tn) paste("http://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF/", tn, sep=""))
  demsDL <- mapply(function(URL, tileName) download.file(URL, paste(tempdir(), tileName, sep="/"), method=method, quiet=quiet), tileURL, tile)
  unzip(paste(tempdir(), tile, sep="/"), exdir=tempdir())
  dems <- lapply(list.files(tempdir(), pattern=".tif", full.names = T), raster)
  # Merge the tiles into a single raster
  if(length(dems)>1){
    DEM <- do.call("merge", dems)}else{DEM <- dems[[1]]}
  # Reduce to the extent of the study area
  DEM <- crop(DEM, bb)
  return(DEM)
}


library(move)
Dem15Leroy <- srtm15p(leroy)
Dem90mLeroy <- srtm90m(leroy, quiet=T)
plot(Dem15Leroy)
lines(leroy)
plot(Dem90mLeroy)
lines(leroy)
plot(extract(Dem15Leroy, leroy, method="bilinear")~extract(Dem90mLeroy, leroy, method="bilinear"), pch=20)
abline(lm(extract(Dem15Leroy, leroy)~extract(Dem90mLeroy, leroy)-1), lwd=2, col="firebrick")
