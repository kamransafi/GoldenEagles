setwd("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/")
library(move)
library(lubridate)
library(doMC)
library(scales)
registerDoMC(cores=4)
load("./EnvData/Items4Plots.Rdata")
DEMlr <- raster("./EnvData/DEMlr.grd")
hillBlurLR <- raster("./EnvData/hillBlurLR.grd")
frameN <- 1:nrow(frame)
lf <- list.files("/mnt/NetBackup/frames", pattern=".png", full.names=T)
if(length(lf)!=0){
  lf <- lf[file.size(lf)!=0]
  frameN <- frameN[! frameN %in% as.numeric(sub(".png","",  (matrix(unlist((strsplit(lf, "/"))), ncol=5, byrow=T)[,5])))]
}

foreach(i=frameN) %dopar% {
  #for(i in frameN){
  dir.create(paste(tempdir(), "/tmp", i, sep=""))
  rasterOptions(tmpdir=paste(tempdir(), "/tmp", i, sep="")) 
  png(paste("/mnt/NetBackup/frames/", sprintf("%04d", i), ".png", sep=""), height=1044, width=3740) #"/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/test.png"
  plot(t(studyArea), type="n", asp=1, xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
  image(DEMlr^(1/1.5), col=alpha(terrain.colors(512), 1), asp=1, add=T)
  plot(hillBlurLR, col=alpha(grey(0:512/521), 0.35), add=T, legend=FALSE, main='')
  plot(wrld, add=T, lwd=10, border="grey70")
  plot(wrld, add=T, lwd=5, border="grey10")
  plot(lakes1, add=T, col="skyblue3", border="skyblue4")
  rect(6.550842, 47.39645, 8.329512, 47.73932, col=alpha("white", 0.7), border="white")
  rect(6.550842, 47.39645, (6.550842+(8.329512-6.550842)*frame[i,2]), 47.73932, col=alpha("grey", 0.7), border="grey")
  text(6.550842+0.5*(8.329512-6.550842), 47.39645+0.5*(47.73932-47.39645), adj=0.5, 
       month(as.POSIXct(frame[i,1], origin="1970-01-01"), label=T, abbr=F),
       cex=7)
  lapply(redEagle, function(x){
    if(length(x[timestamps(x)<=frame[i,1],])>2){
      lines(x[timestamps(x)<=frame[i,1],], lwd=3)
      lines(x[timestamps(x)<=frame[i,1],], lwd=1.5, col="firebrick")
    }})
  lapply(redEagle, function(x){
    points(x=coordinates(x)[1,1], 
           y=coordinates(x)[1,2],
           pch=16, cex=10/3)
    points(x=coordinates(x)[1,1], 
           y=coordinates(x)[1,2],
           pch=20, cex=10/3, col="indianred")
  })
  dev.off()
  unlink(paste(tempdir(), "/tmp", i, sep=""), recursive = TRUE)
  print(paste("Done with frame # ", i, ".", sep=""))
} 


