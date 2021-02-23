setwd("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/")
library(move)
library(lubridate)
load("./EnvData/Items4Plots.Rdata")
DEMlr <- raster("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/DEMlr.grd")
hillBlurLR <- raster("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/hillBlurLR.grd")


jpeg("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/plot_roosts.jpg", height=1044, width=3740)
plot(t(studyArea), type="n", asp=1, xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
image(DEMlr^(1/1.5), col=alpha(terrain.colors(512), 1), asp=1, add=T)
plot(hillBlurLR, col=alpha(grey(0:512/521), 0.35), add=T, legend=FALSE, main='')
plot(wrld, add=T, lwd=10, border="grey70")
plot(wrld, add=T, lwd=5, border="grey10")
plot(lakes1, add=T, col="skyblue3", border="skyblue4")
lapply(Eagles[c(-14,-37,-38,-49)], function(x){
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=16, cex=10/3)
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=20, cex=10/3, col="indianred")
})
lapply(Eagles[c(14,37,38,49)], function(x){
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=16, cex=10/3)
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=20, cex=10/3, col="lavender")
})
dev.off()


jpeg("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/plot_tracks.jpg", height=1044, width=3740)
plot(t(studyArea), type="n", asp=1, xaxt="n", yaxt="n", ylab=NA, xlab=NA, bty="n")
image(DEMlr^(1/1.5), col=alpha(terrain.colors(512), 1), asp=1, add=T)
plot(hillBlurLR, col=alpha(grey(0:512/521), 0.35), add=T, legend=FALSE, main='')
plot(wrld, add=T, lwd=10, border="grey70")
plot(wrld, add=T, lwd=5, border="grey10")
plot(lakes1, add=T, col="skyblue3", border="skyblue4")
lapply(Eagles[c(-14,-37,-38,-49)], function(x){
  lines(x, lwd=3, col="black")
})
lapply(Eagles[c(-14,-37,-38,-49)], function(x){
  lines(x, lwd=1.5, col="firebrick")
})

lapply(Eagles[c(14,37,38,49)], function(x){
  lines(x, lwd=3)
})
lapply(Eagles[c(14,37,38,49)], function(x){
  lines(x, lwd=1.5, col="grey")
})
lapply(Eagles[c(-14,-37,-38,-49)], function(x){
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=16, cex=10/3)
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=20, cex=10/3, col="indianred")
})
lapply(Eagles[c(14,37,38,49)], function(x){
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=16, cex=10/3)
  points(x=coordinates(x)[1,1], 
         y=coordinates(x)[1,2],
         pch=20, cex=10/3, col="lavender")
})
dev.off()

