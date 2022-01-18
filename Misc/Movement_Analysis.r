# Need to add some comments and make it more efficient

###############################
DEM <- raster("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/EnvData/SRTM/DEM.grd")
studyArea <- bbox(extent(c(6.147962, 15.512771, 45.388163, 48.016078)))



IndDist <- unlist(lapply(biRDS, function(x) sum(distance(readRDS(x))))) #total distance
#last positions
lastPos <- lapply(biRDS, function(x) tail(timestamps(readRDS(x)),1))
Cont <- unlist(lapply(lastPos, function(x) difftime(Sys.time(), x, units="days")))
Cont[grep("17", biRDS)]
Cont[grep("18", biRDS)]
Cont[grep("19", biRDS)]
Cont[grep("20", biRDS)]
TrackDays <- unlist(lapply(biRDS, function(x){
  tmp <- readRDS(x)
  ceiling(as.numeric(difftime(tail(timestamps(tmp),1), head(timestamps(tmp),1), units="days")))
}))


NSD <- lapply(biRDS,  function(x) nsd(readRDS(x)))
ex <- lapply(NSD, excur)
classEx <- lapply(ex, function(x) kmeans(sqrt(x[-nrow(x),2]), 2))

###############################
DD <- lapply(biRDS, function(x) dailyDist(readRDS(x))) #
unlist(lapply(DD, function(x) median(x$Distance)))
dat <- NULL
dat$x <- sqrt(unlist(lapply(DD, function(x) mean(x$Distance))))
dat$y <- log(unlist(lapply(DD, function(x) mean(sum(x$n)*nrow(x)))))
mod <- gam(x~s(y), data=dat)
summary(mod)
py <- predict(mod, newdata=data.frame(y=log(seq(30000, 5e8, length.out=1000))))
px <- seq(30000, 5e8, length.out=1000)
pdf("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/DD.pdf", height=5, width=5)
plot(sqrt(unlist(lapply(DD, function(x) mean(x$Distance))))~unlist(lapply(DD, function(x) mean(sum(x$n)*nrow(x)))), log="x", 
     xlab="Number of positions * tracking days", ylab=expression(paste("Daily distance [ ", sqrt(m), " ]", sep="")), pch=20, bty="n", col="grey")
points(sqrt(unlist(lapply(DD, function(x) mean(x$Distance))))~unlist(lapply(DD, function(x) mean(sum(x$n)*nrow(x)))),
       pch=1)
lines(py~px, lwd=1.5)
dev.off()

movePara1Hz <- lapply(biRDS, function(x) moveDat1Hz(readRDS(x), topo=DEM))
moveParaall <- lapply(biRDS, function(x) moveDatall(readRDS(x), topo=DEM))

medGS1 <- (unlist(lapply(movePara1Hz, function(x) x[1,13])))
medGS2 <- (unlist(lapply(moveParaall, function(x) x[1,13])))
medCR1 <- (unlist(lapply(movePara1Hz, function(x) x[2,13])))
medCR2 <- (unlist(lapply(moveParaall, function(x) x[2,13])))
peakGS1 <- (unlist(lapply(movePara1Hz, function(x) x[1,24])))
peakGS2 <- (unlist(lapply(moveParaall, function(x) x[1,24])))
peakCR1 <- (unlist(lapply(movePara1Hz, function(x) x[2,24])))
peakDR1 <- (unlist(lapply(movePara1Hz, function(x) x[2,2])))
pdf("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/Boxplot.pdf", height=7, width=7)
par(mar=c(8, 5, 4, 2))
boxplot(cbind(medGS1, peakGS1, medGS2, peakGS2, medCR1, peakCR1, peakDR1), outline=F, 
        ylab=expression(paste("Speed m" , s^-1,  sep="")), 
        names=c("Median GS 1Hz", "Peak GS 1Hz", "Median GS HQ", "Peak GS HQ", "Median CR 1Hz", "Peak CR 1Hz", "Peak DR 1Hz"), las=2,
        ylim=c(-10,30))
abline(h=0)
dev.off()

asl1 <- (unlist(lapply(movePara1Hz, function(x) x[3,13])))
asl2 <- (unlist(lapply(moveParaall, function(x) x[3,13])))
ag1 <- (unlist(lapply(movePara1Hz, function(x) x[4,13])))
ag2 <- (unlist(lapply(moveParaall, function(x) x[4,13])))
ag2ms1 <- (unlist(lapply(movePara1Hz, function(x) x[5,13])))
ag2ms2 <- (unlist(lapply(moveParaall, function(x) x[5,13])))
agstat1 <- (unlist(lapply(movePara1Hz, function(x) x[6,13])))
agstat2 <- (unlist(lapply(moveParaall, function(x) x[6,13])))

pdf("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/BoxplotHeights.pdf", height=7, width=7)
par(mar=c(9, 5, 4, 2))
par(mfrow=c(1,2))
boxplot(cbind(asl1, asl2), ylab="Altitude a.s.l. in m", outline=F,
        names=c("Median height 1Hz", "Median height HQ"),
        , las=2)

boxplot(cbind(ag1, ag2, ag2ms1, ag2ms2, agstat1, agstat2), outline=F, 
        ylab="Height a.g.l. in m", 
        names=c("Median height 1Hz", 
        "Median height HQ", "Height @motion 1Hz", "Height @motion HQ", 
        "Height stationary 1Hz", "Height stationary HQ"), las=2)

dev.off()

DD.all <- cbind(do.call(rbind, DD), ID=unlist(lapply(1:length(DD), function(x) rep(x, nrow(DD[[x]])))))

  mod <- gamm(log(Distance)~n+s(cumDay), data=DD.all[DD.all$cumDay>199,], random=list(ID=~1))
  print(summary(mod$gam))
  plot(mod$gam, residuals=T, main=i)

################################

mean(unlist(lapply(DD, function(x) median(unlist(x), na.rm=T))))
mean(unlist(lapply(DD, function(x) max(unlist(x), na.rm=T))))
library(adehabitatHR)
HRsize <- lapply(biRDS, function(x){
  tmp <- readRDS(x)
  ts <- floor(difftime(timestamps(tmp), head(timestamps(tmp), 1), units = "days"))
  DDay <- unique(ts)[which(unlist((lapply(unique(ts), function(x) n.locs(tmp[ts==x,]))))>5)]
  tmpB <- lapply(DDay, function(y) {
    br <- tmp[ts==y,]
    if(n.locs(br)>5){
      ma <- cbind(mcp.area(br, percent=95, unout="m2"), y)
    }
    return(ma)})
  arR <- do.call("rbind", tmpB)
  arR$ID <- gsub(" ", "_", sub("./Data/", "", x))
  return(arR)
})

flatHR <- do.call("rbind", HRsize)
pdf("/home/kami/Documents/Projects/Aktiv/GoldenEagles/Analysis/DailyHR.pdf", height=7, width=7)
plot(c(1,1), type="n", ylim=range(flatHR$a), xlim=c(0,180), bty="n", xlab="Days since tagging", ylab="Daily MCP size", log="y")
for(i in unique(flatHR$ID)){
  lines(a~y, data=flatHR[flatHR$ID==i & flatHR$y<180,], lwd=0.75)
}
dev.off()

mcpA <- lapply(1:length(birds), MCParea)
AD <- matrix(unlist(mcpA), ncol=2, byrow=T)
