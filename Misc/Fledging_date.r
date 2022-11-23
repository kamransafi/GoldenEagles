### Testing a fledging date finding algorithm
### Author: Kamran Safi
### Date: 22.11.22

# The idea is to use the GPS locations and based on a simple kmeans clustering 
# on the inter-day distances and the standard variation of the inter-day distances
# identify the day when a bigger than usual step or a larger than usual variation
# in distances was recorded. The inter-day distance is calculated as the distances
# for locations of a given day relative to the mean position in the previous day.

# The code excludes individuals with less than 850 locations and reverts
# to using means of the inter-day distances only, when no standard variance
# could be calculated. This happens when sampling rates were so low that 
# no sd can be calculated.

# Also the location data of only the first year is considered.

# BoxCox transformation is applied for better clustering results.

# Improvement: There could be another argument as distance to and variance
# of distances to the known nest location. Maybe in combination with 
# 3 cluster classes this can be extended to detection of fledging and
# emigration timing.
# The loops can be replaced by apply functions, but given the speed limitation
# being mainly data download this is not super crucial.
# Limiting download to the relevant dates will be a more urgent improvement. Idea
# here would be to use start of deployment and go from there.
# Another idea could be to include HDOP into evaluating the cluster assignment.

library(move)
library(lubridate)
library(adehabitatLT)
library(scales)
#MB login
creds <- movebankLogin()
#get animal names for serial processing
inds <- getMovebankAnimals("LifeTrack Golden Eagle Alps", creds)
#reduce to GPS only
inds <- inds[inds$sensor_type_id==653,]
#reduce to more than 850 locations and the names only
inds <- inds[inds$number_of_events>850,"local_identifier"]
#remove care adult birds
rm_inds <- c("Mellau21 (eobs 6988)", "Memmingen20 (eobs 7507)", "Aosta1_20 (eobs 7511)", "Aosta2_20 (eobs 7558)", "Aosta21 (eobs7590)")
inds <- inds[!inds %in% rm_inds]
#start a pdf for visual
pdf("/home/kami/Documents/Research/Projects/GoldenEagles/Analysis/FledingDates.pdf", height=7, width=7)
#go through the study by individual note: this can be paralellized with foreach
for(ind in inds){
  print(paste("Doing ", ind, ".", sep=""))
  #dld all the data for one individual
  eagleAll <- getMovebankData(study="LifeTrack Golden Eagle Alps", sensorID=653, 
                              animalName=ind, login=creds, 
                              removeDuplicatedTimestamps=T) 
  #restrict to year of tagging
  eagle <- eagleAll[year(timestamps(eagleAll))==min(year(timestamps(eagleAll))),]
  #calculate mean daily position
  Dcoords <- cbind(tapply(coordinates(eagle)[,1], yday(timestamps(eagle)), FUN=mean), tapply(coordinates(eagle)[,2], yday(timestamps(eagle)), FUN=mean))
  #count the number of year days for looping through
  DayCount <- unique(yday(timestamps(eagle)))
  #container for mean inter-day distances
  DM <- NULL
  #contained for sd of inter-day distances
  VarM <- NULL
  #loop through each year day. Can also be made more efficient
  for(i in 1:(length(DayCount)-1)){
    #calculate for each day the mean of inter-day distances
    tmp <- mean(spDistsN1(cbind(coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],1], coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],2]), Dcoords[i,], longlat=T)*1000)
    #calculate for each day the sd of inter-day distances
    tvar <- sd(spDistsN1(cbind(coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],1], coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],2]), Dcoords[i,], longlat=T)*1000)
    DM <- c(DM, tmp)
    VarM <- c(VarM, tvar)
  }
  #BoxCox transform of means
  bcD <- boxcox(DM~DayCount[-1], seq(-2, 2, length=1000), plotit = F)
  tD <- DM^bcD$x[which.max(bcD$y)]
  #apply a kmeans to the transformed mean inter-day distances
  DClus <- kmeans(tD[!is.na(tD)], 2)
  #if there were NAs in the data, put them back
  DClass <- rep(NA, length.out=length(tD))
  DClass[!is.na(tD)] <- DClus$cluster
  
  #The same procedure for the SD of inter-day distances
  bcVar <- boxcox(VarM~DayCount[-1], seq(-2, 2, length=1000), plotit = F)
  tVar <- VarM^bcVar$x[which.max(bcVar$y)]
  VarClus <- kmeans(tVar[!is.na(tVar)], 2)
  VarClass <- rep(NA, length.out=length(tVar))
  VarClass[!is.na(tVar)] <- VarClus$cluster
  
  #choose fledging date as the day either of the measures switched class
  FledgeD <- which(abs(diff(VarClass))+abs(diff(DClass))==1)[1]+1
  
  #if there are many NAs in the VarClass, base the fledging date on mean inter-day distance only
  #remember the criterion in Donly and use this for the plot.
  if(is.na(FledgeD)){
    FledgeD <- which(abs(diff(DClass))==1)[1]+1
    Donly <- TRUE
  }else{
    Donly <- FALSE
  }
  
  plot(eagle, pch=16, col=alpha("grey10", 0.5), type="l", bty="n", xlab="Longitude", ylab="Latitude")
  points(eagle[yday(timestamps(eagle))%in%DayCount[1:FledgeD],], pch=16, col=alpha("firebrick", 0.3), cex=0.3)
  points(eagle[yday(timestamps(eagle))%in%DayCount[1:FledgeD],], pch=1, col=alpha("firebrick", 0.7), cex=0.3)
  legend("topleft", legend=paste("Fledging date: ", as.Date(floor(1/2*(DayCount[FledgeD]+DayCount[FledgeD+1])), origin=paste(year(timestamps(eagle)[1]), "-01-01", sep=""), sep="")), bty="n")
  if(Donly==T){
    legend("bottomright", legend="Used distance only!", bty="n")
  }
  title(ind)
  print("Super success!")
  print("--------------------------------------------------")
}
dev.off()


#funX <- tapply(coordinates(eagle)[,1], yday(timestamps(eagle)), FUN=ecdf)
#funY <- tapply(coordinates(eagle)[,2], yday(timestamps(eagle)), FUN=ecdf)
# dX <- funX[[i]](coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],1])
# dY <- funY[[i]](coordinates(eagle)[yday(timestamps(eagle))==DayCount[i+1],2])
#eagle_trans <- spTransform(eagle, center=T)

# 
# DDist <- spDistsN1(Dcoords, Dcoords[1,])
# 
# seg <- lavielle(VarM, Lmin=7, Kmax=10)
# chooseseg(seg)
# fp <- findpath(seg, 4)
# fp
# 
# DNest <- spDistsN1(coordinates(eagle), coordinates(eagle)[50,])
# plot(distance(eagle))
# ts <- timestamps(eagle)
# CumD <- cumsum(distance(eagle))
# Dratio <- abs(cumsum(DNest)[-1]/cumsum(DNest)[-length(cumsum(DNest))])
# SlopE <- diff(DNest)/as.numeric(diff(timestamps(eagle), units="secs"))
# plot(cumsum(SlopE)~ts[-1], type="l")
# plot(CumD~ts, type="l")
# plot(log(cumsum(DNest))~ts, type="b")#, xlim=c(min(ts), 1.662e9),)
# plot(Dratio^-1~ts[-1], type="l")
# 
# spl <- smooth.spline(ts, CumD, df=3)
# lines(spl, col="red")
# d2 <- function(x) predict(spl, x , deriv=2)$y
# 
# library(adehabitatLT)
# eagleLT <- as(eagle_trans, "ltraj")
# i <- fpt(eagleLT, seq(5, 3000, length=300))
# VarFPT<- varlogfpt(i)
# plot(VarFPT)
# plot(i,125)
# 
# x <- seq(min(ts), max(ts), length=length(ts))
# plot(x, d2(x), type="l")
# abline(h=0, lty=2)
# 
# polyroot(f=d2(x))  # first inflection point
# uniroot(f=d2, interval=c(0.5, 0.9))  # second inflection point
