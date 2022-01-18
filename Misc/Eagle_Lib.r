# Change Projection in CH coordinates
ProjectCH <- function(x)
{
  Bproj <- spTransform(x, CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs"))
  return(Bproj)
}

# A simple outlier removal function based on velocity and an empirical distribution function fitted to the speed values measured.
# This function caculates for each location the distance to the previous and next point divided by the respective time lags.
# It returns the quantile of the realised velocity based on the entire data for each segment speed. Based on these quantiles 
# it is possible to identify which points are characterised by unusual speeds (jump back and forth)
OffIndex <- function(move, lag=1, time.adj=T, longlat=T)
{
  x <- coordinates(move)
  n <- n.locs(move)
  D1 <- matrix(cbind(x[,1][seq(1, n-(2*lag), 1)], 
                     x[,2][seq(1, n-(2*lag), 1)],
                     x[,1][seq((1+2*lag), n, 1)],
                     x[,2][seq((1+2*lag), n, 1)]), ncol=4)
  D2 <- matrix(cbind(x[,1][seq(1, n-(2*lag), 1)], 
                     x[,2][seq(1, n-(2*lag), 1)],
                     x[,1][seq((1+lag), (n-lag), 1)],
                     x[,2][seq((1+lag), (n-lag), 1)]), ncol=4)
  D3 <- matrix(cbind(x[,1][seq((lag+1), (n-lag) , 1)], 
                     x[,2][seq((lag+1), (n-lag), 1)],
                     x[,1][seq((2*lag+1), n, 1)],
                     x[,2][seq((2*lag+1), n, 1)]), ncol=4)
  
  if(longlat==F)
  {
    d13 <- sqrt(((D1[,3]-D1[,1])^2)+((D1[,4]-D1[,2])^2))
    d12 <- sqrt(((D2[,3]-D2[,1])^2)+((D2[,4]-D2[,2])^2))
    d23 <- sqrt(((D3[,3]-D3[,1])^2)+((D3[,4]-D3[,2])^2))
  }else{
    d13 <- apply(D1, 1 , function(x) spDistsN1(matrix(x[1:2], nrow=1), x[3:4], longlat=T)) * 1000
    d12 <- apply(D2, 1 , function(x) spDistsN1(matrix(x[1:2], nrow=1), x[3:4], longlat=T)) * 1000
    d23 <- apply(D3, 1 , function(x) spDistsN1(matrix(x[1:2], nrow=1), x[3:4], longlat=T)) * 1000
  }
  d123 <- d12+d23
  Tx <- timestamps(move)
  dT13 <- as.vector(difftime(Tx[seq((1+2*lag), n, 1)], Tx[seq(1, n-(2*lag), 1)], units="secs"))
  if(time.adj==T)
  {
    di <- ((d123/dT13)-(d13/dT13))
  }else{
    di <- (d123-d13)
  }
  edi <- ecdf(di)
  return(c(rep(NA, lag), edi(di), rep(NA, lag)))
}


# Downloads all movement data for a local copy of a study
localCopy <- function(x, study, login, ExFolder="./Data", NoWarnings=FALSE){
  if(!dir.exists(ExFolder)){
    dir.create(ExFolder)
  }
  if(NoWarnings==TRUE){
    options(warn = -1)
    anim <- suppressMessages(getMovebankData(study, creds, animal=x, removeDuplicatedTimestamps=T))
    options(warn = 1)
  }else{
    anim <- getMovebankData(study, creds, animal=x, removeDuplicatedTimestamps=T)
  }
  ExFile <- paste(ExFolder, "/", gsub(" ", "_", x),"@", gsub("-", "", Sys.Date()), ".rds", sep="")
  saveRDS(anim, file = ExFile)
  message(paste("Saved ", x, " in ", ExFile, ".", sep=""))
}


# Take the first location in every hour of the day (if present) and thin the data to at least 1 pos/h
simple1hThin <- function(x=NULL, study, creds, file.name=NULL){
  if(is.null(x) & is.null(file.name)){stop("Define either an Animal name or an input.")}
  if(!is.null(file)){
    eagle <- readRDS(file.name)
    mess <- file.name
  }else{
    eagle <- getMovebankData(study, creds, animal=x, removeDuplicatedTimestamps=T)
    mess <- x
  }
  eagle1h <- eagle[!duplicated(cbind(year(timestamps(eagle)), month(timestamps(eagle)), day(timestamps(eagle)), hour(timestamps(eagle)))),]
  message(paste("Done with ", mess, ".", sep=""))
  return(eagle1h)
}

# Distance daily accumulated. Returns a dataframe with distance, numbers of observations per day, day of the year, Year, days since start.
dailyDist <- function(x){
  DD <- tapply(distance(x), list(yday(timestamps(x))[-1], year(timestamps(x))[-1]), sum)
  n <- tapply(distance(x), list(yday(timestamps(x))[-1], year(timestamps(x))[-1]), length)
  d <- NULL
  y <- as.numeric(dimnames(DD)[[2]])
  for(i in 1:ncol(DD)){
    d <- rbind(d, data.frame(Distance=DD[,i][!is.na(DD[,i])],
                             n=n[,i][!is.na(DD[,i])],
                             yDay=as.numeric(dimnames(DD)[[1]][!is.na(DD[,i])]),
                             Year=rep(as.numeric(y[i]),length(DD[,i][!is.na(DD[,i])]))))
  }
  d <- d[order(d$Year, d$yDay),]
  d$cumDay <- yday(timestamps(x)[1])+as.numeric(unique(difftime(floor_date(timestamps(x)[-1], unit="days"), floor_date(timestamps(x)[1], unit="days"), units="days")))
  row.names(d) <- 1:nrow(d)
  print(paste("Done with ", as.character(x@idData$individual.local.identifier), "!", sep=""))
  return(d)
  
}


# Function to calculate the net squared distance in map units from a point of reference (defaults to first coordinates).
# Returns a data frame of squared distances and the timestamp.
nsd <- function(x, reference=coordinates(x)[1,]){
  d <- data.frame(N2D=spDistsN1(coordinates(x), reference, longlat = FALSE)^2, 
                  datetime=timestamps(x))
  print(paste("Done with ", as.character(x@idData$individual.local.identifier), "!", sep=""))
  return(d)
}

# Function to determine an excursion outside an area based on the point of reference in the nsd function
# and a distance threshold. The function returns the current state, and the durations of the excursions. 
excur <- function(x, threshold=5000){
  d <- sqrt(x[,1])
  if(max(d)<=threshold)
  {
    print("Tag never moved beyond threshold distance.")
    d <- NULL
  }else{
    d[which(d<=threshold)] <- 0
    d[which(d>threshold)] <- 1
    tStart <- x[which(diff(d)==1),2]
    tEnd <- x[which(diff(d)==-1),2]
    if(length(tStart)-length(tEnd)>0)
    {
      print(paste("Bird outside territory since: ", as.character(tail(tStart, 1)), sep=""))
      d <- data.frame(start=tStart, duration=c(unlist(lapply(1:length(tEnd), function(x) difftime(tEnd[x], tStart[x], units="days"))), NA))
    }
    if(length(tStart)-length(tEnd)==0)
    {
      print(paste("Still in birth territory as of: ", as.character(tail(x[,2], 1)), sep=""))
      d <- data.frame(start=tStart, duration=unlist(lapply(1:length(tEnd), function(x) difftime(tEnd[x], tStart[x], units="days"))))
    }
    if(length(tStart)-length(tEnd)<0)
    {
      message("Something is wrong!")
      d <- NULL
    }
  }
  return(d)
}

# Function to return 99% MCP area and the time difference between the first and last location.
# The projection needs to be in Euclidian space
MCParea <- function(x){
  A <- mcp(x, percent=99, unin="m", unout = "km2")
  dT <- difftime(tail(timestamps(x), 1), head(timestamps(x), 1), units="days")
  return(c(A$area, dT))
}

######### CONTINUE HERE ################
# function to select 1Hz sections and calculate metrics
moveDat1Hz <- function(x, topo=DEM){
  gs <- quantile(distance(x)[timeLag(x)==1], sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  crate <- quantile(diff(x@data$height.above.ellipsoid, lag=1)[timeLag(x)==1], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  HaG<- x@data$height.above.ellipsoid - extract(DEM, x, method="bilinear")
  hag <- quantile(HaG, probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  hasl <- quantile(x@data$height.above.ellipsoid[-1][timeLag(x)==1], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  hag <- quantile(HaG[-1][timeLag(x)==1], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  hag.plus2ms <- quantile(HaG[-1][timeLag(x)==1 & distance(x)>=2], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  hag.less1ms <- quantile(HaG[-1][timeLag(x)==1 & distance(x)<1], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
  df<- data.frame(rbind(as.vector(gs), as.vector(crate), as.vector(hasl), as.vector(hag), as.vector(hag.plus2ms), as.vector(hag.less1ms)))
  names(df) <- names(hag)
  row.names(df) <- c("Ground_speed", "Climb_rate", "Height_asl", "Height_ground", "Height_ground2ms", "Height_groundless1ms")
  return(df)
}

moveDatall <- function(x, topo=DEM){
  x <- x[!is.na(x@data$gps.satellite.count),]
  x <- x[x@data$gps.satellite.count >= 6,]
  if(n.locs(x)<5){
    df<- data.frame(matrix(rep(NA, 6*25), nrow=6))
  }else{
    gs <- quantile(distance(x) / timeLag(x), sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    crate <- quantile(diff(x@data$height.above.ellipsoid, lag=1), probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    HaG<- x@data$height.above.ellipsoid - extract(DEM, x, method="bilinear")
    hag <- quantile(HaG, probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    hasl <- quantile(x@data$height.above.ellipsoid, probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    hag <- quantile(HaG, probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    hag.plus2ms <- quantile(HaG[gs>=2], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    hag.less1ms <- quantile(HaG[gs<1], probs=sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99)), na.rm=T)
    df<- data.frame(rbind(as.vector(gs), as.vector(crate), as.vector(hasl), as.vector(hag), as.vector(hag.plus2ms), as.vector(hag.less1ms)))
  }
  names(df) <- paste0(sort(c(seq(0,1, length.out=21), 0.025, 0.975, 0.01, 0.99))*100, "%")
  row.names(df) <- c("Ground_speed", "Climb_rate", "Height_asl", "Height_ground", "Height_ground2ms", "Height_groundless1ms")
  return(df)
}
