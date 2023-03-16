
# for both functions, acc timestamp and gps timestamp should be in the same time zone and both in UTC time

#Here GPS is the reference, the function returns the GPS data with associated, when available, the ID of the closest ACC event
ACCtoGPS <- function(ACCdata,GPSdata,
                     timeTolerance,
                     accEventCol=NULL,
                     ColsToAssociate=NULL,
                     ACCtimeCol="timestamp", GPStimeCol="timestamp"){
  require(data.table)
  if(is.null(accEventCol)){accEventCol <- grep("event_id|event.id", names(ACCdata), value=T)}
  if(all(class(ACCdata[,ACCtimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F) | 
     all(class(GPSdata[,GPStimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F)){
    stop("ACC or GPS time column is not in POSIXct format.")}
  #create the empty columns that I want to fill in during the loop
  GPSdata$acc_event_id <- NA
  GPSdata$diff_acc_time_s <- NA
  GPSdata$acc_closest_timestamp <- NA
  ACCdata <- ACCdata[ACCdata[,ACCtimeCol] > min(GPSdata[,GPStimeCol]) & ACCdata[,ACCtimeCol] < max(GPSdata[,GPStimeCol]),]
  if(nrow(ACCdata)==0){stop("There are no ACC data available in the GPS time range. Consider chacking that the two datasets are in the same time zone.")}
  GPSdata <- rbindlist(lapply(1:nrow(GPSdata), function(h){
    #create a subset of acc that occurr within a certain time interval from each gps point (+/- timeTolerance)
    gps.time <- GPSdata[h,GPStimeCol]
    acc.sub <- ACCdata[ACCdata[,ACCtimeCol] > gps.time-timeTolerance & ACCdata[,ACCtimeCol] < gps.time+timeTolerance,]
    #there could be no acc data in that close interval, but if there are some (nrow > 1) then we can associate them to the gps info:
    if(nrow(acc.sub) >= 1){
      timeDiff <- abs(difftime(acc.sub[,ACCtimeCol], GPSdata[h,GPStimeCol], units="secs")) #calculates the time difference
      min.diff <- which.min(timeDiff) #selects the row of the nearest point in time of the acc to each gps point (h)
      ## take the acc event id from the point corresponding to the minimum time difference and associate it to the gps point h
      GPSdata$acc_event_id[h] <-  acc.sub[min.diff, accEventCol]
      GPSdata$diff_acc_time_s[h] <- round(min(as.difftime(timeDiff, units="secs")), digits=3)
      GPSdata$acc_closest_timestamp[h] <- as.character(acc.sub[min.diff,ACCtimeCol])
    }
    return(GPSdata[h,])
  }))
  if(all(is.na(GPSdata$acc_closest_timestamp))){warning("No ACC data to associate to any of the GPS burst given this time tolerance. All associated ACC information set to NA. Consider increasing the time tolerance.")}
  GPSdata$acc_closest_timestamp <- as.POSIXct(GPSdata$acc_closest_timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  if(!is.null(ColsToAssociate)){
    GPSdata <- merge(GPSdata, ACCdata[,c(accEventCol,ColsToAssociate)], by.x="acc_event_id", by.y=accEventCol, all.x=T)
    GPSdata <- GPSdata[order(GPSdata$timestamp),]
  }
  return(as.data.frame(GPSdata))  #the function returns the GPS dataset with the additional ACC columns + a time.diff column
}


#Here ACC is the reference, the function returns the ACC data with associated, when available, the ID of the closest GPS event

GPStoACC <- function(GPSdata,ACCdata,
                     timeTolerance,
                     gpsEventCol=NULL,
                     ColsToAssociate=NULL,
                     ACCtimeCol="timestamp", GPStimeCol="timestamp"){
  require(data.table)
  if(is.null(gpsEventCol)){gpsEventCol <- grep("event_id|event.id", names(GPSdata), value=T)}
  if(all(class(ACCdata[,ACCtimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F) | 
     all(class(GPSdata[,GPStimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F)){
    stop("ACC or GPS time column is not in POSIXct format.")}
  #create the empty columns that I want to fill in during the loop
  ACCdata$gps_event_id <- NA
  ACCdata$diff_gps_time_s <- NA
  ACCdata$gps_closest_timestamp <- NA
  GPSdata <- GPSdata[GPSdata[,GPStimeCol] > min(ACCdata[,ACCtimeCol]) & GPSdata[,GPStimeCol] < max(ACCdata[,ACCtimeCol]),]
  if(nrow(GPSdata)==0){stop("There are no GPS data available in the ACC time range. Consider chacking that the two datasets are in the same time zone.")}
  ACCdata <- rbindlist(lapply(1:nrow(ACCdata), function(h){
    #create a subset of acc that occurr within a certain time interval from each gps point (+/- timeTolerance)
    acc.time <- ACCdata[h,ACCtimeCol]
    gps.sub <- GPSdata[GPSdata[,GPStimeCol] > acc.time-timeTolerance & GPSdata[,GPStimeCol] < acc.time+timeTolerance,]
    #there could be no acc data in that close interval, but if there are some (nrow > 1) then we can associate them to the gps info:
    if(nrow(gps.sub) >= 1){
      timeDiff <- abs(difftime(gps.sub[,GPStimeCol], ACCdata[h,ACCtimeCol], units="secs")) #calculates the time difference
      min.diff <- which.min(timeDiff) #selects the row of the nearest point in time of the acc to each gps point (h)
      ## take the gps event id from the point corresponding to the minimum time difference and associate it to the gps point h
      ACCdata$gps_event_id[h] <-  gps.sub[min.diff, gpsEventCol]
      ACCdata$diff_gps_time_s[h] <- round(min(as.difftime(timeDiff, units="secs")), digits=3)
      ACCdata$gps_closest_timestamp[h] <- as.character(gps.sub[min.diff,GPStimeCol])
    }
    return(ACCdata[h,])
  }))
  if(all(is.na(ACCdata$gps_closest_timestamp))){warning("No GPS data to associate to any of the ACC burst given this time tolerance. All associated GPS information set to NA. Consider increasing the time tolerance.")}
    ACCdata$gps_closest_timestamp <- as.POSIXct(ACCdata$gps_closest_timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
if(!is.null(ColsToAssociate)){
  ACCdata <- merge(ACCdata, GPSdata[,c(gpsEventCol,ColsToAssociate)], by.x="gps_event_id", by.y=gpsEventCol, all.x=T)
  ACCdata <- ACCdata[order(ACCdata$timestamp),]
}    
  return(as.data.frame(ACCdata))  #the function returns the ACC dataset with the additional GPS columns + a time.diff column
}
