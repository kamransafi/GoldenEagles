simpleThin <- function(x){
  eagle <- getMovebankData(study, creds, animal=x, removeDuplicatedTimestamps=T)
  eagle1h <- eagle[!duplicated(cbind(year(timestamps(eagle)), month(timestamps(eagle)), day(timestamps(eagle)), hour(timestamps(eagle)))),]
  message(paste("Done with ", x, ".", sep=""))
  return(eagle1h)
}
