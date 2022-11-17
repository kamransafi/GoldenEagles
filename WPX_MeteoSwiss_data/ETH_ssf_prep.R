#script for preparing the dataset to be sent to the ETH for annotation with COSMO-1 data.
#Sep 14. 2022. Konstanz, DE.
#Elham Nourani, PhD. enourani@ab.mpg.de

library(sf)
library(sp)
library(tidyverse)
library(move)
library(CircStats)
library(circular)
library(fitdistrplus)
library(parallel)
library(lubridate)
library(terra)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84") #replace this with utm N32

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
source("/home/enourani/ownCloud/Work/Projects/functions.R")


#initial work on March 2021 to develop the workflow
data <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/LifeTrack Golden Eagle Alps_march2021_gps.csv") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))


#STEP 1: reduce resolution to hr using dplyr
data_15min <- data %>% 
  drop_na(location.lat) %>% 
  mutate(yday = yday(timestamp),
         date_time = lubridate::round_date(timestamp, "15 minutes")) %>% 
  group_by(individual.local.identifier, yday, date_time) %>%  #the data is from one month in the same year, so only group by ind, day
  slice(1) %>% 
  ungroup()


# STEP 2: generate alternative steps (in the same height) also do in different heights? need a distribution of altitudes
#remove duplicated timestamps
#rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = data),"[",-1)) #get all but the first row of each set of duplicate rows

#data <- data[-rows_to_delete,]


mv <- move(x = data_15min$location.long, y = data_15min$location.lat, time = data_15min$timestamp, proj = wgs, 
           data = data_15min, animal = data_15min$individual.local.identifier)

hr <- 15 #minutes; determine the sub-sampling interval
tolerance <- 5 #minutes; tolerance for sub-sampling
n_alt <- 50 #number of alternative steps.

#prepare cluster for parallel computation
mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)

clusterExport(mycl, c("mv", "hr", "tolerance", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
})

(b <- Sys.time()) 

#used_av_ls <- parLapply(mycl, move_ls, function(group){ # for each group (ie. unique species-flyway combination)

#sp_obj_ls <- parLapply(mycl, split(mv), function(track){ #for each individual (i.e. track),
  
sp_obj_ls <- lapply(split(mv), function(track){ 
  
  #--STEP 1: thin the data 
  track_th <- track %>%
    thinTrackTime(interval = as.difftime(hr, units = 'mins'),
                  tolerance = as.difftime(tolerance, units = 'mins')) #the unselected bursts are the large gaps between the selected ones
  
  #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
  track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
  track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
  
  if(nrow(track_th@data) == 1){
    track_th@data$burst_id <- track_th$burst_id
  } else {for(i in 2:nrow(track_th@data)){
    #if(i == nrow(track_th@data)){
    #  track_th@data$burst_id[i] <- NA #why?!
    #} else
    if(track_th@data[i-1,"selected"] == "selected"){
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"]
    } else {
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"] + 1
    }
  }
  }
  
  #convert back to a move object (from move burst)
  track_th <- as(track_th,"Move")
  
  #--STEP 3: calculate step lengths and turning angles 
  #sl_ and ta_ calculations should be done for each burst.
  burst_ls <- split(track_th, track_th$burst_id)
  burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations
  
  burst_ls <- lapply(burst_ls, function(burst){
    burst$step_length <- c(move::distance(burst),NA) 
    burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
    burst
  })
  
  #put burst_ls into one dataframe
  bursted_sp <- do.call(rbind, burst_ls)
  
  #reassign values
  if(length(bursted_sp) >= 1){
    bursted_sp$track <- track@idData$track
    bursted_sp$individual.taxon.canonical.name <- track@idData$individual.taxon.canonical.name
    bursted_sp$individual.local.identifier <- track@idData$individual.local.identifier
  }
  
  bursted_sp
  
}) %>% 
  Filter(function(x) length(x) > 1, .) #remove segments with no observation

Sys.time() - b 
stopCluster(mycl) #n = 53

#STEP 3: convert height above ellipsoid to m asl
dem <- rast("Alps_east_west_dem.tif")

dem_wgs <- dem %>% 
  project("+proj=longlat +datum=WGS84 +no_defs")

data_elev <- lapply(sp_obj_ls, function(x){
  x_sp <- as(x, "SpatVector")
  x_sp$dem_alt <- terra::extract(x = dem_wgs, y = x_sp, method = "bilinear")[,2] 
  x_sp <- x_sp %>% 
    as("data.frame") %>% 
    drop_na("height.above.ellipsoid") %>% 
    mutate(flight_h = height.above.ellipsoid - dem_alt) %>% 
    filter(flight_h >= 0)
  x_sp
})

save(data_elev, file = paste0("march_2021_sl_", hr, ".RData"))


#--STEP 4: estimate step length and turning angle distributions
#put everything in one df
bursted_df <- data_elev %>%  
  reduce(rbind) %>% 
  as.data.frame()

#estimate von Mises parameters for turning angles
#calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))

#estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")

#plot turning angle and step length distributions
#pdf(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ssf_figs/ta_sl_", hr, "_", tolerance, ".pdf"))
par(mfrow=c(1,2))
hist(sl, freq=F, main="", xlab = "Step length (km)")
plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
#dev.off()

#diagnostic plots for step length distribution
#pdf(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ssf_figs/diagnostics_", hr, "_", tolerance, ".pdf"))
plot(fit.gamma1)
#dev.off()




#--STEP 5: produce alternative steps
#prepare cluster for parallel computation

#remove tracks with fewer than 3

mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)

clusterExport(mycl, c("mv", "hr", "mu", "kappa", "fit.gamma1", "tolerance", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
})

(b <- Sys.time()) 
used_av_track <- parLapply(mycl, data_elev, function(track){ #for each track
  
  #remove bursts that are shorter than 3 points
  burst_ls <- split(track,track$burst_id) %>%  Filter(function(x) nrow(x) >= 3, .) 
  
  lapply(burst_ls,function(burst){ #for each burst,
    #assign unique step id
    burst$step_id <- 1:nrow(burst)
    
    
    lapply(c(2:(nrow(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
      
      current_point<- burst[this_point,]
      previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
      used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
      
      #calculate bearing of previous point
      prev_bearing <- NCEP.loxodrome.na(previous_point$location.lat, current_point$location.lat,
                                        previous_point$location.long, current_point$location.long)
      
    
      coordinates(current_point) <- ~ location.long + location.lat
      proj4string(current_point) <- wgs
      
      current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
      
      #randomly generate n alternative points
      rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                        step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
        #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
        mutate(lon = current_point_m@coords[,1] + step_length*cos(turning_angle),
               lat = current_point_m@coords[,2] + step_length*sin(turning_angle))
      
      
      #covnert back to lat-lon proj
      rnd_sp <- rnd
      coordinates(rnd_sp) <- ~lon+lat
      proj4string(rnd_sp) <- meters_proj
      rnd_sp <- spTransform(rnd_sp, wgs)
      
      #check visually
      # mapview(current_point, color = "red") + mapview(previous_point, color = "orange") + mapview(used_point, color = "yellow") + mapview(rnd_sp, color = "black", cex = 0.5)
      
      #put used and available points together
      df <- used_point %>%  
        slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
        mutate(location.long = c(head(location.long,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
               location.lat = c(head(location.lat,1),rnd_sp@coords[,2]),
               turning_angle = c(head(turning_angle,1),deg(rnd_sp$turning_angle)),
               step_length = c(head(step_length,1),rnd_sp$step_length),
               used = c(1,rep(0,n_alt)))  %>%
        #dplyr::select(-c("coords.x1","coords.x2")) %>% 
        rowwise() %>% 
        mutate(heading = NCEP.loxodrome.na(lat1 = current_point@coords[,2], lat2 = location.lat, lon1 = current_point@coords[,1], lon2 = location.long)) %>% 
        as.data.frame()
      
      df
      
    }) %>% 
      reduce(rbind)
    
  }) %>% 
    reduce(rbind)
  
}) %>% 
  reduce(rbind)

Sys.time() - b 
stopCluster(mycl) 


used_av_track <- used_av_track %>% 
  mutate(stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))

saveRDS(used_av_track, file = "march_2021_used_av.rds") 

#sanity check
dd <- used_av_track [1:300,]
coordinates(dd) <-~ location.long + location.lat
proj4string(dd) <- wgs
mapview(dd, zcol = "used")

#Change format based on Michael's requirements:
#{year} {month}  {day}  {hour}  {min}  {sec}  {longitude}  {latitude}  {height above sea level} ...


used_av_track <- readRDS("march_2021_used_av.rds")

data_to_send <- used_av_track %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(year = year(timestamp),
         month = month(timestamp),
         day = day(timestamp),
         hour = hour(timestamp),
         min = minute(timestamp),
         sec = second(timestamp)) %>% 
  rename(latitude = location.lat,
         longitude = location.long,
         height.above.sea.level = flight_h) %>% 
  dplyr::select(c("year", "month", "day", "hour", "min", "sec", "longitude", "latitude", "height.above.sea.level", 35:37, 39, 43:44))


#save as ascii  with whitespaces as column separator (no commas!)

write.table(data_to_send, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/COSMO_wind/march_2021/golden_ealge_march_2021_tracks.ascii",
            row.names = F, col.names = T, sep = " ", append = T, quote = F)


