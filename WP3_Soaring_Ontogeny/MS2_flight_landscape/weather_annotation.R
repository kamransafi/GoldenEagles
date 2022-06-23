#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: annotation with atmospheric variables
#follows on from data_processing_&_annotation.R
#the above code used a time stamp of 20 minutes, because the res of topographic variables is high. Here, redo random step generation with hourly intervals.
#Mar 10. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(move)
library(parallel)
library(CircStats)
library(circular)
library(fitdistrplus)
library(raster)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
source("/home/enourani/ownCloud/Work/Projects/functions.R")

# STEP 1: step selection prep- generate alternative steps ----------------------------------------------------------------

load("flight_only_70ind.RData") #flight; from embc_segmentation.R #n = 56

#remove problematic individuals (i.e. adults)
load("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/em_fl_dt_recurse_70ind.RData") #emig_fledg_dates
inds_to_remove <- emig_fledg_dates %>% 
  mutate(post_fledging_duration = difftime(emigration_dt,fledging_dt, units = "days")) %>% 
  filter(post_fledging_duration <= 30 )

flight <- flight %>%
  filter(!(individual.local.identifier %in% inds_to_remove$individual.local.identifier)) #n_distinct = 49

#create move object
mv <- move(x = flight$location.long, y = flight$location.lat, time = flight$timestamp, proj = wgs, data = flight, animal = flight$individual.local.identifier)

hr <- 60 #minutes; determine the sub-sampling interval
tolerance <- 10 #minutes; tolerance for sub-sampling
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

sp_obj_ls <- parLapply(mycl, split(mv), function(track){ #for each individual (i.e. track),
  
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

Sys.time() - b #2.4 mins
stopCluster(mycl) 

save(sp_obj_ls, file = paste0("sl_", hr, "_min_70_ind.RData"))

#--STEP 4: estimate step length and turning angle distributions
#put everything in one df
bursted_df <- sp_obj_ls %>%  
  reduce(rbind) %>% 
  as.data.frame() %>% 
  dplyr::select(-c("coords.x1","coords.x2"))

#estimate von Mises parameters for turning angles
#calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))

#estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")

#plot turning angle and step length distributions
pdf(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ssf_figs/ta_sl_", hr, "_", tolerance, ".pdf"))
par(mfrow=c(1,2))
hist(sl, freq=F, main="", xlab = "Step length (km)")
plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
dev.off()

#diagnostic plots for step length distribution
pdf(paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ssf_figs/diagnostics_", hr, "_", tolerance, ".pdf"))
plot(fit.gamma1)
dev.off()

#--STEP 5: produce alternative steps
#prepare cluster for parallel computation
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
used_av_track <- parLapply(mycl, sp_obj_ls, function(track){ #for each track
  
  lapply(split(track,track$burst_id),function(burst){ #for each burst,
    
    #assign unique step id
    burst$step_id <- 1:nrow(burst)
    
    
    lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
      
      current_point<- burst[this_point,]
      previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
      used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
      
      #calculate bearing of previous point
      #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
      prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                        previous_point@coords[,1], current_point@coords[,1])
      
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
      df <- used_point@data %>%  
        slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
        mutate(location.long = c(head(coords.x1,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
               location.lat = c(head(coords.x2,1),rnd_sp@coords[,2]),
               turning_angle = c(head(turning_angle,1),deg(rnd_sp$turning_angle)),
               step_length = c(head(step_length,1),rnd_sp$step_length),
               used = c(1,rep(0,n_alt)))  %>%
        dplyr::select(-c("coords.x1","coords.x2")) %>% 
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

Sys.time() - b #40 sec
stopCluster(mycl) 


used_av_track <- used_av_track %>% 
  mutate(stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))

save(used_av_track, file = paste0("alt_", n_alt, "_", hr, "_min_78_ind.RData"))

# STEP 2: prep for movebank annotation ----------------------------------------------------------------

#create a dataframe with movebank specs
df <- used_av_track %>% 
    mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
    as.data.frame()

#rename columns
colnames(df)[c(2,3)] <- c("location-long","location-lat")

write.csv(df, paste0("alt_", n_alt, "_", hr, "_min_78_ind.csv")) #n = 48

#annotation request fails, remove height above ellipsoid to see if it fixes it

df <- read.csv(paste0("alt_", n_alt, "_", hr, "_min_78_ind.csv"))

df_no_na <- df %>% 
  drop_na()

colnames(df_no_na)[c(3,4)] <- c("location-long","location-lat")

write.csv(df_no_na, paste0("alt_", n_alt, "_", hr, "_min_78_ind_no_na.csv"), row.names = F) #n = 48

# STEP 3: open and process annotated file ----------------------------------------------------------------

ann <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/annotations/Apr_18_22/alt_50_60_min_25_ind.csv-3274898828578094969.csv") %>% 
  rename(u_wind = ECMWF.ERA5.PL.U.Wind,
         v_wind = ECMWF.ERA5.PL.V.Wind) %>% 
  mutate(wind_support= wind_support(u = u_wind, v = v_wind, heading = heading),
         cross_wind= cross_wind(u = u_wind, v = v_wind, heading = heading),
         wind_speed = sqrt(u_wind^2 + v_wind^2),
         abs_cross_wind = abs(cross_wind(u = u_wind, v = v_wind, heading = heading)))

# STEP 4: annotate with day since emigration ----------------------------------------------------------------
#code copied from data_processing_&_annotation.R

#open emigration information
## Hester's files on matching the names
load("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/from_Hester/eagle_names.RData") #eagle_names

#open file with info on fledging and emigration timing (from Svea)
dates <- read.csv("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/data/Goldeneagles_emigration_time_10_2021.csv",
                  stringsAsFactors = F, fileEncoding = "latin1") %>%  
  rowwise() %>% 
  mutate(fledging_timestamp = paste(paste(strsplit(date_fledging, "\\.") %>% map_chr(., 3), #yr 
                                          strsplit(date_fledging, "\\.") %>% map_chr(., 2), #mnth
                                          strsplit(date_fledging, "\\.") %>% map_chr(., 1), sep = "-"),  #day
                                    time_fledging, sep = " "),
         emigration_timestamp = ifelse(is.na(date_emigration), NA , 
                                       paste(paste(strsplit(date_emigration, "\\.") %>% map_chr(., 3), #yr 
                                                   strsplit(date_emigration, "\\.") %>% map_chr(., 2), #mnth
                                                   strsplit(date_emigration, "\\.") %>% map_chr(., 1), sep = "-"),  #day
                                             time_emigration, sep = " "))) %>% 
  ungroup() %>% 
  mutate(fledging_timestamp = as.POSIXct(strptime(fledging_timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         emigration_timestamp = as.POSIXct(strptime(emigration_timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  full_join(eagle_names, by = c("id" = "age_name")) %>% 
  as.data.frame()

dyn_ann <- lapply(split(ann, ann$individual.local.identifier), function(x){
  
  ind_dates <- dates %>% 
    filter(local_identifier == unique(x$individual.local.identifier))
  
  x <- x %>% 
    mutate(days_since_emig = difftime(timestamp,ind_dates$emigration_timestamp, units = c("days")),
           days_since_fled = difftime(timestamp,ind_dates$fledging_timestamp, units = c("days")))
  
  x
}) %>% 
  reduce(rbind)

saveRDS(dyn_ann, file =  "alt_50_60_min_25_ind_dyn_time_ann.RData")

