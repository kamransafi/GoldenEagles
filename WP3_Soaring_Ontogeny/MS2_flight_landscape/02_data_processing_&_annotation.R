#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: random step generation and annotation
#follows on from embc_segmentation.R
#the first attempt will only include static variables. res of static variables is higher
#Jan 28. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move)
#library(ctmm)
library(sf)
library(mapview)
library(parallel)
library(CircStats)
library(circular)
library(fitdistrplus)
library(terra)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84") #replace this with utm N32

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
source("/home/enourani/ownCloud/Work/Projects/functions.R")


# STEP 1: open data and filter out non-commuting flights ----------------------------------------------------------------

bc_output <- readRDS("embc_output_20min_61ind.rds") #n = 61.. some of these have very few data points

flight <- bc_output %>% 
  filter(embc_clst_smth %in% c(3,4)) #use the smoothed values of clustering. include both levels 3 and 4

saveRDS(flight, file = "flight_only_20min_56ind.rds") #n = 56

#did I filter for post dispersal!??? yes!!

# STEP 2: variogram to decide on data resolution ----------------------------------------------------------------

flight <- readRDS("flight_only_20min_56ind.rds")

tel <- as.telemetry(flight)
sv <- lapply(tel, variogram)

xlim <- c(0,60 %#% "minute") # 0-12 hour window

var <- sv[[23]]
plot(var, CTMM = variogram.fit(var, interactive = F),xlim = xlim)
plot(var, CTMM = variogram.fit(var, interactive = F), level = 0.5)

level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,1 %#% "hour") # 0-12 hour window
plot(var,xlim=xlim,level=level)

#conclusion: the variogram looks weird... just go with 20 min for now....

# STEP 3: step selection prep- generate alternative steps ----------------------------------------------------------------

flight <- readRDS("flight_only_20min_56ind.rds")

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
      burst$step_length <- c(NA, move::distance(burst)) #the NA used to be at the end, but it makes more sense for it to be at the beginning. it matches how the random sls are allocated
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
  
  Sys.time() - b #1.87 mins
  stopCluster(mycl) #n = 55
  
  saveRDS(sp_obj_ls, file = paste0("sl_", hr, "_min_55_ind2.rds")) #tolerance of 10 min
  
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
  
  sp_obj_ls <- readRDS("sl_60_min_55_ind2.rds")
  
  #prepare cluster for parallel computation
  mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)
  
  clusterExport(mycl, c("sp_obj_ls", "hr", "mu", "kappa", "fit.gamma1", "tolerance", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call
  
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
          mutate(location.long = c(head(coords.x1,1), rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
                 location.lat = c(head(coords.x2,1), rnd_sp@coords[,2]),
                 turning_angle = c(head(turning_angle,1), deg(rnd_sp$turning_angle)),
                 step_length = c(head(step_length,1), rnd_sp$step_length),
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

  Sys.time() - b # 
  stopCluster(mycl) 

  
used_av_track <- used_av_track %>% 
  mutate(stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))
  
saveRDS(used_av_track, file = paste0("alt_", n_alt, "_", hr, "_min_55_ind2.rds")) #the diff between this round and the previous is the assignment of NA when calculating step lengths. used to be end of burst, now is start of burst 

# STEP 4: summary stats ----------------------------------------------------------------

load("alt_50_60_min_55_ind2.rds") #used_av_track

used_av_track %>% 
  mutate(yr_mn = paste(year(timestamp), month(timestamp), sep = "_"),
         yr_day = paste(year(timestamp), yday(timestamp), sep = "_"),
         yr = year(timestamp)) %>% 
  group_by(individual.local.identifier) %>%
  summarise(n_d = n_distinct(yr_day),
            n_m = n_distinct(yr_mn),
            n_yr = n_distinct(yr),
            min_yr = min(yr)) %>% 
  ungroup() %>% 
  arrange(n_d) %>% 
  as.data.frame()
  
#data quantity per month
mn_summary <- used_av_track %>% 
  mutate(mn = month(timestamp)) %>% 
  group_by(individual.local.identifier,mn) %>% 
  summarise(data = n())


barplot(names.arg = mn_summary$mn, height = mn_summary$data, col = as.factor(mn_summary$individual.local.identifier), beside = F)

# STEP 5: annotation: life stages ----------------------------------------------------------------

used_av_track <- readRDS("alt_50_60_min_55_ind2.rds") #used_av_track

emig_dates <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/fleding_emigration_timing_Mar2023.rds")

used_av_track <- used_av_track %>% 
  left_join(emig_dates %>% dplyr::select(c("individual.local.identifier", "emigration_dt")), by = "individual.local.identifier") %>% 
  mutate(days_since_emig = difftime(timestamp,emigration_dt, units = c("days")) %>% as.numeric() %>%  ceiling(),
         weeks_since_emig = difftime(timestamp,emigration_dt, units = c("weeks")) %>% as.numeric() %>%  ceiling())


# STEP 6: annotation: static ----------------------------------------------------------------

#manually annotate with static variables: elevation, terrain ruggedness index (difference between the maximum and the minimum value of a cell and its 8 surrounding cells),
#and distance to ridge. TRI and distance to ridge were prepared by Louise.
#try with 100 m resolution

#100m layers were built in the earlier version of this script
TRI_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/TRI_100_LF.tif")
ridge_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100.tif")
dem_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/dem_100_LF.tif")

#create a stack using raster paths
topo <- c(dem_100, TRI_100, ridge_100)
names(topo) <- c("dem_100", "TRI_100", "ridge_100")

#reproject tracking data to match topo, extract values from topo, convert back to wgs and save as a dataframe
topo_ann_df <- used_av_track %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_transform(crs = crs(topo)) %>% 
  extract(x = topo, y = ., method = "simple", bind = T) %>%
  terra::project(wgs) %>% 
  data.frame(., geom(.)) %>% 
  dplyr::select(-c("geom", "part", "hole")) %>% 
  rename(location.long = x,
         location.lat = y)
  
saveRDS(topo_ann_df, file = "alt_50_60_min_55_ind_static_r_100.rds")

