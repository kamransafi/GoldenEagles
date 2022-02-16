#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: energy landscape construction
#follows on from embc_segmentation.R
#the first attempt will only include static variables. reasons: 1)movebank annotation isn't working and 2)res of static variables is higher
#Jan 28. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move)
library(ctmm)
library(sf)
library(mapview)
library(parallel)
library(CircStats)
library(circular)
library(fitdistrplus)
library(raster)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/global_seascape/functions.R")

# STEP 1: open data and filter out non-commuting flights ----------------------------------------------------------------

load("embc_output.RData") #input

flight <- input %>% 
  filter(embc_clst_smth == 4) #use the smoothed values of clustering

save(flight, file = "flight_only.RData")

# #plot and see
# flight_sf <- flight %>% 
#   st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)
# mapview(flight_sf, zcol = "individual.local.identifier")


# STEP 2: variogram to decide on data resolution ----------------------------------------------------------------

#create move object
load("flight_only.RData") #flight
mv <- move(x = flight$location.long, y = flight$location.lat, time = flight$timestamp, proj = wgs, data = flight, animal = flight$individual.local.identifier)
tel <- as.telemetry(mv)
sv <- lapply(tel, variogram)

var <- sv[[23]]
plot(var, CTMM = variogram.fit(var, interactive = F),xlim=xlim)
plot(var, CTMM = variogram.fit(var, interactive = F), level = 0.5)

level <- c(0.5,0.95) # 50% and 95% CIs
xlim <- c(0,1 %#% "hour") # 0-12 hour window
plot(var,xlim=xlim,level=level)

#conclusion: the variogram looks weird... just go with 20 min for now....

# STEP 3: step selection prep- generate alternative steps ----------------------------------------------------------------

load("flight_only.RData") #flight
#create move object
mv <- move(x = flight$location.long, y = flight$location.lat, time = flight$timestamp, proj = wgs, data = flight, animal = flight$individual.local.identifier)


hr <- 20 #minutes; determine the sub-sampling interval
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
  
  Sys.time() - b 
  stopCluster(mycl) 
  
  save(sp_obj_ls, file = paste0("sl_", hr, "_min_25_ind.RData"))
  
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

  Sys.time() - b 
  stopCluster(mycl) 

save(used_av_track, file = paste0("alt_", n_alt, "_", hr, "_min_25_ind.RData")) #315843

# STEP 4: summary stats ----------------------------------------------------------------

load("alt_50_20_min_25_ind.RData") #used_av_track

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

# STEP 5: annotation: static ----------------------------------------------------------------

load("alt_50_20_min_25_ind.RData") #used_av_track

#manually annotate with static variables: elevation, slope, ruggedness (difference between the maximum and the minimum value of a cell and its 8 surrounding cells),
#unevenness in slope, aspect and elevation (TPI). (difference between the value of a cell and the mean value of its 8 surrounding cells)
#martina aggregates all to 100 m resolution

load("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/dem_wgs") #dem_wgs spatial res: 25 m

dem <- raster("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF")
slope <- terrain(dem, opt = "slope", unit = "degrees", filename = "slope.tif")
aspect <- terrain(dem, opt = "aspect", unit = "degrees", filename = "aspect.tif")
TRI <- terrain(dem, opt = "tri", unit = "degrees", filename = "TRI.tif")
TPI <- terrain(dem, opt = "TPI", unit = "degrees", filename = "TPI.tif")

#estimate slope and aspect unevenness
sl_uneven <- terrain(slope, opt = "TPI", unit = "degrees", filename = "slope_TPI.tif")
as_uneven <- terrain(aspect, opt = "TPI", unit = "degrees", filename = "aspect_TPI.tif")

#aggregate raster cells to 100 m resolution
dem_100 <- raster::aggregate(x = dem, fact = 4, filename = "dem_100.tif")
slope_100 <- raster::aggregate(x = slope, fact = 4, filename = "slope_100.tif")
aspect_100 <- raster::aggregate(x = aspect, fact = 4, filename = "aspect_100.tif")
sl_uneven_100 <- raster::aggregate(x = sl_uneven, fact = 4, filename = "slope_TPI_100.tif")
as_uneven_100 <- raster::aggregate(x = as_uneven, fact = 4, filename = "aspect_TPI_100.tif")
TRI_100 <- raster::aggregate(x = TRI, fact = 4, filename = "TRI_100.tif")
TPI_100 <- raster::aggregate(x = TPI, fact = 4, filename = "TPI_100.tif")

#create a stack
topo <- stack(list(dem_100, sl_uneven_100, as_uneven_100, TRI_100, TPI_100, slope_100, aspect_100))
topo_wgs <- projectRaster(topo, crs = wgs)

#extract values
st_ann <- cbind(used_av_track, extract(x = TPI_100, y = used_av_track[,c("location.long","location.lat")], method = "bilinear"))

save(st_ann, file = "")

########### old
s <- raster("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/slope_aspect.tif")

slope_wgs <-  projectRaster(slope, crs = wgs) 

save(slope_wgs, file = "/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/slope_wgs.RData")

#extract elevation values
post_em$dem_alt <- extract(x = dem_wgs, y = post_em[,c("location.long","location.lat")], method = "bilinear")



