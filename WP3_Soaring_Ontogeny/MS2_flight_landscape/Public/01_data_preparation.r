#Script for preparing golden eagle tracking data for step-selection analysis as reported in Nourani et al. 2024, eLife (https://doi.org/10.7554/eLife.98818.1)
#script 1/2
# Elham Nourani, PhD. 12.07.2024
# enourani@ab.mpg.de

library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(move2)
library(EMbC)
library(CircStats)
library(circular)
library(fitdistrplus)
library(parallel)
library(mapview)
library(units)

##### STEP 0: open tracking data, functions and define variables #####
data <- read.csv("GPS_data.csv") %>%  #this is the post-dispersal data. weeks_since_emig is the column representing weeks since start of dispersal
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

#source the functions from the github repository
source("functions.r")

meters_proj <- st_crs("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +type=crs")
wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")

##### STEP 1: EmbC segmentation #####

# 1.1: estimate flight altitude -----

#open geoid and dem layers
geo <- rast("EGM96_us_nga_egm96_15.tif")
dem <- rast("dem_Alpine_region.tif")

#extract elevation and geoid values for each tracking point 
data_h <- data %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = "+proj=longlat +datum=WGS84 +no_defs") %>% #create a spatial object
  extract(x = dem, y = ., method = "simple", bind = T) %>% 
  extract(x = geo, y = ., method = "simple", bind = T) %>% 
  data.frame(., geom(.)) %>% 
  dplyr::select(-c("geom", "part", "hole")) %>% 
  rename(dem = eu_dem_v11_E40N20) %>% 
  #calculate flight altitude as height above mean sea level - dem. height above sea level is height above ellipsoid - geoid
  mutate(height_msl = height.above.ellipsoid - geoid_undulation) %>% 
  mutate(height_ground = height_msl - dem)

# 1.2: calculate ground speed -----
#create a move2 object and calculate speed using mt_speed()
mv <- mt_as_move2(data_h %>% arrange(individual.local.identifier, timestamp), time_column = "timestamp", 
                  coords = c("x", "y"), crs = "+proj=longlat +datum=WGS84 +no_defs",
                  track_id_column = "individual.local.identifier", track_attributes = c("tag.local.identifier", "individual.taxon.canonical.name")) %>% 
  mutate(ground_speed = mt_speed(.)) %>% 
  drop_na(ground_speed, height_ground) %>% #remove NA values for speed
  filter(between(height_ground, quantile(height_ground, 0.005), quantile(height_ground, 0.999)) & ground_speed < quantile(ground_speed, 0.999))
  
# 1.3: EmbC segmentation -----
#create a matrix of flight height and speed
m <- data.matrix(mv %>% as.data.frame() %>% select(ground_speed, height_ground))

#call embc
bc <- embc(m)

#investigate the bc (Garriga et al 2016; S2)
X11();sctr(bc)

bc_smth <- smth(bc, dlta = .7)
X11();sctr(bc_smth)

#append cluster labels (1:LL, 2:LH, 3:HL, and 4:HH) to original data
mv$embc_clst <- bc@A
mv$embc_clst_smth <- bc_smth@A

#select a sample track and visualize
smpl <- mv %>% 
  filter(individual.local.identifier == "Viluoch17 (eobs 4570)") %>% 
  slice(1:500)

lines <- smpl %>% 
  mt_track_lines() 

mapview(lines, color = "gray") + mapview(as(smpl, "sf"), zcol = "embc_clst_smth")

#only look at levels 3 and 4
smpl2 <- smpl %>% 
  filter(embc_clst_smth %in% c(3,4))

lines2 <- smpl2 %>% 
  mt_track_lines() 

mapview(lines2, color = "gray") + mapview(as(smpl2, "sf"), zcol = "embc_clst_smth") 

#only keep levels 3 and 4
flight_mv <- mv %>%
  filter(embc_clst_smth %in% c(3,4))

##### STEP 2: step-selection prep - generate alternative steps #####

# 2.1: hourly subset and calculate turning angles and step lengths -----
# Create a list of move2 objects
mv_ls <- split(flight_mv, mt_track_id(flight_mv))

sp_obj_ls <- lapply(mv_ls, function(track){ #for each individual (i.e. track),
  
  # hourly subset and create a burst_id column
  #each sequence of rows with a time difference less than step duration will have the same burst_id.
  bursted_mv <- track %>%
    #hourly subset
    mutate(dt_hr = floor_date(timestamp, "1 hour")) %>% 
    group_by(individual.local.identifier,dt_hr) %>% 
    slice(1) %>% 
    ungroup() %>% 
    arrange(timestamp) %>% 
    #assign burst IDs to groups of rows with less than one hour difference
    mutate(time_lag = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "mins") %>%  as.numeric() %>% round(2))) %>% 
    mutate(burst_id = cumsum(time_lag >= 70) + 1) %>% #one hour + 10 min tolerance
    #remove bursts with less than 3 points
    group_by(burst_id) %>%
    filter(n() >= 3) %>% 
    ungroup()
  
    # Apply the function calculate_metrics to each element of the list and combine the results into a single data frame
    result_df <- bursted_mv %>%
      arrange(timestamp) %>%
      split(.$burst_id) %>% 
      map_df(calculate_metrics)
    
    result_df
    
}) 

# 2.2: estimate turning angles and step length distributions -----

#put everything in one data frame
bursted_df <- sp_obj_ls %>%  
  reduce(rbind) %>% 
  mutate(location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2]) %>% 
  as.data.frame()

#estimate von Mises parameters for turning angles
#calculate the averages (mu).steps: 1) make sure units are in radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
mu <- bursted_df %>% 
  drop_na(turning_angle_rad) %>% 
  pull(turning_angle_rad) %>% 
  mean.circular()

kappa <- bursted_df %>% 
  drop_na(turning_angle_rad) %>% 
  pull(turning_angle_rad) %>% 
  as.numeric() %>% 
  est.kappa()

#estimate gamma distribution for step lengths and convert to km
sl <- bursted_df %>% 
  filter(!is.na(step_length) & step_length >  set_units(0, "m")) %>%   #remove 0s and NAs
  pull(step_length) %>% 
  set_units("km") %>% 
    as.numeric()

fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")

#plot turning angle and step length distributions
par(mfrow=c(1,2))
hist(sl, freq=F, main="", xlab = "Step length (km)")
plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
hist(bursted_df %>% drop_na(turning_angle_rad) %>% pull(turning_angle_rad),freq = F,main = "", xlab = "Turning angles (radians)")
plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")

#diagnostic plots for step length distribution
plot(fit.gamma1)

# 2.2: generate alternative steps -----

n_alt <- 50
#prepare cluster for parallel computation
mycl <- makeCluster(5) #the number of CPUs to use (adjust this based on your machine)

clusterExport(mycl, c("sp_obj_ls", "mu", "kappa", "fit.gamma1", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(tidyverse)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
  library(units)
})

(start_time <- Sys.time()) 
used_av_track <- parLapply(mycl, sp_obj_ls, function(track){ #for each track
  
  lapply(split(track,track$burst_id),function(burst){ #for each burst,
    
    #assign unique step id
    burst$step_id <- 1:nrow(burst)
    
    lapply(c(2:(nrow(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
      
      current_point<- burst[this_point,]
      previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
      used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
      
      #calculate bearing of previous point
      prev_bearing <- NCEP.loxodrome.na(st_coordinates(previous_point)[,2], st_coordinates(current_point)[,2],
                                        st_coordinates(previous_point)[,1], st_coordinates(current_point)[,1])
      
      current_point_m <- st_transform(current_point, meters_proj) #convert to meters proj
      
      #randomly generate n alternative points
      rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                        step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
        #find the geographic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
        mutate(lon = st_coordinates(current_point_m)[,1] + step_length*cos(turning_angle),
               lat = st_coordinates(current_point_m)[,2] + step_length*sin(turning_angle))
      
      rnd_sf <- rnd %>% 
        st_as_sf(coords = c("lon", "lat"), crs = meters_proj) %>% 
        st_transform(wgs)

      #check visually
      # mapview(current_point, color = "red") + mapview(previous_point, color = "orange") + mapview(used_point, color = "yellow") + mapview(rnd_sf, color = "black", cex = 0.5)
      
      #put used and available points together
      df <- used_point %>%  
        slice(rep(row_number(), n_alt + 1)) %>% #paste each row n_alt times for the used and alternative steps
        mutate(location.long = c(head(st_coordinates(.)[,1],1), st_coordinates(rnd_sf)[,1]), 
               location.lat = c(head(st_coordinates(.)[,2],1), st_coordinates(rnd_sf)[,2]),
               turning_angle_deg = c(head(turning_angle_deg, 1), set_units(deg(rnd_sf$turning_angle), "degrees")),
               step_length = c(head(step_length,1), set_units(rnd_sf$step_length, "m")),
               used = c(1, rep(0, n_alt)))  %>%
        st_drop_geometry() %>% 
        dplyr::select(-c("turning_angle_rad","dt_hr", "time_lag", "ground_speed", "eobs.horizontal.accuracy.estimate")) %>% 
        rowwise() %>% 
        mutate(heading = NCEP.loxodrome.na(lat1 = st_coordinates(current_point)[,2], lat2 = location.lat, lon1 = st_coordinates(current_point)[,1], lon2 = location.long)) %>% 
        as.data.frame()
      
      df
      
    }) %>% 
      bind_rows()
    
  }) %>% 
    bind_rows()
  
}) %>% 
  bind_rows()

Sys.time() - start_time # Check how long this step took. On my machine this took 
stopCluster(mycl) 

#add stratum ID
used_av_track <- used_av_track %>% 
  mutate(stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))

##### STEP 3: step-selection prep - annotation with topographic info #####

#open topographic layers
TRI_100 <- rast("TRI_100_LF.tif")
ridge_100 <- rast("ridge_100_LF.tif")

#reproject tracking data to match topography, extract values, and convert back to wgs and save as a dataframe
topo_ann_df <- used_av_track %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_transform(crs = crs(TRI_100)) %>% 
  extract(x = TRI_100, y = ., method = "simple", bind = T) %>%
  extract(x = ridge_100, y = ., method = "simple", bind = T) %>% 
  terra::project(wgs) %>% 
  data.frame(., geom(.)) %>% 
  dplyr::select(-c("geom", "part", "hole")) %>% 
  rename(location.long = x,
         location.lat = y,
         ridge_100 = distance_to_ridge_line_mask,
         TRI_100 = TRI)

#The output of this step is stored as "SSF_input_data.csv" in the Edmond repository: https://doi.org/10.17617/3.FM4EJC
#For the next steps, including data analysis and plotting, see "02_data_analysis.r" in the Github repository