#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: data clean up and segmentation
#full data set downloaded from Movebank on Jan. 13. 2022
#Jan 11. 2022. Elham Nourani. Konstanz, DE
#Update June 15. 20222: use 70 individuals (emigration timing was estimated using the recurse package)

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(EMbC)
library(mapview)
library(terra)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")


# STEP 1: open data ----------------------------------------------------------------

#open file with all data (until jan 13. 2022)
data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))

#1 min subset (get rid of super bursts)
data_1min <- data %>% 
  mutate(dt_1min = round_date(timestamp, "1 minute")) %>% 
  group_by(individual.local.identifier,dt_1min) %>% 
  slice(1) %>%
  as.data.frame()

saveRDS(data_1min, file = "/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps_1min.rds")

#download the new eagle data (start from jan 13 2022)
load("mbnklogin.rdata") #Elhamlogin
data_new <- getMovebankData(study = "LifeTrack Golden Eagle Alps", 
                              removeDuplicatedTimestamps = T, timestamp_start = "202201130000000",
                              login = Elhamlogin)

time_lag <- unlist(lapply(timeLag(data_new, units = "mins"),  c, NA))

#1 min subset (get rid of super bursts)
data_new_1min <- as.data.frame(data_new) %>% 
  mutate(dt_1min = round_date(timestamp, "1 minute")) %>% 
  group_by(local_identifier,dt_1min) %>% 
  slice(1) %>%
  as.data.frame()

saveRDS(data_new_1min, file = "/home/enourani/Desktop/Golden_Eagle_data/GPS_Jan13_Jun16_2022/LifeTrack Golden Eagle Alps_1min.rds")

###### bind data together and remove duplicates

cols <- c("timestamp","location.long", "location.lat", "individual.local.identifier",
"tag.local.identifier","heading" , "height.above.ellipsoid")

data_new_1min <- readRDS("/home/enourani/Desktop/Golden_Eagle_data/GPS_Jan13_Jun16_2022/LifeTrack Golden Eagle Alps_1min.rds") %>%
rename(location.lat = location_lat,
location.long = location_long,
individual.local.identifier = local_identifier,
height.above.ellipsoid= height_above_ellipsoid,
tag.local.identifier = tag_local_identifier) %>%
dplyr::select(cols)

data_1min <- readRDS("/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps_1min.rds") %>%
dplyr::select(cols)

all_data <- bind_rows(data_1min2, data_new_1min2)

#remove duplicated rows (mostly from data overlap on Jan 13 2022)
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = all_data$individual.local.identifier, timestamps = all_data$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate


all_data <- all_data[-rows_to_delete,] 

saveRDS(all_data, file = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/GPSdata_78ind_UntilJune162022_1min.rds")

### only keep data with emigration info
all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/GPSdata_78ind_UntilJune162022_1min.rds")

#open file containing emigration dates
load("/home/enourani/Desktop/Hester_GE/recurse/em_fl_dt_recurse_70ind.RData") #emig_fledg_dates
          
#only keep individuals that we have emigration info
data_w_info <- all_data %>% 
  filter(individual.local.identifier %in% emig_fledg_dates$individual.local.identifier) #n = 70

# STEP 2: assign life stages: filter post-emigration ----------------------------------------------------------------

ind_ls <- split(data_w_info, data_w_info$individual.local.identifier)

data_stage <- lapply(ind_ls, function(x){
  d <- emig_fledg_dates %>% 
    filter(individual.local.identifier ==  unique(x$individual.local.identifier))
  
  x <- x %>% 
    mutate(stage = ifelse(timestamp < d$fledging_dt, "pre-fledging",
                          ifelse(
                            between(timestamp, d$fledging_dt, d$emigration_dt), "post_fledging", #this category includes the first day of fledging
                            ifelse(
                              timestamp >= d$emigration_dt,"post_emigration", NA))))
  
  x
}) %>% 
  reduce(rbind)

save(data_stage, file = "data_w_lifestage_1min_70n.RData") #n = 70

#extract post-emigration data
post_em <- data_stage %>% 
  filter(stage == "post_emigration")

save(post_em, file = "post_em_df_1min_70n.RData")

#some summary stats
data_stage %>% 
  arrange(individual.local.identifier, stage, timestamp) %>% 
  group_by(individual.local.identifier, stage) %>% 
  summarise(start = head(timestamp,1), end = tail(timestamp,1)) %>%
  mutate(duration = as.numeric(end - start))



# STEP 3: estimate flight height ----------------------------------------------------------------

load( "post_em_df_15min_70.RData") #post_em

#prepare DEM (merge east and west Alps)

dem <- rast("/home/enourani/Desktop/golden_eagle_static_layers/French_Alps/dem_100.tif")
dem2 <- rast("/home/enourani/Desktop/golden_eagle_static_layers/dem_100.tif")

Alps <- st_read("/home/enourani/ownCloud/Work/GIS_files/Alpine_perimeter/Alpine_Convention_Perimeter_2018_v2.shp") %>% 
  st_transform(crs(dem)) %>% 
  as("SpatVector")
  
dem_east <- mask(dem2,Alps)
dem_west <- mask(dem, Alps)

Alps_topo <- merge(dem_east,dem_west)

writeRaster(Alps_topo, "Alps_east_west_dem.tif")
  
#reproject
alps_topo_wgs <- Alps_topo %>% 
  project("+proj=longlat +datum=WGS84 +no_defs")

#open EU-DEM
#dem <- raster("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF")

#dem_wgs <- projectRaster(dem, crs = wgs) 

#save(dem_wgs, file = "/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/dem_wgs.RData")

#extract elevation values
#load("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/dem_wgs.RData")
post_em$dem_alt <- extract(x = alps_topo_wgs, y = post_em[,c("location.long","location.lat")], method = "bilinear")[,2] 


#calculate flight height as ellipsoid-dem
post_em$flight_h <- post_em$height.above.ellipsoid - post_em$dem_alt

save(post_em, file = "post_em_df_dem_70ind.RData")

# STEP 4: subset to 1 min intervals and estimate ground speed ----------------------------------------------------------------

load("post_em_df_dem_70ind.RData") #post_em

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = post_em$individual.local.identifier, timestamps = post_em$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows


post_em <- post_em[-rows_to_delete,] 

post_em <- post_em %>% 
  arrange(individual.local.identifier,timestamp) %>% 
  as.data.frame()

#convert to a move object.
mv <- move(x = post_em$location.long, y = post_em$location.lat, time = post_em$timestamp, proj = wgs, data = post_em, animal = post_em$individual.local.identifier)
mv$speed <- unlist(lapply(speed(mv),c, NA ))

save(mv, file = "post_em_mv_70_15min.RData") #move object with 15 min frequency and NA values

#remove NA values of flight_h
mv_no_na <- mv[!is.na(mv$flight_h),]

#calculate ground speed and time interval between points
mv_no_na$speed <- unlist(lapply(speed(mv_no_na),c, NA ))
mv_no_na$time_lag <- unlist(lapply(timeLag(mv_no_na, units = "mins"),  c, NA))

# STEP 5: EMbC segmentation ----------------------------------------------------------------

#prep for embc (for all individuals pooled. no reliability function will be used for speed, because data collection frequency is pretty uniform)

#deal with outliers and NAs
# embc_input <- as.data.frame(mv_mnt) %>% 
#   drop_na(speed) %>% #remove NA values for speed
#   mutate(flight_h = ifelse(flight_h > quantile(flight_h, 0.999), quantile(flight_h, 0.999), #replace flight height values higher than the 90% quantile with the 90% quantile value
#                               ifelse(flight_h < quantile(flight_h, 0.001), quantile(flight_h, 0.001), #replace flight height values lower than 10% quantile with the 10% quantile value
#                                      flight_h)),
#          speed = ifelse(speed > quantile(speed, 0.999), quantile(speed, 0.999), #replace flight height values than the 90% quantile with the 90% quantile value
#                                   speed))
# 
# embc_input <- as.data.frame(mv_mnt) %>%  #the medians dont change, so is there any point in reassigning outlier values??
#   drop_na(speed) %>% #remove NA values for speed
#   mutate(flight_h = ifelse(flight_h > 3000, 3000, #replace flight height values higher than 3000 with 3000
#                            ifelse(flight_h < -100, -100, #replace flight height values lower than -100 with -100
#                                   flight_h)),
#          speed = ifelse(speed > 50, 50, #replace flight height values higher than 10000 with 10000
#                         speed))
# 
# #create a matrix of flight height and speed
# m_spd <- data.matrix(embc_input[,c("speed","flight_h")])
# 
# #call embc
# bc_spd <- embc(m_spd)
# 
# #investigate the bc (Garriga et al 2016; S2)
# sctr(bc_spd)
# 
# #smooth the labeling (deals with single assignments that are surrounded by other labels)
# bc_smth <- smth(bc_spd,dlta = 0.5)
# 
# X11()
# par(mfrow = c(1,2))
# sctr(bc_spd)
# sctr(bc_smth)
# 
# 
# #investigate the bc (Garriga et al 2016; S2)
# sctr(bc_spd)
# 
# lkhp(bc_spd)
# 
# stts(bc_spd)
# 
# #distribution of variables in each category
# hist(bc_spd@X[which(bc_spd@A%in%c(1,3)),1],breaks=seq(0,max(bc_spd@X[,1]),
#                                                     max(bc_spd@X[,1])/50),include.lowest=TRUE,xlim=range(bc_spd@X[,1]),
#      xlab='velocity (m/s)',main='velocity distribution for LOW turns',cex.main=0.8)
# abline(v=bc_spd@R[1,3],col=1,lwd=1.5,lty='dashed')
# hist(bc_spd@X[which(bc_spd@A%in%c(2,4)),1],breaks=seq(0,max(bc_spd@X[,1]),
#                                                     max(bc_spd@X[,1])/50),include.lowest=TRUE,xlim=range(bc_spd@X[,1]),
#      xlab='velocity (m/s)',main='velocity distribution for HIGH turns',cex.main=0.8)
# abline(v=bc_spd@R[4,2],col=1,lwd=1.5,lty='dotdash')
# 
# table(bc_spd@A)




###what if i remove the outliers?

input <- as.data.frame(mv_no_na) %>%
  drop_na(speed) %>% #remove NA values for speed
  filter(between(flight_h, quantile(flight_h, 0.005), quantile(flight_h, 0.999)) & speed < quantile(speed, 0.999)) %>% 
  as.data.frame()

#create a matrix of flight height and speed
m <- data.matrix(input[,c("speed","flight_h")])

#call embc
bc <- embc(m)

#investigate the bc (Garriga et al 2016; S2)
X11();sctr(bc)

bc_smth <- smth(bc,dlta = 0.6)
X11();sctr(bc_smth)

#append cluster labels (1:LL, 2:LH, 3:HL, and4:HH) to original data
input$embc_clst <- bc@A
input$embc_clst_smth <- bc_smth@A

saveRDS(input, file = "embc_output_70ind.rds")

#select a sample track and visualize
smpl <- input %>% 
  filter(individual.local.identifier == "Trimmis20 (eobs 7041)")

coordinates(smpl) <- ~ location.long + location.lat
proj4string(smpl) <- wgs
ln <- SpatialLines(list(Lines(list(Line(smpl)), "line1")))
proj4string(ln) <- wgs

mapview(ln, color = "gray") + mapview(smpl, zcol = "embc_clst")

