#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript
#full data set downloaded from Movebank on Jan. 13. 2022
#Jan 11. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(EMbC)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# STEP 1: open data and resolve mismatches in individual names ----------------------------------------------------------------

#open file with all data
data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  mutate(id = strsplit(individual.local.identifier, " ") %>% map_chr(., 1),
         timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))


#n_distinct(data$id) is 78


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


#n_distinct(dates$id) is 32
          
#only keep individuals that we have emigration info
emig_dates <- dates %>% 
  drop_na(emigration_timestamp)

#n_distinct(emig_dates$id) is 25

data_w_info <- data %>% 
  filter(data$individual.local.identifier %in% emig_dates$local_identifier) #n = 25

save(data_w_info, file = "all_data_w_emig.RData")


# STEP 2: assign life stages: filter post-emigration ----------------------------------------------------------------

ind_ls <- split(data_w_info, data_w_info$individual.local.identifier)

data_stage <- lapply(ind_ls, function(x){
  d <- emig_dates %>% 
    filter(local_identifier ==  unique(x$individual.local.identifier))
  
  x <- x %>% 
    mutate(stage = ifelse(timestamp < d$fledging_timestamp, "pre-fledging",
                          ifelse(
                            between(timestamp, d$fledging_timestamp, d$emigration_timestamp), "post_fledging", #this category includes the first day of fleding
                            ifelse(
                              timestamp >= d$emigration_timestamp,"post_emigration", NA))))
  
  x
}) %>% 
  reduce(rbind)


#extract post-emigration data
post_em <- data_stage %>% 
  filter(stage == "post_emigration")

save(post_em, file = "post_em_df.RData")


# STEP 3: estimate flight height ----------------------------------------------------------------

load( "post_em_df.RData")

#open EU-DEM
dem <- raster("/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF")

dem_wgs <- projectRaster(dem, crs = wgs) 

save(dem_wgs, file = "/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/dem_wgs")

#extract elevation values
post_em$dem_alt <- extract(x = dem_wgs, y = post_em[,c("location.long","location.lat")], method = "bilinear")

#calculate flight height as ellipsoid-dem
post_em$flight_h <- post_em$height.above.ellipsoid - post_em$dem_alt

save(post_em, file = "post_em_df_dem.RData")


# STEP 4: subset to 1 min intervals and estimate ground speed ----------------------------------------------------------------

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = post_em$individual.local.identifier, timestamps = post_em$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows


post_em <- post_em[-rows_to_delete,] 

post_em <- post_em[order(c(post_em$individual.local.identifier, post_em$timestamp)),]

#convert to a move object.
mv <- move(x = post_em$location.long, y = post_em$location.lat, time = post_em$timestamp, proj = wgs, data = post_em, animal = post_em$individual.local.identifier)
mv$speed <- unlist(lapply(speed(mv),c, NA ))

save(mv, file = "post_em_mv.RData") #move object with original data frequency and NA values

#remove NA values of flight_h
mv_no_na <- mv[!is.na(mv$flight_h),]

#keep the first instance of every minute.
df_mnt <- as.data.frame(mv_no_na) %>%
  group_by(individual.local.identifier, year(timestamp), month(timestamp), day(timestamp), hour(timestamp), minute(timestamp)) %>%
  slice(1) %>% 
  dplyr::select(!c("coords.x2" ,"coords.x1", "year(timestamp)","month(timestamp)", "day(timestamp)","hour(timestamp)", "minute(timestamp)")) %>% #remove columns that I don't need
  as.data.frame()

#convert back to a move object
mv_mnt <- move(x = df_mnt$location.long, y = df_mnt$location.lat, time = df_mnt$timestamp, proj = wgs, data = df_mnt, animal = df_mnt$individual.local.identifier)

#calculate ground speed and time interval between points
mv_mnt$speed <- unlist(lapply(speed(mv_mnt),c, NA ))
mv_mnt$time_lag <- unlist(lapply(timeLag(mv_mnt, units = "mins"),  c, NA))

save(mv_mnt, file = "post_em_mv_minutely.RData") #move object with minutely data frequency and no NA values for flight height

# STEP 5: EMbC segmentation ----------------------------------------------------------------

load("post_em_mv_minutely.RData") #mv_mnt

#prep for embc (for all individuals pooled. no reliability function will be used for speed, because data collection frequency is pretty uniform)

#deal with outliers and NAs
embc_input <- as.data.frame(mv_mnt) %>% 
  drop_na(speed) %>% #remove NA values for speed
  mutate(flight_h = ifelse(flight_h > quantile(mv_mnt$flight_h, 0.9), quantile(mv_mnt$flight_h, 0.9), #replace flight height values higher than the 90% quantile with the 90% quantile value
                              ifelse(flight_h < quantile(mv_mnt$flight_h, 0.1), quantile(mv_mnt$flight_h, 0.1), #replace flight height values lower than 10% quantile with the 10% quantile value
                                     flight_h)),
         speed = ifelse(speed > quantile(mv$speed, 0.9, na.rm = T), quantile(mv$speed, 0.9, na.rm = T), #replace flight height values than the 90% quantile with the 90% quantile value
                                  speed))

#create a matrix of flight height and speed
m_spd <- data.matrix(embc_input[,c("speed","flight_h")])

#call embc
bc_spd <- embc(m_spd)

#investigate the bc
sctr(bc_spd)


