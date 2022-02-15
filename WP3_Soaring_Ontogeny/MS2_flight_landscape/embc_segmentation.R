#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: data clean up and segmentation
#full data set downloaded from Movebank on Jan. 13. 2022
#Jan 11. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(EMbC)
library(mapview)

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

save(dem_wgs, file = "/home/enourani/ownCloud/Work/GIS_files/EU_DEM/eu_dem_v11_E40N20/dem_wgs.RData")

#extract elevation values
post_em$dem_alt <- extract(x = dem_wgs, y = post_em[,c("location.long","location.lat")], method = "bilinear")

#calculate flight height as ellipsoid-dem
post_em$flight_h <- post_em$height.above.ellipsoid - post_em$dem_alt

save(post_em, file = "post_em_df_dem.RData")


# STEP 4: subset to 1 min intervals and estimate ground speed ----------------------------------------------------------------

load("post_em_df_dem.RData") #post_em

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

input <- as.data.frame(mv_mnt) %>%
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
input$embc_clst_smth <- bc_smth2@A

save(input, file = "embc_output.RData")

#select a sample track and visualize
smpl <- input %>% 
  filter(individual.local.identifier == "Trimmis20 (eobs 7041)")

coordinates(smpl) <- ~ location.long + location.lat
proj4string(smpl) <- wgs
ln <- SpatialLines(list(Lines(list(Line(smpl)), "line1")))
proj4string(ln) <- wgs

mapview(ln, color = "gray") + mapview(smpl, zcol = "embc_clst")

mapview(dd, color = "gray") + mapview(smpl, zcol = "embc_clst_smth")
