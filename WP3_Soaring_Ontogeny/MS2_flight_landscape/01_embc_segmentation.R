#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: data clean up and segmentation
#full data set downloaded from Movebank on Jan. 13. 2022
#Jan 11. 2022. Elham Nourani. Konstanz, DE

###use the new dates for the two individuals:
#Matsch19 (eobs 7035) 2020-03-08 12:17:00
#Güstizia18 (eobs 5942) 2019-03-24 11:18:31

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(EMbC)
library(mapview)
library(terra)
library(data.table); setDTthreads(percent = 65) #set this so getMovebankData works!

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")


# STEP 1: Download data ----------------------------------------------------------------

#open file with all data (March 15. 2023)
  data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp, format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))

#1 min subset to get rid of super bursts!
data_1min <- data %>% 
  mutate(dt_1min = round_date(timestamp, "1 minute")) %>% 
  group_by(individual.local.identifier,dt_1min) %>% 
  slice(1) 

#remove data and clean up 
rm(data); gc(gc())

saveRDS(data_1min, file = "/home/enourani/Desktop/Golden_Eagle_data/All_gps_mar23/LifeTrack_Golden_Eagle_Alps_1min.rds")

### only keep data with emigration info

#open file containing emigration dates
emig_dates <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/fleding_emigration_timing_Mar2023.rds") #this file only has the juvies.
emig_dates[emig_dates$individual.local.identifier == "Johnsbach21-1 (eobs 7582)", "individual.local.identifier"] <- "Johnsbach1_21 (eobs 7582)"          
emig_dates[emig_dates$individual.local.identifier == "Johnsbach21-2 (eobs 7585)", "individual.local.identifier"] <- "Johnsbach2_21 (eobs 7585)"

#only keep individuals that we have emigration info
data_w_info <- data_1min %>% 
  filter(individual.local.identifier %in% emig_dates$individual.local.identifier) #n = 69

# STEP 2: assign life stages: filter post-emigration ----------------------------------------------------------------

ind_ls <- split(data_w_info, data_w_info$individual.local.identifier)

data_stage <- lapply(ind_ls, function(x){
  
  d <- emig_dates %>% 
    filter(individual.local.identifier ==  unique(x$individual.local.identifier))
  
  x <- x %>% 
    mutate(stage = ifelse(timestamp >= d$emigration_dt, "post_emigration", "pre_emigration")) #don't use the fledging dates. it might not be accurate
  
  x
}) %>% 
  reduce(rbind)

save(data_stage, file = "data_w_lifestage_1min_69n.RData")

#reduce frequency to 20 minutes already. better to have a rather uniform res when segmenting using embc

#extract post-emigration data
post_em_20m <- data_stage %>% 
  mutate(dt_20min = round_date(timestamp, "20 minutes")) %>% 
  group_by(individual.local.identifier,dt_20min) %>% 
  slice(1) %>% 
  filter(stage == "post_emigration")

saveRDS(post_em_20m, file = "post_em_df_20min_68n.rds")

# STEP 3: estimate flight height ----------------------------------------------------------------

post_em <- readRDS("post_em_df_20min_68n.rds")

#open geoid and dem layer
geo <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/EGM96_us_nga_egm96_15.tif")
dem <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/Alps_east_west_dem.tif") %>% 
  project("+proj=longlat +datum=WGS84 +no_defs")


#extract elevation values
post_em$dem <- extract(x = dem, y = post_em[,c("location.long","location.lat")], method = "bilinear")[,2] 
post_em$geoid <- extract(x = geo, y = post_em[,c("location.long","location.lat")], method = "bilinear")[,2] 

#calculate flight height as height above mean sea level - dem. height above sea level is height above ellipsoid - geoid
post_em_h <- post_em %>% 
  mutate(height_msl = height.above.ellipsoid - geoid) %>% 
  mutate(height_ground = height_msl - dem)

saveRDS(post_em_h, file = "post_em_df_h_20min_68ind.rds") #n = 68

# STEP 4: estimate ground speed ----------------------------------------------------------------

post_em_h <- readRDS("post_em_df_h_20min_68ind.rds") 

post_em_h <- post_em_h %>% 
  arrange(individual.local.identifier,timestamp) %>% 
  as.data.frame()

#convert to a move object.
mv <- move(x = post_em_h$location.long, y = post_em_h$location.lat, time = post_em_h$timestamp, proj = wgs, data = post_em_h, animal = post_em_h$individual.local.identifier)
mv$ground_speed <- unlist(lapply(speed(mv),c, NA ))
mv$time_lag <- unlist(lapply(timeLag(mv, units = "mins"),  c, NA))

save(mv, file = "post_em_df_h_20min_68ind_mv.rds") 

# STEP 5: EMbC segmentation ----------------------------------------------------------------

# remove the outliers

mv <- readRDS("post_em_df_h_20min_68ind_mv.rds") #after removing the NAs below, we loose 7 individuals

input <- as.data.frame(mv) %>%
  drop_na(ground_speed, height_ground) %>% #remove NA values for speed
  filter(between(height_ground, quantile(height_ground, 0.005), quantile(height_ground, 0.999)) & ground_speed < quantile(ground_speed, 0.999)) %>% #remove outliers
  as.data.frame() #n = 61 :(

#create a matrix of flight height and speed
m <- data.matrix(input[,c("ground_speed","height_ground")])

#call embc
(b <- Sys.time())
bc <- embc(m)
Sys.time() -b #8 min

#investigate the bc (Garriga et al 2016; S2)
X11();sctr(bc)

bc_smth <- smth(bc,dlta = 0.7)
X11();sctr(bc_smth)

#append cluster labels (1:LL, 2:LH, 3:HL, and4:HH) to original data
input$embc_clst <- bc@A
input$embc_clst_smth <- bc_smth@A

saveRDS(input, file = "embc_output_20min_61ind.rds") #some individuals have very few rows of data left.

#select a sample track and visualize
smpl <- input %>% 
  filter(individual.local.identifier == "Viluoch17 (eobs 4570)")

coordinates(smpl) <- ~ location.long + location.lat
proj4string(smpl) <- wgs
ln <- SpatialLines(list(Lines(list(Line(smpl)), "line1")))
proj4string(ln) <- wgs

mapview(ln, color = "gray") + mapview(smpl, zcol = "embc_clst")

##################
#only look at levels 3 and 4

#select a sample track and visualize
sml <- input %>% 
  filter( embc_clst_smth %in% c(3,4) & individual.local.identifier == "Viluoch17 (eobs 4570)")

coordinates(sml) <- ~ location.long + location.lat
proj4string(sml) <- wgs
ln <- SpatialLines(list(Lines(list(Line(sml)), "line1")))
proj4string(ln) <- wgs

mapview(ln, color = "gray") + mapview(sml, zcol = "embc_clst")
