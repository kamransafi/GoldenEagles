#Script for preparing golden eagle tracking data for step-selection analysis as reported in Nourani et al. 2023
# Elham Nourani, PhD. 25.07.2023
# enourani@ab.mpg.de

library(tidyverse)
library(lubridate)
library(EMbC)

##### STEP 0: open tracking data #####

#filter for post-dispersal period

##### STEP 1: EmbC segmentation #####

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

##### STEP 2: step-selection prep - generate alternative steps #####

##### STEP 3: step-selection prep - annotation with topographic info #####