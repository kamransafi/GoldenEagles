### download and clean eagle data

library(move)
library(tidyverse)
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")

# the names of the birds that were tagged in the nest
birds <- getMovebankReferenceTable(study = 282734839, login = loginStored) %>%
  drop_na(animal_id) %>%
  filter(sensor_type_id == 653 & number_of_location_events > 0) %>%
  # remove the adults that were released for monitoring
  filter(!animal_id %in% c(1286236005, 1083984854, 1193204056, 1453194460, 
                              555579518, 1167138388, 1604141579, 1433272424, 
                              2313422481, 3088777090)) %>% 
  rename(individual.id = animal_id) %>% 
  dplyr::select(individual.id) %>% 
  deframe()

# get & clean the data
start_time <- Sys.time()
clean_data <- lapply(birds[2:length(birds)], function(id){
  print(paste0("Cleaning data for bird ", id, ": ", which(birds == id), 
               " of ", length(birds), "."), quote = F)
  # pull down data
  ind <- getMovebankLocationData(282734839, 653, id, loginStored)

  animalName <- unique(ind$individual.id)
  
  locs_df <- ind %>% 
    drop_na(location.long) %>% 
    mutate(index = row_number())
  # remove duplicated locations because they prevent accurate calculations of distance and speed
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),] %>% 
    # almost all the duplicates have missing values in height, heading, and eobs status
    filter(eobs.status == "")
  
  locs_df <- locs_df %>% 
    filter(!index %in% doubles$index)
  # if there are additional duplicates, the location will differ and DOP will be high for one
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),]
  if(nrow(doubles) > 0){
   doubles <- doubles %>% 
      group_by(timestamp) %>% 
      filter(gps.dop == max(gps.dop))
    locs_df <- locs_df %>% 
      filter(!index %in% doubles$index) %>% 
      dplyr::select(-index)
  }
  # filter out burst data
  locs_df <- locs_df %>% 
    mutate(timelag.sec = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>% 
    filter(timelag.sec > 2)
  # calculate track variables
  mv <- move(x=locs_df$location.long, y=locs_df$location.lat, 
             time=locs_df$timestamp,
             proj=crs("+proj=longlat +ellps=WGS84"),
             animal=locs_df$individual.local.identifier,
             data=locs_df)
  mv$timelag.sec <- c(NA,timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA,(mv$height.above.ellipsoid[-1] - mv$height.above.ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA,move::distance(mv))
  mv$gr.speed <- c(NA, speed(mv))
  # remove ridiculous speeds
  mv <- mv[which(mv$gr.speed <= 50),]

  # save on the hard drive
  saveRDS(mv, file = paste0("/home/hbronnvik/Documents/GE_data/no_burst/",
                         animalName,
                         "_gpsNoDup_moveObj.rds"))
})
Sys.time()-start_time # Time difference of 2.855802 hours with bursts 58.23823 mins without
