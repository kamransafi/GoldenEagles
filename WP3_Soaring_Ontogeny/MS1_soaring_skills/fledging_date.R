#R script to extract the fledging date and time from the golden eagle tracking data
#Elham Nourani, PhD. Konstanz, Germany.
#Apr 22. 2022


setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

library(tidyverse)
library(lubridate)
library(sf)
library(amt) #to estimate mcp using a sf object

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")


#STEP 1: open data and match with Svea's timestamps --------------------------------------------------------

# Hester's files on matching the names
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

#open file with all data
data <- read.csv("/home/enourani/Desktop/Golden_Eagle_data/all_GPS_jan13_22/LifeTrack Golden Eagle Alps.csv", encoding = "UTF-8") %>% 
  mutate(#id = strsplit(individual.local.identifier, " ") %>% map_chr(., 1),
         timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  full_join(dates[,c("fledging_timestamp", "emigration_timestamp", "local_identifier")], by = c( "individual.local.identifier" = "local_identifier"))

#n_distinct(data$individual.local.identifier) is 78

saveRDS(data, file = "GPS_data_dates_matched.rds")


#STEP 2: ID fledging date --------------------------------------------------------

#use a sample of individuals first

set.seed(777)
sample <- data %>% 
  filter(individual.local.identifier %in% sample(unique(individual.local.identifier),3, replace = F)) %>%  #randomly select 3 individuals
  group_by(individual.local.identifier) %>% 
  arrange(timestamp) %>% 
  #rowwise() %>% 
  mutate(data_days = difftime(timestamp, head(timestamp, 1), unit = "days")) %>% 
  filter(data_days <= 50) %>%  #only keep the first 50 days of data
  mutate(data_days_rnd = as.numeric(ceiling(data_days))) %>% #round up the day values
  drop_na("location.long") %>% 
  ungroup()

saveRDS(sample, file = "sample_for_fledging_code.rds")

#use parallel package.
lapply(split(sample, sample$individual.local.identifier), function(ind){ #for each individual,
  
  potential_nest <- ind %>% #extract first point recorded for the individual as the potential nest site
    arrange(timestamp) %>% 
    slice(1)
  
  ind <- ind %>% 
    filter(data_days_rnd >= 1) #get rid of day 0
  
  lapply(split(ind, ind$data_days_rnd), function(day){ #for each day
    
    #estimate MCP
    day_mcp <- day %>% 
      st_as_sf(coods = c("location.long", "location.lat"), crs = 4326) %>% #convert data into sf object
      hr_mcp(percentile = 95)

    #check for overlap between the MCP and the nest
    day_overlap <- day %>% 
      
    
    overlap <- potential_nest %>% 
      st_intersection(day_mcp) #extract a T or F
    
    return(data.frame(ind = ind,
                      day = day,
               overlap = overlap))
    
    
  })
  
} )