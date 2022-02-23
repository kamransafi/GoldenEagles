## Task 1: ratio of soaring modes
## 21.02.22


## 1. Identify all thermal soaring events
## 2. Identify all slope soaring events
## 3. Measure how much time is spent in each
## 4. Sum these over each day
## 5. Calculate the ratio
## 6. Plot the ratio over days since fledging


setwd("C:/Users/Tess Bronnvik/Desktop/Improvement_and_Golden_Eagles")

library(tidyverse)
library(move)
library(lubridate)

### Assign life stages to the birds for which we have gps and acc data associated
###########################################################################################
# the file with the names of the birds in different formats
load("eagle_names.RData")

# the birds that have ACC data associated
objs <- list.files("C:/Users/Tess Bronnvik/Desktop/Improvement_and_Golden_Eagles/associated")#, full.names = T)
str_sub(objs, -18,-1) <- ""
objs <- sub(" ", "\\.", objs)
objs <- sub("Art.San ", "", objs)

# the file with the fledging and emigration dates
stages <- read.csv("Goldeneagles10_2021.csv", stringsAsFactors = F)
stages <- stages[order(stages$id),]
stages <- stages[which(stages$id %in% objs),]
# get the fledging date column re-arranged
stages$date_fledging <- paste0(stages$date_fledging, " ", stages$time_fledging)
stages$date_fledging <- as.POSIXct(gsub("\\.", "\\-", stages$date_fledging), format = "%d-%m-%Y %H:%M:%OS")
# get the emigration date column re-arranged
stages$date_emigration <- paste0(stages$date_emigration, " ", stages$time_emigration)
stages$date_emigration <- as.POSIXct(gsub("\\.", "\\-", stages$date_emigration), format = "%d-%m-%Y %H:%M:%OS")
# set all the missing emigration dates to today so that it is impossible for a bird to have timestamps > emigration date
# to assign life stage, we look at timestamps < emigration date, and therefore cannot use NA
stages$date_emigration[which(is.na(stages$date_emigration))] <- Sys.time()

# open the segmented and associated files
# give each a life_stage column
# either add data or leave NA
gpsAcc_fls <- list.files("gps_acc_age")
gpsAcc_fls <- paste0("gps_acc_age/", gpsAcc_fls)

gpsAcc_ls <- list()

for (i in gpsAcc_fls) {
  load(i)
  if(unique(gaa_df$individual_local_identifier) %in% eagle_names$local_identifier){
    gpsAcc_ls[[length(gpsAcc_ls) + 1]] <- gaa_df
    print(paste0("Successfully loaded ", sub("_.RData", "", sub("\\/", "", gsub("gps_acc_age", "", i))), "."), quote = F)
  } else{rm(i);rm(gaa_df)}
}


# loop through each individual and label its data with life stage
for (i in 1:length(gpsAcc_ls)) {
  # get the associated and segmented ACC & GPS data from a single id
  gaa_df <- gpsAcc_ls[[i]]
  # add the name for life stage data to the data frame with gps and acc
  gaa_df$sname <- NA
  gaa_df$sname <- eagle_names$age_name[which(eagle_names$local_identifier == unique(gaa_df$individual_local_identifier))]
  # for each id, match stages to move, 
  gaa_df$stage <- NA
  # then say before fledge, after fledge, and after emigration
  if (unique(gaa_df$individual_local_identifier) %in% eagle_names$local_identifier){ #our_Inds$eobs
    # timestamps less than date of fledging are classified as 1 (pre-fledging)
    gaa_df$stage[which(gaa_df$timestamp < stages$date_fledging[which(stages$id == unique(gaa_df$sname))])] <- 1
    # timestamps greater than date of fledging & less than date of emigration are classfied as 2 (fledgling)
    gaa_df$stage[which(gaa_df$timestamp > stages$date_fledging[which(stages$id == unique(gaa_df$sname))] & 
                         gaa_df$timestamp < stages$date_emigration[which(stages$id == unique(gaa_df$sname))])] <- 2
    # timesamps greater than date of emigration are classfied as 3 (emigrant)
    gaa_df$stage[which(gaa_df$timestamp > stages$date_emigration[which(stages$id == unique(gaa_df$sname))])] <- 3
  }
  # save the file again
  save(gaa_df, file = paste0("gps_acc_age/", unique(gaa_df$individual_local_identifier), "_gps_acc_age.RData"))
  # signal
  print(paste0("Added life stage information to ", unique(gaa_df$individual_local_identifier), "."), quote = F)
}

# finally, check that there are life stage data associated with each individual
for (i in 1:length(gpsAcc_ls)) {
  df <- gpsAcc_ls[[i]]
  if(is.na(df$stage)){
    print(paste0(unique(gpsAcc_ls[[i]]$individual_local_identifier), " has no stages."))
  }
}

###########################################################################################


### Assign IDs to events, count the time spent, and plot it
###########################################################################################
cfls_ls <- list.files("gps_acc_age")
cfls_ls <- paste0("gps_acc_age/", cfls_ls)

cfls <- data.frame()

for (i in cfls_ls) {
  load(i)
  gaa_df <- gaa_df[which(gaa_df$stage == 2),]
  cfls <- rbind(cfls,gaa_df)
  rm(gaa_df)
  print(paste0("Successfully loaded ", sub("_.RData", "", sub("\\/", "", gsub("gps_acc_age", "", i))), "."), quote = F)
}
ori_cfls <- cfls
#save(ori_cfls, file = "cfls.RData")
load("cfls.RData"); cfls <- ori_cfls

cfls <- cfls %>% 
  # only post-fledging and pre-emigration soaring events
  filter(thermalClust != "other") %>% 
  mutate(year = year(timestamp),
         month = month(timestamp),
         day = day(timestamp), 
         ymd = paste(year, month, day, sep = "-"))

## 2. give IDs to each unique linear event

# where the thermal cluster is linear soaring, True
cfls$linear[which(cfls$thermalClust == "linear")] <- T
# where it is circular or "other", False
cfls$linear[which(cfls$thermalClust != "linear")] <- F
# assign a new number to each slope soaring event 
cfls$linear_event <- inverse.rle(within.list(rle(cfls$linear), 
                                 values[values] <- seq_along(values[values])))
cfls$linear_event[which(cfls$linear_event == 0)] <- NA
# give each a linear ID
cfls$linearID <- paste0(cfls$burstIDcorrect, "_", cfls$linear_event)
cfls$linearID[grep("NA", cfls$linearID)] <- NA



## 3. measure the time spent in each "event"

cfls <- cfls %>%
  # work within slope soaring events
  group_by(individual_local_identifier, year, month, day, burstIDcorrect, linearID) %>%
  arrange(timestamp) %>% 
  # calculate the time difference between the last and first observations
  mutate(line_dt = difftime(.$timestamp[n()],.$timestamp[1])) %>% 
  ungroup() %>% 
  # work within thermal soaring events
  group_by(individual_local_identifier, year, month, day, burstIDcorrect, thermalID) %>%
  arrange(timestamp) %>% 
  # calculate the time difference between the last and first observations
  mutate(circle_dt = difftime(.$timestamp[n()],.$timestamp[1])) %>% 
  ungroup()
# these columns now hold the time spent in an event and the time spent out of it
cfls$line_dt[which(cfls$thermalClust != "linear")] <- NA
cfls$circle_dt[which(cfls$thermalClust != "circular")] <- NA
# these columns now hold only the time spent in an event

# 4. sum the time spent per day
cfls_ls <- split(cfls, cfls$individual_local_identifier)

lapply(cfls_ls, sum())

Lcfls <- cfls %>% 
  # within each slope soaring event 
  group_by(linearID) %>% 
  # select the first observation
  slice(1) %>% 
  # in all the data
  ungroup() %>% 
  # per day
  group_by(individual_local_identifier, ymd) %>% 
  # sum the seconds in linear events
  mutate(day_lines = sum(as.numeric(line_dt), na.rm = T),
         day_circles = NA) %>% 
  # one value per day
  slice(1) %>% 
  ungroup()
Ccfls <- cfls %>% 
  # within each thermal soaring event 
  group_by(thermalID) %>% 
  # select the first observation
  slice(1) %>% 
  # in all the data
  ungroup() %>% 
  # per day
  group_by(individual_local_identifier, ymd) %>% 
  # sum the seconds in linear events
  mutate(day_lines = NA,
         day_circles = sum(as.numeric(circle_dt), na.rm = T)) %>% 
  # one value per day
  slice(1) %>% 
  ungroup()


## 5. Calculate the ratio of linear to circular soaring
# only one observation of each soaring event, both linear and circular
#Ccfls$day_lines <- Lcfls$day_lines[which(Lcfls$ymd == Ccfls$ymd)]
Rcfls <- Ccfls %>% 
  rowwise() %>% 
  mutate(ratiosoar = day_lines/day_circles)


plot(Rcfls$day, Rcfls$ratiosoar)

