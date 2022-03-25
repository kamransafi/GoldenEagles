## Task 1: ratio of soaring modes
## Hester Brønnvik
## hbronnvik@ab.mpg.de
## 21.02.22


## 1. Calculate the total time spent flying in these data (gliding, soaring, flapping)
## 2. Identify all soaring events
## 3. Measure how much time is spent in each
## 4. Sum these over each day
## 5. Calculate the ratio of soaring to flight
## 6. Plot the ratio over days since fledging


setwd("C:/Users/Tess Bronnvik/Desktop/Improvement_and_Golden_Eagles")


library(move)
library(lubridate)
library(tidyverse)
library(ggpubr)

### Assign life stages to the birds for which we have gps data
###########################################################################################
# the file with the names of the birds in different formats
load("eagle_names.RData")

# the birds that have ACC data associated
objs <- list.files("C:/Users/Tess Bronnvik/Desktop/Improvement_and_Golden_Eagles/thermals")#, full.names = T)
objs <- sub("\\) ", "\\)", objs) # correct Stürfis20
objs <- str_sub(objs, 1,-48)
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

# open the segmented files
# give each a life_stage column
# either add data or leave NA
gps_fls <- list.files("thermals")
gps_fls <- paste0("thermals/", gps_fls)

gps_ls <- list()

for (i in gps_fls) {
  load(i)
  if(unique(HRdf$local_identifier) %in% eagle_names$local_identifier){
    gps_ls[[length(gps_ls) + 1]] <- HRdf
    print(paste0("Successfully loaded ", unique(HRdf$local_identifier), "."), quote = F)
  } else{rm(i);rm(HRdf)}
}


# loop through each individual and label its data with life stage
for (i in 1:length(gps_ls)) {
  # get the segmented GPS data from a single id
  gps_df <- gps_ls[[i]]
  # add the name for life stage data to the data frame with gps
  gps_df$sname <- NA
  gps_df$sname <- eagle_names$age_name[which(eagle_names$local_identifier == unique(gps_df$local_identifier))]
  # for each id, match stages to move, 
  gps_df$stage <- NA
  # then say before fledge, after fledge, and after emigration
  if (unique(gps_df$local_identifier) %in% eagle_names$local_identifier){ #our_Inds$eobs
    # timestamps less than date of fledging are classified as 1 (pre-fledging)
    gps_df$stage[which(gps_df$timestamp < stages$date_fledging[which(stages$id == unique(gps_df$sname))])] <- 1
    # timestamps greater than date of fledging & less than date of emigration are classfied as 2 (fledgling)
    gps_df$stage[which(gps_df$timestamp > stages$date_fledging[which(stages$id == unique(gps_df$sname))] & 
                         gps_df$timestamp < stages$date_emigration[which(stages$id == unique(gps_df$sname))])] <- 2
    # timesamps greater than date of emigration are classfied as 3 (emigrant)
    gps_df$stage[which(gps_df$timestamp > stages$date_emigration[which(stages$id == unique(gps_df$sname))])] <- 3
  }
  # save the file again
  save(gps_df, file = paste0("gps_age/", unique(gps_df$local_identifier), gsub("-",".", Sys.Date()), "_gps_age.RData"))
  # signal
  print(paste0("Added life stage information to ", unique(gps_df$local_identifier), "."), quote = F)
}

# finally, check that there are life stage data associated with each individual
for (i in 1:length(gps_ls)) {
  df <- gps_ls[[i]]
  if(NA %in% unique(df$stage)){
    print(paste0(unique(gps_ls[[i]]$local_identifier), " has no stages."))
  }
}

###########################################################################################

### Retrieve GPS data for the birds with life stage data
###########################################################################################
load("loginStored.RData")
eagleStudyId <- 282734839


# get individuals

inds <- vector()
for (i in 1:length(gps_fls)) {
  id <- gps_fls[[i]]@idData$local_identifier
  inds <- c(inds,id)
}

inds_todo <- eagle_names$local_identifier[which(!eagle_names$local_identifier %in% inds)]

gps_fls <- list()

for (i in inds_todo) {
  gps_fls[[length(gps_fls) + 1]] <- getMovebankData(study="LifeTrack Golden Eagle Alps", animalName = i,
                                                    removeDuplicatedTimestamps=T, login=loginStored)
  print(paste0("Retrieved location data for bird ", which(inds_todo == i), ", ", i, "."), quote = F)
}
# memory error for 
# "Tuors1 19 (eobs 7010)" "Tuors2 19 (eobs 7011)" "Viluoch17 (eobs 4570)"
#save(gps_fls, file = "gps_data_29_23.02.22.RData")


preps <- list.files("prepped/")
# preps <- str_sub(preps, 1, -35)
preps <- paste0("prepped/", preps)
load(preps[1])

###########################################################################################


### Assign IDs to events, count the time spent, and plot it
###########################################################################################
cfls_ls <- list.files("gps_age")
cfls_ls <- paste0("gps_age/", cfls_ls)

cfls <- data.frame()

for (i in cfls_ls) {
  load(i)
  gps_df <- gps_df[which(gps_df$stage == 2),]
  cfls <- rbind(cfls,gps_df)
  rm(gps_df)
  print(paste0("Successfully loaded ", str_sub(sub("_.RData", "", sub("\\/", "", gsub("gps_age", "", i))), 1, -11), "."), quote = F)
}
ori_cfls <- cfls
#save(ori_cfls, file = "cfls.RData")
#load("cfls.RData"); cfls <- ori_cfls

# reclassify the thermal soaring events so that there is no circular gliding
cfls$thermalClust <- "other"
cfls$thermalClust[which(cfls$gr.speed >= 2 & cfls$soarClust == "soar" & cfls$turnAngle_smooth >= 300)] <- "circular"
cfls$thermalClust[which(cfls$gr.speed >= 2 & cfls$soarClust=="soar" & cfls$thermalClust != "circular")] <- "linear"
cfls$thermalClust <- factor(cfls$thermalClust, levels=c("circular","linear","other"))

# remove unnecessary and unwieldy data
{
cfls$nick_name <- NULL
cfls$earliest_date_born <- NULL
cfls$latest_date_born <- NULL
cfls$comments <- NULL
cfls$death_comments <- NULL
cfls$barometric_pressure <- NULL
cfls$eobs_status <- NULL
cfls$eobs_temperature <- NULL
cfls$eobs_type_of_fix <- NULL
cfls$eobs_used_time_to_get_fix <-NULL
cfls$sensor <- NULL
cfls$sensor_type_id <- NULL
cfls$sensor_type_id <- NULL
cfls$sensor_type_ids <- NULL
cfls$timestamp_start <- NULL
cfls$timestamp_end <- NULL
}


# append the fledging dates for each individual
cfls$fledging_date <- NA
tcfls <- cfls
cfls <- data.frame()

inds <- unique(tcfls$local_identifier)
for (i in unique(tcfls$local_identifier)) {
  df <- tcfls[which(tcfls$local_identifier == i),]
  df$fledging_date <- stages$date_fledging[which(stages$id == 
                                                   eagle_names$age_name[which(eagle_names$local_identifier == i)])]
  cfls <- rbind(cfls, df)
  
}

## 1. calculate the time spent in flight per day
cfls <- cfls[which(cfls$gr.speed >= 2),] # only flight data
cfls$ymd <- format(cfls$timestamp, format = "%Y-%m-%d")
cfls$id_date <- paste0(cfls$local_identifier," ", cfls$ymd)

id_date <- unique(cfls$id_date)

cfls_t <- data.frame()

for (i in id_date) {
  df <- cfls[which(cfls$id_date == i),]
  df <- df[order(df$timestamp),]
  df$flight_time <- as.numeric(difftime(df$timestamp[nrow(df)], df$timestamp[1], units = "sec"))
  cfls_t <- rbind(cfls_t, df)
  rm(df)
  print(paste0("Calculated flight time per day for ", i, "."), quote = F)
}

cfls <- cfls_t; rm(cfls_t)

## 2. give IDs to each unique soaring event

# a data frame to hold the IDs
cfls2 <- data.frame()
# a vector to separate the individuals and days so that there is no chance that the data
# for one ends on a thermal and the next starts on a thermal, which are then counted as a single
# event



for (i in id_date) {
  df <- cfls[which(cfls$id_date == i), ]
  df <- df[order(df$timestamp),]
  df$linear <- NA
  # where the thermal cluster is linear soaring, True
  df$linear[which(df$thermalClust == "linear")] <- T
  # where it is circular or "other", False
  df$linear[which(df$thermalClust != "linear")] <- F
  # assign a new number to each slope soaring event 
  df$linearID <- inverse.rle(within.list(rle(df$linear), 
                                           values[values] <- seq_along(values[values])))
  df$linearID[which(df$linearID == 0)] <- NA
  cfls2 <- rbind(cfls2, df)
  rm(df)
  print(paste0("Assigned orographic event IDs for ", i, "."), quote = F)
}
cfls <- cfls2; rm(cfls2)

# a data frame to hold the IDs
cfls3 <- data.frame()


# and to each thermal event
for (i in unique(cfls$id_date)) {
  df <- cfls[which(cfls$id_date == i), ]
  df <- df[order(df$timestamp),]
  df$thermal <- NA
  # where the thermal cluster is thermal soaring, True
  df$thermal[which(df$thermalClust == "circular")] <- T
  # where it is linear or "other", False
  df$thermal[which(df$thermalClust != "circular")] <- F
  # assign a new number to each slope soaring event 
  df$thermalID <- inverse.rle(within.list(rle(df$thermal), 
                                         values[values] <- seq_along(values[values])))
  df$thermalID[which(df$thermalID == 0)] <- NA
  cfls3 <- rbind(cfls3, df)
  rm(df)
  print(paste0("Assigned thermal event IDs for ", i, "."), quote = F)
}
cfls <- cfls3; rm(cfls3)

## 3. measure the time spent in each soaring "event"

# pare down the data to relieve computing stress
# using only the data for linear or thermal soaring
cfls <- cfls[which(cfls$thermalClust != "other"),]

# create a column to store the linear event difference in time 
cfls$line_dt <- NA

templ <- data.frame()
tempt <- data.frame()

for (i in id_date) {
  df <- cfls[which(cfls$id_date == i),]
  line <- unique(df$linearID)
  for (j in line) {
    df2 <- df[which(df$linearID == j),]
    df2$line_dt <- as.numeric(difftime(df2$timestamp[nrow(df2)], df2$timestamp[1], units = "sec"))
    templ <- rbind(templ, df2)
    rm(df2)
  }
  print(paste0("Calculated time spent in each linear soaring event for ", i, "."), quote = F)
  circle <- unique(df$thermalID)
  for (k in circle){
    df2 <- df[which(df$thermalID == k),]
    df2$circle_dt <- as.numeric(difftime(df2$timestamp[nrow(df2)], df2$timestamp[1], units = "sec"))
    tempt <- rbind(tempt, df2)
    rm(df2)
  }
  print(paste0("Calculated time spent in each thermal soaring event for ", i, "."), quote = F)
  rm(df)
}

#save(templ, file = "templ_1Mar.RData")
#save(tempt, file = "tempt_1Mar.RData")

# 4. sum the time spent per day

thermal_times <- data.frame()

# the thermal events in a day for an individual
for (i in unique(tempt$id_date)) {
  # single individual and day
  df <- tempt[which(tempt$id_date == i),]
  # its thermal events
  circle <- unique(df$thermalID)
  all_obs <- data.frame()
  for (j in circle) {
    # single thermal
    event <- df[which(df$thermalID == j),]
    # the first observation of it (all observations hold the same amount of time)
    obs <- event[1,]
    # added to the first observation of all other thermals on that day
    all_obs <- rbind(all_obs, obs)
  }
  # sum the time spent on thermals per day and add it to the id_date df
  df$tpd <- sum(all_obs$circle_dt)
  thermal_times <- rbind(thermal_times, df)
  print(paste0("Calculated time spent in thermals for ", i, "."), quote = F)
}

linear_times <- data.frame()

for (i in unique(templ$id_date)) {
  # single individual and day
  df <- templ[which(templ$id_date == i),]
  # its thermal events
  line <- unique(df$linearID)
  all_obs <- data.frame()
  for (j in line) {
    # single thermal
    event <- df[which(df$linearID == j),]
    # the first observation of it (all observations hold the same amount of time)
    obs <- event[1,]
    # added to the first observation of all other thermals on that day
    all_obs <- rbind(all_obs, obs)
  }
  # sum the time spent on thermals per day and add it to the id_date df
  df$tpd <- sum(all_obs$line_dt)
  linear_times <- rbind(linear_times, df)
  print(paste0("Calculated time spent linear soaring for ", i, "."), quote = F)
}

rm(templ);rm(tempt)


thermal_times <- thermal_times %>% 
  group_by(id_date) %>% 
  slice(1)
linear_times <- linear_times %>% 
  group_by(id_date) %>% 
  slice(1)

## 5. Calculate the ratio of soaring time to flight time


# calculate the number of days between each observation and fledging
thermal_times$dsf <- as.numeric(as.Date(thermal_times$ymd) - as.Date(paste(year(thermal_times$fledging_date), month(thermal_times$fledging_date), day(thermal_times$fledging_date), sep = "-")))
linear_times$dsf <- as.numeric(as.Date(linear_times$ymd) - as.Date(paste(year(linear_times$fledging_date), month(linear_times$fledging_date), day(linear_times$fledging_date), sep = "-")))

# calculate the ratio of soaring to total time spent flying per day
thermal_times$ttf <- thermal_times$tpd/thermal_times$flight_time
linear_times$ltf <- linear_times$tpd/linear_times$flight_time

#save(rcfls, file = "rcfls_01.3.22.RData")

## 6. Plot the ratio of soaring:flight over days since fledging

thermal_plot <- ggplot(thermal_times, aes(x = dsf, y = ttf)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Relative time spent thermal soaring (s)") +
  scale_x_continuous(breaks = seq.int(0, 500, by = 50)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
linear_plot <- ggplot(linear_times, aes(x = dsf, y = ltf)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Relative time spent slope soaring (s)") + 
  scale_x_continuous(breaks = seq.int(0, 500, by = 50)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
ggarrange(thermal_plot,linear_plot, common.legend = F)

###########################################################################################

