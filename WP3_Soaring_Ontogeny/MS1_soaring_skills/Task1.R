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

cfls <- cfls %>% 
  # only post-fledging and pre-emigration soaring events
  filter(thermalClust != "other") %>% 
  mutate(year = year(timestamp),
         month = month(timestamp),
         day = day(timestamp), 
         ymd = paste(year, month, day, sep = "-"))

## 2. give IDs to each unique linear event

# a data frame to hold the IDs
cfls2 <- data.frame()
# a vector to separate the individuals and days so that there is no chance that the data
# for one ends on a thermal and the next starts on a thermal, which are then counted as a single
# event
cfls$id_date <- paste0(cfls$local_identifier," ", cfls$ymd)
id_date <- unique(cfls$id_date)

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
  print(paste0("Assigned orographic event IDs for ", i, "."), quote = F)
}

# a data frame to hold the IDs
cfls3 <- data.frame()


# and to each thermal event
for (i in unique(cfls$id_date)) {
  df <- cfls2[which(cfls2$id_date == i), ]
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
  print(paste0("Assigned thermal event IDs for ", i, "."), quote = F)
}


## 3. measure the time spent in each "event"

cfls <- cfls3 %>%
  # work within slope soaring events
  group_by(id_date, linearID) %>%
  # calculate the time difference between the last and first observations
  mutate(line_dt = difftime(.$timestamp[n()],.$timestamp[1])) %>% 
  ungroup() %>% 
  # work within thermal soaring events
  group_by(id_date, thermalID) %>%
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
  group_by(id_date, linearID) %>% 
  # select the first observation
  slice(1) %>% 
  # in all the data
  ungroup() %>% 
  # per day
  group_by(id_date) %>% 
  # sum the seconds in linear events
  mutate(day_lines = sum(as.numeric(line_dt), na.rm = T)) %>% 
  # one value per day
  slice(1) %>% 
  ungroup()
Ccfls <- cfls %>% 
  # within each thermal soaring event 
  group_by(id_date, thermalID) %>% 
  # select the first observation
  slice(1) %>% 
  # in all the data
  ungroup() %>% 
  # per day
  group_by(id_date) %>% 
  # sum the seconds in linear events
  mutate(day_circles = sum(as.numeric(circle_dt), na.rm = T)) %>% 
  # one value per day
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(thermalClust, timestamp, local_identifier, ymd, id_date, linear, 
         linearID, thermal, thermalID, circle_dt, day_circles)


## 5. Calculate the ratio of linear to circular soaring

# merge the data by id_date so that a single frame holds both linear and thermal soaring
Rcfls <- merge(Lcfls, Ccfls, by = "id_date")

# append the fledging dates for each individual
Rcfls$fledging_date <- NA
rcfls <- data.frame()

inds <- unique(Rcfls$local_identifier.x)
for (i in inds) {
  df <- Rcfls[which(Rcfls$local_identifier.x == i),]
  df$fledging_date <- stages$date_fledging[which(stages$id == 
            eagle_names$age_name[which(eagle_names$local_identifier == i)])]
  rcfls <- rbind(rcfls, df)
  
}

# calculate the number of days between each observation and fledging
rcfls$dsf <- as.numeric(as.Date(rcfls$ymd.x) - as.Date(paste(year(rcfls$fledging_date), month(rcfls$fledging_date), day(rcfls$fledging_date), sep = "-")))

rcfls <- rcfls %>% 
  rowwise() %>% 
  mutate(ratiosoar = day_lines/day_circles)

ggplot(rcfls[which(rcfls$ratiosoar != Inf),], aes(x = dsf, y = ratiosoar, color = local_identifier.x)) + 
  geom_point() +
  #geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Ratio of slope to thermal soaring (s)") + 
  ylim(c(0, 150)) +
  scale_x_continuous(breaks = seq.int(0, 500, by = 50)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))


