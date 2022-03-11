## Task 2: wind speed exposure within thermals
## Hester Brønnvik
## hbronnvik@ab.mpg.de
## 10.03.2022

## 1. Identify the maximum wind speed experienced each day
## 2. Plot against days since fledging

library(lubridate)
library(ggplot2) 
library(ggpubr)

cfls_ls <- list.files("gps_acc_age") 
cfls_ls <- paste0("gps_acc_age/", cfls_ls)

cfls <- data.frame()

for (i in cfls_ls) {
  load(i)
  gaa_df <- gaa_df[which(gaa_df$stage == 2),]
  cfls <- rbind(cfls, gaa_df)
  rm(gaa_df)
  print(paste0("Successfully loaded ", str_sub(sub("_.RData", "", sub("\\/", "", gsub("gps_acc_age", "", i))), 1, -13), "."), quote = F)
}

# only the thermal soaring data
cfls <- cfls[which(cfls$thermalClust == "circular" & cfls$soarClust == "soar"),]

# append the fledging dates
cfls$fledging_date <- NA
tcfls <- cfls
cfls <- data.frame()

inds <- unique(tcfls$individual_local_identifier)
for (i in unique(tcfls$individual_local_identifier)) {
  df <- tcfls[which(tcfls$individual_local_identifier == i),]
  df$fledging_date <- stages$date_fledging[which(stages$id == 
                                                   eagle_names$age_name[which(eagle_names$local_identifier == i)])]
  cfls <- rbind(cfls, df)
  
}

# calculate the time difference between fledging and the data
cfls$dsf <- as.numeric(cfls$timestamp - cfls$fledging_date)


cfls$id_date <- paste(str_sub(cfls$individual_local_identifier, 1, -13), paste(year(cfls$timestamp), month(cfls$timestamp), day(cfls$timestamp), sep = "-"), sep = "_")

# calculate the wind speed at each location
# mod(180 + 180/pi * atan2(u,v), 360)
# sqrt(u^2 + v^2)
cfls$wind_speed <- sqrt((cfls$windX)^2 + (cfls$windY)^2)


# extract the maximum wind speed on each day
pd_ws <- data.frame()

for (i in unique(cfls$id_date)) {
  # for each unique individual on each day
  df <- cfls[which(cfls$id_date == i),]
  # identify the maximum wind speed encountered
  df2 <- df[which.max(df$wind_speed),]
  # save the date and maximum for plotting
  pd_ws <- rbind(pd_ws, df2)
  # signal
  print(paste0("Separated wind speed for ", i, "."), quote = F)
  rm(df);rm(df2);rm(i)
}

## 2. Plot
ggplot(pd_ws, aes(x = dsf, y = wind_speed, color = individual_local_identifier)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Days since fledging", y = "Maximum wind speed") +
  scale_x_continuous(breaks = seq.int(0, max(pd_ws$dsf), by = 40)) +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))



