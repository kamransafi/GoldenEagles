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

cfls_ls <- list.files("gps_acc_age")
cfls_ls <- paste0("gps_acc_age/", cfls_ls)

cfls <- data.frame()

for (i in cfls_ls) {
  load(i)
  cfls <- rbind(cfls,gaa_df)
  inds <- unique(cfls$individual_local_identifier[which(cfls$stage == 2)])
  cfls <- cfls[which(cfls$individual_local_identifier %in% inds),]
  print("Successfully loaded a file.")
}
ori_cfls <- cfls
#save(ori_cfls, file = "cfls.RData")
load("cfls.RData"); cfls <- ori_cfls

## give IDs to each unique linear event

# where the thermal cluster is linear soaring, True
cfls$linear[which(cfls$thermalClust == "linear")] <- T
# where it is circular or other, False
cfls$linear[which(cfls$thermalClust != "linear")] <- F
# assign a new number to each slope soaring event 
cfls$linear_event <- inverse.rle(within.list(rle(cfls$linear), 
                                             values[values] <- seq_along(values[values])))
cfls$linear_event[which(cfls$linear_event == 0)] <- NA
# give each a linear ID
cfls$linearID <- paste0(cfls$burstIDcorrect, "_", cfls$linear_event)
cfls$linearID[grep("NA", cfls$linearID)] <- NA