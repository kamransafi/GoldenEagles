#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: energy landscape construction
#follows on from embc_segmentation.R
#Jan 28. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(EMbC)
library(mapview)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")