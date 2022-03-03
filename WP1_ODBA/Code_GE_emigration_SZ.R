## DATA CLEANING AND PREPARATION - PAPER GOLDEN EAGLE ACTIVITY
# Svea-Sophie Zimmermann
# 06.05.2021


# necessary libraries

library(tidyverse)
library(dplyr)
library(magrittr)
library(readr)
library(lubridate)
library(geosphere)
library(move)
library(stringr)
library(shiny)
library(viridis)
library(calibrate)
library(spatstat)
library(raster)
library(swfscMisc)
library(rgdal)
library(data.table)


setwd("~/Vogelwarte/Data/GPS")

## DATA CLEANUP ---------------------------------------------------
# (following Kami's book) 

setwd("~/Vogelwarte/Data/GPS/Raw data")
GE <- list.files(pattern = "*.csv", full.names = F)
GE <- sapply(GE, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id")
setwd("~/Vogelwarte/Data/GPS")

# check timestamps and remove timestamp duplicates
# then check for location duplicates

GE_red <- GE %>%
  mutate(timestamp = as.POSIXct(timestamp, format="%F %T ", tz="UTC")) %>%
  distinct(timestamp, `location-long`, `location-lat`, `individual-local-identifier`, .keep_all= TRUE) %>%
  arrange(`individual-local-identifier`, timestamp)

dup <- anyDuplicated(GE_red[,c('timestamp', 'individual-local-identifier')], 
                     fromLast = TRUE)
table(dup)  # row 1553899 seems to have a location duplicate


# remove duplicate locations

while (anyDuplicated(GE_red[,c('timestamp', 'individual-local-identifier')], 
                     incomparables = FALSE, fromLast = TRUE)) {
  dist1 <- sum(distHaversine(GE_red[c(dup-2, dup+1), c("location-long", "location-lat")],
                             GE_red[c(dup-1), c("location-long", "location-lat")]))
  dist2 <- sum(distHaversine(GE_red[c(dup-2, dup+1), c("location-long", "location-lat")], 
                             GE_red[dup, c("location-long", "location-lat")]))
  if(dist1 < dist2)
  {
    GE_red <- GE_red[-dup,]
  }else{
    GE_red <- GE_red[-(dup-1),]
  }
}

#GE_red <- GE_red[GE_red$`event-id` !=5318992185,] # I removed a second duplicate manually, because it was overlooked for some reason


# thin data to one point per hour

GE_red <- GE_red %>%
  drop_na(`location-long`, `location-lat`) # drop NA locations now, otherwise we might retain empty points for first point per hour

GE_red <- GE_red %>%  # extract year, month, day and hour from timestamps in new columns
  mutate(year = format(ymd_hms(as.character(GE_red$timestamp, tz = "UTC")),'%Y')) %>%
  mutate(month = format(ymd_hms(as.character(GE_red$timestamp, tz = "UTC")),'%m')) %>%
  mutate(day = format(ymd_hms(as.character(GE_red$timestamp, tz = "UTC")),'%d')) %>%
  mutate(hour = format(ymd_hms(as.character(GE_red$timestamp, tz = "UTC")),'%H'))

dup <- duplicated(GE_red[, c('hour', 'day', 'month', 'year', 'individual-local-identifier')], 
                  fromLast = TRUE)  # checking to see how many points we will retain
table(dup)

GE1h <- GE_red[!duplicated(GE_red[,c("hour", "day", "month", "year", "individual-local-identifier")]),] # retains first data point per hour by day, month, year and bird


# 9 of the 2019 birds sent locations throughout the night in July -> drop those

GE1h <- GE1h[!(GE1h$hour > '19' | GE1h$hour < '04'),] # retains data between 4am and 7pm (summer schedule)


# save file

GE1h$id <- str_replace_all(GE1h$id, ".csv", '') # remove .csv from id
names(GE1h) <- gsub(x = names(GE1h), pattern = "-", replacement = "_") # replace "-" in all column headers for easier selection
GE1h <- as.data.frame(GE1h) # as dataframe so it can easily be converted into a move object

write.table(GE1h, file="eagle_1h.csv", sep=",", row.names = FALSE)



## SOME CODE FOR INSPECTION ------------------------------------------

dat <- read_csv("eagle_1h.csv", 
                locale = locale(encoding = "latin1"))
dat <- as.data.frame(dat)


# turn dataframe into move object

mALL <- move(x=dat$location_long, y=dat$location_lat,
             time=as.POSIXct(dat$timestamp,
                             format="%Y-%m-%d %H:%M:%S", tz="UTC"), 
             data=dat, 
             proj=CRS("+proj=longlat"), 
             animal=dat$id
)

# see names for selecting individuals
namesIndiv(mALL)

# check how many locations per animal
n.locs(mALL) # select individual with [["id"]]

# check time distribution
timeLags <- timeLag(mALL, units='hours')
timeLagsVec <- unlist(timeLags)
summary(timeLagsVec)

max(timeLagsVec)

hist(timeLagsVec, breaks=50, main=NA, xlab="Time lag in hours")

hist(timeLagsVec[timeLagsVec<15], breaks="FD",
     main=NA, xlab="Time lag in hours")
summary(as.factor(round(timeLagsVec, 1)), maxsum=12)

# check speed
speeds <- unlist(speed(mALL))
summary(speeds)  # maximum speed seems to be within limit of GE (up to 35 m/s)



## DATA PREPARATION FOR DISPERSAL VISUALIZATION --------------------------

dat <- read_csv("eagle_1h.csv", locale = locale(encoding = "latin1"))
dat <- subset(dat, select = c('id', 'event_id', 'timestamp', 'location_long', 'location_lat',
                              'individual_local_identifier', 'sensor_type'))  # subset for cleaner table

# unlist individual data sets

my_list <- split(dat, dat$id)
list2env(my_list, envir = .GlobalEnv)


# create new columns containing nest locations

Drosloeng17 <- Drosloeng17 %>%
  mutate(nest_long="9.97340", .after = 5) %>%
  mutate(nest_lat="46.63350", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Viluoch17 <- Viluoch17 %>%
  mutate(nest_long="9.93533", .after = 5) %>%
  mutate(nest_lat="46.64785", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Almen18 <- Almen18 %>%
  mutate(nest_long="9.47662", .after = 5) %>%
  mutate(nest_lat="46.74778", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Romerio18 <- Romerio18 %>%
  mutate(nest_long="10.11425", .after = 5) %>%
  mutate(nest_lat="46.28694", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Guestizia18 <- Guestizia18 %>%
  mutate(nest_long="10.08859", .after = 5) %>%
  mutate(nest_lat="46.72670", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Mals18 <- Mals18 %>%
  mutate(nest_long="10.57476", .after = 5) %>%
  mutate(nest_lat="46.73462", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Nalps18 <- Nalps18 %>%
  mutate(nest_long="8.77537", .after = 5) %>%
  mutate(nest_lat="46.64801", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Schlanders18 <- Schlanders18 %>%
  mutate(nest_long="10.78324", .after = 5) %>%
  mutate(nest_lat="46.67331", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Schlappin1.18 <- Schlappin1.18 %>%
  mutate(nest_long="9.93796", .after = 5) %>%
  mutate(nest_lat="46.90393", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Schlappin2.18 <- Schlappin2.18 %>%
  mutate(nest_long="9.93796", .after = 5) %>%
  mutate(nest_lat="46.90393", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Tasna18 <- Tasna18 %>%
  mutate(nest_long="10.20200", .after = 5) %>%
  mutate(nest_lat="46.79734", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Umbrail18 <- Umbrail18 %>%
  mutate(nest_long="10.44419", .after = 5) %>%
  mutate(nest_lat="46.58266", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Almen19 <- Almen19 %>%
  mutate(nest_long="9.47662", .after = 5) %>%
  mutate(nest_lat="46.74778", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Dischma1.19 <- Dischma1.19 %>%
  mutate(nest_long="9.87436", .after = 5) %>%
  mutate(nest_lat="46.76993", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Dischma2.19 <- Dischma2.19 %>%
  mutate(nest_long="9.87436", .after = 5) %>%
  mutate(nest_lat="46.76993", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Fahrntal19 <- Fahrntal19 %>%
  mutate(nest_long="11.393766", .after = 5) %>%
  mutate(nest_lat="46.755829", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Flueela19 <- Flueela19 %>%
  mutate(nest_long="9.88206", .after = 5) %>%
  mutate(nest_lat="46.80506", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Grosio19 <- Grosio19 %>%
  mutate(nest_long="10.26932", .after = 5) %>%
  mutate(nest_lat="46.31300", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Kastelbell19 <- Kastelbell19 %>%
  mutate(nest_long="10.9301952", .after = 5) %>%
  mutate(nest_lat="46.6578762", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Matsch19 <- Matsch19 %>%
  mutate(nest_long="10.70951", .after = 5) %>%
  mutate(nest_lat="46.744396", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Nalps19 <- Nalps19 %>%
  mutate(nest_long="8.77537", .after = 5) %>%
  mutate(nest_lat="46.64801", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Sampuoir1.19 <- Sampuoir1.19 %>%
  mutate(nest_long="10.20885", .after = 5) %>%
  mutate(nest_lat="46.74940", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Sampuoir2.19 <- Sampuoir2.19 %>%
  mutate(nest_long="10.20885", .after = 5) %>%
  mutate(nest_lat="46.74940", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Seta19 <- Seta19 %>%
  mutate(nest_long="9.72978", .after = 5) %>%
  mutate(nest_lat="46.83071", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Siat19 <- Siat19 %>%
  mutate(nest_long="9.16043", .after = 5) %>%
  mutate(nest_lat="46.80863", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Sinestra1.19 <- Sinestra1.19 %>%
  mutate(nest_long="10.32794", .after = 5) %>%
  mutate(nest_lat="46.84676", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Sinestra2.19 <- Sinestra2.19 %>%
  mutate(nest_long="10.32794", .after = 5) %>%
  mutate(nest_lat="46.84676", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Torta19 <- Torta19 %>%
  mutate(nest_long="10.04729", .after = 5) %>%
  mutate(nest_lat="46.64272", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Trenzeira19 <- Trenzeira19 %>%
  mutate(nest_long="10.17477", .after = 5) %>%
  mutate(nest_lat="46.60908", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Tuors1.19 <- Tuors1.19 %>%
  mutate(nest_long="9.77489", .after = 5) %>%
  mutate(nest_lat="46.64178", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Tuors2.19 <- Tuors2.19 %>%
  mutate(nest_long="9.77489", .after = 5) %>%
  mutate(nest_lat="46.64178", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

ValGrande19 <- ValGrande19 %>%
  mutate(nest_long="10.387419", .after = 5) %>%
  mutate(nest_lat="46.280201", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Windlahn19 <- Windlahn19 %>%
  mutate(nest_long="11.4121", .after = 5) %>%
  mutate(nest_lat="46.628117", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Avers20 <- Avers20 %>%
  mutate(nest_long="9.49018", .after = 5) %>%
  mutate(nest_lat="46.47156", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Cornasc20 <- Cornasc20 %>%
  mutate(nest_long="10.11553", .after = 5) %>%
  mutate(nest_lat="46.29438", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Flueela20 <- Flueela20 %>%
  mutate(nest_long="9.88206", .after = 5) %>%
  mutate(nest_lat="46.80506", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Punteglias20 <- Punteglias20 %>%
  mutate(nest_long="8.97243", .after = 5) %>%
  mutate(nest_lat="46.76106", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Siat20 <- Siat20 %>%
  mutate(nest_long="9.16043", .after = 5) %>%
  mutate(nest_lat="46.80863", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Sils20 <- Sils20 %>%
  mutate(nest_long="9.76896", .after = 5) %>%
  mutate(nest_lat="46.44910", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Stuerfis20 <- Stuerfis20 %>%
  mutate(nest_long="9.62078", .after = 5) %>%
  mutate(nest_lat="47.04375", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)

Trimmis20 <- Trimmis20 %>%
  mutate(nest_long="9.58673", .after = 5) %>%
  mutate(nest_lat="46.87403", .after = 6) %>%
  mutate_at("nest_long", as.numeric) %>%
  mutate_at("nest_lat", as.numeric)


# add distance to nest with distHaversine function (calculates dist between two longlat points in m)

GE1h_disp <- rbind(Almen18, Almen19, Avers20, Cornasc20, Dischma1.19, Dischma2.19,
                   Drosloeng17, Fahrntal19, Flueela19, Flueela20, Grosio19, Guestizia18,
                   Kastelbell19, Mals18, Matsch19, Nalps18, Nalps19, Punteglias20, Romerio18,
                   Sampuoir1.19, Sampuoir2.19, Schlanders18, Schlappin1.18, Schlappin2.18, 
                   Seta19, Siat19, Siat20, Sils20, Sinestra1.19, Sinestra2.19, Stuerfis20, 
                   Tasna18, Torta19, Trenzeira19, Trimmis20, Tuors1.19, Tuors2.19, Umbrail18,
                   ValGrande19, Viluoch17, Windlahn19)


GE1h_disp$dist_to_nest <- distHaversine(GE1h_disp[, c('location_long', 'location_lat')], 
                                        GE1h_disp[, c('nest_long', 'nest_lat')])


write.table(GE1h_disp, file="eagle1h_disp.csv", sep=",", row.names = FALSE)


## DISPERSAL VISUALIZER For PRELIM EMIGRATION -----------------------------------
# adapted from Dispersal Visualizer v.0.2 Vogelwarte ugk, March 18th 2019

dat <- read_csv("eagle1h_natal_phases.csv", locale = locale(encoding = "latin1"))
sdat <-  dat[dat$id =="Trenzeira19", ]  # select each bird here

sdat <- sdat[,c(4,5,6,7,3,13,8,9)]  # select only relevant columns (long lat, id, dist to nest etc)
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat


# interactive plot to determine preliminary date of emigration with dist to nest based on visual pattern

ui <- basicPage(
  plotOutput("plot1", click = "plot_click"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    par(mfrow=c(1,2))
    
    mycol = plasma(length(sdat$timestamp))
    
    
    plot(sdat$x1, sdat$y1, pch = 20, cex=1, col = "white", 
         xlab = "Longitude", ylab = "Latitude", main = " ")
    lines(sdat$location_long, sdat$location_lat, col = "black")
    
    points(sdat$x1, sdat$y1, pch = 20, cex=1, col = mycol, 
           xlab = "x", ylab = "y")
    #textxy(sdat$x1, sdat$y1, sdat$timestamp, cex = 0.7)
    
    #Neststandort
    points(sdat[1,]$nest_long, sdat[1,]$nest_lat, pch = 23, cex = 1, col = "black", bg="green")
    points(sdat[1,]$cent_long, sdat[1,]$cent_lat, pch = 20, cex = 1, col = "black", bg="black")
    
    plot(sdat$timestamp, sdat$dist_to_cent, pch = 20, abline(h=8000), col = mycol, 
         xlab = "Year 2019-2020",  cex = 1.2, ylab = "Distance to centroid [m]")
    lines(sdat$timestamp, sdat$dist_to_cent, col = "black")
    #textxy(sdat$timestamp, sdat$dist_to_nest, sdat$timestamp)
    
  })  
  
  output$info <- renderPrint({
    # With base graphics, need to tell it what the x and y variables are.
    nearPoints(sdat[,1:10], input$plot_click, xvar = "timestamp", yvar = "dist_to_cent")
    # nearPoints() also works with hover and dblclick events
  })
}

shinyApp(ui, server)


## EXTRACT CENTROIDS OF NATAL TERRITORIES --------------------------------

dat <- read_csv("eagle1h_disp.csv", locale = locale(encoding = "latin1"))

# Droslöng17
sdat <-  dat[dat$id =="Drosloeng17", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2017-08-06 17:00:00") & timestamp <= as.POSIXct("2018-03-20 09:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Viluoch17
sdat <-  dat[dat$id =="Viluoch17", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2017-07-27 09:00:00") & timestamp <= as.POSIXct("2018-03-07 13:30:09"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # add individual lat long range
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Almen18
sdat <-  dat[dat$id =="Almen18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-13 15:20:00") & timestamp <= as.POSIXct("2019-01-30 13:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# ArtSanRomerio18
sdat <-  dat[dat$id =="Romerio18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-04 13:30:00") & timestamp <= as.POSIXct("2019-03-16 09:00:10"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Güstizia18
sdat <-  dat[dat$id =="Guestizia18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-07-28 11:00:00") & timestamp <= as.POSIXct("2019-02-23 10:00:53"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Mals18
sdat <-  dat[dat$id =="Mals18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-09-02 16:00:12"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Nalps18
sdat <-  dat[dat$id =="Nalps18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-13 11:00:00") & timestamp <= as.POSIXct("2018-10-01 13:00:38"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Schlanders18
sdat <-  dat[dat$id =="Schlanders18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-27 13:00:41") & timestamp <= as.POSIXct("2018-09-12 11:00:10"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Schlappin1.18
sdat <-  dat[dat$id =="Schlappin1.18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-02 10:00:00") & timestamp <= as.POSIXct("2019-03-26 12:00:09"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Schlappin2.18
sdat <-  dat[dat$id =="Schlappin2.18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-07-30 13:30:00"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Tasna18
sdat <-  dat[dat$id =="Tasna18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-07-28 09:20:00") & timestamp <= as.POSIXct("2019-04-06 13:00:10"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Umbrail18
sdat <-  dat[dat$id =="Umbrail18", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2018-08-07 12:00:00") & timestamp <= as.POSIXct("2019-03-28 08:00:20"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Almen19
sdat <-  dat[dat$id =="Almen19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-22 15:00:00") & timestamp <= as.POSIXct("2020-02-14 15:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Dischma1.19
sdat <-  dat[dat$id =="Dischma1.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-30 10:00:00") & timestamp <= as.POSIXct("2020-03-04 12:00:28"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Dischma2.19
sdat <-  dat[dat$id =="Dischma2.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-01 12:00:00") & timestamp <= as.POSIXct("2020-03-28 14:00:09"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


#Fahrntal19
sdat <-  dat[dat$id =="Fahrntal19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-03 08:00:41") & timestamp <= as.POSIXct("2019-11-14 11:00:08"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Flüela19
sdat <-  dat[dat$id =="Flueela19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-30 10:00:53") & timestamp <= as.POSIXct("2020-01-10 12:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Grosio19
sdat <-  dat[dat$id =="Grosio19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-27 10:30:12") & timestamp <= as.POSIXct("2020-03-07 08:00:08"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Kastelbell19
sdat <-  dat[dat$id =="Kastelbell19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-08 07:40:23") & timestamp <= as.POSIXct("2020-02-05 10:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Matsch19
sdat <-  dat[dat$id =="Matsch19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-24 12:00:23") & timestamp <= as.POSIXct("2020-03-08 09:00:08"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.4, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Nalps19
sdat <-  dat[dat$id =="Nalps19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-01 12:00:00") & timestamp <= as.POSIXct("2020-02-01 19:01:53"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Sampuoir1.19
sdat <-  dat[dat$id =="Sampuoir1.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-29 15:00:00") & timestamp <= as.POSIXct("2020-02-23 15:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Sampuoir2.19
sdat <-  dat[dat$id =="Sampuoir2.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-05 05:20:00") & timestamp <= as.POSIXct("2020-02-29 12:00:09"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Seta19
sdat <-  dat[dat$id =="Seta19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-03 07:00:00") & timestamp <= as.POSIXct("2020-02-22 12:00:10"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Siat19
sdat <-  dat[dat$id =="Siat19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-01 05:20:00"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Sinestra1.19
sdat <-  dat[dat$id =="Sinestra1.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-23 10:00:00") & timestamp <= as.POSIXct("2020-03-18 11:00:00"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Sinestra2.19
sdat <-  dat[dat$id =="Sinestra2.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-27 11:00:00") & timestamp <= as.POSIXct("2020-03-06 09:00:13"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Torta19
sdat <-  dat[dat$id =="Torta19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-20 11:30:00") & timestamp <= as.POSIXct("2020-02-21 07:20:41"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Trenzeira19
sdat <-  dat[dat$id =="Trenzeira19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-08-06 14:00:16") & timestamp <= as.POSIXct("2020-02-12 10:00:11"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Tuors1.19
sdat <-  dat[dat$id =="Tuors1.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-26 08:00:00") & timestamp <= as.POSIXct("2019-12-13 19:00:41"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Tuors2.19
sdat <-  dat[dat$id =="Tuors2.19", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2019-07-25 19:00:00") & timestamp <= as.POSIXct("2019-12-13 15:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.08, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Avers20
sdat <-  dat[dat$id =="Avers20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-31 07:00:00") & timestamp <= as.POSIXct("2021-03-10 11:00:00"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.1, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Cornasc20
sdat <-  dat[dat$id =="Cornasc20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-28 10:30:00") & timestamp <= as.POSIXct("2021-04-08 14:00:08"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Flüela20
sdat <-  dat[dat$id =="Flüela20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-08-06 10:30:00") & timestamp <= as.POSIXct("2020-09-09 13:00:00"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Punteglias20
sdat <-  dat[dat$id =="Punteglias20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-23 19:00:41") & timestamp <= as.POSIXct("2020-11-13 13:00:10"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Siat20
sdat <-  dat[dat$id =="Siat20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-16 04:40:00") & timestamp <= as.POSIXct("2021-02-28 13:00:23"))

# create density plot for time after fledging and before emigration
# # and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr))) # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Sils20
sdat <-  dat[dat$id =="Sils20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-29 10:30:00") & timestamp <= as.POSIXct("2021-04-01 14:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Stürfis20
sdat <-  dat[dat$id =="Stuerfis20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-30 10:00:00") & timestamp <= as.POSIXct("2021-03-20 10:00:23"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# Trimmis20
sdat <-  dat[dat$id =="Trimmis20", ]

sdat <- sdat[,c(4,5,6,7,3,10)]
sdat$x1 <- sdat$location_long
sdat$y1 <- sdat$location_lat

# subset with period between fledging and emigration
sdat_natal <- sdat %>%
  filter(timestamp >= as.POSIXct("2020-07-31 8:00:00") & timestamp <= as.POSIXct("2020-12-06 13:00:41"))

# create density plot for time after fledging and before emigration
# and calculate centroid coordinates for natal territory

rlon <- range(sdat_natal$location_long) # range is needed for creating a window in next step
rlat <- range(sdat_natal$location_lat)

w <- owin(rlon, rlat)    # range of lat lon sets window for next step
natal <- ppp(sdat_natal$location_long, sdat_natal$location_lat, window=w)
z <- density.ppp(natal, sigma = 0.05, edge = FALSE)
dr = raster(z)
centr <- as.data.frame(xyFromCell(dr, which.max(dr)))  # calculates centroid coordinates for "brightest pixel" which represents the middle
plot(z)
points(natal)
points(centr, col="red")


# create new columns containing centroid locations

GE_disp <- read_csv("eagle1h_disp.csv", 
                    locale = locale(encoding = "latin1"))

# unlist individual data sets
my_list <- split(GE_disp, GE_disp$id)
list2env(my_list, envir = .GlobalEnv)


Drosloeng17 <- Drosloeng17 %>%
  as.data.frame() %>%
  mutate(cent_long="9.974545", .after = 7) %>%
  mutate(cent_lat="46.64245", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Viluoch17 <- Viluoch17 %>%
  as.data.frame() %>%
  mutate(cent_long="9.935285", .after = 7) %>%
  mutate(cent_lat="46.65311", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Almen18 <- Almen18 %>%
  as.data.frame() %>%
  mutate(cent_long="9.485651", .after = 7) %>%
  mutate(cent_lat="46.74960", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Romerio18 <- Romerio18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.11959", .after = 7) %>%
  mutate(cent_lat="46.28770", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Guestizia18 <- Guestizia18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.05765", .after = 7) %>%
  mutate(cent_lat="46.72517", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Mals18 <- Mals18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.57402", .after = 7) %>%
  mutate(cent_lat="46.73634", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Nalps18 <- Nalps18 %>%
  as.data.frame() %>%
  mutate(cent_long="8.762783", .after = 7) %>%
  mutate(cent_lat="46.63818", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Schlanders18 <- Schlanders18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.77231", .after = 7) %>%
  mutate(cent_lat="46.66919", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Schlappin1.18 <- Schlappin1.18 %>%
  as.data.frame() %>%
  mutate(cent_long="9.932979", .after = 7) %>%
  mutate(cent_lat="46.90496", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Schlappin2.18 <- Schlappin2.18 %>%
  as.data.frame() %>%
  mutate(cent_long="9.941998", .after = 7) %>%
  mutate(cent_lat="46.90391", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Tasna18 <- Tasna18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.20496", .after = 7) %>%
  mutate(cent_lat="46.80255", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Umbrail18 <- Umbrail18 %>%
  as.data.frame() %>%
  mutate(cent_long="10.43750", .after = 7) %>%
  mutate(cent_lat="46.57738", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Almen19 <- Almen19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.481770", .after = 7) %>%
  mutate(cent_lat="46.74760", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Dischma1.19 <- Dischma1.19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.889641", .after = 7) %>%
  mutate(cent_lat="46.76480", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Dischma2.19 <- Dischma2.19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.893993", .after = 7) %>%
  mutate(cent_lat="46.76714", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Fahrntal19 <- Fahrntal19 %>%
  as.data.frame() %>%
  mutate(cent_long="11.40867", .after = 7) %>%
  mutate(cent_lat="46.77367", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Flueela19 <- Flueela19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.883983", .after = 7) %>%
  mutate(cent_lat="46.80626", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Grosio19 <- Grosio19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.28703", .after = 7) %>%
  mutate(cent_lat="46.33298", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Kastelbell19 <- Kastelbell19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.92447", .after = 7) %>%
  mutate(cent_lat="46.65614", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Matsch19 <- Matsch19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.66710", .after = 7) %>%
  mutate(cent_lat="46.72800", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Nalps19 <- Nalps19 %>%
  as.data.frame() %>%
  mutate(cent_long="8.761305", .after = 7) %>%
  mutate(cent_lat="46.63917", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Sampuoir1.19 <- Sampuoir1.19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.20102", .after = 7) %>%
  mutate(cent_lat="46.74399", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Sampuoir2.19 <- Sampuoir2.19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.20010", .after = 7) %>%
  mutate(cent_lat="46.74454", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Seta19 <- Seta19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.729777", .after = 7) %>%
  mutate(cent_lat="46.82978", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Siat19 <- Siat19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.159059", .after = 7) %>%
  mutate(cent_lat="46.81337", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Sinestra1.19 <- Sinestra1.19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.327010", .after = 7) %>%
  mutate(cent_lat="46.86307", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Sinestra2.19 <- Sinestra2.19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.327680", .after = 7) %>%
  mutate(cent_lat="46.86840", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Torta19 <- Torta19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.030510", .after = 7) %>%
  mutate(cent_lat="46.64838", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Trenzeira19 <- Trenzeira19 %>%
  as.data.frame() %>%
  mutate(cent_long="10.185560", .after = 7) %>%
  mutate(cent_lat="46.61689", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Tuors1.19 <- Tuors1.19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.786865", .after = 7) %>%
  mutate(cent_lat="46.64471", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Tuors2.19 <- Tuors2.19 %>%
  as.data.frame() %>%
  mutate(cent_long="9.791679", .after = 7) %>%
  mutate(cent_lat="46.64550", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Avers20 <- Avers20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.492579", .after = 7) %>%
  mutate(cent_lat="46.44250", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Cornasc20 <- Cornasc20 %>%
  as.data.frame() %>%
  mutate(cent_long="10.11040", .after = 7) %>%
  mutate(cent_lat="46.28576", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Flueela20 <- Flueela20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.892892", .after = 7) %>%
  mutate(cent_lat="46.79224", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Punteglias20 <- Punteglias20 %>%
  as.data.frame() %>%
  mutate(cent_long="8.973435", .after = 7) %>%
  mutate(cent_lat="46.76727", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Siat20 <- Siat20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.161629", .after = 7) %>%
  mutate(cent_lat="46.81464", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Sils20 <- Sils20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.768729", .after = 7) %>%
  mutate(cent_lat="46.45514", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Stuerfis20 <- Stuerfis20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.624991", .after = 7) %>%
  mutate(cent_lat="47.04147", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)

Trimmis20 <- Trimmis20 %>%
  as.data.frame() %>%
  mutate(cent_long="9.585814", .after = 7) %>%
  mutate(cent_lat="46.87831", .after = 8) %>%
  mutate_at("cent_long", as.numeric) %>%
  mutate_at("cent_lat", as.numeric)


# add distance to centroid with distHaversine function (calculates dist between two latlong points in m)
library(geosphere)

GE_disp <- rbind(Almen18, Almen19, Avers20, Cornasc20, Dischma1.19, Dischma2.19, Drosloeng17,
                 Fahrntal19, Flueela19, Flueela20, Grosio19, Guestizia18, Kastelbell19, 
                 Mals18, Matsch19, Nalps18, Nalps19, Punteglias20, Romerio18, Sampuoir1.19,
                 Sampuoir2.19, Schlanders18, Schlappin1.18, Schlappin2.18, Seta19, Siat19,
                 Siat20, Sils20, Sinestra1.19, Sinestra2.19, Stuerfis20, Tasna18, Torta19,
                 Trenzeira19, Trimmis20, Tuors1.19, Tuors2.19, Umbrail18, Viluoch17)


GE_disp$dist_to_cent <- distHaversine(GE_disp[, c('location_long', 'location_lat')], 
                                      GE_disp[, c('cent_long', 'cent_lat')])

write.table(GE_disp, file="eagle1h_disp_nc.csv", sep=",", row.names = FALSE)



