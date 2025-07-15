### Polygons
# access the login to Movebank
# library(rgeos) 
library(move)
library(sf)
library(tidyverse)
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 12), text = element_text(size = 15)))

load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")
birds <- getMovebankReferenceTable(study = 282734839, login = loginStored) %>%
  drop_na(animal_id) %>%
  filter(sensor_type_id == 653)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs") # map projection
wgs_sf <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")

# the data pulled off Movebank
files <- list.files("/home/hbronnvik/Documents/chapter4/clean_data_24", full.names = T)

# check one bird at a time
# build <- getMovebankData(282734839, 4081264449, loginStored, removeDuplicatedTimestamps = T)
build <- readRDS(files[grep(1621089215, files)])
build@idData$individual.local.identifier
# downsample to reduce burden and clarify early area use
build <- build[timestamps(build) <= build$timestamp[1]+days(365)]
build <- build[!month(build$timestamp) %in% c(6, 7)]
build <- build[which(build$timelag.sec > 1),]
# visualize
build$obs <- 1:nrow(build)
mapview::mapview(build, zcol = "obs")

# hand selected points for a polygon
xs <- c(12.79452,12.79255,12.78397,12.77967,12.82465,12.84491,12.86104,12.83023)
ys <- c(47.55684,47.56581,47.56894,47.57936,47.61201,47.62184,47.59088,47.56240)
poly_xy <- data.frame(x = xs, y = ys, id = build@idData$individual.id)
# saveRDS(poly_xy, file = paste0("/home/hbronnvik/Documents/chapter4/polygons/", build@idData$individual.id, "_natal_est.rds"))
# visualize
home <- poly_xy %>%
  st_as_sf(coords = c("x", "y"), crs = wgs_sf) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")
mapview::mapview(build, zcol = "obs") + mapview::mapview(home, color = "black", alpha.regions = 0, lwd = 2)

###------------------------------
birds <- readRDS("/home/hbronnvik/Documents/chapter4/golden_eagle_ids.rds")
# the point data for the birds
files <- list.files("/home/hbronnvik/Documents/GE_data/no_burst", full.names = T)
# the polygons drawn above
polys <- lapply(list.files("/home/hbronnvik/Documents/chapter4/polygons", full.names = T), readRDS) %>% reduce(rbind)
# map projections
wgs_sf <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")
m_sf <- sf::st_crs("+proj=moll +ellps=WGS84")
# the dates from the Sempach work for comparison
sempach_dates <- read.csv("/home/hbronnvik/Documents/chapter4/expert_dates.csv") %>% 
  mutate(vis_date = as.Date(paste(expert_disp_date_year, expert_disp_date_month, expert_disp_date_day, sep = "-"))) %>% 
  dplyr::select(individual.local.identifier, vis_date)
# organize by name
file_ids <- unlist(lapply(files, function(x){sub("burst/", "", str_split(x, "_")[[1]][3])}))

birds <- data.frame(file = files) %>% 
  rowwise() %>% 
  mutate(individual.id = as.numeric(sub("burst/", "", str_split(file, "_")[[1]][3]))) %>% 
  ungroup() %>% 
  left_join(birds) %>% 
  arrange(animal_local_identifier)

# a threshold for defining what counts as having dispersed
at_home <- 25 # hours inside the territory

dispersals <- lapply(1:nrow(birds), function(n){
  # collect the data
  mv <- readRDS(birds$file[n])
  animal_name <- mv@idData$individual.id
  print(n);print(mv@idData$individual.local.identifier)
  if(animal_name %in% unique(polys$id)){
    # filter to the time range before the parents nest again
    mv <- mv[timestamps(mv) <= mv$timestamp[1]+days(365)]
    mv <- mv[!month(mv$timestamp) %in% c(6,7)]
    # filter out the bursts
    mv <- mv[which(mv$timelag.sec > 1),]
    if(nrow(mv) > 0){
      # make a polygon from the x,y coordinates recorded for the bird
      # measure the distance from each point to the polygon
      bird <- mv %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c("location.long", "location.lat"), crs = wgs_sf) %>% 
        sf::st_transform(crs = m_sf)
      home_p <- polys %>%
        filter(id == animal_name) %>% # i <- unique(unique_polygons$id)[36]
        st_as_sf(coords = c("x", "y"), crs = wgs_sf) %>%
        summarize(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON") %>% 
        sf::st_transform(crs = m_sf)
      # intersect the movement and the territory
      inside <- st_distance(bird, home_p) %>% 
        as.numeric()
      ind <- mv %>% 
        as.data.frame() %>% 
        mutate(d = inside/1000) # convert to km
      
      # ggplot(ind, aes(location.long, location.lat, color = d==0)) +
      #   geom_point()
      
      ind <- ind %>% 
        mutate(outside = d > 0,
               new_state = ifelse(outside == lag(outside), F, T))
      ind$new_state[1] <- F
      ind <- ind %>% 
        mutate(state = cumsum(new_state)) %>% 
        group_by(state) %>% 
        mutate(time_in_state = as.numeric(difftime(timestamp[n()], timestamp[1], units = "hours"))) %>% 
        ungroup()
      
      dispersal <- ind %>% 
        filter(outside == F & time_in_state > at_home) %>% 
        slice(n()) %>% 
        dplyr::select(timestamp) %>% 
        deframe()
      
      if(mv@idData$individual.local.identifier %in% sempach_dates$individual.local.identifier){
        og_date <- sempach_dates %>% 
          filter(individual.local.identifier == mv@idData$individual.local.identifier) %>% 
          dplyr::select(vis_date) %>% 
          deframe()
      }else{og_date <- NA}
      # the number of days that exist after dispersal and within a year of tagging
      # this allows us to guess whether the bird really dispersed or died on an excursion
      if(length(dispersal) > 0){
        if(max(mv$timestamp) == dispersal){
          did_disperse <- "firebrick"
        }else{
          post_dispersal <- ind %>% 
            filter(timestamp > dispersal) %>% 
            summarize(range = as.numeric(difftime(timestamp[n()], timestamp[1], units = "days"))) %>% 
            deframe()
          did_disperse <- ifelse(post_dispersal > 10, "black", "firebrick") 
        }
      }else{
        did_disperse <- "firebrick"
      }
      
      p <- ind %>% 
        ggplot(aes(timestamp, d, color = time_in_state)) +
        geom_path() +
        geom_vline(xintercept = dispersal, lty = 2, color = did_disperse) +
        # geom_vline(aes(xintercept = dispersal+days(14)), lty = 2, color = "firebrick") +
        geom_vline(xintercept = as.POSIXct(og_date), col = "gray70") +
        labs(x = "Time", y = "Distance from natal territory", 
             color = "Hours without \nborder crossings",
             title = paste0("Individual: ", mv@idData$individual.local.identifier, "\n Date:", dispersal),
        )
      return(p)
    } 
  }
})
# remove entries from birds that did not disperse (leaving 79)
dispersals <- dispersals[!sapply(dispersals, function(x) is.null(x))]

# pdf("/home/hbronnvik/Documents/chapter4/figures/dispersal_plots_0417.pdf",
#     width = 11.7, height = 8.5, onefile = TRUE)
# for (i in 1:length(dispersals)) {
#   print(dispersals[[i]])
# }
# dev.off()

dispersal_dates <- lapply(1:nrow(birds), function(n){
  # if(birds$individual.id[n] %in% file_ids){
    print(n)
    # collect the data
    mv <- readRDS(birds$file[n])
    animal_name <- mv@idData$individual.id
    print(mv@idData$individual.local.identifier)
    if(animal_name %in% unique(polys$id)){
      # filter to the time range before the parents nest again
      mv <- mv[timestamps(mv) <= mv$timestamp[1]+days(365)]
      mv <- mv[!month(mv$timestamp) %in% c(6,7)]
      # filter out the bursts
      mv <- mv[which(mv$timelag.sec > 1),]
      # only consider birds that actually have data
      if(animal_name == 4228812291){ # Fliegenluckenmauer24 (eobs 8729) has post-mortem data (?)
        mv <- mv[date(mv$timestamp) < "2024-08-29"]
      }
      if(animal_name == 523652947){ # Mals18 (eobs 6225) was recaptured and driven
        mv <- mv[date(mv$timestamp) < "2018-11-03"]
      }
      if(nrow(mv) > 40){
        # measure the distance from each point to the polygon
        bird <- mv %>% 
          as.data.frame() %>% 
          sf::st_as_sf(coords = c("location.long", "location.lat"), crs = wgs_sf) %>% 
          sf::st_transform(crs = m_sf)
        home_p <- polys %>%
          filter(id == animal_name) %>% # i <- unique(unique_polygons$id)[36]
          st_as_sf(coords = c("x", "y"), crs = wgs_sf) %>%
          summarize(geometry = st_combine(geometry)) %>%
          st_cast("POLYGON") %>% 
          sf::st_transform(crs = m_sf)
        # intersect the movement and the territory
        inside <- st_distance(bird, home_p) %>% 
          as.numeric()
        ind <- mv %>% 
          as.data.frame() %>% 
          mutate(d = inside)
        
        # sp
        # make a polygon from the x,y coordinates recorded for the bird
        # poly_xy <- polys %>% 
        #   filter(id == animal_name)
        # xym <- cbind(poly_xy$x, poly_xy$y)
        # p <- Polygon(xym)
        # ps <- Polygons(list(p), unique(poly_xy$id))
        # sps <- SpatialPolygons(list(ps))
        # proj4string(sps) <- wgs
        # 
        # ind_sp <- mv %>%
        #   as.data.frame()
        # coordinates(ind_sp) <- ~location.long+location.lat
        # proj4string(ind_sp) <- proj4string(sps)
        # dists <- apply(gDistance(ind_sp, sps, byid=TRUE),2,min) %>%
        #   as.data.frame() %>%
        #   rownames_to_column()
        # colnames(dists)[2] <- "d"
        # ind <- ind_sp %>%
        #   as.data.frame() %>%
        #   rownames_to_column() %>%
        #   full_join(dists, by = "rowname")
        
        # ggplot(ind, aes(location.long, location.lat, color = d==0)) +
        #   geom_polygon(data = poly_xy, aes(x = x, y = y), color = "black") +
        #   geom_point()
        
        ind <- ind %>% 
          mutate(outside = d > 0,
                 new_state = ifelse(outside == lag(outside), F, T))
        ind$new_state[1] <- F
        ind <- ind %>% 
          mutate(state = cumsum(new_state)) %>% 
          group_by(state) %>% 
          mutate(time_in_state = as.numeric(difftime(timestamp[n()], timestamp[1], units = "hours"))) %>% 
          ungroup()
        
        dispersal <- ind %>% 
          filter(outside == F & time_in_state > at_home) %>% 
          slice(n()) %>% 
          dplyr::select(timestamp) %>% 
          deframe()
        
        # the number of days that exist after dispersal and within a year of tagging
        # this allows us to guess whether the bird really dispersed or died on an excursion
        post_dispersal <- ind %>% 
          filter(timestamp > dispersal)
        if(nrow(post_dispersal) > 0){
         post_dispersal <- post_dispersal %>% 
            summarize(range = as.numeric(difftime(timestamp[n()], timestamp[1], units = "days"))) %>% 
            deframe()
         did_disperse <- ifelse(post_dispersal > 10, T, F)
         }else{
           did_disperse <- F 
           post_dispersal <- 0
           }

        d <- data.frame(individual.id = animal_name,
                        individual.local.identifier = mv@idData$individual.local.identifier,
                        dispersal_date = dispersal,
                        did_disperse = did_disperse, 
                        post_dispersal = post_dispersal)
        return(d)
      }else{
        d <- data.frame(individual.id = animal_name,
                        individual.local.identifier = mv@idData$individual.local.identifier,
                        dispersal_date = NA,
                        did_disperse = F, 
                        post_dispersal = 0)
        return(d)
      }
    }else{
      d <- data.frame(individual.id = animal_name,
                      individual.local.identifier = mv@idData$individual.local.identifier,
                      dispersal_date = NA,
                      did_disperse = F, 
                      post_dispersal = 0)
      return(d)
    }
  # }else{
  #   d <- data.frame(individual.id = birds$animal_id[n],
  #                   individual.local.identifier = birds$animal_local_identifier[n],
  #                   dispersal_date = NA,
  #                   did_disperse = F)
  # }
}) %>% 
  reduce(rbind)%>% 
  arrange(individual.local.identifier)

dates <- dispersal_dates %>% 
  # add on the expert assessments for the birds that have them
  left_join(sempach_dates, by = join_by(individual.local.identifier))%>%
  # calculate how much difference there is between these and those
  mutate(change = as.numeric(difftime(date(dispersal_date), vis_date, units = "days"))) %>% 
  # add on tagging dates
  left_join(birds %>% dplyr::select(individual.id, deploy_on_timestamp))  %>% 
  # remove birds that actually did not disperse
  # filter(did_disperse == T) %>% 
  # hand-filter the birds that dispersed, but during a gap so we don't know when
  filter(!individual.id %in% c(#2191867481, 908235585, 1193193426, 1193192892,
                               # and the Swedish birds
                               2831595227, 2831607583, 4016104249, 4016082269, 4016095602))
# write.csv(dates, "/home/hbronnvik/Documents/chapter4/dispersal_dates_20231204.csv", row.names = F)
# saveRDS(dates, "/home/hbronnvik/Documents/chapter4/dispersal_dates_20250417.rds")

# png(filename = paste0("/home/hbronnvik/Documents/chapter4/figures/dates_difference_0724.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
dates %>% 
  mutate(stay = as.numeric(difftime(vis_date, deploy_on_timestamp, units = "days"))) %>% 
  ggplot(aes(stay, change)) +
  geom_hline(yintercept = c(-7,7), lty = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_point() +
  labs(x = "Duration of stay using the visual estimate (days)", y = "Change in estimate (days)") +
  scale_y_continuous(n.breaks = 50) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))
# dev.off()
