### find excursions and their qualities for Golden Eagles
### Hester Bronnvik
### 2024-10-12

library(move)
# library(rgeos)
library(sf)
library(tidyverse)
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 12), text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#6F101E", "#8F1D28", "#BF413E", "#FACEB6", "white", 
                              "#A4CCE3", "#3D84B4", "#1D5497", "#0A3162"))
###------------------------------
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")
birds <- getMovebankReferenceTable(study = 282734839, login = loginStored) %>%
  drop_na(animal_id) %>%
  filter(sensor_type_id == 653) %>% 
  dplyr::select(animal_local_identifier, animal_id, deploy_on_timestamp)

files <- list.files("/home/hbronnvik/Documents/GE_data/no_burst", full.names = T)
polys <- lapply(list.files("/home/hbronnvik/Documents/chapter4/polygons", full.names = T), readRDS) %>% 
  reduce(rbind)

# wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs") # map projection
m_sf <- sf::st_crs("+proj=moll +ellps=WGS84")
wgs_sf <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")

sempach_dates <- read.csv("/home/hbronnvik/Documents/chapter4/expert_dates.csv") %>% 
  mutate(vis_date = as.Date(paste(expert_disp_date_year, expert_disp_date_month, expert_disp_date_day, sep = "-"))) %>% 
  dplyr::select(individual.local.identifier, vis_date)

dates <- readRDS("/home/hbronnvik/Documents/chapter4/dispersal_dates_20250127.rds") %>% 
  mutate(pfdp = as.numeric(difftime(dispersal_date, deploy_on_timestamp, units = "days"))) 

file_ids <- unlist(lapply(files, function(x){sub("burst/", "", str_split(x, "_")[[1]][3])}))

files <- data.frame(file = list.files("/home/hbronnvik/Documents/GE_data/no_burst", full.names = T)) %>% 
  rowwise() %>% 
  mutate(file_id = sub("burst/", "", str_split(file, "_")[[1]][3])) %>% 
  ungroup() #%>% 
  # hand-filter the birds that dispersed, but during a gap so we don't know when
  # filter(!file_id %in% c(2191867481, 908235585, 1193193426, 1193192892,
  #                        # and the Swedish birds
  #                        2831595227, 2831607583,
  #                        # bird+tag recapped and driven
  #                        523652947) &
  #          # only consider birds that dispersed and thus for which we have all the excursions they did
  #          file_id %in% dates$individual.id)

birds <- birds %>% 
  filter(animal_id %in% files$file_id) %>% 
  arrange(animal_local_identifier)

dates$individual.local.identifier[!dates$individual.id %in% polys$id]

excursion_info <- lapply(1:nrow(files), function(n){
  print(n)
  # collect the data
  mv <- readRDS(files$file[n])
  animal_name <- mv@idData$individual.id
  if(animal_name %in% unique(polys$id)){
    # for birds that dispersed, filter to their time inside the natal territory
    dispersal <- dates %>% 
      drop_na(dispersal_date) %>% 
      filter(individual.id == animal_name) %>% 
      dplyr::select(dispersal_date) %>% 
      deframe()
    # for birds that did not disperse, use the final time
    if(length(dispersal) == 0){
      dispersal <- max(mv$timestamp[mv$timestamp < (mv$timestamp[1]+days(600))]) # ID 1600584021 has three data from 2 years after it died
    }
    # filter to the time range before dispersal/death
    mv <- mv[timestamps(mv) <= dispersal]
    # filter out the bursts
    # mv <- mv[which(mv$timelag.sec > 1),]
    if(nrow(mv) > 40){
      # make a polygon from the x,y coordinates recorded for the bird
      # make the polygon in sf and project to meters
      home_p <- polys %>%
        filter(id == animal_name) %>%
        st_as_sf(coords = c("x", "y"), crs = wgs_sf) %>%
        summarize(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON") %>% 
        sf::st_transform(crs = m_sf)
      # transform the data to meters projection
      bird <- mv %>% 
        as.data.frame() %>% 
        sf::st_as_sf(coords = c("location.long", "location.lat"), crs = wgs_sf) %>% 
        sf::st_transform(crs = m_sf)
      # measure the distance from each point to the MCP
      ind <- mv %>%
        as.data.frame() %>% 
        mutate(d = as.numeric(st_distance(bird, home_p, sparse = T)))
      
      if(!"heading" %in% colnames(ind)){
        ind$heading <- NA
        ind$study.name.1 <- NA
      }
      
      ind <- ind %>% 
        mutate(outside = d > 0,
               new_state = ifelse(outside == lag(outside), F, T))
      ind$new_state[1] <- F
      ind <- ind %>% 
        mutate(state = cumsum(new_state)) %>% 
        group_by(state) %>% 
        mutate(time_in_state = as.numeric(difftime(timestamp[n()], timestamp[1], units = "hours"))) %>% 
        ungroup()
      
      # remove unrealistically early data (only error could mean a bird was outside a polygon on day 1)
      # see Nalps19 (eobs 5861) whose tag was on for the drive from Siat19
      tagging <- dates %>% 
        filter(individual.id == animal_name) %>% 
        dplyr::select(deploy_on_timestamp) %>% 
        deframe()
      
      excursions <- ind %>% 
        filter(outside == T & date(timestamp) > tagging) 
      
      if(nrow(excursions) > 0){
        excursions <- excursions %>% 
          group_by(state) %>% 
          mutate(excursion_id = paste0(individual.id, "_", gsub("-", "_", date(min(timestamp))), "_", state),
                 excursion_duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "hours")),
                 # the highest distance achieved
                 excursion_dist = max(d)/1000,
                 # the distances between points
                 sl = distVincentyEllipsoid(cbind(location.long, location.lat), 
                                            cbind(lag(location.long), lag(location.lat))),
                 # all the distances between points + the distance from the polygon to the first point
                 # and + the distance from the last point back to the polygon
                 excursion_dist_sum = (sum(sl, na.rm = T)+d[1]+d[n()])/1000,
                 # how many points were captured
                 excursion_observations = n()) %>% 
          ungroup() %>% 
          mutate(poly_area = as.numeric(st_area(home_p)),
                 n_excursions = length(unique(excursion_id)))
        
        # ex_oi <- excursions %>% 
        #   filter(excursion_id == "2186854726_2023_03_27_109") %>% 
        #   mutate(sl = distVincentyEllipsoid(cbind(location.long, location.lat), 
        #                                     cbind(lag(location.long), lag(location.lat))))
        
        return(excursions)
      }
    }else{
      
      return(NULL)
    }
  }else{
    return(NULL)
  }
})
# leaving 84
excursion_info <- excursion_info[!sapply(excursion_info, function(x) is.null(x))]
excursion_info <- data.table::rbindlist(excursion_info, fill = T)

excursion_info <- excursion_info %>%
  arrange(individual.local.identifier)

# saveRDS(excursion_info, file = "/home/hbronnvik/Documents/chapter4/excursions_250127.rds")
excursion_info <- readRDS("/home/hbronnvik/Documents/chapter4/excursions_250127.rds")
