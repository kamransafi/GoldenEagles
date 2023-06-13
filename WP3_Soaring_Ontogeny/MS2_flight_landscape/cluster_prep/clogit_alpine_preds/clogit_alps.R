#R script to predict the flyability of the alpine region for golden eagle juvies
#Elham Nourani, PhD
#Jan 31, 2023. Konstanz, DE


library(dplyr)
library(stringr)
library(survival)
library(purrr)


#STEP 1: open data, define variables and formula -----------------------------------------------------------------

data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds")
topo_df <- readRDS("topo_df_100_LF.rds")
ssf <- readRDS("ssf_model_for_raven.rds")

#STEP 2: make predictions for each week since dispersal
wks_ls <- split(data, data$weeks_since_emig)

n <- 10 #make n predictions for each week

areas <- lapply(wks_ls, function(x){
  
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  areas_wk <- data.frame()
  
  #make predictions 10 times for each week. select step lengths randomly for each trial
  
  for(i in 1:n){
    
    #generate a new dataset
    new_data <- topo_df %>%
      mutate(step_length = sample(x$step_length, nrow(topo_df), replace = T),
             step_length_z = (step_length - attr(data[,colnames(data) == "step_length_z"],'scaled:center'))/attr(data[,colnames(data) == "step_length_z"],'scaled:scale'),
             weeks_since_emig = one_week,
             weeks_since_emig_z =  (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
             stratum = sample(x$stratum, nrow(topo_df), replace = T))
    
    #predict using the model
    preds <- predict(ssf, newdata = new_data, type = "risk")
    preds_prob <- preds/(preds+1)
    
    preds_pr <- new_data %>%
      mutate(preds = preds,
             probs = preds_prob)
    
    #skip saving this file if on local machine. it's 1 GB
    saveRDS(preds_pr, paste0("alpine_preds_wk_", week_i, "_trial_", i, ".rds"))
    
    #calculate suitable areas
    area_.7 <- preds_pr %>%
      filter(probs >= 0.7) %>% 
      summarize(pixels = n()) %>% #count the 
      mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
             area_km2 = round(area_m2/1e6,3),
             week_since_dispersal = week_i,
             trial = i)
    
    
    areas_wk <- rbind(areas_wk, area_.7)
  }
  
  areas_wk
  
}) %>% 
  reduce(rbind)

saveRDS(areas, file = "alpine_preds_areas.rds")



