#code for using the SSF built in energy_landscape_modeling_method1.R to make flight suitability maps for regular peirods of time for the entire Alpine region
#May 16. 2022. Konstanz, DE.
#Elham Nourani, PhD.

library(tidyverse)


#I need to make predictions with the INLA model for the entire Alpine region. i.e. I need to append enough rows with NA values to the dataset 
#options: 1) do it in batches, 2) only make predictions for unique combos of dem and tri, then associate the prediction to the rest of the cells
#with the same values. 3) run on the cluster.
#AND, let's not forget that I need to do this for different timestamps (weeks)!! :p


# STEP 1: open topo and tracking data ----------------------------------------------------------------
#open topo stack for the Alps (layers for east and west Alps were merged in Alps_topo_layers_prep.R)
alps_topo_wgs <- readRDS("Alps_dem_tri_wgs.rds")

#how many rows do we have: 17,539,565
alps_topo_df <- alps_topo_wgs %>% 
  as.data.frame(xy = T)

saveRDS(alps_topo_df, file = "Alps_dem_tri_wgs_df.rds")

#unique dem-tri combos: 17,435,266 lol
alps_topo_df %>% 
  group_by(dem_100,TRI_100)

#open eagle data
load("alt_50_20_min_25_ind_static_time_ann_weeks.RData") #ann_cmpl

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig)), #120 unique weeks
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))  

# STEP 2: create new datasets: one for each interval ----------------------------------------------------------------

#what intervals? there are 120 weeks. maybe do week 1, week2, month1, month2, month 6, yr 1? (in total, one year post fledging)... it may have plateaued by then. 
#depending on how cooperative the cluster is, try more regular intervals (esp. early on)



#next steps: 1) create the new datasets to get an idea of what I'm working with. 2) figure out the cluster 
#alos, at some point, add temperature to the original model