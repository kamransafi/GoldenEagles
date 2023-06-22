#script for making Alpine prediction maps for golden eagles
#21.06.2023. Canberra, AU.
#Elham Nourani, PhD.

library(foreach)
library(doParallel)

# Set the number of cores to be used
num_cores <- 8
registerDoParallel(cores = num_cores)

# Create a list of unique week values
week_list <- unique(data$weeks_since_emig)

# Empty dataframe to save the areas
areas_df <- data.frame(week_since_dispersal = integer(156),
                       area_m2 = numeric(156),
                       area_km2 = numeric(156))

# Create a foreach loop to iterate over the week values in parallel
foreach(week_i = week_list, .packages = c("ggplot2", "oce", "dplyr", "gtools", "glmmTMB")) %dopar% {
  # Generate a new dataset for the current week
  topo_df <- topo_df %>%
    mutate(step_length = mean(x$step_length),
           step_length_z = 0,
           weeks_since_emig = week_i,
           weeks_since_emig_z = (week_i - attr(data[, colnames(data) == "weeks_since_emig_z"], "scaled:center")) / attr(data[, colnames(data) == "weeks_since_emig_z"], "scaled:scale"),
           stratum_ID = NA,
           animal_ID = sample(x$animal_ID, nrow(topo_df), replace = TRUE))
  
  # Predict using the model
  topo_df$preds <- predict(TMB_M, newdata = topo_df, type = "link", newparams = F, allow.new.levels = F, se.fit = F)
  #Sys.time() - b #24 min for one week
  topo_df$probs <- gtools::inv.logit(preds)
  
  gc()
  
  # Extract flyable areas
  r_.6 <- topo_df %>% 
    filter(probs >= 0.6) %>% 
    dplyr::select(c("location.long", "location.lat", "probs"))
  
  # Calculate suitable areas
  area_.6 <- r_.6 %>% 
    summarize(pixels = n()) %>% 
    mutate(area_m2 = pixels * 100 * 100,
           area_km2 = round(area_m2 / 1e6, 3),
           week_since_dispersal = week_i)
  
  # Combine the results with the global data frame (areas_df)
  areas_df[one_week,] <- area_.6
  
  # Plot
  p <- ggplot() +
    geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = r_.6, aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = 0.2), guide = "none") +
    labs(title = paste0("Week ", as.numeric(week_i), " since dispersal"), x = "", y = "") +
    theme_void()
  
  ggsave(plot = p, filename = paste0(week_i, "_alpine_pred.png"), dev = "cairo_pdf", width = 7, height = 5, dpi = 300)
  
  #clean up
  rm(p, r_0.6); gc(); gc()

}

# Stop the parallel processing
stopImplicitCluster()
