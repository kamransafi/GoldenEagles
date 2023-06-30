#script for constructing the energy landscape for golden eagles
#using glmmTMB, which allows for random effects to be included
#follows from 03_energy_landscape_modeling_method1.R, 04_INLA_output_plots.R; workflow similar to 03_04_clogit_workflow.R
#21.06.2023. Canberra, AU.
#Elham Nourani, PhD.

library(tidyverse)
library(glmmTMB)
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(corrr)
library(terra)
library(ggnewscale)


setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")
setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#colors
clr <- oce::oceColorsPalette(100)[9] #was [2] before
clr_light <- oce::oceColorsPalette(100)[10]
clr2 <- oce::oceColorsPalette(100)[80]

#open data. data was prepared in 03_energy_landscape_modeling_method1.R. 
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  mutate(animal_ID = as.numeric(as.factor(ind1)), #animal ID and stratum ID should be numeric
         stratum_ID = as.numeric(as.factor(stratum)))

#look at correlation
data %>% 
dplyr::select(c("TRI_100", "step_length", "weeks_since_emig", "ridge_100", "dem_100")) %>% 
  correlate() #TRI and ridge = 0.44

# STEP 1: run the model ------------------------------------------------------------------ 
#this is based on Muff et al:
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

TMB_struc <- glmmTMB(used ~ -1 + TRI_100_z * step_length_z * weeks_since_emig_z + 
                     ridge_100_z * step_length_z * weeks_since_emig_z + (1|stratum_ID) + 
                     (0 + ridge_100_z | animal_ID) + 
                     (0 + TRI_100_z | animal_ID), 
                   family = poisson, data = data, doFit = FALSE,
                   #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                   map = list(theta = factor(c(NA,1:2))), #2 is the n of random slopes
                   #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                   start = list(theta = c(log(1e3),0,0))) #add a 0 for each random slope. in this case, 2


TMB_m <- glmmTMB:::fitTMB(TMB_struc)
summary(TMB_m)

saveRDS(TMB_m, file = "TMB_model.rds")

#extract coefficient estimates and confidence intervals
confint(TMB_m)

#extract individual-specific random effects:
ranef(TMB_m)[[1]]$animal_ID

# STEP 2: model validation ------------------------------------------------------------------ 

#calculate the RMSE
performance_rmse(TMB_m)

# STEP 3: PLOT coefficient estimates ------------------------------------------------------------------ 

graph <- confint(TMB_m) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!(Factor %in% c("weeks_since_emig_z", "Std.Dev.ridge_100_z|animal_ID", "Std.Dev.TRI_100_z|animal_ID.1")))
  
colnames(graph)[c(2,3)] <- c("Lower", "Upper") 

         labels <- rev(c("TRI", "Step length", "Distance to ridge", "TRI: Step length", "TRI: Week",
               "Step length: Week", "Step length: Distance to ridge", "Distance to ridge: Week",
               "TRI: Step length: Week", "Distance to ridge: Step length: Week"))

VarOrder <- rev(unique(graph$Factor))
graph$Factor <- factor(graph$Factor, levels = VarOrder)

#plot
X11(width = 7, height = 2) 

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(color = clr, size = 1.7)  +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = clr, linewidth = .7) +
  labs(x = "Estimate", y = "") +
  scale_y_discrete(labels = labels) +
  xlim(-.71, .25) +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

ggsave(coefs, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/coeffs.pdf", 
       width = 7, height = 2, dpi = 300)

# STEP 4: PLOT individual-specific coefficients ------------------------------------------------------------------ 

rnd <- ranef(TMB_m) %>% 
  as.data.frame() %>% 
  filter(grpvar == "animal_ID") %>% 
  mutate(ID = rep(levels(as.factor(data$ind1)),2),
         variable = c(rep("Distance_to_ridge", 55), rep("TRI", 55)),
         Lower = condval - condsd,
         Upper = condval + condsd)  #assign individuals' names
  
order <- rnd %>% 
  filter(variable == "Distance_to_ridge") %>% 
  arrange(desc(condval)) %>% 
  pull(ID)

rnd$ID <- factor(rnd$ID, levels = order)

#cols <- c(TRI = "lightcoral",
#          Distance_to_ridge = "cornflowerblue")

cols <- c(TRI = clr2,
          Distance_to_ridge = clr)

X11(width = 7, height = 9)
coefs_inds <- ggplot(rnd, aes(x = condval, y = ID, color = variable)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(size = 2, position = position_dodge(width = .7))  +
  geom_linerange(aes(xmin = Lower, xmax = Upper), size = 0.8, position = position_dodge(width = .7)) +
  scale_color_manual(values = cols) + 
  scale_y_discrete(labels = order) +
  labs(x = "Difference from fixed effect estimate", y = "") +
  theme_minimal() +
  theme(text = element_text(size = 8), #font size should be between 6-8
        axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
        axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6)))

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/ind_coefs.pdf", 
       width = 7, height = 9, dpi = 300)


# STEP 5: predictions for interaction terms + PLOT ------------------------------------------------------------------ 

#new dataset created in 03_energy_landscape_modeling_method1.r ... this file was updated to have stratum_ID and animal_ID as numeric variables
#accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
#(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
new_data <- readRDS("new_data_only_ssf_preds.rds") #%>% 
  mutate(stratum_ID = NA)
  #       animal_ID = NA) #set these to NA for population-level predictions


#predict using the model
preds <- predict(TMB_m, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

#prepare for plotting
y_axis_var <- c("TRI_100", "ridge_100")
x_axis_var <- "weeks_since_emig"
labels <- data.frame(vars = y_axis_var, 
                     label = c("Terrain Ruggedness Index", "Distance to ridge (km)"))

for (i in y_axis_var){
  
  label <- labels %>% 
    filter(vars == i) %>% 
    pull(label)
  
  #interaction to be plotted
  interaction_term <- paste0("wk_", i)
  pred_r <- preds_pr %>% 
    filter(interaction == interaction_term) %>%  #only keep the rows that contain data for this interaction term
    dplyr::select(c(which(names(.) %in% c(x_axis_var, i)), "probs")) %>% 
    terra::rast(type = "xyz") %>%
    #focal(w = 3, fun = mean, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) #%>% 
  #rename(probs = focal_mean)
  
  #X11(width = 6.9, height = 3.5)
  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = probs)) +
    scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                         na.value = "white", name = "Intensity of use")+
    #guides(fill = guide_colourbar(title.position = "top")) +
    guides(fill = guide_colourbar(title.vjust = .95)) + #the legend title needs to move up a bit
    labs(x = "Week since dispersal", y = label) + #add label for GRC plot
    theme_minimal() +
    theme(plot.margin = margin(0, 15, 0, 0, "pt"),
          legend.direction="horizontal",
          legend.position = "bottom",
          legend.key.width=unit(.7,"cm"),
          legend.key.height=unit(.25,"cm"),
          text = element_text(size = 8), #font size should be between 6-8
          axis.title.x = element_text(hjust = 1, margin = margin(t=6)), #align the axis labels
          axis.title.y = element_text(angle = 90, hjust = 1, margin=margin(r=6))) 
  
  assign(paste0(i, "_p"), pred_p)
}


#plot all interaction plots together
X11(width = 4.5, height = 4.8)
combined <- TRI_100_p + ridge_100_p & theme(legend.position = "bottom")
p <- combined + plot_layout( guides = "collect", nrow = 2)

ggsave(p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/interactions.pdf", 
       width = 4.5, height = 4.6, dpi = 300)

# STEP 6: predictions for the Alps + PLOTS ------------------------------------------------------------------ 

#read in the ssf model
TMB_M <- readRDS( "TMB_model.rds")

#Alpine df prepared in 03_04_clogit_workflow.R
topo_df <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/topo_df_100_LF.rds")

#prepare dist to ridge to be used as the base layer for plotting
ridge_100 <- rast("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/ridge_100_LF.tif") %>% 
  as.data.frame(xy = T) %>% 
  drop_na(distance_to_ridge_line_mask)
ridge_100[ridge_100$distance_to_ridge_line_mask > 5000, "distance_to_ridge_line_mask"] <- NA #do this for the plot to look nicer... NA values will be white

#areas_df <- data.frame(pixels = numeric(156),
#                       area_m2 = numeric(156),
#                       area_km2 = numeric(156),
#                       week_since_dispersal = integer(156))

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/tmb_figs/alpine_preds/")

wks_ls <- split(data, data$weeks_since_emig)

gc()

(start <- Sys.time())
for(x in wks_ls[c(75:156)]){
  
  #extract week number
  one_week <- x %>% 
    distinct(weeks_since_emig) %>% 
    pull(weeks_since_emig)
  
  #create week ID to be used for naming files and plots
  week_i <- one_week %>% 
    str_pad(3,"left","0")
  
  #generate a new dataset
  topo_df <- topo_df %>%
    mutate(step_length = mean(x$step_length),
           step_length_z = 0,
           weeks_since_emig = one_week,
           weeks_since_emig_z = (one_week - attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:center'))/attr(data[,colnames(data) == "weeks_since_emig_z"],'scaled:scale'),
           stratum_ID = NA,
           animal_ID = sample(x$animal_ID, nrow(topo_df), replace = T)) #set the grouping variables to NA as per glmm.TMB recommendation
  
  #predict using the model
  #(b <- Sys.time())
  topo_df$preds <- predict(TMB_M, newdata = topo_df, type = "link")
  #Sys.time() - b #24 min for one week
  
  topo_df$probs <- gtools::inv.logit(topo_df$preds)
  
  gc()
  
  #skip saving this file. it's 1 GB. just get out of it what I need
  #saveRDS(preds_pr, paste0("/home/enourani/Documents/alpine_preds/alpine_preds_wk_", week_i, ".rds"))
  
  #calculate suitable areas
  area_.6 <- topo_df %>% 
    filter(probs >= .6) %>% 
    summarize(pixels = n()) %>% #count the 
    mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
           area_km2 = round(area_m2/1e6,3),
           week_since_dispersal = week_i)
  
  #save area as a file
  save(area_.6, file = paste0("area_alps_wk_", week_i, ".r"))
  
  
  #plot the raw map for a few weeks
  if(one_week %in% c(1,24,52,104,156)){
   # tiff(paste0("alps_wk_", week_i, ".tiff"), width = 7, height = 5, units = "in", res = 400)
    t <- topo_df %>% 
            ggplot() +
            geom_tile(aes(x = location.long, y = location.lat, fill = probs)) +
            scale_fill_gradientn(colours = oce::oceColorsPalette(100), limits = c(0,1),
                                 na.value = "white", name = "Intensity of use") +
            labs(x = "", y = "", title = paste0("Week ", one_week, " since dispersal")) +
            theme_void()
    
    ggsave(plot = t, filename = paste0("alps_wk_", week_i, ".tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400)
   # dev.off()
    rm(t)
    gc()
  }

  
  #plot
  #tiff(paste0(week_i, "_alpine_pred.tiff"), width = 7, height = 5, units = "in", res = 400)
  #p <- ggplot() +
  p <- ggplot() +  
  geom_tile(data = ridge_100, aes(x = x, y = y, fill = scale(distance_to_ridge_line_mask))) +
    scale_fill_gradientn(colors = grey.colors(100), guide = "none", na.value = "white") +
    new_scale_fill() +
    stat_density_2d(data = topo_df %>% filter(probs >= .6) %>% dplyr::select("location.lat", "location.long", "probs"), 
                    aes(x = location.long, y = location.lat, fill = after_stat(level)), geom = "polygon") +
    scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(100)[51:100], alpha = .2)) +
    labs(title = paste0("Week ", one_week, " since dispersal"), x = "", y = "") +
    theme_void()
  
  #dev.off()
  ggsave(plot = p, filename = paste0(week_i, "_alpine_pred.tiff"), device = "tiff", width = 7, height = 4.5, dpi = 400) #3min
  
  rm(p)
  gc()
  print(paste0("week ", week_i, " of 156 done!"))

}

Sys.time() - start #30 min per week


# ## add areas of weeks 1 - 41 to the dataframe
# files <- list.files("/home/enourani/Documents/alpine_preds", full.names = T, pattern = ".rds")[1:41]
# 
# areas_1_41 <- data.frame()
# 
# for (i in files){
#   preds_pr <- readRDS(i)
#   
#   week_i <- str_sub(i, -7, -5)
#   
#   area_.7 <- preds_pr %>%
#     filter(probs >= 0.7) %>% 
#     summarize(pixels = n()) %>% #count the 
#     mutate(area_m2 = pixels * 100 * 100, #the resolution of the cell size
#            area_km2 = round(area_m2/1e6,3),
#            week_since_dispersal = week_i)
#   
#   areas_1_41 <- rbind(areas_1_41, area_.7)
#   
# }
# 
# areas_df <- areas %>% 
#   reduce(rbind) %>% 
#   bind_rows(areas_1_41,.) %>% 
#   mutate(week_since_dispersal = as.numeric(week_since_dispersal))
# 
# saveRDS(areas_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/xxx")
# 



# STEP 9: xxx ------------------------------------------------------------------ 
