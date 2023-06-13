# create figures for INLA model for landscape level changes in GE flight. follows up on energy_landscape_modeling_method1.R and /home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/
#Jul 19. 2022. Konstanz, DE
#Elham Nourani

#load packages
library(tidyverse)
library(terra)
library(sf)
library(scales)
library(patchwork) #patching up interaction plots

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

result_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/cluster_computing/GE_inla_static/results/Jun1_2023/"

# PLOT 1: coefficient plots ----------------------------------------------------------------------------------------------------

graph <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/cluster_computing/GE_inla_static/results/Jun1_2023/graph_M_main100_hrly_nopred.rds")

#remove weeks since dispersal
graph <- graph[graph$Factor != "weeks_since_emig_z",]
#droplevels(graph$Factor)
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

graph$Factor_n <- as.numeric(graph$Factor)

#plot
X11(width = 8, height = 6)

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", linewidth = 0.5) +
  geom_point(color = "cornflowerblue", size = 2)  +
  labs(x = "Estimate", y = "") +
  #scale_y_discrete(labels = rev(c("TRI", "Step length", "Distance to ridge", "TRI: Step length", "TRI: Week",
  #                                "Step length: Week", "Step length: Distance to ridge", "Distance to ridge: Week", "TRI: Step length: Week", "Distance to ridge: Step length: Week"))) +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", linewidth = 1) +
  theme_classic()
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_coefs_May23.png", 
       width = 8, height = 6, dpi = 300)


# PLOT 2: interaction plots ----------------------------------------------------------------------------------------------------

preds <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/May1223_slrm/preds_M_main100_hrly_try2.rds")


#prepare for plotting
y_axis_var <- c("TRI_100", "ridge_100", "step_length")
x_axis_var <- "weeks_since_emig"

labels <- data.frame(vars = y_axis_var, 
                     label = c("Terrain Ruggedness Index", "Distance to ridgeline (km)", "Step length (m)"))

#older versions of this script contained backtransforming the z-scores in a for loop.

for (i in y_axis_var){
  
  label <- labels %>% 
    filter(vars == i) %>% 
    pull(label)
  
  #interaction to be plotted
  interaction_term <- paste0("wk_", i)
  
  #create a raster and run a moving window to make the pattern easier to see.
  pred_r <- preds %>% 
    filter(interaction == interaction_term) %>%  #only keep the rows that contain data for this interaction term
    dplyr::select(c(which(names(.) %in% c(x_axis_var, i)), "prob_pres")) %>% 
    terra::rast(type = "xyz") %>%
    focal(w = 3, fun = median, na.policy = "all", na.rm = T) %>%
    as.data.frame(xy = T) %>% 
    rename(prob_pres = focal_median)
  
  #plot
  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = prob_pres)) +
    scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                         na.value = "white", name = "Intensity of use") +
    labs(x = "", y = i) +
    theme_classic()
  
  #save the plot
  ggsave(plot = pred_p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_", interaction_term,".png"), 
         width = 9, height = 4, dpi = 500)
}


#old:
y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

for (i in y_axis_var){
  
  i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
  i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')
  
  #summarize values, so each (x,y) combo has one probability value
  avg_pred <- preds %>% 
    group_by_at(c(which(names(preds) == i),which(names(preds) == x_axis_var))) %>%  #group by weeks since emigration and i
    summarise(avg_pres = mean(prob_pres)) %>% 
    ungroup() %>% 
    mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
           dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
    rename(y = which(names(.) == i), #y axis variables
           x = which(names(.) == x_axis_var)) %>% #weeks since emig
    as.data.frame()
  
  saveRDS(avg_pred, file = paste0("inla_pred_Mar23_100_", i,".rds"))
  
}



#######use ggplot

#Fill in the NAs to get rid of white spaces without altering the existing values
#create rasters
r_dem <- readRDS("inla_pred_Mar23_100_dem_100_z.rds") %>% 
  dplyr::select(x,y,avg_pres) %>% 
  rast(type = "xyz") %>% 
  focal(w = 7, fun = mean, na.policy = "only", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(avg_pres = focal_mean)
  
r_rug <- readRDS("inla_pred_Mar23_100_TRI_100_z.rds") %>% 
  dplyr::select(x,y,avg_pres) %>% 
  rast(type = "xyz") %>% 
  focal(w = 7, fun = mean, na.policy = "only", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(avg_pres = focal_mean)

#create plots
p_dem <- ggplot(data = r_dem) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "lightslateblue", mid = "seashell2", high = "firebrick1",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "", y = "Elevation \n (m)") +
  theme_classic()

p_rugg <- ggplot(data = r_rug) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "lightslateblue", mid = "seashell2", high = "firebrick1",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "Weeks since dispersal", y = "Terrain Ruggedness \n Index") +
  theme_classic()


#put both plots in one device
X11(width = 9, height = 4)
combined <- p_dem + p_rugg & theme(legend.position = "right")
(p_2  <- combined + plot_layout(guides = "collect", nrow = 2))

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_interaction_FEB23_200m_NAsfilled.png", 
       width = 9, height = 4, dpi = 300)
#the plot ws so much more clear for the 100m res :( consider trying with 100m again! tile the alpine region for predictions!

# PLOT 3: individual variation plots ----------------------------------------------------------------------------------------------------

data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/new_data_ssf_inla_preds.rds")
rnd <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/May1223_slrm/rnd_coeff_M_main100_hrly_try2.rds")
graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/May1223_slrm/graph_M_main100_hrly_try2.rds")


#the ind names are messed up. fix it using the emigration date dataframe
emig_dates <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/fleding_emigration_timing_Mar2023.rds")


#plot the effect of the grouped effect

#Table with summary of random effects; ID is for the unique individuals
tab_tri <- rnd$ind2
tab_ridge <- rnd$ind3

ind_IDs <- unique(data$ind1)

tri_vars <- rnd[[2]] %>% 
  mutate(coef = mean + graph %>% filter(Factor == "TRI_100_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "TRI_100_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "TRI_100_z") %>% pull(Upper),
         variable = "TRI")


two_vars <- rnd[[3]] %>% 
  mutate(coef = mean + graph %>% filter(Factor == "ridge_100_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "ridge_100_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "ridge_100_z") %>% pull(Upper),
         variable = "Distance_to_ridge") %>% 
  bind_rows(tri_vars)

cols <- c(TRI = "lightcoral",
          Distance_to_ridge = "cornflowerblue")

X11(width = 7, height = 9)
(coefs_inds <- ggplot(two_vars, aes(x = coef, y = ID, color = variable)) +
    geom_vline(xintercept = graph %>% filter(Factor == "TRI_100_z") %>% pull(Estimate), linetype="dashed", 
               color = "lightcoral", size = 0.5) +
    geom_vline(xintercept = graph %>% filter(Factor == "ridge_100_z") %>% pull(Estimate), linetype="dashed", 
               color = "cornflowerblue", size = 0.5) +
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = cols) + 
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "") +
    theme_classic()) 

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_coefs_inds_May23.png", 
       width = 7.5, height = 10, dpi = 300)
