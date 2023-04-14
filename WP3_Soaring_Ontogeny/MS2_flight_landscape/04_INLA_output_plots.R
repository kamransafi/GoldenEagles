# create figures for INLA model for landscape level changes in GE flight. follows up on energy_landscape_modeling_method1.R and /home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/
#Jul 19. 2022. Konstanz, DE
#Elham Nourani

#load packages
library(tidyverse)
library(terra)
library(sf)
library(scales)
library(patchwork) #patching up interaction plots

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

# PLOT 1: coefficient plots ----------------------------------------------------------------------------------------------------

graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Apr23_seasonality_100/graph_M_main100_hrly.rds")

#remove weeks since dispersal
graph <- graph[graph$Factor != "weeks_since_emig_z",]
#droplevels(graph$Factor)
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

graph$Factor_n <- as.numeric(graph$Factor)

#plot in ggplot2... later on reorder the variables and make the names more cohesive
X11(width = 8, height = 6)

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", size = 0.5) +
  geom_point(color = "cornflowerblue", size = 2)  +
  #xlim(-0.1,0.6) +
  #scale_y_discrete(name = "",
  #                 labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM")) +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", size = 1) +
  theme_classic()
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_coefs_static100_sl.png", 
       width = 4.7, height = 2.7, dpi = 300)


# PLOT 2: interaction plots ----------------------------------------------------------------------------------------------------

all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_100_daytemp_inlaready_wks.rds")
preds <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Mar23_temp_100/preds_M_main100.rds")

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

data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/all_inds_annotated_static_apr23.rds")
rnd <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Apr23_seasonality_100/rnd_coeff_M_main100_hrly.rds")
graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Apr23_seasonality_100/graph_M_main100_hrly.rds")



#plot the effect of the grouped effect

#Table with summary of random effects; ID is for the unique individuals
tab_dem <-  rnd$ind1
tab_tri <-  rnd$ind2
tab_slope_tpi <-  rnd$ind3

ind_IDs <- unique(data$ind1)

dem_coefs <- tab_dem %>% 
  mutate(ind_ID = rep(ID[1:length(ind_IDs)], 12), #the tab datasets only have the ind id for the first season and then just numbers
         #month = factor(rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), each = length(ind_IDs)),
         month = factor(rep(1:12, each = length(ind_IDs)),
                        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

X11(width = 10, height = 10)
(dem_inds <- ggplot(dem_coefs, aes(x = mean, y = ind_ID)) +
    geom_vline(xintercept = 0, linetype="dashed", size = 0.5) +
    geom_linerange(aes(xmin = dem_coefs$'0.025quant', xmax = dem_coefs$'0.975quant'), size = 0.6, color = "#a9c4f5") +
    geom_point(size = 1.5, color =  "#6495ed") +
    facet_wrap(vars(month), ncol = 6) +
    labs(x = "Estimate", y = "") +
    ggtitle("Elevation") +
    theme_classic() +
    theme(axis.text = element_text(size = 6.5)))

ggsave(plot = dem_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_season_coefs_dem.png", 
       width = 10, height = 10, dpi = 300)


tri_coefs <- tab_tri %>% 
  mutate(ind_ID = rep(ID[1:length(ind_IDs)], 12), #the tab datasets only have the ind id for the first season and then just numbers
         #month = factor(rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), each = length(ind_IDs)),
         month = factor(rep(1:12, each = length(ind_IDs)),
                        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

(tri_inds <- ggplot(tri_coefs, aes(x = mean, y = ind_ID)) +
    geom_vline(xintercept = 0, linetype="dashed", size = 0.5) +
    geom_linerange(aes(xmin = tri_coefs$'0.025quant', xmax = tri_coefs$'0.975quant'), size = 0.6, color = "#a9c4f5") +
    geom_point(size = 1.5, color =  "#6495ed") +
    facet_wrap(vars(month), ncol = 6) +
    labs(x = "Estimate", y = "") +
    ggtitle("Terrain Ruggedness Index") +
    theme_classic() +
    theme(axis.text = element_text(size = 6.5)))

ggsave(plot = tri_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_season_coefs_tri.png", 
       width = 10, height = 10, dpi = 300)


slope_tpi_coefs <- tab_slope_tpi %>% 
  mutate(ind_ID = rep(ID[1:length(ind_IDs)], 12), #the tab datasets only have the ind id for the first season and then just numbers
  #month = factor(rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), each = length(ind_IDs)),
  month = factor(rep(1:12, each = length(ind_IDs)),
                 labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

(tpi_inds <- ggplot(slope_tpi_coefs, aes(x = mean, y = ind_ID)) +
    geom_vline(xintercept = 0, linetype="dashed", size = 0.5) +
    geom_linerange(aes(xmin = slope_tpi_coefs$'0.025quant', xmax = slope_tpi_coefs$'0.975quant'), size = 0.6, color = "#a9c4f5") +
    geom_point(size = 1.5, color =  "#6495ed") +
    facet_wrap(vars(month), ncol = 6) +
    labs(x = "Estimate", y = "") +
    ggtitle("Slope Variation Index") +
    theme_classic() +
    theme(axis.text = element_text(size = 6.5)))

ggsave(plot = tpi_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_season_coefs_slope_tpi.png", 
       width = 10, height = 10, dpi = 300)

