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

graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Feb23_no_temp_200m/graph_M_main200.rds")

#remove weeks since dispersal
graph <- graph[graph$Factor != "weeks_since_emig_n_z",]
#droplevels(graph$Factor)
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

graph$Factor_n <- as.numeric(graph$Factor)

#plot in ggplot2
X11(width = 4.7, height = 2.7)

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", size = 0.5) +
  geom_point(color = "cornflowerblue", size = 2)  +
  xlim(-0.1,0.6) +
  scale_y_discrete(name = "",
                   labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM")) +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", size = 1) +
  theme_classic()
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_coefs_static200_FEB23.png", 
       width = 4.7, height = 2.7, dpi = 300)


# PLOT 2: interaction plots ----------------------------------------------------------------------------------------------------

all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_200_inlaready_wks.rds")
preds <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Feb23_no_temp_200m/preds_M_main200.rds")

y_axis_var <- c("dem_200_z", "TRI_200_z")
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
  
  saveRDS(avg_pred, file = paste0("inla_pred_FEB23_200m_", i,".rds"))
  
}

#######use ggplot

#Fill in the NAs to get rid of white spaces without altering the existing values
#create rasters
r_dem <- readRDS("inla_pred_FEB23_200m_dem_200_z.rds") %>% 
  dplyr::select(x,y,avg_pres) %>% 
  rast(type = "xyz") %>% 
  focal(w = 7, fun = mean, na.policy = "only", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(avg_pres = focal_mean)
  
r_rug <- readRDS("inla_pred_FEB23_200m_TRI_200_z.rds") %>% 
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

rnd <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Feb23_no_temp_200m/rnd_coeff_M_main200.rds")

graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/Feb23_no_temp_200m/graph_M_main200.rds")

#!!!!!!!!make sure to add the coefficient to these.
names <- rnd[[2]]$ID

dem <- rnd[[2]] %>% 
  mutate(coef = mean + graph %>% filter(Factor == "dem_200_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "dem_200_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "dem_200_z") %>% pull(Upper),
         variable = "DEM")

two_vars <- rnd[[3]] %>% 
  mutate(coef = mean + graph %>% filter(Factor == "TRI_200_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "TRI_200_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "TRI_200_z") %>% pull(Upper),
         variable = "TRI") %>% 
  bind_rows(dem)

cols <- c(DEM = "lightcoral",
          TRI = "cornflowerblue")

#plot two_vars
X11(width = 7, height = 9)
(coefs_inds <- ggplot(two_vars, aes(x = coef, y = ID, color = variable)) +
    geom_vline(xintercept = graph %>% filter(Factor == "dem_200_z") %>% pull(Estimate), linetype="dashed", 
               color = "lightcoral", size = 0.5) +
    geom_vline(xintercept = graph %>% filter(Factor == "TRI_200_z") %>% pull(Estimate), linetype="dashed", 
               color = "cornflowerblue", size = 0.5) +
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = cols) + 
    labs(x = "Estimate", y = "") +
    theme_classic()) 

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_coefs_FEB23_200m.png", 
       width = 7.5, height = 10, dpi = 300)


