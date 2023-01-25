# create figures for INLA model for landscape level changes in GE flight. follows up on energy_landscape_modeling_method1.R and /home/mahle68/ownCloud/Work/cluster_computing/GE_inla_static/
#Jul 19. 2022. Konstanz, DE
#Elham Nourani


#load packages
library(tidyverse)
library(terra)
library(sf)
library(scales)
library(patchwork) #patching up interaction plots

setwd("/home/mahle68/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")


# PLOT 1: coefficient plots ----------------------------------------------------------------------------------------------------

graph <- readRDS("/home/mahle68/ownCloud/Work/cluster_computing/GE_inla_static/results/graph_M_pred.rds")

#remove weeks since dispersal
graph <- graph[graph$Factor != "weeks_since_emig_n_z",]
droplevels(graph$Factor)
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#min <- min(graph$Lower,na.rm = T)
#max <- max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

#plot in ggplot2
X11(width = 4.7, height = 2.7)

coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray", size = 0.5) +
  geom_point(color = "cornflowerblue", size = 2)  +
  xlim(-0.1,0.5) +
  scale_y_discrete(name = "",
                   labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM")) +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", size = 1) +
  theme_classic()
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_coefs_static_w_rnd_wk_48minimal.png", 
       width = 4.7, height = 2.7, dpi = 300)




# PLOT 2: interaction plots ----------------------------------------------------------------------------------------------------

all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_inlaready_wks.rds")
#new_data <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/alt_50_20_min_48_ind_static_inlaready_wmissing_wks.rds")
preds <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/preds_M_preds.rds")

y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

labels <- data.frame(var = c("dem_100_z", "TRI_100_z","weeks_since_emig_n_z"),
                     label = c("Altitude (m.asl)", "Terrain ruggedness index", "Weeks since dispersal"))


#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

for (i in y_axis_var){
  
  #fitted_df <- preds %>% 
  #  mutate(new_data[is.na(new_data$used), i])
  #extract scale and center values needed to back transform the scaled values
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
  
  saveRDS(avg_pred, file = paste0("inla_pred_w_random_n48_df_", i,".rds"))
  
  #create a raster
  #r <- rast(avg_pred[,c("x", "y", "avg_pres")], crs = "+proj=longlat +datum=WGS84") #use the ggplot code at the end to plot it. reorder the columns
  
  
  #saveRDS(r, file = paste0("inla_pred_w_random_n48_", i,".rds"))
  
}

#######use ggplot

r_dem <- readRDS("inla_pred_w_random_n48_df_dem_100_z.rds")
r_rug <- readRDS("inla_pred_w_random_n48_df_TRI_100_z.rds")

X11(width = 6, height = 4)

p_rugg <- ggplot(data = r_rug) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  #scale_fill_gradientn(colours = oce::oceColorsPalette(80), limits = c(0,1), 
  scale_fill_gradient2(low = "lightslateblue", mid = "white", high = "firebrick1",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  theme_minimal() +
  labs(x = "Weeks since dispersal", y = "Ruggedness")

ggsave(plot = p_rugg, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_48n_",i, ".png"), 
       width = 6, height = 4, dpi = 300)



X11(width = 6, height = 4)

p_dem <- ggplot(data = r_dem) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  #scale_fill_gradientn(colours = oce::oceColorsPalette(80), limits = c(0,1), 
  scale_fill_gradient2(low = "lightslateblue", mid = "white", high = "firebrick1",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  theme_minimal() +
  labs(x = "Weeks since dispersal", y = "Elevation (m)")

ggsave(plot = p_dem, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_48n_",i, ".png"), 
       width = 6, height = 4, dpi = 300)

#put both plots in one device
X11(width = 10, height = 4)
combined <- p_dem + p_rugg & theme(legend.position = "right")
p_2 <- combined + plot_layout(guides = "collect")

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_n48.png", 
       width = 10, height = 4, dpi = 300)


# PLOT 3: individual variation plots ----------------------------------------------------------------------------------------------------

rnd <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/rnd_coeff_M_preds.rds")

#!!!!!!!!make sure to add the coefficient to these.


names <- rnd[[2]]$ID

dem <- rnd[[2]] %>% 
  mutate(variable = "DEM")

two_vars <- rnd[[3]] %>% 
  mutate(variable = "TRI") %>% 
  bind_rows(dem)

cols <- c(DEM = "lightcoral",
          TRI = "cornflowerblue")

#plot two_vars
X11(width = 7, height = 9)
(coefs_inds <- ggplot(two_vars, aes(x = mean, y = ID, color = variable)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray25", size = 0.5) +
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = two_vars[, 4], xmax = two_vars[, 6]), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = cols) + 
    theme_minimal()) +
  labs(x = "Estimate", y = "")


ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/ind_coefs_static_w_rnd_wk_48.png", 
       width = 7.5, height = 10, dpi = 300)


