# create figures for INLA model for landscpe level changes in GE flight. follows up on energy_landscape_modeling_method1.R and /home/mahle68/ownCloud/Work/cluster_computing/GE_inla_static/
#Jul 19. 2022. Konstanz, DE
#Elham Nourani


#load packages
library(tidyverse)
library(terra)
library(sf)
library(scales)


rnd <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/rnd_coeff_M_preds.rds")


# PLOT 1: coefficient plots ----------------------------------------------------------------------------------------------------

graph <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/results/graph_M_pred.rds")

#remove weeks since dispersal
graph <- graph[graph$Factor != "weeks_since_emig_n_z",]
droplevels(graph$Factor)
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames




min <- min(graph$Lower,na.rm = T)
max <- max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

#remove weeks since emigration. because it is NA
#plot
X11(width = 4.7, height = 2.7)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/static_landscape_figs/coeffs_48n_wpred.png", 
    width = 4.7, height = 2.7, units = "in", res = 300)

par(cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 8.5, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-0.1,0.5), ylim = c(0.7,4.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph$Estimate, graph$Factor_n, col = "cornflowerblue", pch = 20, cex = 2)
arrows(graph$Lower, graph$Factor_n,
       graph$Upper, graph$Factor_n,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line


#add axes
axis(side= 1, at = seq(-0.1, 0.5, by =  0.1), labels = seq(-0.1, 0.5, by =  0.1), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at = c(1:4),
     labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM"),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

dev.off()

# PLOT 2: interaction plots ----------------------------------------------------------------------------------------------------

all_data <- readRDS("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/alt_50_20_min_48_ind_static_inlaready_wks.rds")
new_data <- readRDS("/home/enourani/ownCloud/Work/cluster_computing/GE_inla_static/alt_50_20_min_48_ind_static_inlaready_wmissing_wks.rds")
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

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_2vars.png", 
       width = 10, height = 4, dpi = 300)