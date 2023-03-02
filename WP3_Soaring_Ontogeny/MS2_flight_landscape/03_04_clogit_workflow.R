#script for constructing the energy landscape for golden eagles
#trial with clogit. INLA predictions are all close to 0
#follows from 03_energy_landscape_modeling_method1.R and 04_INLA_output_plots.R
#02.03.2023. Konstanz, DE.
#Elham Nourani, PhD.


library(tidyverse)
library(magrittr)
library(corrr)
library(jtools) #to plot coeffs of clogit
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot)
library(terra)
library(gstat) #for interpolations
library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots
library(patchwork) #patching up interaction plots


setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

all_data <- readRDS("alt_50_20_min_48_ind_static_100_daytemp_inlaready_wks.rds") #this has the limit on TRI range

# STEP 1: ssf modeling ----------------------------------------------------------------

## mid_step: investigate using clogit
# control for monthly temperature....

form1a <- used ~ dem_100_z * weeks_since_emig_n_z *t2m_z + 
  TRI_100_z * weeks_since_emig_n_z + 
  strata(stratum)

ssf <- clogit(form1a, data = all_data)
summary(ssf)

X11(width = 14, height = 12)
coefs <- plot_summs(ssf, point.size = 7, point.shape = 16) +
  labs(y = NULL) +
  labs(x = NULL) +
  scale_y_discrete( labels = c("Elevation:Weeks since dispersal:Temperature", "Weeks since dispersal:Ruggedness", "Week since dispersal: Temperature","Elevation:Temperature",
                               "Elevation:Weeks since dispersal","Ruggedness",  "Temperature",  "Elevation")) +
  theme_classic() +
  theme(axis.text.y = element_text(face="bold", size = 20),
                axis.text.x = element_text(face="bold", size = 20))
  

ggsave(plot = coefs, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/clogit_coeffs.png", 
       width = 14, height = 12, dpi = 500)


# STEP 2: predict using the ssf  ----------------------------------------------------------------

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations
grd_dem <- expand.grid(x = seq(from = min(all_data$weeks_since_emig_n_z), to = max(all_data$weeks_since_emig_n_z), by = 0.15),
                   y = seq(from = min(all_data$dem_100_z), to = max(all_data$dem_100_z), by = 0.15))

grd_tri <- expand.grid(x = seq(from = min(all_data$weeks_since_emig_n_z), to = max(all_data$weeks_since_emig_n_z), by = 0.15),
                       y = seq(from = min(all_data$TRI_100_z), to = max(all_data$TRI_100_z), by = 0.15))

grd_temp <- expand.grid(x = seq(from = min(all_data$t2m_z), to = max(all_data$t2m_z), by = 0.15),
                       y = seq(from = min(all_data$dem_100_z), to = max(all_data$dem_100_z), by = 0.15))

n <- nrow(grd_dem) #the dem grid has more number of rows, so use it

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 20), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 20), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 20), n, replace = T),
         t2m_z = sample(seq(min(all_data$t2m_z),max(all_data$t2m_z), length.out = 20), n, replace = T))


#predict using the model
preds <- predict(ssf, newdata = new_data, type = "risk")
preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise %>% 
  mutate(probs = preds/(preds+1))

#prepare for plotting
y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

for (i in y_axis_var){
  
  i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
  i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')
  
  #summarize values, so each (x,y) combo has one probability value
  avg_pred <- preds_pr %>% 
    group_by_at(c(which(names(preds_pr) == i),which(names(preds_pr) == x_axis_var))) %>%  #group by weeks since emigration and i
    summarise(avg_pres = mean(probs)) %>% 
    ungroup() %>% 
    mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
           dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
    rename(x = which(names(.) == x_axis_var), #weeks since emig
           y = which(names(.) == i)) %>% #y axis variables
    dplyr::select(c("x","y","avg_pres")) %>%  #reorder the columns for making xyz raster
    as.data.frame()
  
  #saveRDS(avg_pred, file = paste0("inla_pred_clogit_", i,".rds"))
  
  r <- avg_pred %>% 
    rast(type = "xyz") %>% 
    focal(w = 7, fun = median, na.policy = "only", na.rm = T) %>% 
    as.data.frame(xy = T) %>%
    rename(avg_pres = focal_median)
  
  saveRDS(r, file = paste0("inla_pred_clogit_", i,".rds"))

}


#create plots
p_dem <- readRDS("inla_pred_clogit_dem_100_z.rds") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "", y = "Elevation \n (m)") +
  theme_classic()

p_rugg <- readRDS("inla_pred_clogit_TRI_100_z.rds") %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                       na.value = "white", name = "Intensity of use") +
  labs(x = "Weeks since dispersal", y = "Terrain Ruggedness \n Index") +
  theme_classic()


#put both plots in one device
X11(width = 9, height = 4)
combined <- p_dem + p_rugg & theme(legend.position = "right")
(p_2  <- combined + plot_layout(guides = "collect", nrow = 2))




#interpolate. for visualization purposes
dem_r <- readRDS("inla_pred_clogit_dem_100_z.rds")
dem_rast <- dem_r %>% 
  rast(type = "xyz") 

surf.1 <- Tps(as.matrix(dem_r[, c(1,2)], col = 2), dem_r[,3])

grd <- expand.grid(x = seq(from = ext(dem_rast)[1],to = ext(dem_rast)[2],by = 2),
                   y = seq(from = ext(dem_rast)[3],to = ext(dem_rast)[4],by = 2))

grd$coords <- matrix(c(grd$x,grd$y),ncol=2)

surf.1.pred <- predict.Krig(surf.1,grd$coords)
interpdf <- data.frame(grd$coords, surf.1.pred)

colnames(interpdf) <- c("weeks","dem","prob_pres")

coordinates(interpdf) <- ~ weeks + dem
gridded(interpdf) <- TRUE
interpr <- raster(interpdf)
