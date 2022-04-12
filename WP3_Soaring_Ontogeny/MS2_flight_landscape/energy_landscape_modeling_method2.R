#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#This version tries the method of fitting one model per individual. using the iSSA framework (selection-free movement kernel + habitat selection kernel)
#Apr 11. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot) #to plot coeffs of INLA
library(jtools) #to plot coeffs of clogit

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")

#open annotated data (static variables and time since fledging and emigration)
load("alt_50_20_min_25_ind_static_time_ann.RData") #cmpl_ann

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)), #round up
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_")) #forgot to include stratum id in the previous code ) %>%  


Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# STEP 1: data exploration ----------------------------------------------------------------

#number of days available for each individual
cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% #round up days since emigration. we don't want zeros
  ggplot(aes(x = individual.local.identifier, y = max_day)) +
  geom_col()

cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_day = max(ceiling(as.numeric(days_since_emig)))) %>% 
  summarize(mode = Mode(max_day), # 537
            mean = mean(max_day)) #444


#terrain ~ days since emigration ..... no patterns in the plots

plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "dem_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "slope_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "aspect_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "slope_TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "TRI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("days_since_emig", "aspect_TPI_100")])


ggplot(cmpl_ann[cmpl_ann$used == 1,], aes(as.numeric(days_since_emig), aspect_TPI_100)) +
  geom_point() +
  stat_smooth(aes(group = individual.local.identifier), method = "lm") +
  theme_minimal() +
  theme(legend.position = "none")


# STEP 2: summary plots ----------------------------------------------------------------

#one set of boxplots for select days: 10, 50 , 100, 500

data_int <- cmpl_ann %>%
  filter(days_since_emig_n %in% c(1, 10, 30 , 100, 300)) %>% 
  mutate(days_f = as.factor(days_since_emig_n))


variables <- c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100" )

X11(width = 6, height = 9)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/box_plots.png", width = 6, height = 9, units = "in", res = 300)

par(mfrow= c(4,2), 
    oma = c(0,0,3,0), 
    las = 1,
    mgp = c(0,1,0))
for(i in 1:length(variables)){
  boxplot(data_int[,variables[i]] ~ data_int[,"days_f"], data = data_int, boxfill = NA, border = NA, main = variables[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(data_int[data_int$used == 1, variables[i]] ~ data_int[data_int$used == 1,"days_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9),  lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1, "days_f"])) - 0.15)
  boxplot(data_int[data_int$used == 0, variables[i]] ~ data_int[data_int$used == 0, "days_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1 , "days_f"])) + 0.15)
  
}
dev.off()


#also consider making the plots for periods of time and not only one day

# STEP 3: check for collinearity ----------------------------------------------------------------

cmpl_ann %>% 
  dplyr::select(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "days_since_emig_n")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #slope and TRI are correlated (.98)

# STEP 4: standardize variables ----------------------------------------------------------------

all_data <- cmpl_ann %>% 
  mutate_at(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "days_since_emig_n"),
            list(z = ~(scale(.)))) 

# STEP 5: ssf modeling ----------------------------------------------------------------

#select one individual to get the workflow straight: "Almen19 (eobs 7001)" or "Droslöng17 (eobs 5704)"

one_ind <- all_data %>% 
  filter(individual.local.identifier == "Droslöng17 (eobs 5704)")


#model
#slope on its own is a percentage, dont use it. use unevenness in the slope. aspect on its own is also weird. use unevenness in aspect. treat days since as a continous variable
form1a <- used ~  cos(turning_angle) * days_since_emig_n_z + 
  log(step_length) * days_since_emig_n_z + 
  dem_100_z * days_since_emig_n_z + 
  slope_TPI_100_z * days_since_emig_n_z + 
  aspect_TPI_100_z * days_since_emig_n_z  + 
  TRI_100_z * days_since_emig_n_z + 
  TPI_100_z * days_since_emig_n_z + 
  strata(stratum)
ssf_full <- clogit(form1a, data = one_ind) #in the clogit model, scaling the days makes the results easier to understand 
summary(ssf_full) 
plot_summs(ssf_full)

ssf_all_inds <- clogit(form1a, data = all_data) 
plot_summs(ssf_all_inds)

#model without interaction term for movement variables
form1aa <- used ~  cos(turning_angle) + log(step_length) + 
  dem_100_z * days_since_emig_n_z + 
  slope_TPI_100_z * days_since_emig_n_z + 
  aspect_TPI_100_z * days_since_emig_n_z  + 
  TRI_100_z * days_since_emig_n_z + 
  TPI_100_z * days_since_emig_n_z + 
  strata(stratum)

ssf_all_inds2 <- clogit(form1aa, data = all_data) 
plot_summs(ssf_all_inds2)

#model with only habitat selection 
form1b <- used ~ dem_100_z * days_since_emig_n_z + 
  slope_TPI_100_z * days_since_emig_n_z + 
  aspect_TPI_100_z * days_since_emig_n_z  + 
  TRI_100_z * days_since_emig_n_z + 
  TPI_100_z * days_since_emig_n_z + 
  strata(stratum)

ssf_HS <- clogit(form1b, data = one_ind) 
plot_summs(ssf_HS)

HS_all_inds <- clogit(form1b, data = all_data) 
plot_summs(HS_all_inds)

plot_summs(ssf_all_inds, ssf_all_inds2, HS_all_inds)
plot_summs(ssf_full,ssf_HS)


#create new data for predictions
n <- 100
new_data <- one_ind %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,  #regular intervals for all variables, so we can make a raster later on
         days_since_emig_n = sample(seq(min(one_ind$days_since_emig_n),max(one_ind$days_since_emig_n), length.out = 20), n, replace = T), 
         days_since_emig_n_z = sample(seq(min(one_ind$days_since_emig_n_z),max(one_ind$days_since_emig_n_z), length.out = 20), n, replace = T), #more levels for time than the other variables
         dem_100_z = sample(seq(min(one_ind$dem_100_z),max(one_ind$dem_100_z), length.out = 20), n, replace = T),
         slope_TPI_100_z = sample(seq(min(one_ind$slope_TPI_100_z),max(one_ind$slope_TPI_100_z), length.out = 20), n, replace = T),
         aspect_TPI_100_z = sample(seq(min(one_ind$aspect_TPI_100_z),max(one_ind$aspect_TPI_100_z), length.out = 20), n, replace = T),
         TRI_100_z = sample(seq(min(one_ind$TRI_100_z),max(one_ind$TRI_100_z), length.out = 20), n, replace = T),
         TPI_100_z = sample(seq(min(one_ind$TPI_100_z),max(one_ind$TPI_100_z), length.out = 20), n, replace = T)) #%>% 
  #full_join(one_ind)

#predictions for only the habitat selection model
preds <- predict(ssf_HS, newdata = new_data, type = "risk")
preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise %>% 
  mutate(probs = preds/(preds+1))

# STEP 6: ssf output plots: the interaction plots ----------------------------------------------------------------
#create a raster of the prediction


y_axis_var <- c("dem_100_z", "slope_TPI_100_z", "aspect_TPI_100_z", "TRI_100_z", "TPI_100_z")
x_axis_var <- "days_since_emig_n_z"

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

lapply(y_axis_var, function(i){
  
#extract scale and center values needed to back transform the scaled values. do this using the whole dataset.
i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')

#summarize values, so each (x,y) combo has one probability value
avg_pred <- preds_pr %>% 
  group_by_at(c(i,x_axis_var)) %>%  #group by days since emigration and i
  summarise(avg_pres = mean(probs)) %>% 
  ungroup() %>% 
  mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
         dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
  #rename(x_backtr = 1, #days since fledging
  #       i_backtr = 2) %>%  #y axis variables
  as.data.frame()

#create a raster
coordinates(avg_pred) <- c(x_axis_var, i)
gridded(avg_pred) <- TRUE
r <- raster(avg_pred)


#interpolate. for visualization purposes
surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2), as.data.frame(r,xy = T)[,3])

grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 2),
                   y = seq(from = extent(r)[3],to = extent(r)[4],by = 2))

grd$coords <- matrix(c(grd$x,grd$y), ncol = 2)

surf.1.pred <- predict.Krig(surf.1,grd$coords)
interpdf <- data.frame(grd$coords, surf.1.pred)

colnames(interpdf) <- c("x_backtr","i_backtr","prob_pres")

coordinates(interpdf) <- ~ x_backtr + i_backtr
gridded(interpdf) <- TRUE
interpr <- raster(interpdf)

#save the raster


})

X11()

spplot(interpr)


#create a color palette
cuts <- c(0, 0.25,0.5,0.75,1) #set breaks
pal <- colorRampPalette(c("aliceblue", "lightskyblue1", "khaki2", "sandybrown", "salmon2","tomato"))
colpal <- pal(100)


#manually determine the range of y axis for each variable
if(i == "dem_100_z"){
  y_axis_r <- c(68,3977)
  y_axis_lab <- c(100, seq(500, 3500, 1000)) #keep all labels at n = 5 to keep everything neat
} else if(i == "slope_100_z"){
  y_axis_r <- c(0,70)
  y_axis_lab <- seq(10, 50, 10)
} else if(i == "aspect_100_z"){
  y_axis_r <- c(4,355)
  y_axis_lab <- c(5, 90, 180, 270, 350)
}


#range of x axis
x_axis_r <- c(1, 829)
#labels of x axis
x_axis_lab <- seq(100,800, 100)


#plot
X11(width = 5, height = 4)

par(cex = 0.7,
    oma = c(1,3.5,1,1),
    mar = c(1, 1, 1, 1.5),
    bty = "n",
    mgp = c(1,0.5,0)
)


raster::plot(interpr, col = colpal, axes = F, box = F, legend = F, ext = extent(c(x_axis_r[1], x_axis_r[2], y_axis_r[1], y_axis_r[2]))) #crop to the extent of observed data

#add axes
axis(side = 1, at = x_axis_lab, #x_axis
     labels = x_axis_lab,
     tick = T , col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck = -.015, line = -2.7, cex.axis = 0.7) #tick marks smaller than default by this proportion

axis(side = 2, at = y_axis_lab, labels = y_axis_lab, 
     tick = T ,col = NA, col.ticks = 1, tck = -.015, las = 2, cex.axis = 0.7)

lines(x = c(-21.9, -21.9), y = c(-9.9,13.9)) #x and y axis lines
#abline(h =-10)

#axis titles
mtext(strsplit(x_axis_var, split = "z")[[1]], 1, line = -4, cex = 0.9, font = 3)
mtext(strsplit(i, split = "z")[[1]], 2, line = 1.2, cex = 0.9)

#add legend
plot(interpr, legend.only = T, horizontal = F, col = colpal, legend.args = list("Probability of use", side = 4, font = 1, line = 1.5, cex = 0.7),
     legend.shrink = 0.4,
     #smallplot= c(0.12,0.7, 0.06,0.09),
     axis.args = list(at = seq(0,1,0.25), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = seq(0,1,0.25), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -0.8, cex.axis = 0.7))

