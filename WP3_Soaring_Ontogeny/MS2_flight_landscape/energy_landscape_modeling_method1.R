#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#Feb 22. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(corrr)
library(INLA)
library(fields)
library(raster)
library(survival)
library(ggregplot)
library(terra)
library(gstat) #for interpolations

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

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

# quick and dirty using clogit. for what time period???

form1a <- used ~ dem_100_z * days_since_emig_n_z + 
  #aspect_100_z * days_since_emig_n_z + #AIC with aspect = 47442.4; without aspect:  47444.29 (both without individual!!)
  TRI_100_z * days_since_emig_n_z + 
  TPI_100_z * days_since_emig_n_z + 
  slope_TPI_100_z * days_since_emig_n_z + 
  aspect_TPI_100_z * days_since_emig_n_z  + 
  strata(stratum)

ssf_full <- clogit(form1a, data = all_data)
summary(ssf_full) #so, include at least dem, TRI, TPI, slope_TPI and aspect_TPI in the model

#just as a try, include a 3-way interaction
form1b <- used ~ dem_100_z * days_since_emig_n_z * as.factor(ind1) +
  TPI_100_z * days_since_emig_n_z  * as.factor(ind1) + 
  strata(stratum)

ssf_test <- clogit(form1b, data = all_data)
summary(ssf_test) #lol. it works, but too much hassle to do for all variables and get anything useful out of

# INLA formula using interaction terms for time and predictor variables.
# control for month of year or temperature....

#repeat the random effect
all_data <- all_data %>% 
  mutate(ind1 = factor(individual.local.identifier),
         ind2 = factor(individual.local.identifier),
         ind3 = factor(individual.local.identifier),
         ind4 = factor(individual.local.identifier),
         days_f1 = factor(days_since_emig_n),
         days_f2 = factor(days_since_emig_n),
         days_f3 = factor(days_since_emig_n),
         days_f4 = factor(days_since_emig_n))

save(all_data, file = "alt_50_20_min_25_ind_static_inlaready.RData")


load("alt_50_20_min_25_ind_static_inlaready.RData")

mean.beta <- 0
prec.beta <- 1e-4 

#add one new row to unique strata instead of entire empty copies of strata. assign day since emigration and terrain values on a regular grid
set.seed(500)

n <- 1000
new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,
         days_since_emig_n = sample(seq(min(all_data$days_since_emig_n),max(all_data$days_since_emig_n), length.out = 10), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T), #regular intervals for wind support and delta t, so we can make a raster later on
         slope_100_z = sample(seq(min(all_data$slope_100_z),max(all_data$slope_100_z), length.out = 10), n, replace = T),
         slope_TPI_100_z = sample(seq(min(all_data$slope_TPI_100_z),max(all_data$slope_TPI_100_z), length.out = 10), n, replace = T),
         aspect_100_z = sample(seq(min(all_data$aspect_100_z),max(all_data$aspect_100_z), length.out = 10), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 10), n, replace = T),
         TPI_100_z = sample(seq(min(all_data$TPI_100_z),max(all_data$TPI_100_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)

saveRDS(new_data,"alt_50_20_min_25_ind_static_inlaready_wmissing.rds")


#model formula. slope and TRI are correlated. This version includes on TRI and dem. Martina's preprint suggests TRI, dem and slope tpi, but the latter was insig
formulaM <- used ~ -1 + 
  dem_100_z * days_since_emig_n_z +
  TRI_100_z * days_since_emig_n_z + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, TRI_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
  
  
#model
(b <- Sys.time())
M_marti_c <- inla(formulaM, family = "Poisson", 
          control.fixed = list(
            mean = mean.beta,
            prec = list(default = prec.beta)),
          data = all_data, 
          num.threads = 10,
          control.predictor = list(compute = TRUE, link = 1), 
          control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b  #without random effects: 1.76 min. with random effects: 16.7 min

save(M_marti_c, file = "inla_model_w_random.RData")



Efxplot(list(M_marti, M_marti_b, M_marti_c)) # all of them are very similar in terms of cpo and Mlik

#Model for predictions
(b <- Sys.time())
M_pred3 <- inla(formulaM, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #500 missing values: 7.47236 mins

save(M_pred3, file = "inla_model_predw_random.RData")



# #try link = 1
# (b <- Sys.time())
# M_pred2 <- inla(formulaM, family = "Poisson", 
#                control.fixed = list(
#                  mean = mean.beta,
#                  prec = list(default = prec.beta)),
#                data = new_data, 
#                num.threads = 10,
#                control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted.
#                control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
# Sys.time() - b #1.069171 mins without random effects. 7.47236 mins
# 
# save(M_pred2, file = "inla_model_pred2_norandom.RData")


# STEP 6: ssf output plots: the interaction plots ----------------------------------------------------------------

load("inla_model_predw_random.RData") #M_pred3
new_data <- readRDS("alt_50_20_min_25_ind_static_inlaready_wmissing.rds")
load("alt_50_20_min_25_ind_static_inlaready.RData")

#extract predicted values
used_na <- which(is.na(new_data$used))

y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "days_since_emig_n_z"

labels <- data.frame(var = c("dem_100_z", "TRI_100_z","days_since_emig_n_z"),
                     label = c("Altitude (m.asl)", "Terrain ruggedness index", "Days since dispersal"))

#extract probability of presence for missing values
fitted_values <- data.frame(days_since_emig_n_z = new_data[is.na(new_data$used) ,"days_since_emig_n_z"],
                            preds = M_pred3$summary.fitted.values[used_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds))) #Inverse-logit transformation to get the probability (between 0 and 1)

#extract center and scale values for time variable, to be used for back transformation. The y-axis attributes will be extracted in the for loop
x_axis_attr_scale <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:scale')
x_axis_attr_center <- attr(all_data[,colnames(all_data) == x_axis_var],'scaled:center')

for (i in y_axis_var){
  
  fitted_df <- fitted_values %>% 
    mutate(new_data[is.na(new_data$used), i])
  #extract scale and center values needed to back transform the scaled values
  i_scale <- attr(all_data[,colnames(all_data) == i],'scaled:scale')
  i_center <- attr(all_data[,colnames(all_data) == i],'scaled:center')
  
  #summarize values, so each (x,y) combo has one probability value
  avg_pred <- fitted_df %>% 
    group_by_at(c(1,4)) %>%  #group by days since emigration and i
    summarise(avg_pres = mean(prob_pres)) %>% 
    ungroup() %>% 
    mutate(dplyr::select(.,all_of(i)) * i_scale + i_center, #back-transform the values for plotting
           dplyr::select(.,all_of(x_axis_var)) * x_axis_attr_scale + x_axis_attr_center ) %>% #these columns replace the original columns 1 and 2
    rename(x = 1, #days since fledging
           y = 2) %>%  #y axis variables
    as.data.frame()
  
  #create a raster
  r <- rast(avg_pred, crs = "+proj=longlat +datum=WGS84")
  
  # #interpolation using empty raster
  # empty_r <- rast(extent = ext(r), ncol = 191, nrow = 1086, crs = "+proj=longlat +datum=WGS84") #based on the interpolation using tps
  # r_interp <- terra::interpolate(empty_r, r)
  # 
  # 
  #interpolation using the nearest neighbor method. Tps got rid of rare large values!
 
  gs <- gstat(formula = avg_pres ~ 1, locations =~ x + y, data = avg_pred, nmax = 5, set = list(idp = 0)) #nmax is the number of neighbors; idp = all neighbors equally weighted
  nn <- interpolate(r, gs, debug.level = 0)
  
  
  #create empty raster for interpolation
  # surface <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2), as.data.frame(r,xy = T)[,3])
  # grd <- expand.grid(x = seq(from = ext(r)[1],to = ext(r)[2],by = 4),
  #                    y = seq(from = ext(r)[3],to = ext(r)[4],by = 4))
  # 
  # #make sure elevation starts at 0. this messes up the proportions for some reason.
  # # if(i == "dem_100_z"){
  # #   grd <- expand.grid(x = seq(from = ext(r)[1],to = ext(r)[2],by = 2),
  # #                      y = seq(from = 0,to = ext(r)[4],by = 2))
  # # }
  # 
  # grd$coords <- matrix(c(grd$x,grd$y), ncol=2)
  # 
  # preds <- predict.Krig(surface,grd$coords)
  # interpdf <- data.frame(grd$coords, preds)
  # 
  # colnames(interpdf) <- c("x_backtr","y_backtr","prob_pres")
  # 
  # interpr <- rast(interpdf, crs = "+proj=longlat +datum=WGS84")
  # 
  # 
  #create a color palette
  cuts <- c(0, 0.25,0.5,0.75,1) #set breaks
  pal <- colorRampPalette(c("aliceblue", "lightskyblue1", "khaki2", "sandybrown", "salmon2","tomato"))
  colpal <- pal(80)
  
  #plot
  X11(width = 5, height = 4)
  
  par(cex = 0.7,
      oma = c(1,2,0, 4),
      mar = c(1, 2, 2, 5),
      tcl = -0.5,
      bty = "n"
  )
  
  plot(nn,1, col = colpal, axes = F, box = F, legend = T, plg = list(at = c(0.25,0.35,0.45), labels = c(0.25,0.35,0.45))) #plg includes legend options
  
  #add axes
  axis(side = 1, at = x_axis_lab, #x_axis
       labels = x_axis_lab,
       tick = T , col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
       tck = -.015, line = 0, cex.axis = 0.7) #tick marks smaller than default by this proportion
  
  axis(side = 2, at = y_axis_lab, labels = y_axis_lab, 
       tick = T ,col = NA, col.ticks = 1, tck = -.015, las = 2, cex.axis = 0.7)
  
  #x and y axis lines
  abline(v = ext(nn)[1])
  abline(h = ext(nn)[3])
  
  #axis titles
  mtext(labels %>% filter(var == i) %>% pull(label), 2, line = 3, cex = 0.9, font = 3)
  mtext(labels %>% filter(var == x_axis_var) %>% pull(label), 1, line = 2.5, cex = 0.9, font = 3)
  
  #add legend title
  graphics::text(labels = "Probability of use", x = , line = 2.5, cex = 0.9, font = 3)
}


saveRDS(nn, file = "inla_pred_w_random_tri.rds")
saveRDS(nn, file = "inla_pred_w_random_dem.rds")


#try the plot with ggplot
 ggplot() +
  geom_raster(data = wind_df %>% filter(unique_hour == i), aes(x = lon, y = lat, fill = wind_speed))


 gplot(nn) +
   geom_tile(aes(fill = value)) +
   scale_fill_gradientn(colours = colpal, limits = c(0.2,0.7), name = "Intensity of use")
 
 
 
