#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#Feb 22. 2022. Elham Nourani. Konstanz, DE

library(tidyverse)
library(corrr)
library(INLA)
library(fields)
library(raster)

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

# quick and dirty using clogit? for what time period???



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
set.seed(200)

n <- 500
new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,
         days_since_emig_n = sample(seq(min(all_data$days_since_emig_n),max(all_data$days_since_emig_n), length.out = 10), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T), #regular intervals for wind support and delta t, so we can make a raster later on
         slope_100_z = sample(seq(min(all_data$slope_100_z),max(all_data$slope_100_z), length.out = 10), n, replace = T),
         aspect_100_z = sample(seq(min(all_data$aspect_100_z),max(all_data$aspect_100_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)


#model formula. slope and TRI are correlated
formulaM <- used ~ -1 + dem_100_z * days_since_emig_n_z + slope_100_z * days_since_emig_n_z + aspect_100_z * days_since_emig_n_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, slope_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, aspect_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

#model
(b <- Sys.time())
M <- inla(formulaM, family = "Poisson", 
          control.fixed = list(
            mean = mean.beta,
            prec = list(default = prec.beta)),
          data = all_data, 
          num.threads = 10,
          control.predictor = list(compute = TRUE, link = 1), 
          control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b  #5.515833 mins

save(M, file = "inla_model_1.RData")

#Model for predictions
(b <- Sys.time())
M_pred <- inla(formulaM, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #500 missing values: 7.47236 mins

save(M_pred, file = "inla_model_pred1.RData")

#try link = 1
(b <- Sys.time())
M_pred2 <- inla(formulaM, family = "Poisson", 
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data, 
               num.threads = 10,
               control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
Sys.time() - b #1.069171 mins without random effects. 7.47236 mins


# STEP 6: ssf output plots ----------------------------------------------------------------

#interaction plots
#extract predicted values
used_na <- which(is.na(new_data$used))

y_axis_var <- c("dem_100_z", "slope_100_z", "aspect_100_z")
x_axis_var <- "days_since_emig_n_z"


#extract information for rows that had NAs as response variables
preds <- data.frame(dem_100_z = new_data[is.na(new_data$used) ,"dem_100_z"],
                    slope_100_z = new_data[is.na(new_data$used) ,"slope_100_z"],
                    aspect_100_z = new_data[is.na(new_data$used) ,"aspect_100_z"],
                    days_since_emig_n_z = new_data[is.na(new_data$used) ,"days_since_emig_n_z"],
                    preds = M_pred$summary.fitted.values[used_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds))) #this should be between 0-1

#create a raster of predictions for each variable

for (i in y_axis_var){
  y <- preds
  
  avg_pred <- preds %>% 
    group_by(x_axis_var, i) %>%  
    summarise(avg_pres = mean(prob_pres)) %>% 
    ungroup() %>% 
    mutate(wspt_backtr = slope_100_z * attr(all_data$slope_100_z, 'scaled:scale') + attr(all_data$slope_100_z, 'scaled:center'),
           dt_backtr = days_since_emig_n_z * attr(all_data$days_since_emig_n_z, 'scaled:scale') + attr(all_data$days_since_emig_n_z, 'scaled:center')) %>% 
    dplyr::select(-c("days_since_emig_n_z","slope_100_z")) %>% 
    as.data.frame()
  
}
avg_preds_slope <- preds %>% 
  group_by(days_since_emig_n_z, slope_100_z) %>%  
  summarise(avg_pres = mean(prob_pres)) %>% 
  ungroup() %>% 
  mutate(wspt_backtr = slope_100_z * attr(all_data$slope_100_z, 'scaled:scale') + attr(all_data$slope_100_z, 'scaled:center'),
         dt_backtr = days_since_emig_n_z * attr(all_data$days_since_emig_n_z, 'scaled:scale') + attr(all_data$days_since_emig_n_z, 'scaled:center')) %>% 
  dplyr::select(-c("days_since_emig_n_z","slope_100_z")) %>% 
  as.data.frame()


coordinates(avg_preds_slope) <-~ wspt_backtr + dt_backtr 
gridded(avg_preds_slope) <- TRUE
r <- raster(avg_preds_slope)


#interpolate. for visualization purposes
surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])

grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 2),
                   y = seq(from = extent(r)[3],to = extent(r)[4],by = 2))

grd$coords <- matrix(c(grd$x,grd$y),ncol=2)

surf.1.pred <- predict.Krig(surf.1,grd$coords)
interpdf <- data.frame(grd$coords,surf.1.pred)

colnames(interpdf) <- c("slope","days_since_emig","prob_pres")

coordinates(interpdf) <- ~ slope + days_since_emig
gridded(interpdf) <- TRUE
interpr <- raster(interpdf)


#create a color palette
cuts <- c(0, 0.25,0.5,0.75,1) #set breaks
pal <- colorRampPalette(c("aliceblue", "lightskyblue1", "khaki2", "sandybrown", "salmon2","tomato"))
colpal <- pal(200)



#plot
X11(width = 5, height = 4)

par(cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(0, 0, 0, 1.5),
    bty = "n",
    mgp = c(1,0.5,0)
)


raster::plot(interpr, col = colpal, axes = F, box = F, legend = F, ext = extent(c(-22, 28.9, -9.7, 14))) #crop to the extent of observed data

#add axes
axis(side = 1, at = seq(-20,30,10),
     labels = seq(-20,30,10),
     tick = T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck = -.015, line = -5.75, cex.axis = 0.7) #tick marks smaller than default by this proportion

axis(side = 2, at = c(-5, 0,5, 10), labels = c(-5, 0,5, 10), 
     tick = T ,col = NA, col.ticks = 1, tck = -.015, las = 2, cex.axis = 0.7)

lines(x = c(-21.9, -21.9), y = c(-9.9,13.9))
abline(h =-10)

#axis titles
mtext("Wind support (m/s)", 1, line = -4, cex = 0.9, font = 3)
mtext(expression(italic(paste(Delta,"T", "(°C)"))), 2, line = 1.2, cex = 0.9)

#add legend
plot(interpr, legend.only = T, horizontal = F, col = colpal, legend.args = list("Probability of use", side = 4, font = 1, line = 1.5, cex = 0.7),
     legend.shrink = 0.4,
     #smallplot= c(0.12,0.7, 0.06,0.09),
     axis.args = list(at = seq(0,1,0.25), #same arguments as any axis, to determine the length of the bar and tick marks and labels
                      labels = seq(0,1,0.25), 
                      col = NA, #make sure box type in par is set to n, otherwise axes will be drawn on the legend :p
                      col.ticks = NA,
                      line = -0.8, cex.axis = 0.7))