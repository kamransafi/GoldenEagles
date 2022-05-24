#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R
#the first attempt will only include static variables. reasons: 1) movebank annotation isn't working and 2) res of static variables is higher
#Feb 22. 2022. Elham Nourani. Konstanz, DE

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

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#open annotated data (static variables and time since fledging and emigration)
load("alt_50_20_min_25_ind_static_time_ann_weeks.RData") #ann_cmpl

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig)), #120 unique weeks
         stratum = paste(individual.local.identifier, burst_id, step_id, sep = "_"))

 
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
  mutate_at(c("dem_100", "TRI_100", "days_since_emig_n", "weeks_since_emig_n"),
            list(z = ~(scale(.)))) 

#see previous versions for exploration using clogit
# STEP 5: ssf modeling ----------------------------------------------------------------

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

saveRDS(all_data, file = "alt_50_20_min_25_ind_static_inlaready_wks.rds")


all_data <- readRDS("alt_50_20_min_25_ind_static_inlaready_wks.rds")



#add one new row to unique strata instead of entire empty copies of strata. assign day since emigration and terrain values on a regular grid
set.seed(500)

n <- 1000
new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,
        # weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 10), n, replace = T), 
         weeks_since_emig_n_z = sample(seq(min(all_data$weeks_since_emig_n_z),max(all_data$weeks_since_emig_n_z), length.out = 10), n, replace = T),
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T), #regular intervals, so we can make a raster later on
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)

saveRDS(new_data,"alt_50_20_min_25_ind_static_inlaready_wmissing_wks.rds")

new_data <- readRDS("alt_50_20_min_25_ind_static_inlaready_wmissing_wks.rds")

#model formula. slope and TRI are correlated. This version includes on TRI and dem. Martina's preprint suggests TRI, dem and slope tpi, but the latter was insig
formulaM <- used ~ -1 + 
  dem_100_z * weeks_since_emig_n_z +
  TRI_100_z * weeks_since_emig_n_z + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, TRI_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
  
mean.beta <- 0
prec.beta <- 1e-4 


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

save(M_marti_c, file = "inla_model_w_random_wks.RData")

Efxplot(M_marti_c) # all of them are very similar in terms of cpo and Mlik

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
Sys.time() - b #1000 missing values: 17 mins

saveRDS(M_pred3, file = "inla_model_predw_random_wks.rds")

#make plot for coefficients ---------------------- 

# posterior means of coefficients
graph <- as.data.frame(summary(M_marti_c)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

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

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/static_landscape_figs/coeffs.png", 
    width = 4.7, height = 2.7, units = "in", res = 300)

par(cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 8.5, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-0.05,0.05), ylim = c(0.7,5.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph$Estimate, graph$Factor_n, col = "cornflowerblue", pch = 20, cex = 2)
arrows(graph$Lower, graph$Factor_n,
       graph$Upper, graph$Factor_n,
       col = "cornflowerblue", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line


#add axes
axis(side= 1, at = seq(-0.05, 0.05, by =  0.01), labels = seq(-0.05, 0.05, by =  0.01), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at = c(1:5),
     labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "Weeks since dispersal","DEM"),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

dev.off()


# STEP 6: ssf output plots: the interaction plots ----------------------------------------------------------------

M_pred3 <- readRDS("inla_model_predw_random_wks.rds") 
new_data <- readRDS("alt_50_20_min_25_ind_static_inlaready_wmissing.rds")
all_data <- readRDS("alt_50_20_min_25_ind_static_inlaready_wks.rds")


#extract predicted values
used_na <- which(is.na(new_data$used))

y_axis_var <- c("dem_100_z", "TRI_100_z")
x_axis_var <- "weeks_since_emig_n_z"

labels <- data.frame(var = c("dem_100_z", "TRI_100_z","weeks_since_emig_n_z"),
                     label = c("Altitude (m.asl)", "Terrain ruggedness index", "Weeks since dispersal"))

#extract probability of presence for missing values
fitted_values <- data.frame(weeks_since_emig_n_z = new_data[is.na(new_data$used) ,"weeks_since_emig_n_z"],
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
  r <- rast(avg_pred, crs = "+proj=longlat +datum=WGS84") #use the ggplot code at the end to plot it
 
  
  saveRDS(r, file = paste0("inla_pred_w_random_", i,".rds"))
  
}

#######use ggplot

r_dem <- readRDS("inla_pred_w_random_dem_100_z.rds")
r_rug <- readRDS("inla_pred_w_random_TRI_100_z.rds")

X11(width = 6, height = 4)

p_rugg <- ggplot(data = as.data.frame(r_rug, xy = T)) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradientn(colours = oce::oceColorsPalette(80), limits = c(0,0.5), 
                       na.value = "white", name = "Intensity of use") +
  theme_minimal() +
  labs(x = "Weeks since dispersal", y = "Ruggedness")

ggsave(plot = p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_",i, ".png"), 
       width = 6, height = 4, dpi = 300)



X11(width = 6, height = 4)

p_dem <- ggplot(data = as.data.frame(r_dem, xy = T)) +
  geom_tile(aes(x = x, y = y, fill = avg_pres)) +
  scale_fill_gradientn(colours = oce::oceColorsPalette(80), limits = c(0,0.5), 
                       na.value = "white", name = "Intensity of use") +
  theme_minimal() +
  labs(x = "Weeks since dispersal", y = "Elevation (m)")

ggsave(plot = p, filename = paste0("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_",i, ".png"), 
    width = 6, height = 4, dpi = 300)


#put both plots in one device
X11(width = 10, height = 4)
combined <- p_dem + p_rugg & theme(legend.position = "right")
p_2 <- combined + plot_layout(guides = "collect")

ggsave(plot = p_2, filename = "/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/inla_preds_static_w_rnd_wk_2vars.png", 
       width = 10, height = 4, dpi = 300)

