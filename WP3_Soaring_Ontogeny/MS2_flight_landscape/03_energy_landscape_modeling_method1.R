#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R (and temp_download_&prep.R)
#only includes static variables
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
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#open annotated data
cmpl_ann <- readRDS("alt_50_20_min_70_ind_static_time_ann_temp.rds") #data with montly temperature added in

cmpl_ann <- cmpl_ann %>% 
  mutate(days_since_emig_n = ceiling(as.numeric(days_since_emig)),#round up
         weeks_since_emig_n = ceiling(as.numeric(weeks_since_emig))) #135 unique weeks. median =  30.00, 3rd quart. = 61

# STEP 1: data exploration ----------------------------------------------------------------

#number of days available for each individual
cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_week = max(ceiling(as.numeric(weeks_since_emig)))) %>% #round up days since emigration. we don't want zeros
  ggplot(aes(x = individual.local.identifier, y = max_week)) +
  geom_col()

cmpl_ann %>% 
  group_by(individual.local.identifier) %>% 
  summarize(max_week = max(ceiling(as.numeric(weeks_since_emig)))) %>% 
  summarize(mode = Mode(max_week), #  22
            median = median(max_week), # 62.5
            mean = mean(max_week)) #  63


#terrain ~ days since emigration ..... no patterns in the plots

plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "dem_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "slope_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "aspect_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "slope_TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "TPI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "TRI_100")])
plot(cmpl_ann[cmpl_ann$used == 1, c("weeks_since_emig_n", "aspect_TPI_100")])


ggplot(cmpl_ann[cmpl_ann$used == 1 & cmpl_ann$weeks_since_emig_n <= 60,], aes(as.numeric(weeks_since_emig_n), slope_TPI_100)) +
  geom_point() +
  stat_smooth(aes(group = individual.local.identifier), method = "lm") +
  #geom_smooth() +
  theme_minimal() +
  theme(legend.position = "none")


# STEP 2: summary plots ----------------------------------------------------------------

#one set of boxplots for select days: 10, 50 , 100, 500

data_int <- cmpl_ann %>%
  filter(weeks_since_emig_n %in% c(1, 2, 4, 10, 50, 100)) %>% 
  mutate(weeks_f = as.factor(weeks_since_emig_n))


variables <- c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100" )

X11(width = 6, height = 9)

png("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/paper_prep/initial_figs/box_plots_n48.png", width = 6, height = 9, units = "in", res = 300)

par(mfrow= c(4,2), 
    oma = c(0,0,3,0), 
    las = 1,
    mgp = c(0,1,0))
for(i in 1:length(variables)){
  boxplot(data_int[,variables[i]] ~ data_int[,"weeks_f"], data = data_int, boxfill = NA, border = NA, main = variables[i], xlab = "", ylab = "")
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c(alpha("darkgoldenrod1", 0.9),"gray"), bty = "n", cex = 0.8)
  }
  boxplot(data_int[data_int$used == 1, variables[i]] ~ data_int[data_int$used == 1,"weeks_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = alpha("darkgoldenrod1", 0.9),  lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1, "weeks_f"])) - 0.15)
  boxplot(data_int[data_int$used == 0, variables[i]] ~ data_int[data_int$used == 0, "weeks_f"], outcol = alpha("black", 0.2),
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(data_int[data_int$used == 1 , "weeks_f"])) + 0.15)
  
}
dev.off()


# STEP 3: check for collinearity ----------------------------------------------------------------

cmpl_ann %>% 
  dplyr::select(c("dem_100", "slope_100", "aspect_100", "TRI_100", "TPI_100", "slope_TPI_100", "aspect_TPI_100", "month_temp" , "weeks_since_emig_n")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.6) #slope and TRI are correlated (.98)

# STEP 4: standardize variables ----------------------------------------------------------------

all_data <- cmpl_ann %>% 
  mutate_at(c("dem_100", "TRI_100", "slope_TPI_100", "weeks_since_emig_n","month_temp" , "days_since_emig_n"),
            list(z = ~(scale(.)))) 

# STEP 5: ssf modeling ----------------------------------------------------------------

## mid_step: investigate using clogit
# control for monthly temperature....

form1a <- used ~ dem_100_z * weeks_since_emig_n_z + 
  TRI_100_z * weeks_since_emig_n_z + 
  month_temp_z +
  strata(stratum)

ssf <- clogit(form1a, data = all_data)
summary(ssf)
plot_summs(ssf)

# INLA formula using interaction terms for time and predictor variables.

#repeat the random effect
all_data <- all_data %>% 
  mutate(ind1 = factor(individual.local.identifier),
         ind2 = factor(individual.local.identifier),
         ind3 = factor(individual.local.identifier),
         ind4 = factor(individual.local.identifier))

saveRDS(all_data, file = "alt_50_20_min_48_ind_static_temp_inlaready_wks.rds")


all_data <- readRDS("alt_50_20_min_48_ind_static_temp_inlaready_wks.rds")

#add one new row to unique strata instead of entire empty copies of strata. assign week since emigration and terrain values on a regular grid, so we can make a raster later on
set.seed(500)

#n needs to be large enough to cover the whole range of 
n <- 1500

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% 
  mutate(used = NA,
         month_temp_z = all_data %>%  filter(month_temp == median(month_temp)) %>%  slice(1) %>%  pull(month_temp_z), #assign mean temp for all rows
         weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 10), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 10), n, replace = T)) %>% 
  full_join(all_data)

saveRDS(new_data,"alt_50_20_min_48_ind_static_temp_inlaready_wmissing_wks_n1500.rds")


#the model will be run on the cluster. see cluster_prep/order_of_business_main_model.txt
