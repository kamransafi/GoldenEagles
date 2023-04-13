#script for analysis of golden eagle data for the dynamics of the energy landscape manuscript: modeling the energy landscape
#follows on from embc_segmentation.R and data_processing_&_annotation.R (and temp_download_&prep.R)
#only includes static variables. 
#try adding daily temp to see. seems like dem and temp are interacting!.. from 03_method2
#Feb 22. 2022. Elham Nourani. Konstanz, DE
#revisiting this code on April 11. 2023:
# - use hourly steps instead of 20 minutes
# - include step length in the analysis
# - consider using the 25 resolution for the terrain


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
#library(rastervis) #remotes::install_github("oscarperpinan/rastervis")
library(patchwork) #patching up interaction plots
library(oce) #color palette for interaction plots

setwd("/home/enourani/ownCloud/Work/Projects/GE_ontogeny_of_soaring/R_files/")
#setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


#open annotated data. prepared in 03_04b_clogit_randomization.R
data <-  readRDS("alt_50_60_min_55_ind_static_100.rds") %>% 
  mutate_at(c("step_length", "dem_100", "TRI_100", "slope_TPI_100", "weeks_since_emig"), list(z = ~(scale(.)))) %>%   #calculate the z scores
  mutate(ind1 = individual.local.identifier,
         ind2 = individual.local.identifier,
         ind3 = individual.local.identifier,
         month = month(timestamp))

saveRDS(data, file = "all_inds_annotated_static_apr23.rds") #for a version with dynamic variables, see all_inds_annotated_apr23.rds

# STEP 2: CLOGIT- ssf modeling exploration  ----------------------------------------------------------------

## mid_step: investigate using clogit
# control for monthly temperature....


form1 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100_z * step_length_z * weeks_since_emig_z 
  slope_TPI_100_z * step_length_z * weeks_since_emig_z + 
  strata(stratum)

ssf1 <- clogit(form1, data = data)
summary(ssf1)
plot_summs(ssf1)
modelsummary(ssf1) #AIC = 196770.1

form2 <- used ~ dem_100_z * step_length_z * weeks_since_emig_z + 
  TRI_100_z * step_length_z * weeks_since_emig_z + 
  slope_TPI_100_z * step_length_z * weeks_since_emig_z + 
  t2m +
  strata(stratum)

ssf2 <- clogit(form2, data = data)
summary(ssf2)
plot_summs(ssf2)
modelsummary(ssf2) #AIC = 189219.6


form1a <- used ~ dem_200_z * TRI_200_z + 
  strata(stratum)

form1a <- used ~ dem_200_z * TRI_200_z * weeks_since_emig_n + 
  strata(stratum)

form1a <- used ~ weeks_since_emig_n + 
  strata(stratum)

#I go more into detail of using this model in 03_04_clogit_workflow.R

# STEP 3: INLA parameterization ----------------------------------------------------------------
#the aim is to take seasonality into account

#do a test with a small sample of the data
sample <- data %>% 
  filter(stratum %in% sample(unique(data$stratum), 300, replace = F))


mean.beta <- 0
prec.beta <- 1e-4

#define variables and formula. use one variable for now to test the model 

#model without seasonality
F1 <- used ~ -1 + dem_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(3,0.05))))


F1 <- used ~ -1 + dem_100_z * step_length_z  +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6), fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(3,0.05))))



M2 <- inla(F1, family = "Poisson", 
          control.fixed = list(
            mean = mean.beta,
            prec = list(default = prec.beta)),
          data = sample, 
          num.threads = 7,
          control.predictor = list(compute = TRUE, link = 1), 
          control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))


M1 <- inla(F1, family = "Poisson",
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = sample,
           num.threads = 5, 
           control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
           control.compute = list(config = TRUE, cpo = F)) #deactivate cpo to save computing power
           #control.inla(strategy = "adaptive", int.strategy = "eb"),
           #inla.mode = "experimental", verbose = F)

Efxplot(M1)

##model with seasonality
F2 <- used ~ -1 +
  dem_100_z * weeks_since_emig_z * step_length_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE))

f(year.num, model = "iid",
  group = month.num, control.group = list(model = "ar1", 
                                          scale.model = TRUE)) #but this will go over the intercept, which we dont have. unless i use it when specifying the hyperparameter



M2 <- inla(F2, family = "Poisson",
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = sample,
           num.threads = 5)



Efxplot(list(M1,M2))

##full model with seasonality (to run on raven)
F_full <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE)) +
  f(ind2, TRI_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE)) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = TRUE))



#plot the effect of the grouped effect
X11()
par(mfrow = c(4, 3))
par(mar = c(3, 1.25, 1.5, 1.5))
#Table with summary of random effects; ID is for the YEAR
tab <-  M2$summary.random$ind1
inds <- unique(sample$ind1)

for(month in 1:12) {
  aux <- tab[(month-1) * n_distinct(sample$ind1) + 1:n_distinct(sample$ind1), ] #extract data for all individuals for month 1
  
  plot(as.factor(inds), aux[, "mean"], type = "l", xlab = "", ylab = "",
       main = paste0("Month ", month), ylim = c(-2, 2), las = 2, cex.axis = 0.75)
  lines(as.factor(inds), aux[, "0.025quant"], lty = 2)
  lines(as.factor(inds), aux[, "0.975quant"], lty = 2)
}

x.years <- 1949:1960

for(month in 1:12) {
  aux <- tab[(month-1) * 12 + 1:12, ]
  
  plot(x.years, aux[, "mean"], type = "l", xlab = "", ylab = "",
       main = paste0("Month ", month), ylim = c(4, 6.5), las = 2, cex.axis = 0.75)
  lines(x.years, aux[, "0.025quant"], lty = 2)
  lines(x.years, aux[, "0.975quant"], lty = 2)
}



#plot the output
# posterior means of coefficients
graph <- as.data.frame(summary(M1)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

graph <- graph %>%
  mutate(Factor = rownames(graph))

graph <- graph[graph$Factor != "weeks_since_emig_z",]
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
  #xlim(-0.1,0.6) +
  #scale_y_discrete(name = "",
  #                 labels = c("Weeks since dispersal * TRI","Weeks since dispersal * DEM", "TRI", "DEM")) +
  geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", size = 1) +
  theme_classic()


################################################################

#to take the seasonality into account, look at this. section 3.5.5.
#I can do dem, tri and slope tpi as functions of month. can i do it at the same time as defining the random slope
airp.data$month.num <- as.numeric(airp.data$month)
airp.data$year.num <- as.numeric(airp.data$year)
airp.iid.ar1 <- inla(log(airp) ~ 0 +  f(year.num, model = "iid",
                                        group = month.num, control.group = list(model = "ar1", 
                                                                                scale.model = TRUE)),
                     data = airp.data)
summary(airp.iid.ar1)

# STEP 3: create new data for prediction ----------------------------------------------------------------

#add one new row to unique strata instead of entire empty copies of strata. assign week since emigration and terrain values on a regular grid, so we can make a raster later on
set.seed(7777)

#n needs to be large enough to cover the whole range of 
n <- 1000

new_data <- all_data %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = T) %>% 
  mutate(used = NA,
         weeks_since_emig_n = sample(seq(min(all_data$weeks_since_emig_n),max(all_data$weeks_since_emig_n), length.out = 10), n, replace = T), 
         dem_100_z = sample(seq(min(all_data$dem_100_z),max(all_data$dem_100_z), length.out = 10), n, replace = T),
         TRI_100_z = sample(seq(min(all_data$TRI_100_z),max(all_data$TRI_100_z), length.out = 10), n, replace = T),
         t2m_z = sample(seq(min(all_data$t2m_z),max(all_data$t2m_z), length.out = 10), n, replace = T)) %>% #set to the mean of temperature. it is zero, because we are working with a z-score 
  full_join(all_data)

saveRDS(new_data,"alt_50_20_min_48_ind_static_daytemp_100_inlaready_wmissing_wks_n1000.rds")


#the model will be run on the cluster. see cluster_prep/order_of_business_main_model.txt
