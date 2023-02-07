#R script to predict the flyability of the alpine region for golden eagle juvies
#Elham Nourani, PhD
#Jan 31, 2023. Konstanz, DE


#setwd("GE_ALPS")



library(tidyverse)
library(INLA)
inla.setOption(pardiso.license = "pardiso.lic")

#STEP 1: open data and assign variables -----------------------------------------------------------------

weeks <- c(1, 2, 3, 4, 24, 48, 104) 

alps_data <- readRDS("alps_alt_50_20_min_48_ind_wmissing_Jun_temp.rds") %>% 
  as.data.frame() 

#try with a sample to test the cluster
#alps_data <- alps_data %>% 
#  group_by(used) %>% 
#  slice(1:1000)

weeks_z_info <- readRDS("weeks_since_z_info.rds")

#define variables and formula
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

#STEP 2: build one model per week ------------------------------------------------------------------------

results_list <- list()

for(x in weeks){ 
  #create the new data
  new_data <- alps_data %>%
    mutate(weeks_since_emig_n = replace_na(weeks_since_emig_n, x)) %>%
    rowwise() %>%
    mutate(weeks_since_emig_n_z = (weeks_since_emig_n - weeks_z_info$mean_wks)/weeks_z_info$sd_wks) #estimate z-score based on original data
  
  #run model
  
  M_pred <- inla(formulaM, family = "Poisson", 
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = new_data, 
                 num.threads = 20,
                 control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                 control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = F), #deactivate cpo to save computing power
                 inla.mode="experimental")
                 
                 # posterior means of coefficients
                 graph <- as.data.frame(summary(M_pred)$fixed)
                 colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
                 colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
                 colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
                 
                 graph <- graph %>% 
                   mutate(Factor = rownames(graph),
                          week = x)

                 
                 #extract info to make prediction plots--------------------------------------------------------------------------
                 used_na <- which(is.na(new_data$used))
                 
                 #extract information for rows that had NAs as response variables
                 preds <- data.frame(location.long = new_data[is.na(new_data$used) ,"location.long"],
                                     location.lat = new_data[is.na(new_data$used) ,"location.lat"],
                                     dem_100_z = new_data[is.na(new_data$used) ,"dem_100_z"],
                                     TRI_100_z = new_data[is.na(new_data$used) ,"TRI_100_z"],
                                     weeks_since_emig_n_z = new_data[is.na(new_data$used) ,"weeks_since_emig_n_z"],
                                     preds = M_pred$summary.fitted.values[used_na,"mean"], 
                                     preds_sd = M_pred$summary.fitted.values[used_na,"sd"]) %>%  
                   mutate(prob_pres = exp(preds)/(1+exp(preds)),#this should be between 0-1
                          week = x) 
                 
                 week_results <- list(graph, preds)
                 
                 results_list <- append(results_list, week_results)
                 
                 rm(M_pred,preds, graph, new_data)
                 gc(gc())
}
                   
saveRDS(results_list, file = "results_list.rds")
