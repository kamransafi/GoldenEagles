#R script to predict the flyability of the alpine region for golden eagle juvies
#Elham Nourani, PhD
#Jan 31, 2023. Konstanz, DE


#setwd("slrm_jobs")

library(dplyr)
library(INLA)
inla.setOption(pardiso.license = "pardiso.lic")

#STEP 1: open data and assign variables -----------------------------------------------------------------

weeks_z_info <- readRDS("weeks_since_z_info.rds")

weeks <- data.frame(weeks_since_emig_n = c(1, 2, 3, 4, 24, 48, 104)) %>% 
  rowwise() %>% 
  mutate(weeks_since_emig_n_z = (weeks_since_emig_n - weeks_z_info$mean_wks)/weeks_z_info$sd_wks) %>% 
  ungroup() %>% 
  as.data.frame()

alps_data <- readRDS("alps_alt_50_20_min_48_ind_static200_wmissing.rds")

#define variables and formula
formulaM <- used ~ -1 +
  dem_200_z * weeks_since_emig_n_z +
  TRI_200_z * weeks_since_emig_n_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_200_z, model = "iid",
    hyper = list(theta=list(initial = log(1), fixed = F, prior="pc.prec", param = c(3,0.05)))) +
  f(ind2, TRI_200_z,  model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(3,0.05))))

f(year.num, model = "iid",
  group = month.num, control.group = list(model = "ar1", 
                                          scale.model = TRUE)) #but this will go over the intercept, which we dont have. unless i use it when specifying the hyperparameter

mean.beta <- 0
prec.beta <- 1e-4

#STEP 2: build one model per week ------------------------------------------------------------------------

results_list <- list()

for(x in 1:nrow(weeks)){
  
  new_data <- alps_data
  new_data[is.na(new_data$weeks_since_emig_n),"weeks_since_emig_n"] <- weeks[x,"weeks_since_emig_n"]
  new_data[is.na(new_data$weeks_since_emig_n_z),"weeks_since_emig_n_z"] <- weeks[x,"weeks_since_emig_n_z"]
  
  #run the model
  #b <- Sys.time()
  M_pred <- inla(formulaM, family = "Poisson",
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = new_data,
                 num.threads = 1, # be careful of memory. it was on 20 and 40 before. the job finished, but there was a memory error in the error log.
                 control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                 control.compute = list(openmp.strategy =  "pardiso", config = TRUE, cpo = F), #deactivate cpo to save computing power
                 control.inla(strategy = "adaptive", int.strategy = "eb"),
                 inla.mode="experimental", verbose = F)
  #Sys.time()- b #25 min for one model
  
  # posterior means of coefficients
  graph <- as.data.frame(summary(M_pred)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph <- graph %>%
    mutate(Factor = rownames(graph),
           week = weeks[x,"weeks_since_emig_n"])
  
  
  #extract info to make prediction plots--------------------------------------------------------------------------
  used_na <- which(is.na(new_data$used))
  
  #extract information for rows that had NAs as response variables. append to the original new_data (with NAs as used)
  preds <- new_data %>% 
    filter(is.na(used)) %>% 
    dplyr::select(c("location.long", "location.lat", "dem_200", "TRI_200", "dem_200_z", "TRI_200_z", "weeks_since_emig_n_z", "weeks_since_emig_n")) %>% 
    mutate(preds = M_pred$summary.fitted.values[used_na,"mean"],
           preds_sd = M_pred$summary.fitted.values[used_na,"sd"]) %>% 
    mutate(prob_pres = exp(preds)/(1+exp(preds)))
  
  week_results <- list(graph, preds)
  
  results_list <- append(results_list, week_results)
  
  rm(M_pred, preds, graph, new_data)
  gc(gc())
}

saveRDS(results_list, file = "results_list.rds")
