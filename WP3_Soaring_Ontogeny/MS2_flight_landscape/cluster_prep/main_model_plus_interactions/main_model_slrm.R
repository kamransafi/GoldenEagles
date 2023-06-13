#R script to predict the flyability of the alpine region for golden eagle juvies
#Elham Nourani, PhD
#Jan 31, 2023. Konstanz, DE


#setwd("GE_interactions/")

library(dplyr)
library(INLA)
inla.setOption(pardiso.license = "pardiso.lic")

#STEP 1: open data, define variables and formula -----------------------------------------------------------------

data <- readRDS("new_data_ssf_inla_preds.rds") %>% filter(!(is.na(used)))


#define variables and formula
F_full <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  ridge_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", cyclic = T, scale.model = TRUE)) +
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, ridge_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


mean.beta <- 0
prec.beta <- 1e-4

#STEP 2: build the model ------------------------------------------------------------------------

M_main <-  inla(F_full, family = "Poisson",
                control.fixed = list(
                  mean = mean.beta,
                  prec = list(default = prec.beta)),
                data = data,
                num.threads = 1, 
                control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted.
                control.compute = list(openmp.strategy =  "pardiso", config = TRUE, cpo = T))


#extract model validation results
eval <- data.frame(CPO = mean(M_main$cpo$cpo, na.rm = T),
                   Mlik = as.numeric(M_main$mlik[,1][2])) 

saveRDS(eval, file = "eval_M_main100_hrly_nopred.rds")

# extract info for coefficient plots

# posterior means of coefficients
graph <- as.data.frame(summary(M_main)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

graph$Factor <- rownames(graph)

saveRDS(graph, file = "graph_M_main100_hrly_nopred.rds")

#-----------------------------------------------------------------------------------------------
# STEP 2: EXTRACT & SAVE PREDICTIONS
# 
# #extract predictions for interaction terms
# used_na <- which(is.na(data$used))
# 
# preds <- data %>% 
#   slice(used_na) %>% #extract information for rows that had NAs as response variables
#   mutate(preds = M_main$summary.fitted.values[used_na,"mean"]) %>% 
#   mutate(prob_pres = exp(preds)/(1 + exp(preds))) #this should be between 0-1
# 
# saveRDS(preds, file = "preds_M_main100_hrly.rds")


#extract info for ind- variation

summ_rndm <- M_main$summary.random

saveRDS(summ_rndm, file = "rnd_coeff_M_main100_hrly_nopred.rds")


