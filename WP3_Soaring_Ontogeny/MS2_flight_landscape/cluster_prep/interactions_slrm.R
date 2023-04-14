#R script to predict the flyability of the alpine region for golden eagle juvies
#Elham Nourani, PhD
#Jan 31, 2023. Konstanz, DE


#setwd("GE_interactions/")

library(dplyr)
library(INLA)
inla.setOption(pardiso.license = "pardiso.lic")

#STEP 1: open data, define variables and formula -----------------------------------------------------------------

file_ls <- list.files(pattern = "inla_preds_input", full.names = T)

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

mean.beta <- 0
prec.beta <- 1e-4

#STEP 2: build one model for each interaction term ------------------------------------------------------------------------

results_list <- list()

for(i in file_ls){

new_data <- readRDS(i)

#run model
(b <- Sys.time())
M_pred <- inla(F_full, family = "Poisson",
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data,
               num.threads = 10, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy =  "pardiso", config = TRUE, cpo = T), #deactivate cpo #to save computing power
               control.inla(strategy = "adaptive", int.strategy = "eb"),
               inla.mode = "experimental", verbose = F)
Sys.time()-b #22 minutes

# posterior means of coefficients
graph <- as.data.frame(summary(M_pred)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

graph <- graph %>%
  mutate(Factor = rownames(graph))

#extract info to make prediction plots--------------------------------------------------------------------------
used_na <- which(is.na(new_data$used))

#extract information for rows that had NAs as response variables. append to the original new_data (with NAs as used)
preds <- new_data %>% 
  filter(is.na(used)) %>% 
  dplyr::select(c("location.long", "location.lat", "dem_100_z", "TRI_100_z", "weeks_since_emig_z")) %>% 
  mutate(preds = M_pred$summary.fitted.values[used_na,"mean"],
         preds_sd = M_pred$summary.fitted.values[used_na,"sd"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds)))

results <- list(graph, preds)

results_list <- append(results_list, results)

rm(M_pred, preds, graph, new_data)

gc(gc())

}

saveRDS(results_list, file = "results_list.rds")
