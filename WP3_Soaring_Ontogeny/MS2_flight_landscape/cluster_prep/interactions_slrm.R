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
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))), #seasonality only for dem
    group = month, control.group = list(model = "ar1", scale.model = F, cyclic = T)) + #scale was true before
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

mean.beta <- 0
prec.beta <- 1e-4

#STEP 2: build one model for each interaction term ------------------------------------------------------------------------

results_list <- list()

for(i in file_ls){

new_data <- readRDS(i)

#work on a sample
#sample <- new_data %>% 
#  filter(!(is.na(used))) %>% 
#  filter(stratum %in% sample(unique(new_data$stratum), 300, replace = F)) %>%  #sample from original data
#  full_join(new_data[1:200,]) #add some missing data


#run model
(b <- Sys.time())
M_pred <- inla(F_full, family = "Poisson",
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data,
               num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
               control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted using the link function.
               control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
               #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
               inla.mode = "experimental")
Sys.time()-b #22 min

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
  dplyr::select(c("location.long", "location.lat", "dem_100", "TRI_100", "step_length", "weeks_since_emig",
                  "dem_100_z", "TRI_100_z", "step_length_z", "weeks_since_emig_z")) %>% 
  mutate(preds = M_pred$summary.fitted.values[used_na,"mean"],
         preds_sd = M_pred$summary.fitted.values[used_na,"sd"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds)))

results <- list(graph, preds)

results_list <- append(results_list, results)

rm(M_pred, preds, graph, new_data)

gc(gc())

}

saveRDS(results_list, file = "results_list.rds")

#-------------------------------------------
#try predictions without seasonality. the above model produces really large predictions. predictions are pretty much the same with this model
F_OG <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


(b <- Sys.time())
M_OG <- inla(F_OG, family = "Poisson",
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data,
               num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
               control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted using the link function.
               control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
               #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
               inla.mode = "experimental")
Sys.time()-b


#try predictions without NA values in step length. there are 20,000 NAs!! predictions are still quite large, but also there are many Inf predictions!!
new_data_cmpl <- new_data %>% 
  filter(!is.na(step_length))

(b <- Sys.time())
M_cmpl <- inla(F_full, family = "Poisson",
             control.fixed = list(
               mean = mean.beta,
               prec = list(default = prec.beta)),
             data = new_data_cmpl,
             num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
             control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted using the link function.
             control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
             #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
             inla.mode = "experimental")
Sys.time()-b

#try predictions without step length.. still very large
F_ns <-used ~ -1 +
  dem_100_z * weeks_since_emig_z +
  TRI_100_z * weeks_since_emig_z +
  slope_TPI_100_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))), #seasonality only for dem
    group = month, control.group = list(model = "ar1", scale.model = F, cyclic = T)) + #scale was true before
  f(ind2, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, slope_TPI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


(b <- Sys.time())
M_ns <- inla(F_ns, family = "Poisson",
             control.fixed = list(
               mean = mean.beta,
               prec = list(default = prec.beta)),
             data = new_data,
             num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
             control.predictor = list(compute = TRUE, link = 1), #this means that NA values will be predicted using the link function.
             control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
             #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
             inla.mode = "experimental")
Sys.time() - b

#make predictions, change inla call arguments
(b <- Sys.time())
M_pred2 <- inla(F_full, family = "Poisson",
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = new_data,
               num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted using the link function.
               #control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
               #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
               inla.mode = "experimental")
Sys.time()-b #22 min

used_na <- which(is.na(new_data$used))

#extract information for rows that had NAs as response variables. append to the original new_data (with NAs as used)
preds <- new_data %>% 
  filter(is.na(used)) %>% 
  dplyr::select(c("location.long", "location.lat", "dem_100", "TRI_100", "step_length", "weeks_since_emig",
                  "dem_100_z", "TRI_100_z", "step_length_z", "weeks_since_emig_z")) %>% 
  mutate(preds = M_pred2$summary.fitted.values[used_na,"mean"],
         preds_sd = M_pred2$summary.fitted.values[used_na,"sd"]) %>% 
  mutate(prob_pres = exp(preds)/(1+exp(preds)))


############ no predictions right now. just seeeing whether distance to ridge has seasonality or not

F_ridge <- used ~ -1 +
  dem_100_z * step_length_z * weeks_since_emig_z +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  distance_ridge_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))),
    group = month, control.group = list(model = "ar1", scale.model = F, cyclic = T)) + #scale was true before
  f(ind2, distance_ridge_z, model = "iid",
    hyper = list(theta=list(initial = log(1), prior = "pc.prec", param = c(3,0.05))), 
    group = month, control.group = list(model = "ar1", scale.model = F, cyclic = T)) #scale was true before

(b <- Sys.time())
M_r <- inla(F_ridge, family = "Poisson",
                control.fixed = list(
                  mean = mean.beta,
                  prec = list(default = prec.beta)),
                data = data,
                num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
                control.predictor = list(compute = TRUE), #this means that NA values will be predicted using the link function.
                #control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
                #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
                inla.mode = "experimental")
Sys.time()-b #22 min

############ no predictions right now. so from above it seems that dist to ridge is not corr with seasonality. so remove dem and seasonality!

F_ridge2 <- used ~ -1 +
  TRI_100_z * step_length_z * weeks_since_emig_z +
  slope_TPI_100_z * step_length_z * weeks_since_emig_z +
  distance_ridge_z * step_length_z * weeks_since_emig_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, TRI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind2, slope_TPI_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, distance_ridge_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
M_r2 <- inla(F_ridge2, family = "Poisson",
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = data,
            num.threads = 20, # be careful of memory. it was on 20 and 40 before. the job finished, but #there was a memory error in the error log.
            control.predictor = list(compute = TRUE), #this means that NA values will be predicted using the link function.
            #control.compute = list(openmp.strategy =  "huge", config = TRUE, cpo = T), #deactivate cpo #to save computing power
            #control.inla = list(strategy = "adaptive", int.strategy = "eb"), #"eb": results at the posterior mode for the hyperparameters
            inla.mode = "experimental")
Sys.time()-b #22 min
