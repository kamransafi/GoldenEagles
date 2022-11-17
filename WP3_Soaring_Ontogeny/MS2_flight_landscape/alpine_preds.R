#run INLA models on the cluster with the Alpine region data

#see updated version here: /home/mahle68/ownCloud/Work/cluster_computing/GE_inla_static/alps_preds

weeks <- c(1, 2, 3, 4, 24, 48, 104) 

alps_data <- readRDS( "inla_preds_for_cluster/generic_alt_50_20_min_48_ind_wmissing.rds")

#define variables and formula
formula_pred <- used ~ -1 + 
  dem_100_z * weeks_since_emig_n_z +
  TRI_100_z * weeks_since_emig_n_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(ind1, dem_100_z, model = "iid",
    hyper=list(theta=list(initial=log(1)ss,fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, TRI_100_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

mean.beta <- 0
prec.beta <- 1e-4 

#alternatively, this can be done on the cluster:
lapply(weeks, function(x){
  
  #create the new data
  new_data <- alps_data %>%
    mutate(weeks_since_emig_n = i) #what about the z?

  #run model
  
  M_pred <- inla(formula_pred, family = "Poisson", 
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = new_data, 
                 num.threads = 20,
                 control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
                 control.compute = list(openmp.strategy = "huge", config = TRUE), #deactivate CPO
                 inla.mode="experimental")

  #save results
  eval <- data.frame(CPO = mean(M_pred$cpo$cpo, na.rm = T), #0.95 
                     Mlik = as.numeric(M_pred$mlik[,1][2])) #-161036.6
  
  saveRDS(eval, file = "eval_M_pred.rds")
  
  # extract info for coefficient plots
  
  # posterior means of coefficients
  graph <- as.data.frame(summary(M_pred)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  #graph$Model<-i
  graph$Factor <- rownames(graph)
  
  saveRDS(graph, file = "graph_M_pred.rds")
  
  #extract info to make prediction plots
  
  used_na <- which(is.na(new_data$used))
  
  #extract information for rows that had NAs as response variables
  preds <- data.frame(dem_100_z = new_data[is.na(new_data$used) ,"dem_100_z"],
                      TRI_100_z = new_data[is.na(new_data$used) ,"TRI_100_z"],
                      weeks_since_emig_n_z = new_data[is.na(new_data$used) ,"weeks_since_emig_n_z"],
                      preds = M_pred$summary.fitted.values[used_na,"mean"]) %>% 
    mutate(prob_pres = exp(preds)/(1+exp(preds))) #this should be between 0-1
  
  
  saveRDS(preds, file = "preds_M_preds.rds")
  
})