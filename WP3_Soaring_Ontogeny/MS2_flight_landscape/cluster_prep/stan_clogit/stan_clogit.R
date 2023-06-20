#R script to run the clogit using stan. this will enable me to add individual random effects
#Elham Nourani, PhD
#Jun 12, 2023. Canberra, AU.


library(dplyr)
library(rstanarm)



#STEP 1: open data, define variables and formula -----------------------------------------------------------------
data <- readRDS("all_inds_annotated_static_3yrs_apr23.rds") %>% 
  mutate(stratum_IDs = stratum,
         stratum = as.numeric(as.factor(stratum))) %>% 
  arrange(stratum) #to use stan_clogit, data should be ordered by startum


#STEP 2: run the model -----------------------------------------------------------------
f <- used ~ TRI_100_z * step_length_z * weeks_since_emig_z + 
  ridge_100_z * step_length_z * weeks_since_emig_z + 
  (TRI_100_z | ind1) + (ridge_100_z | ind2)

(b <- Sys.time())
ssf <- stan_clogit(f, 
                   strata = stratum,
                   data = data,
                   QR = TRUE,
                   chains = 5,
                   cores = 5,
                   warmup = 100,
                   iter = 200)
Sys.time() - b

saveRDS(ssf, file = "ssf_stan_clogit_model.rds")



