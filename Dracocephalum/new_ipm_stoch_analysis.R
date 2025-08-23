## This script calculates lambda for different environmental variable levels
## under asymptotic behaviour (i.e. long run simulations)
rm(list=ls())

library(dplyr)
library(lme4)
library(ipmr)
library(mgcv)
library(parallel)
library(forecast)

setwd("/Users/esinickin/Desktop/Comparative-demography-project-main/Dracocephalum")

source("R/functions_ipmr.R")
source("R/functions_GAM.R")

result_dir = "results/"

VR_FLM <- readRDS("results/rds/VR_FLM.rds")
state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")
climate_models <- readRDS("results/rds/ARIMA_clim_mods.rds")
CHELSA_data = "data/CHELSA_data.csv"

lag = 24
n_it = 1 #5000

# param/model list 
params <- list(
  surv_mod = VR_FLM$surv,
  s_int = coef(VR_FLM$surv)[1],
  s_stems = coef(VR_FLM$surv)[2],
  s_site_CR = 0,
  s_site_HK = coef(VR_FLM$surv)[3],
  s_site_KS = coef(VR_FLM$surv)[4],
  s_site_RU = coef(VR_FLM$surv)[5],
  s_slope = coef(VR_FLM$surv)[6],
  s_rock = coef(VR_FLM$surv)[7],
  s_soil_depth = coef(VR_FLM$surv)[8],
  
  grow_mod = VR_FLM$growth,
  grow_sd = sd(resid(VR_FLM$growth)),
  g_int = coef(VR_FLM$growth)[1],
  g_stems = coef(VR_FLM$growth)[2],
  g_site_CR = 0,
  g_site_HK = coef(VR_FLM$growth)[3],
  g_site_KS = coef(VR_FLM$growth)[4],
  g_site_RU = coef(VR_FLM$growth)[5],
  g_slope = coef(VR_FLM$growth)[6],
  g_rock = coef(VR_FLM$growth)[7],
  g_soil_depth = coef(VR_FLM$growth)[8],
  
  pflower_mod = VR_FLM$flower_p,
  fp_int = coef(VR_FLM$flower_p)[1],
  fp_stems = coef(VR_FLM$flower_p)[2],
  fp_site_CR = 0,
  fp_site_HK = coef(VR_FLM$flower_p)[3],
  fp_site_KS = coef(VR_FLM$flower_p)[4],
  fp_site_RU = coef(VR_FLM$flower_p)[5],
  fp_slope = coef(VR_FLM$flower_p)[6],
  fp_rock = coef(VR_FLM$flower_p)[7],
  fp_soil_depth = coef(VR_FLM$flower_p)[8],
  
  seedp_mod = VR_FLM$seedp,
  sp_int = coef(VR_FLM$seedp)[1],
  sp_stems = coef(VR_FLM$seedp)[2],
  sp_site_CR = 0,
  sp_site_HK = coef(VR_FLM$seedp)[3],
  sp_site_KS = coef(VR_FLM$seedp)[4],
  sp_site_RU = coef(VR_FLM$seedp)[5],
  sp_slope = coef(VR_FLM$seedp)[6],
  sp_rock = coef(VR_FLM$seedp)[7],
  sp_soil_depth = coef(VR_FLM$seedp)[8],
  
  seedn_mod = VR_FLM$seedn,
  sn_int = coef(VR_FLM$seedn)[1],
  sn_stems = coef(VR_FLM$seedn)[2],
  sn_site_CR = 0,
  sn_site_HK = coef(VR_FLM$seedn)[3],
  sn_site_KS = coef(VR_FLM$seedn)[4],
  sn_site_RU = coef(VR_FLM$seedn)[5],
  sn_slope = coef(VR_FLM$seedn)[6],
  sn_rock = coef(VR_FLM$seedn)[7],
  sn_soil_depth = coef(VR_FLM$seedn)[8],
  
  seed_surv1 = 0.45,  ## Probability of seed being viable at the next census 
  seed_surv2 = 0.089,  ## Probability of viable seed surviving first year in seed bank. 
  ## (i.e. produced in t, surv to census t+1 (with % above), not germinated in t+1 and 
  seed_surv3 = 0.663,  ## Probability of viable seed in seedbank (yr 1) surviving 2nd year
  
  germ_mean = mean(state_independent_variables$est_germination_rate$germ),     
  
  sdl_surv_mod = state_independent_variables$sdl_surv,
  sdl_s_int = lme4::fixef(state_independent_variables$sdl_surv)[1],
  
  sdl_d_int = lme4::fixef(state_independent_variables$sdl_size_d)[1],
  sdl_size_d_sd = sd(resid(state_independent_variables$sdl_size_d))
  )


## Set integration params
L <- 0 
U <- 4.672829 * 1.1   # Max observed size + 10%
n = 100


##------------------------------------------------------------------------------
## Calculate lambda at different co-variate levels
##------------------------------------------------------------------------------

### Loop through different populations and env_param levels 
localities <- c("Cr", "Hk", "Ks", "Ru")
shading <- seq(0,6, length.out = 4)
slope <- seq(0, 50, length.out = 6)
rock <- seq(0, 80, length.out = 6)
soil_depth <- seq(0,10, length.out = 6)

model <- c("ACCESS1", "CESM1", "CMCC", "MIROC5")
scenario <- c("rcp45", "rcp85")

hist <- expand.grid(localities = localities, 
                           shading = shading,
                    slope = 14,
                    rock = 30,
                    soil_depth = 5,
                           time = "hist",
                           scenario = NA,
                           model = NA
) 
fut <- expand.grid(localities = localities, 
                          shading = shading, 
                          slope = 14,
                          rock = 30,
                          soil_depth = 5,
                          time = "future",
                          scenario = scenario,
                          model = model
) 

slope <- expand.grid(localities = "CR", 
                   shading = 3, 
                   slope = slope,
                   rock = 30,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

rock <- expand.grid(localities = localities, 
                   shading = 3, 
                   slope = 14,
                   rock = rock,
                   soil_depth = 5,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

soil <- expand.grid(localities = localities, 
                   shading = 3, 
                   slope = 14,
                   rock = 30,
                   soil_depth = soil_depth,
                   time = "hist",
                   scenario = NA,
                   model = NA
) 

df_env <- rbind(hist, fut, rock, slope, soil) %>% 
  mutate(localities = as.character(localities))

rep <- rep(c(1:nrow(df_env)))


# ### Set up parallel
# cl <- makeCluster(detectCores() - 1)
# clusterExport(cl=cl, c("df_env", "ipm_loop", "run_ipm",
#                        "params", "climate_models",
#                        "U", "L", "n", "n_it", "lag",
#                        "sampling_env", "FLM_clim_predict"))
# 
# clusterEvalQ(cl, c(library("ipmr"), library("dplyr"), library("forecast")))
# 
# df <- parLapply(cl,
#                   as.list(rep),
#                   function(x) tryCatch(ipm_loop(i = x, df_env = df_env,
#                                                 params = params,
#                                                 climate_models = climate_models,
#                                                 n_it = n_it,
#                                                 U = U, L = L, n = n,
#                                                 save = T),
#                                        error = function(e) NULL)) %>%
#   bind_rows()
# 
# stopCluster(cl)
# 
# write.csv(df, file = "results/overview_lambda_env_levels.csv",
#           row.names = F)


i = 1 
# not sure what i is in the code above ====> Sanne: i is used to 
# loop through every entry of "rep", so that we loop through all the rows 
# of the df_env (in other words, all the different combinations of 
# environmental variables)

# df <- ipm_loop(i = i, df_env = df_env,
#          params = params,
#          climate_models = climate_models,
#          n_it = n_it,
#          U = U, L = L, n = n,
#          save = T)


###################################################################
# SENSITIVTIY ANALYSES ####################################

# We want to calculate the sensitivities of lambda to climatic drivers precip and pet
# We use the method adopted from Morris et al. 2020 where they calculated the "scaled" sensitivities
# which is calculated like that: Sens_to_precip = abs((lambda_max - lambda_min) / ((max.precip - min.precip)/sd.precip))
# for that we need to get lambda_max and lambda_min, where the climatic driver in question (in this example precip) is first set to its maximum value and then lambda is calculated, and then set to its minimum value and then get lambda
# But we want 2 types of sensitivities, once without covariation and once with
# WITHOUT COVARIATION:
# means that the other climatic driver (in this example pet) is set to its mean value
# WITH COVARIATION:
# means that the other climatic driver is set to its observed value when precip was at its maximum (then get lambda_max) and when it was at its minimum (to get lambda_min)

## NOTE FROM SANNE ---------------------------------------------------------------
# Lets say you'll take simulations based on historical variables with the
# following "other" environmental variables:
df_env[i,]
# Then the climate formatting becomes quite easy, as long as we go straight to the
# ipm_run, instead of the ipm_loop function:
clim_mod <- list(pet = climate_models$clim_hist_model$Cr.pet_scaled,
                 pr = climate_models$clim_hist_model$Cr.pr_scaled,
                 tas = climate_models$clim_hist_model$Cr.tas_scaled)

clim_sim <- lapply(clim_mod, function(x)
  simulate(x, nsim = ((n_it * 12) + (3*lag))) %>%
    ts(., start= c(2023,1), frequency = 12))

# I know from your comments below you do want to do this for all populations and then get averages
# but with changing "i" and selecting the appropriate models for clim_mod that should be easy
# END NOTE ------------------------------------------------------------------------------
# 

# 1) Sens to precipitation covariation ########
# here we want to calculate the sensitivities assuming covariation
# which means we want to set the other climatic driver pet to its observed value when precipitation was at its maximum and then later at its minimum
# but here I thought I could do it with the ipm_loop() function where I set the year to the year where precip was max or min

clim_data <- read.csv(CHELSA_data) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled", "tas_scaled"), 
                        response_t1 = T,
                        lag = 24)

clim_df <- clim_data %>% select(year_t0, population, pr_scaledcovar, pet_scaledcovar, tas_scaledcovar) %>%
  rename(localities = population) %>% 
  rowwise() %>%
  mutate(pr_mean = mean(pr_scaledcovar[,1:12]),
         pet_mean = mean(pet_scaledcovar[,1:12]),
         tas_mean = mean(tas_scaledcovar[,1:12]))


##NOTE FROM SANNE: --------------------------------------------------------------
## if you only want a single max/min year, you will have to change the 
## mutate() function above so that the mean is only taken from the first 12 columns
## of the "_scaledcovar" variables
##END NOTE ---------------------------------------------------------------------


# Find the average max and min year for pr
print(paste0("averages maximum pr = ", clim_df$year_t0[which.max(clim_df$pr_mean)]))
print(paste0("averages minimum pr = ", clim_df$year_t0[which.min(clim_df$pr_mean)]))

pr_max <- which.max(clim_df$pr_mean)
pr_min <- which.min(clim_df$pr_mean)

clim_pr_max <- clim_sim
clim_pr_max$pr <- ts(rep(clim_df$pr_scaledcovar[pr_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_pr_max$pet <- ts(rep(clim_df$pet_scaledcovar[pr_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_pr_max$tas <- ts(rep(clim_df$tas_scaledcovar[pr_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)

clim_pr_min <- clim_sim
clim_pr_min$pr <- ts(rep(clim_df$pr_scaledcovar[pr_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_pr_min$pet <- ts(rep(clim_df$pet_scaledcovar[pr_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_pr_min$tas <- ts(rep(clim_df$tas_scaledcovar[pr_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)


# environmental params
env_params_pr_max <- append(
  clim_pr_max,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  ))

env_params_pr_min <- append(
  clim_pr_min,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  )) 

# IPM loop
# resample 100 times to get uncertainties around lambda
# first get sensitivities for each population separately and then at the end average across populations

# empty list for sensitivities for each pop
CR.sens.precip.cov=NULL
HK.sens.precip.cov=NULL
KS.sens.precip.cov=NULL
RU.sens.precip.cov=NULL

# also calculate log ratios
#abs(log(lambda(mpm.max)/lambda(mpm.min)))
CR.sens.precip.cov_l_ratio=NULL
HK.sens.precip.cov_l_ratio=NULL
KS.sens.precip.cov_l_ratio=NULL
RU.sens.precip.cov_l_ratio=NULL

for(u in 1:100){
  
  # pop CR
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  CR.sens.precip.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  CR.sens.precip.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop HK
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  HK.sens.precip.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  HK.sens.precip.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop KS
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  KS.sens.precip.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  KS.sens.precip.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop RU
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  RU.sens.precip.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  RU.sens.precip.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
}

sens.precip.cov=data.frame(CR=CR.sens.precip.cov,
                           HK=HK.sens.precip.cov,
                           KS=KS.sens.precip.cov,
                           RU=RU.sens.precip.cov)

sens.precip.cov$mean.sens=rowMeans(sens.precip.cov)


sens.precip.cov_l_ratio=data.frame(CR=CR.sens.precip.cov_l_ratio,
                           HK=HK.sens.precip.cov_l_ratio,
                           KS=KS.sens.precip.cov_l_ratio,
                           RU=RU.sens.precip.cov_l_ratio)

sens.precip.cov_l_ratio$l_ratio=rowMeans(sens.precip.cov_l_ratio)


# 2) Sens to precipitation no covariation ###############################
# meaning that we assume that the other climatic driver is constant at its mean

clim_pr_max <- clim_sim
clim_pr_max$pr <- ts(rep(max(clim_pr_max$pr), length(clim_pr_max$pr)), start= c(2023,1), frequency = 12)
clim_pr_max$pet <- ts(rep(0.000001, length(clim_pr_max$pet)), start= c(2023,1), frequency = 12)
clim_pr_max$tas <- ts(rep(0.000001, length(clim_pr_max$tas)), start= c(2023,1), frequency = 12)

clim_pr_min <- clim_sim
clim_pr_min$pr <- ts(rep(min(clim_pr_min$pr), length(clim_pr_min$pr)), start= c(2023,1), frequency = 12)
clim_pr_min$pet <- ts(rep(0.000001, length(clim_pr_min$pet)), start= c(2023,1), frequency = 12)
clim_pr_min$tas <- ts(rep(0.000001, length(clim_pr_min$tas)), start= c(2023,1), frequency = 12)


# environmental params
env_params_pr_max <- append(
  clim_pr_max,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  ))

env_params_pr_min <- append(
  clim_pr_min,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  )) 


# IPM loop

# empty list for sensitivities for each pop
CR.sens.precip.no.cov=NULL
HK.sens.precip.no.cov=NULL
KS.sens.precip.no.cov=NULL
RU.sens.precip.no.cov=NULL

# log ratios
CR.sens.precip.no.cov_l_ratio=NULL
HK.sens.precip.no.cov_l_ratio=NULL
KS.sens.precip.no.cov_l_ratio=NULL
RU.sens.precip.no.cov_l_ratio=NULL

for(u in 1:100){
  # pop CR
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  CR.sens.precip.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  CR.sens.precip.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop HK
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  HK.sens.precip.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  HK.sens.precip.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop KS
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  KS.sens.precip.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  KS.sens.precip.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  

  # pop RU
  ipm_max <- run_ipm(params = params, env_params = env_params_pr_max, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_pr_min, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  RU.sens.precip.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_pr_max$pr, na.rm = T) - min(clim_pr_min$pr, na.rm = T))/1))
  
  RU.sens.precip.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
}


sens.precip.no.cov=data.frame(CR=CR.sens.precip.no.cov,
                           HK=HK.sens.precip.no.cov,
                           KS=KS.sens.precip.no.cov,
                           RU=RU.sens.precip.no.cov)

sens.precip.no.cov$mean.sens=rowMeans(sens.precip.no.cov)


sens.precip.no.cov_l_ratio=data.frame(CR=CR.sens.precip.no.cov_l_ratio,
                                   HK=HK.sens.precip.no.cov_l_ratio,
                                   KS=KS.sens.precip.no.cov_l_ratio,
                                   RU=RU.sens.precip.no.cov_l_ratio)

sens.precip.no.cov_l_ratio$l_ratio=rowMeans(sens.precip.no.cov_l_ratio)


# 3) Sens to temperature covariation#####################

# Find the average max and min year for pr
print(paste0("averages maximum tas = ", clim_df$year_t0[which.max(clim_df$tas_mean)]))
print(paste0("averages minimum tas = ", clim_df$year_t0[which.min(clim_df$tas_mean)]))

tas_max <- which.max(clim_df$tas_mean)
tas_min <- which.min(clim_df$tas_mean)

clim_tas_max <- clim_sim
clim_tas_max$pr <- ts(rep(clim_df$pr_scaledcovar[tas_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_tas_max$pet <- ts(rep(clim_df$pet_scaledcovar[tas_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_tas_max$tas <- ts(rep(clim_df$tas_scaledcovar[tas_max,-1], n_it*1.5), start= c(2023,6), frequency = 12)

clim_tas_min <- clim_sim
clim_tas_min$pr <- ts(rep(clim_df$pr_scaledcovar[tas_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_tas_min$pet <- ts(rep(clim_df$pet_scaledcovar[tas_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)
clim_tas_min$tas <- ts(rep(clim_df$tas_scaledcovar[tas_min,-1], n_it*1.5), start= c(2023,6), frequency = 12)


# environmental params
env_params_tas_max <- append(
  clim_tas_max,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  ))

env_params_tas_min <- append(
  clim_tas_min,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  )) 


CR.sens.temp.cov=NULL
HK.sens.temp.cov=NULL
KS.sens.temp.cov=NULL
RU.sens.temp.cov=NULL

CR.sens.temp.cov_l_ratio=NULL
HK.sens.temp.cov_l_ratio=NULL
KS.sens.temp.cov_l_ratio=NULL
RU.sens.temp.cov_l_ratio=NULL


for(u in 1:100){
  
  # pop CR
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  CR.sens.temp.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  CR.sens.temp.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop HK
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  HK.sens.temp.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  HK.sens.temp.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop KS
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  KS.sens.temp.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  KS.sens.temp.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop RU
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  RU.sens.temp.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  RU.sens.temp.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
}


sens.temp.cov=data.frame(CR=CR.sens.temp.cov,
                              HK=HK.sens.temp.cov,
                              KS=KS.sens.temp.cov,
                              RU=RU.sens.temp.cov)

sens.temp.cov$mean.sens=rowMeans(sens.temp.cov)

sens.temp.cov_l_ratio=data.frame(CR=CR.sens.temp.cov_l_ratio,
                         HK=HK.sens.temp.cov_l_ratio,
                         KS=KS.sens.temp.cov_l_ratio,
                         RU=RU.sens.temp.cov_l_ratio)

sens.temp.cov_l_ratio$l_ratio=rowMeans(sens.temp.cov_l_ratio)



# 4) Sens to temperature no cov ###########################

clim_tas_max <- clim_sim

clim_tas_max$tas <- ts(rep(max(clim_tas_max$tas), length(clim_tas_max$tas)), start= c(2023,1), frequency = 12)
clim_tas_max$pet <- ts(rep(0.000001, length(clim_tas_max$pet)), start= c(2023,1), frequency = 12)
clim_tas_max$pr <- ts(rep(0.000001, length(clim_tas_max$pr)), start= c(2023,1), frequency = 12)


clim_tas_min <- clim_sim
clim_tas_min$tas <- ts(rep(min(clim_tas_min$tas), length(clim_tas_min$tas)), start= c(2023,1), frequency = 12)
clim_tas_min$pet <- ts(rep(0.000001, length(clim_tas_min$pet)), start= c(2023,1), frequency = 12)
clim_tas_min$pr <- ts(rep(0.000001, length(clim_tas_min$pr)), start= c(2023,1), frequency = 12)


# environmental params
env_params_tas_max <- append(
  clim_tas_max,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  ))

env_params_tas_min <- append(
  clim_tas_min,
  list(
    lags = lag,
    shading = df_env$shading[i],
    slope = df_env$slope[i],
    rock = df_env$rock[i],
    soil_depth = df_env$soil_depth[i]
  )) 


CR.sens.temp.no.cov=NULL
HK.sens.temp.no.cov=NULL
KS.sens.temp.no.cov=NULL
RU.sens.temp.no.cov=NULL

CR.sens.temp.no.cov_l_ratio=NULL
HK.sens.temp.no.cov_l_ratio=NULL
KS.sens.temp.no.cov_l_ratio=NULL
RU.sens.temp.no.cov_l_ratio=NULL

for(u in 1:100){
  
  # pop CR
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "CR", 
                     n_it = n_it, U = U, L = L, n = n)
  
  CR.sens.temp.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  
  CR.sens.temp.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  
  # pop HK
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "HK", 
                     n_it = n_it, U = U, L = L, n = n)
  
  HK.sens.temp.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  HK.sens.temp.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop KS
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "KS", 
                     n_it = n_it, U = U, L = L, n = n)
  
  KS.sens.temp.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  KS.sens.temp.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
  # pop RU
  ipm_max <- run_ipm(params = params, env_params = env_params_tas_max, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  ipm_min <- run_ipm(params = params, env_params = env_params_tas_min, 
                     locality = "RU", 
                     n_it = n_it, U = U, L = L, n = n)
  
  RU.sens.temp.no.cov[u] = abs(lambda(ipm_max, log = F) - lambda(ipm_min, log = F)/((max(clim_tas_max$tas, na.rm = T) - min(clim_tas_min$tas, na.rm = T))/1))
  
  RU.sens.temp.no.cov_l_ratio[u]=abs(log(lambda(ipm_max,log=F)/lambda(ipm_min,log=F)))
  
}


sens.temp.no.cov=data.frame(CR=CR.sens.temp.no.cov,
                         HK=HK.sens.temp.no.cov,
                         KS=KS.sens.temp.no.cov,
                         RU=RU.sens.temp.no.cov)

sens.temp.no.cov$mean.sens=rowMeans(sens.temp.no.cov)

sens.temp.no.cov_l_ratio=data.frame(CR=CR.sens.temp.no.cov_l_ratio,
                            HK=HK.sens.temp.no.cov_l_ratio,
                            KS=KS.sens.temp.no.cov_l_ratio,
                            RU=RU.sens.temp.no.cov_l_ratio)

sens.temp.no.cov_l_ratio$l_ratio=rowMeans(sens.temp.no.cov_l_ratio)


# SAVE OUTPUT #######################
output_df=data.frame(precip.no.cov=sens.precip.no.cov$mean.sens,
                   precip.cov=sens.precip.cov$mean.sens,
                   temp.no.cov=sens.temp.no.cov$mean.sens,
                   temp.cov=sens.temp.cov$mean.sens)

#write.csv(output_df,"SensDraco.csv",row.names = F)

# re-arrange df 


Sens_Dracocephalum=data.frame(species="Dracocephalum austriacum",
                              study.doi="Evers et al. in Prep",
                              year.of.publication="2024",
                              group="Plants",
                              continent="Europe",
                              driver=rep(c("temperature","precipitation"),each=200),
                              driver.type="C",
                              stage.age="all",
                              vital.rates="all",
                              sens=c(output_df$temp.no.cov,output_df$temp.cov,
                                     output_df$precip.no.cov,output_df$precip.cov),
                              cov=rep(c(0,1),each=100),
                              mat=2, # Source: Dostálek & Münzbergová 2013
                              n.vr=6, # number of vital rates models
                              n.pam=59, #length of params df (number of parameters)
                              dens=0,
                              biotic_interactions=0,
                              lambda.sim=0,
                              study.length=16,
                              
                              l_ratio=c(sens.temp.no.cov_l_ratio$l_ratio,
                                        sens.temp.cov_l_ratio$l_ratio,
                                        sens.precip.no.cov_l_ratio$l_ratio,
                                        sens.precip.cov_l_ratio$l_ratio))

w#rite.csv(Sens_Dracocephalum, "Sensitivities_Dracocephalum.csv",row.names = F)

