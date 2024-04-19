## This script calculates lambda for different environmental variable levels
## under asymptotic behaviour (i.e. long run simulations)


# This script was created by Sanne Evers
# and edited by Esin Ickin (30.01.2024)

## after line 142, the sensitivity analyses start



#-------------------------------------
# See my notes below! L171
#-------------------------------------


rm(list=ls())

setwd("/Users/esinickin/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/Dracocephalum")

library(dplyr)
library(lme4)
library(ipmr)
library(mgcv)
library(parallel)
#library(forecast)

source("R/functions_ipmr.R")
source("R/functions_GAM.R")

result_dir = "results/"

VR_FLM <- readRDS("results/rds/VR_FLM.rds")
state_independent_variables <- readRDS("results/rds/state_independent_VR.rds")
CHELSA_data = "data/CHELSA_data.csv"
data_for_modeling = "data/Dracocephalum_with_vital_rates.csv"


n_it = 10 #000  Number of iterations for the IPM.

# param/model list 
params <- list(
  surv_mod = VR_FLM$surv,
  s_int = coef(VR_FLM$surv)[1],
  s_stems = coef(VR_FLM$surv)[2],
  s_site_CR = 0,
  s_site_HK = coef(VR_FLM$surv)[3],
  s_site_KS = coef(VR_FLM$surv)[4],
  s_site_RU = coef(VR_FLM$surv)[5],
  
  grow_mod = VR_FLM$growth,
  g_int = coef(VR_FLM$growth)[1],
  g_stems = coef(VR_FLM$growth)[2],
  g_site_CR = 0,
  g_site_HK = coef(VR_FLM$growth)[3],
  g_site_KS = coef(VR_FLM$growth)[4],
  g_site_RU = coef(VR_FLM$growth)[5],
  grow_sd = sd(resid(VR_FLM$growth)),
  
  pflower_mod = VR_FLM$flower_p,
  fp_int = coef(VR_FLM$flower_p)[1],
  fp_stems = coef(VR_FLM$flower_p)[2],
  fp_site_CR = 0,
  fp_site_HK = coef(VR_FLM$flower_p)[3],
  fp_site_KS = coef(VR_FLM$flower_p)[4],
  fp_site_RU = coef(VR_FLM$flower_p)[5],
  fp_slope = coef(VR_FLM$flower_p)[6],
  
  pabort_mod = VR_FLM$abort_p,
  ab_int = coef(VR_FLM$abort_p)[1],
  ab_stems = coef(VR_FLM$abort_p)[2],
  ab_site_CR = 0,
  ab_site_HK = coef(VR_FLM$abort_p)[3],
  ab_site_KS = coef(VR_FLM$abort_p)[4],
  ab_site_RU = coef(VR_FLM$abort_p)[5],
  
  nseed_mod = VR_FLM$n_seeds,
  ns_int = coef(VR_FLM$n_seeds)[1],
  ns_stems = coef(VR_FLM$n_seeds)[2],
  ns_site_CR = 0,
  ns_site_HK = coef(VR_FLM$n_seeds)[3],
  ns_site_KS = coef(VR_FLM$n_seeds)[4],
  ns_site_RU = coef(VR_FLM$n_seeds)[5],
  
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
L <- min(VR_FLM$growth$model$ln_stems_t0, na.rm = T)
U <- max(VR_FLM$growth$model$ln_stems_t0, na.rm = T) * 1.1
n = 100


##------------------------------------------------------------------------------
## Calculate lambda at different co-variate levels
##------------------------------------------------------------------------------

### Climate dataframe

demo_data <- read.csv(data_for_modeling)
clim_data <- read.csv(CHELSA_data) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = 24)
#rm(clim_df)
clim_df <- clim_data %>% select(year_t0, population, pr_scaledcovar, pet_scaledcovar) %>%
  rename(localities = population)

climate_df = clim_df

### Values to loop through. Populations and shading levels
#localities <- c("CR", "HK", "KS", "RU")
#shading <- seq(0,6, length.out = 4)
#years <- c(min(clim_df$year_t0):max(clim_df$year_t0))

#loop_df <- expand.grid(loc = localities, 
#                       shading = shading, 
#                       year = years)

### Run as array job
#i = 110    # set row from loop_df.

#df <- ipm_loop(yr = loop_df$year[i], 
#               loc = loop_df$loc[i],  
#               shading = loop_df$shading[i], 
#               params = params,
#               climate_df = clim_df,
#               n_it = n_it,
#               U = U, L = L, n = n) 

# write.csv(df, file = file.path(result_dir, paste0("ipm_stoch_", i, ".csv")), row.names = F)



# SENSITIVTIY ANALYSES ####################################

# We want to calculate the sensitivities of lambda to climatic drivers precip and pet, right?
# We use the method adopted from Morris et al. 2020 where they calculated the "scaled" sensitivities
# which is calculated like that: Sens_to_precip = abs((lambda_max - lambda_min) / ((max.precip - min.precip)/sd.precip))
# for that we need to get lambda_max and lambda_min, where the climatic driver in question (in this example precip) is first set to its maximum value and then lambda is calculated, and then set to its minimum value and then get lambda
# But we want 2 types of sensitivities, once without covariation and once with
# WITHOUT COVARIATION:
# means that the other climatic driver (in this example pet) is set to its mean value
# WITH COVARIATION:
# means that the other climatic driver is set to its observed value when precip was at its maximum (then get lambda_max) and when it was at its minimum (to get lambda_min)



# test just for me to see how env_params looks like
#env_params <- list(
# pet = clim_df$pet_scaledcovar[which(clim_df$year_t0 == 2003 & clim_df$localities == "HK"),],
#  pr = clim_df$pr_scaledcovar[which(clim_df$year_t0 == 2003 & clim_df$localities == "HK"),],
#  lags = 24,
#  shading = shading
#)


# 1) Sens to precipitation covariation ########
# here we want to calculate the sensitivities assuming covariation
# which means we want to set the other climatic driver pet to its observed value when precipitation was at its maximum and then later at its minimum
# but here I thought I could do it with the ipm_loop() function where I set the year to the year where precip was max or min

# now first get the covariation
clim_pr_cov <- read.csv(CHELSA_data) 

# pet when pr was at its maximum
pet_when_pr_max=clim_pr_cov$pet_scaled[which(clim_pr_cov$pr_scaled==max(clim_pr_cov$pr_scaled,na.rm=T))][1]

clim_pr_max_cov <- read.csv(CHELSA_data) %>%
  mutate(across(contains("pet"), ~ pet_when_pr_max),
         across(contains("pr"), ~ max(.x, na.rm = T))) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = 24) %>% 
  select(year_t0, population, pr_scaledcovar, pet_scaledcovar) %>%
  rename(localities = population)


# pet when pr was at its minimum
pet_when_pr_min=clim_pr_cov$pet_scaled[which(clim_pr_cov$pr_scaled==min(clim_pr_cov$pr_scaled,na.rm=T))][1]

clim_pr_min_cov <- read.csv(CHELSA_data) %>%
  mutate(across(contains("pet"), ~ pet_when_pr_min),
         across(contains("pr"), ~ min(.x, na.rm = T))) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = 24) %>% 
  select(year_t0, population, pr_scaledcovar, pet_scaledcovar) %>%
  rename(localities = population)

# loop through the years
# and loop through all the locations
# and then average it
loc <- c("CR","HK","KS","RU")
year <- c(2003:2021)

loop_df <- expand.grid(loc = loc, 
                       year = year)

# start sens analysis
sens.precip.cov=NULL

for(u in 1:50){ # to get uncertainties
  lambda_max=NULL
  mean_lambda_max=NULL
  lambda_min=NULL
  mean_lambda_min=NULL
  for(i in 1:nrow(loop_df)){ # to get averages across years and locations
    ipm_max <- ipm_loop(yr = loop_df$year[i], 
                        loc = loop_df$loc[i],  
                        shading = 3, # average shading level
                        params = params,
                        climate_df = clim_pr_max_cov,
                        n_it = n_it,
                        U = U, L = L, n = n)
    
    lambda_max[i] <- exp(ipm_max$lambda)
    
    ipm_min <- ipm_loop(yr = loop_df$year[i], 
                        loc = loop_df$loc[i],  
                        shading = 3, # average shading level
                        params = params,
                        climate_df = clim_pr_min_cov,
                        n_it = n_it,
                        U = U, L = L, n = n)
    
    # lambdas are log transformed, so need to backtransform them
    lambda_min[i] <- exp(ipm_min$lambda)
  }
  
  mean_lambda_max[u] <- mean(lambda_max)
  mean_lambda_min[u] <- mean(lambda_min)
  
  # now we could get lambda and calculate the "scaled sensitivities" according to Morris et al. 2020
  # which goes like this: S = abs((lambda_max - lambda_min) / ((max.precip - min.precip)/sd.precip))
  sens.precip.cov[u] = abs((mean_lambda_max[u] - mean_lambda_min[u]) / ((max(clim_pr_max_cov[,3], na.rm = T) - min(clim_pr_min_cov[,3], na.rm = T))/1))
  
  
}

hist(sens.precip.cov)



# 2) Sens to PET covariation #############
# meaning we set the other climatic driver to its observed value when PET was at its max and later min
# this is going to be the similar to the code below 2)

# pr when pet was at its maximum
pr_when_pet_max=clim_pr_cov$pr_scaled[which(clim_pr_cov$pet_scaled==max(clim_pr_cov$pet_scaled,na.rm=T))][1]

clim_pet_max_cov <- read.csv(CHELSA_data) %>%
  mutate(across(contains("pet"), ~ max(.x, na.rm = T)),
         across(contains("pr"), ~ pr_when_pet_max)) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = 24) %>% 
  select(year_t0, population, pr_scaledcovar, pet_scaledcovar) %>%
  rename(localities = population)


# pr when pet was at its minimum
pr_when_pet_min=clim_pr_cov$pr_scaled[which(clim_pr_cov$pet_scaled==min(clim_pr_cov$pet_scaled,na.rm=T))][1]

clim_pet_min_cov <- read.csv(CHELSA_data) %>%
  mutate(across(contains("pet"), ~ min(.x, na.rm = T)),
         across(contains("pr"), ~ pr_when_pet_min)) %>%
  climate_wider_for_gam(clim_data = ., 
                        variables = c("pr_scaled", "pet_scaled"), 
                        demo_data = demo_data, 
                        response_t1 = T,
                        lag = 24) %>% 
  select(year_t0, population, pr_scaledcovar, pet_scaledcovar) %>%
  rename(localities = population)


# start sens analysis
sens.pet.cov=NULL

for(u in 1:50){ # to get uncertainties
  lambda_max=NULL
  mean_lambda_max=NULL
  lambda_min=NULL
  mean_lambda_min=NULL
  for(i in 1:nrow(loop_df)){ # to get averages across years and locations
    ipm_max <- ipm_loop(yr = loop_df$year[i], 
                        loc = loop_df$loc[i],  
                        shading = 3, # average shading level
                        params = params,
                        climate_df = clim_pet_max_cov,
                        n_it = n_it,
                        U = U, L = L, n = n)
    
    lambda_max[i] <- exp(ipm_max$lambda)
    
    ipm_min <- ipm_loop(yr = loop_df$year[i], 
                        loc = loop_df$loc[i],  
                        shading = 3, # average shading level
                        params = params,
                        climate_df = clim_pet_min_cov,
                        n_it = n_it,
                        U = U, L = L, n = n)
    
    # lambdas are log transformed, so need to backtransform them
    lambda_min[i] <- exp(ipm_min$lambda)
  }
  
  mean_lambda_max[u] <- mean(lambda_max)
  mean_lambda_min[u] <- mean(lambda_min)
  
  # now we could get lambda and calculate the "scaled sensitivities" according to Morris et al. 2020
  # which goes like this: S = abs((lambda_max - lambda_min) / ((max.precip - min.precip)/sd.precip))
  sens.pet.cov[u] = abs((mean_lambda_max[u] - mean_lambda_min[u]) / ((max(clim_pet_max_cov[,3], na.rm = T) - min(clim_pet_min_cov[,3], na.rm = T))/1))
  
}

hist(sens.pet.cov)


# 3) Save Output ####################################

output_df=data.frame(sens.pet.cov=sens.pet.cov,sens.precip.cov=sens.precip.cov)

write.csv(output_df,"SensCov.csv",row.names = F)