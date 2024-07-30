############################################################################
#
# This script calculates the min, max, and mean of each covariate
# and projects the dewy-pine populations using an individual-based
# model (IBM) under two scenarios for each focal covariate related
# to each focal vital rate (reproductive and non-reproductive survival, and
# reproduction probability): 
# 
# (1) Focal covariate = min value, 
#     other covariates of focal vital rate = observed when focal covariate = min,
#     covariates of other vital rates = mean
# (2) Focal covariate = max value, 
#     other covariates of focal vital rate = observed when focal covariate = max,
#     covariates of other vital rates = mean
#
# Author: Eva Conquet
#
###########################################################################

###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())


## 1.2. Loading libraries ----
# -----------------------

library(lubridate)
library(dplyr)
library(metRology)
library(mgcv)
library(crch)
library(snowfall)


## 1.3. Loading data ----
# ------------------

# Dewy-pine data
droso = read.csv("Data/droso_natural.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_NoDormancyLoss.csv")

seeds_per_flower = 9.8


# Correction factors
seedSurv_df = data.frame(seedSurv = c(0.45, 0.45, 0.35, 0.35, 0.33, 
                                      0.84, 0.84, 0.82, 0.8, 0.82),
                         TSF = rep(c("0", "1", "2", "3", "4"), 2),
                         LS = c(rep("HLS", 5), rep("LLS", 5)))
seedSurv_df = seedSurv_df[which(seedSurv_df$LS == "LLS"), ]


# Vital-rate models
load("Data/Survival_GAM_Undisturbed.RData")
load("Data/Growth_GAM_Undisturbed.RData")
load("Data/FloweringProb_GAM_Undisturbed.RData")
load("Data/NbFlowers_GAM_Undisturbed.RData")
load("Data/SeedlingSize_GAM_Undisturbed.RData")


# Average density per square for covariate standardization
nbSquares = aggregate(quadratID ~ time + site, data = droso, FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, data = droso, function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, data = density_per_square, function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Year- and population-specific climatic variables for covariate standardization
summerT_timeSeries = aggregate(summerT_unscaled ~ time + site, 
                               data = droso, mean)
prevwinterT_timeSeries = aggregate(prevwinterT_unscaled ~ time + site, 
                                   data = droso, mean)
fallR_timeSeries = aggregate(fallR_unscaled ~ time + site, 
                             data = droso, mean)
prevfallR_timeSeries = aggregate(prevfallR_unscaled ~ time + site, 
                                 data = droso, mean)




###########################################################################
#
# 2. Building vital-rate functions ----
#
###########################################################################

## 2.1. Survival ----
# --------------

non_repro_survival_function = function(size, abLarge, fallR, summerT, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso$size_unscaled, na.rm = T)) / (2 * sd(droso$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))
  summerT_scaled = (summerT - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Calculate survival
  survival = predict(surv_natural, newdata = data.frame(size = size_scaled,
                                                   abLarge = abLarge_scaled,
                                                   fallR = fallR_scaled,
                                                   summerT = summerT_scaled,
                                                   TSFcont = TSFcont_scaled,
                                                   time = year,
                                                   site = population), type = "response")
  
  return(survival)
}


repro_survival_function = function(size, abLarge, fallR, summerT, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso$size_unscaled, na.rm = T)) / (2 * sd(droso$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))
  summerT_scaled = (summerT - mean(summerT_timeSeries$summerT_unscaled, na.rm = T)) / (2 * sd(summerT_timeSeries$summerT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Calculate survival
  survival = predict(surv_natural, newdata = data.frame(size = size_scaled,
                                                   abLarge = abLarge_scaled,
                                                   fallR = fallR_scaled,
                                                   summerT = summerT_scaled,
                                                   TSFcont = TSFcont_scaled,
                                                   time = year,
                                                   site = population), type = "response")
  
  return(survival)
}


## 3.2. Growth ----
# ------------

growth_function = function(size, abLarge, fallR, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso$size_unscaled, na.rm = T)) / (2 * sd(droso$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  fallR_scaled = (fallR - mean(fallR_timeSeries$fallR_unscaled, na.rm = T)) / (2 * sd(fallR_timeSeries$fallR_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Get parameters from growth model (mean, sd, and degrees of freedom)
  growth_mean = predict(growth_natural, newdata = data.frame(size = size_scaled,
                                                        abLarge = abLarge_scaled,
                                                        fallR = fallR_scaled,
                                                        TSFcont = TSFcont_scaled,
                                                        time = year,
                                                        site = population), type = "response")
  
  growth_sd = family(growth_natural)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  growth_df = family(growth_natural)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  # growth = dnorm(sizeNext, mean = growth_mean,
  #                sd = sd(mgcv::residuals.gam(growth_natural)))
  
  return(list(mean = growth_mean, 
              sd = growth_sd, 
              df = growth_df))
}


## 3.3. Seedling sizes ----
# --------------------

seedling_size_function = function(abLarge, prevwinterT, TSFcont, year, population){
  
  # Standardize covariates
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Get parameters from seedling size model (mean, sd, and degrees of freedom)
  seedling_size_mean = as.numeric(predict(seedlingSize_natural, newdata = data.frame(abLarge = abLarge_scaled,
                                                                                prevwinterT = prevwinterT_scaled,
                                                                                TSFcont = TSFcont_scaled,
                                                                                time = year,
                                                                                site = population), type = "response"))
  
  seedling_size_sd = family(seedlingSize_natural)$getTheta(trans = T)[2] # Using trans = T, the second value of theta is sigma, the standard deviation (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  seedling_size_df = family(seedlingSize_natural)$getTheta(trans = T)[1] # Using trans = T, the first value of theta is nu, the degrees of freedom (https://stats.stackexchange.com/questions/550339/extracting-the-degrees-of-freedom-of-t-distribution-of-a-gam)
  
  # seedling_size = dnorm(sizeNext, mean = seedling_size_mean,
  #                       sd = sd(mgcv::residuals.gam(seedlingSize_natural)))
  
  return(list(mean = seedling_size_mean, 
              sd = seedling_size_sd, 
              df = seedling_size_df))
  
}


## 3.4. Flowering probability ----
# ---------------------------

flowering_function = function(size, abLarge, prevwinterT, prevfallR, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso$size_unscaled, na.rm = T)) / (2 * sd(droso$size_unscaled, na.rm = T))
  abLarge_scaled = (abLarge - mean(yearly_density_per_square$abLarge_unscaled, na.rm = T)) / (2 * sd(yearly_density_per_square$abLarge_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  prevfallR_scaled = (prevfallR - mean(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T)) / (2 * sd(prevfallR_timeSeries$prevfallR_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Calculate flowering probability
  flowering = predict(flowering_natural, newdata = data.frame(size = size_scaled,
                                                         abLarge = abLarge_scaled,
                                                         prevwinterT = prevwinterT_scaled,
                                                         prevfallR = prevfallR_scaled,
                                                         TSFcont = TSFcont_scaled,
                                                         time = year,
                                                         site = population), type = "response")
  
  return(flowering)
}


## 3.5. Number of flowers ----
# -----------------------

nbFlowers_function = function(size, prevwinterT, TSFcont, year, population){
  
  # Standardize covariates
  size_scaled = (size - mean(droso$size_unscaled, na.rm = T)) / (2 * sd(droso$size_unscaled, na.rm = T))
  prevwinterT_scaled = (prevwinterT - mean(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T)) / (2 * sd(prevwinterT_timeSeries$prevwinterT_unscaled, na.rm = T))
  TSFcont_scaled = (TSFcont - mean(droso$TSFcont_unscaled, na.rm = T)) / (2 * sd(droso$TSFcont_unscaled, na.rm = T))
  
  # Calculate number of flowers
  nbFlowers = predict(nbFlow_natural, newdata = data.frame(size = size_scaled,
                                                      prevwinterT = prevwinterT_scaled,
                                                      TSFcont = TSFcont_scaled,
                                                      time = year,
                                                      site = population), type = "response")
  
  return(nbFlowers)
}


## 3.6. Immediate germination (goCont) ----
# ------------------------------------

goCont_function = function(TSF){
  
  if(TSF %in% seq(0, 5)){
    
    goCont = droso_seedbank$value[which(droso_seedbank$vital_rate == "goCont" &
                                          droso_seedbank$TSF == TSF)]
  }
  
  else{
    
    goCont = droso_seedbank$value[which(droso_seedbank$vital_rate == "goCont" &
                                          droso_seedbank$TSF == 5)]
    
  }
  
  if(population == "SierraCarboneraY5") goCont = goCont * 0.4
  else if(population == "SierraRetinY5") goCont = goCont * 0.4
  
  return(goCont)
}


## 3.7. Staying in the seed bank (staySB) ----
# ---------------------------------------

staySB_function = function(TSF){
  
  if(TSF %in% seq(0, 5)){
    
    staySB = droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
                                          droso_seedbank$TSF == TSF)]
    
    # staySB = ifelse(droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
    #                                       droso_seedbank$TSF == TSF)] > 0.7, 0.7, droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
    #                                                                                                            droso_seedbank$TSF == TSF)])
    
    
  }
  
  else{
    
    staySB = droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
                                          droso_seedbank$TSF == 5)]
    
    # staySB = ifelse(droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
    #                                              droso_seedbank$TSF == 5)] > 0.7, 0.7, droso_seedbank$value[which(droso_seedbank$vital_rate == "staySB" &
    #                                                                                                                   droso_seedbank$TSF == 5)])
  }
  
  return(staySB)
}


## 3.8. Germinating out of the seed bank (outSB) ----
# ----------------------------------------------

outSB_function = function(TSF, seedSurv){
  
  if(TSF %in% seq(0, 5)){
    
    outSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "outSB" &
                                         droso_seedbank$TSF == TSF)] *
      seedSurv
    
  }
  
  else{
    
    outSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "outSB" &
                                         droso_seedbank$TSF == 5)] *
      seedSurv
  }
  
  if(population == "SierraCarboneraY5") outSB = outSB * 0.4
  else if(population == "SierraRetinY5") outSB = outSB * 0.4
  
  return(outSB)
}


## 3.9. Going to the seed bank (goSB) ----
# -----------------------------------

goSB_function = function(){
  
  goSB = droso_seedbank$value[which(droso_seedbank$vital_rate == "goSB")]
  
  return(goSB)
}




###########################################################################
#
# 4. Get temperature (next summer and previous winter), ----
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
#
###########################################################################

population = "Vertedero"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

## 4.1. Min, max, standard deviation, and mean of each covariate ----
# --------------------------------------------------------------

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# # Size
# # We first calculate the mean size for each year and take the min and max
# # to avoid having multiple years where the min and max occur.
# 
# yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)
# 
# min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
# max_size = max(yearly_mean_size$size_unscaled, na.rm = T)
# 
# sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)
# 
# 
# # TSF
# min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
# max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)
# 
# sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)


## 4.2. Other covariate values when a given covariate is at its min and max ----
# -------------------------------------------------------------------------

## 4.2.1. Min next summer mean max. daily temperature ----
# ---------------------------------------------------

year_min_summerT = unique(summerT_timeSeries$time[which(summerT_timeSeries$summerT_unscaled == min_summerT &
                                                          summerT_timeSeries$site == population)])

# Previous winter mean max. daily temperature
prevwinterT_value_min_summerT = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_summerT &
                                                                                           prevwinterT_timeSeries$site == population)])  

# Next fall cumulative rainfall
fallR_value_min_summerT = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_summerT &
                                                                         fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_min_summerT = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_summerT &
                                                                                     prevfallR_timeSeries$site == population)])

# Density
dens_value_min_summerT = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_summerT &
                                                                            yearly_density_per_square$site == population)]

# Size
size_value_min_summerT = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_summerT)])

# TSF
TSF_value_min_summerT = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_summerT)])


## 4.2.2. Max next summer mean max. daily temperature ----
# ---------------------------------------------------

year_max_summerT = unique(summerT_timeSeries$time[which(summerT_timeSeries$summerT_unscaled == max_summerT &
                                                          summerT_timeSeries$site == population)])

# Previous winter mean max. daily temperature
prevwinterT_value_max_summerT = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_summerT &
                                                                                           prevwinterT_timeSeries$site == population)])  

# Next fall cumulative rainfall
fallR_value_max_summerT = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_summerT &
                                                                         fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_max_summerT = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_summerT &
                                                                                     prevfallR_timeSeries$site == population)])

# Density
dens_value_max_summerT = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_summerT &
                                                                            yearly_density_per_square$site == population)]

# Size
size_value_max_summerT = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_summerT)])

# TSF
TSF_value_max_summerT = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_summerT)])


## 4.2.3. Min previous winter mean max. daily temperature ----
# -------------------------------------------------------

year_min_prevwinterT = unique(prevwinterT_timeSeries$time[which(prevwinterT_timeSeries$prevwinterT_unscaled == min_prevwinterT &
                                                                  prevwinterT_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_min_prevwinterT = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_prevwinterT &
                                                                                   summerT_timeSeries$site == population)])  

# Next fall cumulative rainfall
fallR_value_min_prevwinterT = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_prevwinterT &
                                                                             fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_min_prevwinterT = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_prevwinterT &
                                                                                         prevfallR_timeSeries$site == population)])

# Density
dens_value_min_prevwinterT = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_prevwinterT &
                                                                                yearly_density_per_square$site == population)]

# Size
size_value_min_prevwinterT = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_prevwinterT)])

# TSF
TSF_value_min_prevwinterT = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_prevwinterT)])


## 4.2.4. Max previous winter mean max. daily temperature ----
# -------------------------------------------------------

year_max_prevwinterT = unique(prevwinterT_timeSeries$time[which(prevwinterT_timeSeries$prevwinterT_unscaled == max_prevwinterT &
                                                                  prevwinterT_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_max_prevwinterT = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_prevwinterT &
                                                                                   summerT_timeSeries$site == population)])  

# Next fall cumulative rainfall
fallR_value_max_prevwinterT = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_prevwinterT &
                                                                             fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_max_prevwinterT = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_prevwinterT &
                                                                                         prevfallR_timeSeries$site == population)])

# Density
dens_value_max_prevwinterT = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_prevwinterT &
                                                                                yearly_density_per_square$site == population)]

# Size
size_value_max_prevwinterT = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_prevwinterT)])

# TSF
TSF_value_max_prevwinterT = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_prevwinterT)])


## 4.2.5. Min next fall cumulative rainfall ----
# -----------------------------------------

year_min_fallR = unique(fallR_timeSeries$time[which(fallR_timeSeries$fallR_unscaled == min_fallR &
                                                      fallR_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_min_fallR = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_fallR &
                                                                             summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_min_fallR = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_fallR &
                                                                                         prevwinterT_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_min_fallR = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_fallR &
                                                                                   prevfallR_timeSeries$site == population)])

# Density
dens_value_min_fallR = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_fallR &
                                                                          yearly_density_per_square$site == population)]

# Size
size_value_min_fallR = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_fallR)])

# TSF
TSF_value_min_fallR = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_fallR)])


## 4.2.6. Max next fall cumulative rainfall ----
# -----------------------------------------

year_max_fallR = unique(fallR_timeSeries$time[which(fallR_timeSeries$fallR_unscaled == max_fallR &
                                                      fallR_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_max_fallR = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_fallR &
                                                                             summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_max_fallR = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_fallR &
                                                                                         prevwinterT_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_max_fallR = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_fallR &
                                                                                   prevfallR_timeSeries$site == population)])

# Density
dens_value_max_fallR = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_fallR &
                                                                          yearly_density_per_square$site == population)]

# Size
size_value_max_fallR = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_fallR)])

# TSF
TSF_value_max_fallR = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_fallR)])


## 4.2.7. Min previous fall cumulative rainfall ----
# ---------------------------------------------

year_min_prevfallR = unique(prevfallR_timeSeries$time[which(prevfallR_timeSeries$prevfallR_unscaled == min_prevfallR &
                                                              prevfallR_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_min_prevfallR = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_prevfallR &
                                                                                 summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_min_prevfallR = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_prevfallR &
                                                                                             prevwinterT_timeSeries$site == population)])

# Next fall cumulative rainfall
fallR_value_min_prevfallR = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_prevfallR &
                                                                           fallR_timeSeries$site == population)])

# Density
dens_value_min_prevfallR = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_prevfallR &
                                                                              yearly_density_per_square$site == population)]

# Size
size_value_min_prevfallR = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_prevfallR)])

# TSF
TSF_value_min_prevfallR = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_prevfallR)])


## 4.2.8. Max previous fall cumulative rainfall ----
# ---------------------------------------------

year_max_prevfallR = unique(prevfallR_timeSeries$time[which(prevfallR_timeSeries$prevfallR_unscaled == max_prevfallR &
                                                              prevfallR_timeSeries$site == population)])

# Next summer mean max. daily temperature
summerT_value_max_prevfallR = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_prevfallR &
                                                                                 summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_max_prevfallR = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_prevfallR &
                                                                                             prevwinterT_timeSeries$site == population)])

# Next fall cumulative rainfall
fallR_value_max_prevfallR = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_prevfallR &
                                                                           fallR_timeSeries$site == population)])

# Density
dens_value_max_prevfallR = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_prevfallR &
                                                                              yearly_density_per_square$site == population)]

# Size
size_value_max_prevfallR = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_prevfallR)])

# TSF
TSF_value_max_prevfallR = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_prevfallR)])


## 4.2.9. Min density ----
# -------------------

year_min_dens = unique(yearly_density_per_square$time[which(yearly_density_per_square$abLarge_unscaled == min_dens &
                                                              yearly_density_per_square$site == population)])

# Next summer mean max. daily temperature
summerT_value_min_dens = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_dens &
                                                                            summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_min_dens = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_dens &
                                                                                        prevwinterT_timeSeries$site == population)])

# Next fall cumulative rainfall
fallR_value_min_dens = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_dens &
                                                                      fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_min_dens = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_dens &
                                                                                  prevfallR_timeSeries$site == population)])

# Size
size_value_min_dens = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_dens)])

# TSF
TSF_value_min_dens = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_dens)])


## 4.2.10. Max density ----
# --------------------

year_max_dens = unique(yearly_density_per_square$time[which(yearly_density_per_square$abLarge_unscaled == max_dens &
                                                              yearly_density_per_square$site == population)])

# Next summer mean max. daily temperature
summerT_value_max_dens = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_dens &
                                                                            summerT_timeSeries$site == population)])  

# Previous winter mean max. daily temperature
prevwinterT_value_max_dens = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_dens &
                                                                                        prevwinterT_timeSeries$site == population)])

# Next fall cumulative rainfall
fallR_value_max_dens = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_dens &
                                                                      fallR_timeSeries$site == population)])

# Previous fall cumulative rainfall
prevfallR_value_max_dens = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_dens &
                                                                                  prevfallR_timeSeries$site == population)])

# Size
size_value_max_dens = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_dens)])

# TSF
TSF_value_max_dens = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_dens)])


# ## 4.2.11. Min size ----
# # -----------------
# 
# year_min_size = unique(yearly_mean_size$time[which(yearly_mean_size$size_unscaled == min_size)])
# 
# # Next summer mean max. daily temperature
# summerT_value_min_size = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_size &
#                                                                             summerT_timeSeries$site == population)])  
# 
# # Previous winter mean max. daily temperature
# prevwinterT_value_min_size = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_size &
#                                                                                         prevwinterT_timeSeries$site == population)])
# 
# # Next fall cumulative rainfall
# fallR_value_min_size = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_size &
#                                                                       fallR_timeSeries$site == population)])
# 
# # Previous fall cumulative rainfall
# prevfallR_value_min_size = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_size &
#                                                                                   prevfallR_timeSeries$site == population)])
# 
# # Density
# dens_value_min_size = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_size &
#                                                                          yearly_density_per_square$site == population)]
# 
# # TSF
# TSF_value_min_size = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_min_size)])
# 
# 
# ## 4.2.12. Max size ----
# # -----------------
# 
# year_max_size = unique(yearly_mean_size$time[which(yearly_mean_size$size_unscaled == max_size)])
# 
# # Next summer mean max. daily temperature
# summerT_value_max_size = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_size &
#                                                                             summerT_timeSeries$site == population)])  
# 
# # Previous winter mean max. daily temperature
# prevwinterT_value_max_size = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_size &
#                                                                                         prevwinterT_timeSeries$site == population)])
# 
# # Next fall cumulative rainfall
# fallR_value_max_size = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_size &
#                                                                       fallR_timeSeries$site == population)])
# 
# # Previous fall cumulative rainfall
# prevfallR_value_max_size = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_size &
#                                                                                   prevfallR_timeSeries$site == population)])
# 
# # Density
# dens_value_max_size = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_size &
#                                                                          yearly_density_per_square$site == population)]
# 
# # TSF
# TSF_value_max_size = mean(droso_pop$TSFcont_unscaled[which(droso_pop$time == year_max_size)])
# 
# 
# ## 4.2.13. Min TSF ----
# # ----------------
# 
# year_min_TSF = unique(droso_pop$time[which(droso_pop$TSFcont_unscaled == min_TSF &
#                                              droso_pop$site == population)])
# 
# # Next summer mean max. daily temperature
# summerT_value_min_TSF = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_min_TSF &
#                                                                            summerT_timeSeries$site == population)])  
# 
# # Previous winter mean max. daily temperature
# prevwinterT_value_min_TSF = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_min_TSF &
#                                                                                        prevwinterT_timeSeries$site == population)])
# 
# # Next fall cumulative rainfall
# fallR_value_min_TSF = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_min_TSF &
#                                                                      fallR_timeSeries$site == population)])
# 
# # Previous fall cumulative rainfall
# prevfallR_value_min_TSF = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_min_TSF &
#                                                                                  prevfallR_timeSeries$site == population)])
# 
# # Density
# dens_value_min_TSF = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_min_TSF &
#                                                                         yearly_density_per_square$site == population)]
# 
# # Size
# size_value_min_TSF = mean(droso_pop$size_unscaled[which(droso_pop$time == year_min_TSF)])
# 
# 
# ## 4.2.14. Max TSF ----
# # ----------------
# 
# year_max_TSF = unique(droso_pop$time[which(droso_pop$TSFcont_unscaled == max_TSF &
#                                              droso_pop$site == population)])
# 
# # Next summer mean max. daily temperature
# summerT_value_max_TSF = unique(summerT_timeSeries$summerT_unscaled[which(summerT_timeSeries$time == year_max_TSF &
#                                                                            summerT_timeSeries$site == population)])  
# 
# # Previous winter mean max. daily temperature
# prevwinterT_value_max_TSF = unique(prevwinterT_timeSeries$prevwinterT_unscaled[which(prevwinterT_timeSeries$time == year_max_TSF &
#                                                                                        prevwinterT_timeSeries$site == population)])
# 
# # Next fall cumulative rainfall
# fallR_value_max_TSF = unique(fallR_timeSeries$fallR_unscaled[which(fallR_timeSeries$time == year_max_TSF &
#                                                                      fallR_timeSeries$site == population)])
# 
# # Previous fall cumulative rainfall
# prevfallR_value_max_TSF = unique(prevfallR_timeSeries$prevfallR_unscaled[which(prevfallR_timeSeries$time == year_max_TSF &
#                                                                                  prevfallR_timeSeries$site == population)])
# 
# # Density
# dens_value_max_TSF = yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$time == year_max_TSF &
#                                                                         yearly_density_per_square$site == population)]
# 
# # Size
# size_value_max_TSF = mean(droso_pop$size_unscaled[which(droso_pop$time == year_max_TSF)])


# Create data frame to store values
covariates_values_1 = expand.grid(focal_cov = "summerT",
                                  focal_min_max = c("min", "max"),
                                  other_cov = c("prevwinterT", "fallR", "prevfallR", 
                                                "dens", "size", "TSF"),
                                  other_mean_obs = c("mean", "obs"))

covariates_values_1$year = NA
covariates_values_1$year[which(covariates_values_1$focal_min_max == "min")] = year_min_summerT
covariates_values_1$year[which(covariates_values_1$focal_min_max == "max")] = year_max_summerT

covariates_values_1$focal_min_max
covariates_values_1$focal_value = rep(c(min_summerT, max_summerT), nrow(covariates_values_1)/2)

covariates_values_1$other_value = NA
covariates_values_1$other_cov[which(covariates_values_1$other_mean_obs == "mean")]
covariates_values_1$other_value[which(covariates_values_1$other_mean_obs == "mean")] = rep(c(mean_prevwinterT,
                                                                                             mean_fallR, 
                                                                                             mean_prevfallR,
                                                                                             mean_dens,
                                                                                             mean_size,
                                                                                             mean_TSF), each = 2)

paste(covariates_values_1$other_cov[which(covariates_values_1$other_mean_obs == "obs")],
      covariates_values_1$focal_min_max[which(covariates_values_1$other_mean_obs == "obs")])
covariates_values_1$other_value[which(covariates_values_1$other_mean_obs == "obs")] = c(prevwinterT_value_min_summerT,
                                                                                        prevwinterT_value_max_summerT,
                                                                                        fallR_value_min_summerT,
                                                                                        fallR_value_max_summerT,
                                                                                        prevfallR_value_min_summerT,
                                                                                        prevfallR_value_max_summerT,
                                                                                        dens_value_min_summerT,
                                                                                        dens_value_max_summerT,
                                                                                        size_value_min_summerT,
                                                                                        size_value_max_summerT,
                                                                                        TSF_value_min_summerT,
                                                                                        TSF_value_max_summerT)


covariates_values_2 = expand.grid(focal_cov = "prevwinterT",
                                  focal_min_max = c("min", "max"),
                                  other_cov = c("summerT", "fallR", "prevfallR",
                                                "dens", "size", "TSF"),
                                  other_mean_obs = c("mean", "obs"))

covariates_values_2$year = NA
covariates_values_2$year[which(covariates_values_2$focal_min_max == "min")] = year_min_prevwinterT
covariates_values_2$year[which(covariates_values_2$focal_min_max == "max")] = year_max_prevwinterT

covariates_values_2$focal_min_max
covariates_values_2$focal_value = rep(c(min_prevwinterT, max_prevwinterT), nrow(covariates_values_2)/2)

covariates_values_2$other_value = NA
covariates_values_2$other_cov[which(covariates_values_2$other_mean_obs == "mean")]
covariates_values_2$other_value[which(covariates_values_2$other_mean_obs == "mean")] = rep(c(mean_summerT,
                                                                                             mean_fallR, 
                                                                                             mean_prevfallR,
                                                                                             mean_dens,
                                                                                             mean_size,
                                                                                             mean_TSF), each = 2)

paste(covariates_values_2$other_cov[which(covariates_values_2$other_mean_obs == "obs")],
      covariates_values_2$focal_min_max[which(covariates_values_2$other_mean_obs == "obs")])
covariates_values_2$other_value[which(covariates_values_2$other_mean_obs == "obs")] = c(summerT_value_min_prevwinterT,
                                                                                        summerT_value_max_prevwinterT,
                                                                                        fallR_value_min_prevwinterT,
                                                                                        fallR_value_max_prevwinterT,
                                                                                        prevfallR_value_min_prevwinterT,
                                                                                        prevfallR_value_max_prevwinterT,
                                                                                        dens_value_min_prevwinterT,
                                                                                        dens_value_max_prevwinterT,
                                                                                        size_value_min_prevwinterT,
                                                                                        size_value_max_prevwinterT,
                                                                                        TSF_value_min_prevwinterT,
                                                                                        TSF_value_max_prevwinterT)


covariates_values_3 = expand.grid(focal_cov = "fallR",
                                  focal_min_max = c("min", "max"),
                                  other_cov = c("summerT", "prevwinterT", "prevfallR",
                                                "dens", "size", "TSF"),
                                  other_mean_obs = c("mean", "obs"))

covariates_values_3$year = NA
covariates_values_3$year[which(covariates_values_3$focal_min_max == "min")] = year_min_fallR
covariates_values_3$year[which(covariates_values_3$focal_min_max == "max")] = year_max_fallR

covariates_values_3$focal_min_max
covariates_values_3$focal_value = rep(c(min_fallR, max_fallR), nrow(covariates_values_3)/2)

covariates_values_3$other_value = NA
covariates_values_3$other_cov[which(covariates_values_3$other_mean_obs == "mean")]
covariates_values_3$other_value[which(covariates_values_3$other_mean_obs == "mean")] = rep(c(mean_summerT,
                                                                                             mean_prevwinterT, 
                                                                                             mean_prevfallR,
                                                                                             mean_dens,
                                                                                             mean_size,
                                                                                             mean_TSF), each = 2)

paste(covariates_values_3$other_cov[which(covariates_values_3$other_mean_obs == "obs")],
      covariates_values_3$focal_min_max[which(covariates_values_3$other_mean_obs == "obs")])
covariates_values_3$other_value[which(covariates_values_3$other_mean_obs == "obs")] = c(summerT_value_min_fallR,
                                                                                        summerT_value_max_fallR,
                                                                                        prevwinterT_value_min_fallR,
                                                                                        prevwinterT_value_max_fallR,
                                                                                        prevfallR_value_min_fallR,
                                                                                        prevfallR_value_max_fallR,
                                                                                        dens_value_min_fallR,
                                                                                        dens_value_max_fallR,
                                                                                        size_value_min_fallR,
                                                                                        size_value_max_fallR,
                                                                                        TSF_value_min_fallR,
                                                                                        TSF_value_max_fallR)


covariates_values_4 = expand.grid(focal_cov = "prevfallR",
                                  focal_min_max = c("min", "max"),
                                  other_cov = c("summerT", "prevwinterT", "fallR",
                                                "dens", "size", "TSF"),
                                  other_mean_obs = c("mean", "obs"))

covariates_values_4$year = NA
covariates_values_4$year[which(covariates_values_4$focal_min_max == "min")] = year_min_prevfallR
covariates_values_4$year[which(covariates_values_4$focal_min_max == "max")] = year_max_prevfallR

covariates_values_4$focal_min_max
covariates_values_4$focal_value = rep(c(min_prevfallR, max_prevfallR), nrow(covariates_values_4)/2)

covariates_values_4$other_value = NA
covariates_values_4$other_cov[which(covariates_values_4$other_mean_obs == "mean")]
covariates_values_4$other_value[which(covariates_values_4$other_mean_obs == "mean")] = rep(c(mean_summerT,
                                                                                             mean_prevwinterT, 
                                                                                             mean_fallR,
                                                                                             mean_dens,
                                                                                             mean_size,
                                                                                             mean_TSF), each = 2)

paste(covariates_values_4$other_cov[which(covariates_values_4$other_mean_obs == "obs")],
      covariates_values_4$focal_min_max[which(covariates_values_4$other_mean_obs == "obs")])
covariates_values_4$other_value[which(covariates_values_4$other_mean_obs == "obs")] = c(summerT_value_min_prevfallR,
                                                                                        summerT_value_max_prevfallR,
                                                                                        prevwinterT_value_min_prevfallR,
                                                                                        prevwinterT_value_max_prevfallR,
                                                                                        fallR_value_min_prevfallR,
                                                                                        fallR_value_max_prevfallR,
                                                                                        dens_value_min_prevfallR,
                                                                                        dens_value_max_prevfallR,
                                                                                        size_value_min_prevfallR,
                                                                                        size_value_max_prevfallR,
                                                                                        TSF_value_min_prevfallR,
                                                                                        TSF_value_max_prevfallR)

covariates_values_5 = expand.grid(focal_cov = "dens",
                                  focal_min_max = c("min", "max"),
                                  other_cov = c("summerT", "prevwinterT", "fallR",
                                                "prevfallR", "size", "TSF"),
                                  other_mean_obs = c("mean", "obs"))

covariates_values_5$year = NA
covariates_values_5$year[which(covariates_values_5$focal_min_max == "min")] = year_min_dens
covariates_values_5$year[which(covariates_values_5$focal_min_max == "max")] = year_max_dens

covariates_values_5$focal_min_max
covariates_values_5$focal_value = rep(c(min_dens, max_dens), nrow(covariates_values_5)/2)

covariates_values_5$other_value = NA
covariates_values_5$other_cov[which(covariates_values_5$other_mean_obs == "mean")]
covariates_values_5$other_value[which(covariates_values_5$other_mean_obs == "mean")] = rep(c(mean_summerT,
                                                                                             mean_prevwinterT, 
                                                                                             mean_fallR,
                                                                                             mean_prevfallR,
                                                                                             mean_size,
                                                                                             mean_TSF), each = 2)

paste(covariates_values_5$other_cov[which(covariates_values_5$other_mean_obs == "obs")],
      covariates_values_5$focal_min_max[which(covariates_values_5$other_mean_obs == "obs")])
covariates_values_5$other_value[which(covariates_values_5$other_mean_obs == "obs")] = c(summerT_value_min_dens,
                                                                                        summerT_value_max_dens,
                                                                                        prevwinterT_value_min_dens,
                                                                                        prevwinterT_value_max_dens,
                                                                                        fallR_value_min_dens,
                                                                                        fallR_value_max_dens,
                                                                                        prevfallR_value_min_dens,
                                                                                        prevfallR_value_max_dens,
                                                                                        size_value_min_dens,
                                                                                        size_value_max_dens,
                                                                                        TSF_value_min_dens,
                                                                                        TSF_value_max_dens)


# covariates_values_6 = expand.grid(focal_cov = "size",
#                                   focal_min_max = c("min", "max"),
#                                   other_cov = c("summerT", "prevwinterT", "fallR",
#                                                 "prevfallR", "dens", "TSF"),
#                                   other_mean_obs = c("mean", "obs"))
# 
# covariates_values_6$year = NA
# covariates_values_6$year[which(covariates_values_6$focal_min_max == "min")] = year_min_size
# covariates_values_6$year[which(covariates_values_6$focal_min_max == "max")] = year_max_size
# 
# covariates_values_6$focal_min_max
# covariates_values_6$focal_value = rep(c(min_size, max_size), nrow(covariates_values_6)/2)
# 
# covariates_values_6$other_value = NA
# covariates_values_6$other_cov[which(covariates_values_6$other_mean_obs == "mean")]
# covariates_values_6$other_value[which(covariates_values_6$other_mean_obs == "mean")] = rep(c(mean_summerT,
#                                                                                              mean_prevwinterT, 
#                                                                                              mean_fallR,
#                                                                                              mean_prevfallR,
#                                                                                              mean_dens,
#                                                                                              mean_TSF), each = 2)
# 
# paste(covariates_values_6$other_cov[which(covariates_values_6$other_mean_obs == "obs")],
#       covariates_values_6$focal_min_max[which(covariates_values_6$other_mean_obs == "obs")])
# covariates_values_6$other_value[which(covariates_values_6$other_mean_obs == "obs")] = c(summerT_value_min_size,
#                                                                                         summerT_value_max_size,
#                                                                                         prevwinterT_value_min_size,
#                                                                                         prevwinterT_value_max_size,
#                                                                                         fallR_value_min_size,
#                                                                                         fallR_value_max_size,
#                                                                                         prevfallR_value_min_size,
#                                                                                         prevfallR_value_max_size,
#                                                                                         dens_value_min_size,
#                                                                                         dens_value_max_size,
#                                                                                         TSF_value_min_size,
#                                                                                         TSF_value_max_size)
# 
# 
# covariates_values_7 = expand.grid(focal_cov = "TSF",
#                                   focal_min_max = c("min", "max"),
#                                   other_cov = c("summerT", "prevwinterT", "fallR",
#                                                 "prevfallR", "dens", "size"),
#                                   other_mean_obs = c("mean", "obs"))
# 
# covariates_values_7$year = NA
# covariates_values_7$year[which(covariates_values_7$focal_min_max == "min")] = year_min_TSF
# covariates_values_7$year[which(covariates_values_7$focal_min_max == "max")] = year_max_TSF
# 
# covariates_values_7$focal_min_max
# covariates_values_7$focal_value = rep(c(min_TSF, max_TSF), nrow(covariates_values_7)/2)
# 
# covariates_values_7$other_value = NA
# covariates_values_7$other_cov[which(covariates_values_7$other_mean_obs == "mean")]
# covariates_values_7$other_value[which(covariates_values_7$other_mean_obs == "mean")] = rep(c(mean_summerT,
#                                                                                              mean_prevwinterT, 
#                                                                                              mean_fallR,
#                                                                                              mean_prevfallR,
#                                                                                              mean_dens,
#                                                                                              mean_size), each = 2)
# 
# paste(covariates_values_7$other_cov[which(covariates_values_7$other_mean_obs == "obs")],
#       covariates_values_7$focal_min_max[which(covariates_values_7$other_mean_obs == "obs")])
# covariates_values_7$other_value[which(covariates_values_7$other_mean_obs == "obs")] = c(summerT_value_min_TSF,
#                                                                                         summerT_value_max_TSF,
#                                                                                         prevwinterT_value_min_TSF,
#                                                                                         prevwinterT_value_max_TSF,
#                                                                                         fallR_value_min_TSF,
#                                                                                         fallR_value_max_TSF,
#                                                                                         prevfallR_value_min_TSF,
#                                                                                         prevfallR_value_max_TSF,
#                                                                                         dens_value_min_TSF,
#                                                                                         dens_value_max_TSF,
#                                                                                         size_value_min_TSF,
#                                                                                         size_value_max_TSF)



covariates_values = rbind(covariates_values_1,
                          covariates_values_2,
                          covariates_values_3,
                          covariates_values_4,
                          covariates_values_5
                          # ,
                          # covariates_values_6,
                          # covariates_values_7
)
rm(list = c("covariates_values_1",
            "covariates_values_2",
            "covariates_values_3",
            "covariates_values_4",
            "covariates_values_5"
            # ,
            # "covariates_values_6",
            # "covariates_values_7"
))




###########################################################################
#
# 5. Individual-based model function ----
#
###########################################################################

## 5.1. Timestep projection function ----
# ----------------------------------

ibm_sim = function(n_sim,              # Number of simulations
                   n_years,            # Number of years per simulation
                   covariates_values,
                   focal_vr, 
                   focal_cov,          # Focal covariate
                   focal_min_max,      # Focal covariate fixed value (min or max)
                   sim = sim,          # Current simulation number
                   first_year,         # Year of start of simulation
                   years_RE,           # Observed years available for random year effect
                   seedbank_size,      # Initial seedbank size
                   recruitCap,         # Maximum number of recruits per quadrat
                   max_nbFlowers,      # Maximum number of flowers per individual
                   population,         # Population ID
                   # ibm_data,           # Previously stored projection data
                   # ibm_data_repro,     # Previously stored projection data, large individuals
                   data_initial,       # Initial dataset
                   seedbank_initial){  # Initial seedbank data
  
  # print("COVARIATES")
  # Focal and observed covariate values for each vital rate
  
  # Focal vital rate = non-reproductive survival
  if(focal_vr == "non_repro_survival"){
    
    # Focal covariate = Next summer temperature
    if(focal_cov == "summerT"){
      
      # Temperature value
      nrSurv_summerT_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for non-reproductive survival = observed 
      # value when summer temperature is the value above
      
      nrSurv_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                        covariates_values$other_mean_obs == "obs" &
                                                                        covariates_values$focal_cov == focal_cov &
                                                                        covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Next fall rainfall
    else if(focal_cov == "fallR"){
      
      # Rainfall value
      nrSurv_fallR_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                        covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for non-reproductive survival = observed 
      # value when fall rainfall is the value above
      
      nrSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Aboveground density of large individuals
    else if(focal_cov == "dens"){
      
      # Density value
      nrSurv_dens_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for non-reproductive survival = observed 
      # value when aboveground density is the value above
      
      nrSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                        covariates_values$other_mean_obs == "obs" &
                                                                        covariates_values$focal_cov == focal_cov &
                                                                        covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      nrSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    
    # Other covariate values for other vital rates = mean value
    
    rSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                       covariates_values$other_mean_obs == "mean")])
    
    flowering_prevwinterT_value = nbFlowers_prevwinterT_value = seedlingSize_prevwinterT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevwinterT" &
                                                                                                                                              covariates_values$other_mean_obs == "mean")])
    
    rSurv_fallR_value = growth_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                                          covariates_values$other_mean_obs == "mean")])
    
    flowering_prevfallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevfallR" &
                                                                             covariates_values$other_mean_obs == "mean")])
    
    rSurv_dens_value = growth_dens_value = flowering_dens_value = seedlingSize_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                                                                                         covariates_values$other_mean_obs == "mean")])
    
    rSurv_size_value = growth_size_value = flowering_size_value = nbFlowers_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                                                                                      covariates_values$other_mean_obs == "mean")])
    
    rSurv_TSF_value = growth_TSF_value = flowering_TSF_value = nbFlowers_TSF_value = seedlingSize_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                                                                                                           covariates_values$other_mean_obs == "mean")])
  }  
  
  # Focal vital rate = reproductive survival
  if(focal_vr == "repro_survival"){
    
    # Focal covariate = Next summer temperature
    if(focal_cov == "summerT"){
      
      # Temperature value
      rSurv_summerT_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for reproductive survival = observed 
      # value when summer temperature is the value above
      
      rSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                          covariates_values$other_mean_obs == "mean")])
      
      rSurv_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                     covariates_values$other_mean_obs == "obs" &
                                                                     covariates_values$focal_cov == focal_cov &
                                                                     covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Next fall rainfall
    else if(focal_cov == "fallR"){
      
      # Rainfall value
      rSurv_fallR_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for reproductive survival = observed 
      # value when fall rainfall is the value above
      
      rSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                         covariates_values$other_mean_obs == "obs" &
                                                                         covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                     covariates_values$other_mean_obs == "obs" &
                                                                     covariates_values$focal_cov == focal_cov &
                                                                     covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Aboveground density of large individuals
    else if(focal_cov == "dens"){
      
      # Density value
      rSurv_dens_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for reproductive survival = observed 
      # value when aboveground density is the value above
      
      rSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                         covariates_values$other_mean_obs == "obs" &
                                                                         covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                       covariates_values$other_mean_obs == "obs" &
                                                                       covariates_values$focal_cov == focal_cov &
                                                                       covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                      covariates_values$other_mean_obs == "obs" &
                                                                      covariates_values$focal_cov == focal_cov &
                                                                      covariates_values$focal_min_max == focal_min_max)])
      
      rSurv_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                     covariates_values$other_mean_obs == "obs" &
                                                                     covariates_values$focal_cov == focal_cov &
                                                                     covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Other covariate values for other vital rates = mean value
    
    nrSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                        covariates_values$other_mean_obs == "mean")])
    
    flowering_prevwinterT_value = nbFlowers_prevwinterT_value = seedlingSize_prevwinterT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevwinterT" &
                                                                                                                                              covariates_values$other_mean_obs == "mean")])
    
    nrSurv_fallR_value = growth_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                                           covariates_values$other_mean_obs == "mean")])
    
    flowering_prevfallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevfallR" &
                                                                             covariates_values$other_mean_obs == "mean")])
    
    nrSurv_dens_value = growth_dens_value = flowering_dens_value = seedlingSize_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                                                                                          covariates_values$other_mean_obs == "mean")])
    
    nrSurv_size_value = growth_size_value = flowering_size_value = nbFlowers_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                                                                                       covariates_values$other_mean_obs == "mean")])
    
    nrSurv_TSF_value = growth_TSF_value = flowering_TSF_value = nbFlowers_TSF_value = seedlingSize_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                                                                                                            covariates_values$other_mean_obs == "mean")])
    
  }  
  
  
  # Focal vital rate = flowering probability
  if(focal_vr == "flowering"){
    
    # Focal covariate = Previous winter temperature
    if(focal_cov == "prevwinterT"){
      
      # Temperature value
      flowering_prevwinterT_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                                 covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for flowering = observed 
      # value when previous winter temperature is the value above
      
      flowering_prevfallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevfallR" &
                                                                               covariates_values$other_mean_obs == "obs" &
                                                                               covariates_values$focal_cov == focal_cov &
                                                                               covariates_values$focal_min_max == focal_min_max)])
      
      flowering_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      flowering_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      flowering_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                         covariates_values$other_mean_obs == "obs" &
                                                                         covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Previous fall rainfall
    else if(focal_cov == "prevfallR"){
      
      # Rainfall value
      flowering_prevfallR_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                               covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for flowering probability = observed 
      # value when previous fall rainfall is the value above
      
      flowering_prevwinterT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevwinterT" &
                                                                                 covariates_values$other_mean_obs == "obs" &
                                                                                 covariates_values$focal_cov == focal_cov &
                                                                                 covariates_values$focal_min_max == focal_min_max)])
      
      flowering_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      flowering_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      flowering_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                         covariates_values$other_mean_obs == "obs" &
                                                                         covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Focal covariate = Aboveground density of large individuals
    else if(focal_cov == "dens"){
      
      # Density value
      flowering_dens_value = unique(covariates_values$focal_value[which(covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      # Other covariate values for flowering probability = observed 
      # value when aboveground density is the value above
      
      flowering_prevwinterT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevwinterT" &
                                                                                 covariates_values$other_mean_obs == "obs" &
                                                                                 covariates_values$focal_cov == focal_cov &
                                                                                 covariates_values$focal_min_max == focal_min_max)])
      
      flowering_prevfallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevfallR" &
                                                                               covariates_values$other_mean_obs == "obs" &
                                                                               covariates_values$focal_cov == focal_cov &
                                                                               covariates_values$focal_min_max == focal_min_max)])
      
      flowering_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                          covariates_values$other_mean_obs == "obs" &
                                                                          covariates_values$focal_cov == focal_cov &
                                                                          covariates_values$focal_min_max == focal_min_max)])
      
      flowering_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                         covariates_values$other_mean_obs == "obs" &
                                                                         covariates_values$focal_cov == focal_cov &
                                                                         covariates_values$focal_min_max == focal_min_max)])
      
    }
    
    # Other covariate values for other vital rates = mean value
    
    nrSurv_summerT_value = rSurv_summerT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "summerT" &
                                                                                              covariates_values$other_mean_obs == "mean")])
    
    nbFlowers_prevwinterT_value = seedlingSize_prevwinterT_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "prevwinterT" &
                                                                                                                covariates_values$other_mean_obs == "mean")])
    
    nrSurv_fallR_value = rSurv_fallR_value = growth_fallR_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "fallR" &
                                                                                                               covariates_values$other_mean_obs == "mean")])
    
    nrSurv_dens_value = rSurv_dens_value = growth_dens_value = flowering_dens_value = seedlingSize_dens_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "dens" &
                                                                                                                                                             covariates_values$other_mean_obs == "mean")])
    
    nrSurv_size_value = rSurv_size_value = growth_size_value = flowering_size_value = nbFlowers_size_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "size" &
                                                                                                                                                          covariates_values$other_mean_obs == "mean")])
    
    nrSurv_TSF_value = rSurv_TSF_value = growth_TSF_value = flowering_TSF_value = nbFlowers_TSF_value = seedlingSize_TSF_value = unique(covariates_values$other_value[which(covariates_values$other_cov == "TSF" &
                                                                                                                                                                              covariates_values$other_mean_obs == "mean")])
  }  
  
  
  # print("SETUP")
  # Year sequence for climatic variable predictions and random effects
  
  # years_obs = c(first_year, sample(years_RE, n_years-1, replace = T)) # Random effects
  years_obs = rep(covariates_values$year[which(covariates_values$focal_cov == focal_cov &
                                                 covariates_values$focal_min_max == focal_min_max)],
                  n_years) # Random effects
  
  # Empty files to hold results 
  
  log_lambda = rep(NA, n_years)
  seedbank_size_vec = rep(NA, n_years)
  seedbank_size_vec[1] = seedbank_size
  extinction = 0
  
  sim_data = NULL # Full individual data across whole simulation
  
  data_droso = data_initial[which(data_initial$time == years_obs[1]), ] # Yearly individual data
  
  # Seedbank seeds data
  data_SB = seedbank_initial
  
  # Assign seeds in initial seedbank to quadrats
  data_SB$quadratID = sample(unique(data_droso$quadratID), size = seedbank_size, 
                             replace = T, 
                             prob = as.numeric(table(data_droso$quadratID)/sum(table(data_droso$quadratID))))
  
  # Highest seed ID (format Seed_XXX) to give names to new seeds 
  max_seed_ID = max(as.numeric(unlist(lapply(strsplit(data_SB$ID, split = "_"), function(x) x[2]))))
  
  time_sim = 1 # Timestep
  sim_data = rbind(sim_data, cbind(data_droso, time_sim)) # Merge full individual 
  # data and yearly individual 
  # data with timestep info
  
  
  # Project the population
  for(i in 2:n_years){
    
    # print(paste("Year", i))
    
    seed_produced = data_droso[0, ] # Data on seeds produced by reproducing plants
    seed_produced[1, ] = NA
    seed_produced$goCont = seed_produced$goSB = seed_produced$staySB = seed_produced$outSB = seed_produced$rownb = NA
    seed_produced = seed_produced[0, ]
    
    
    # Get new years for random effects
    year_obs = years_obs[i-1]
    
    
    ### TSF, CORRECTION FACTORS AND NUMBER OF SQUARES ###
    
    # print("CORRECTION FACTORS AND NUMBER OF SQUARES")
    if(focal_vr == "non_repro_survival") TSF_value = nrSurv_TSF_value
    else if(focal_vr == "repro_survival") TSF_value = rSurv_TSF_value
    else if(focal_vr == "flowering") TSF_value = flowering_TSF_value
    
    TSFcat = ifelse(TSF_value <= 4, TSF_value, 4)
    
    corr_seed_surv = seedSurv_df$seedSurv[which(seedSurv_df$TSF == TSFcat)] # From Paniw et al. 2017 (J. Appl. Ecol.)
    
    
    if(nrow(data_droso) > 0){
      
      ### FLOWERING ###
      
      # print("FLOWERING")
      
      data_droso$flowering = NA
      
      data_droso$flowering = rbinom(n = nrow(data_droso), 
                                    flowering_function(size = flowering_size_value,
                                                       abLarge = flowering_dens_value,
                                                       prevwinterT = flowering_prevwinterT_value,
                                                       prevfallR = flowering_prevfallR_value,
                                                       TSFcont = flowering_TSF_value,
                                                       year = year_obs,
                                                       population = population), 
                                    size = 1)
      
      
      ### RECRUITMENT (number of flowers) ###
      # For the ones that reproduced 
      
      # print("RECRUITMENT")
      
      # Number of flowers per individual reproducing
      
      # print("NUMBER OF FLOWERS")
      
      data_droso$nbFlowers = NA
      
      if(length(data_droso$nbFlowers[which(data_droso$flowering == 1)]) > 0){
        
        data_droso$nbFlowers[which(data_droso$flowering == 1)] = rnbinom(n = nrow(data_droso[which(data_droso$flowering == 1), ]),
                                                                         mu = nbFlowers_function(size = nbFlowers_size_value,
                                                                                                 prevwinterT = nbFlowers_prevwinterT_value,
                                                                                                 TSFcont = nbFlowers_TSF_value,
                                                                                                 year = year_obs,
                                                                                                 population = population),
                                                                         size = 1)
        
        # Cap the number of flowers
        if(any(data_droso$nbFlowers[which(data_droso$flowering == 1)] > max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)], na.rm = T)) data_droso$nbFlowers[which(data_droso$flowering == 1 & data_droso$nbFlowers > max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)])] = max_nbFlowers$nbFlowers[which(max_nbFlowers$TSFcont == TSFcat)]
        
        # Number of seeds per individual reproducing
        
        # print("NUMBER OF SEEDS")
        
        data_droso_sub = data_droso[which(data_droso$flowering == 1 & 
                                            data_droso$nbFlowers > 0), ]
        
        if(any(data_droso_sub$nbFlowers > 0)){
          
          nb_seeds = rpois(nrow(data_droso_sub), lambda = seeds_per_flower) # Number of seeds per flower per individual
          
          seed_produced = data_droso_sub[rep(row.names(data_droso_sub), 
                                             data_droso_sub$nbFlowers * nb_seeds), ] # Dataset of seeds produced
          
          # Format dataset to match other datasets
          seed_produced$size_unscaled = seed_produced$sizeNext = seed_produced$flowering = seed_produced$nbFlowers = seed_produced$survival = NA
          rownames(seed_produced) = seq(1, nrow(seed_produced))
          
          
          # Seeds germinating directly
          
          # print("SEEDS GERMINATING")
          
          seed_produced$goCont = rbinom(n = nrow(seed_produced),
                                        prob = goCont_function(TSF = TSFcat) * corr_seed_surv,
                                        size = 1)
          
          
          # Seeds going to the seedbank
          
          # print("SEEDS TO SEEDBANK")
          
          seed_produced$goSB = NA
          
          if(any(seed_produced$goCont == 0)) seed_produced$goSB[which(seed_produced$goCont == 0)] = rbinom(n = nrow(seed_produced[which(seed_produced$goCont == 0), ]),
                                                                                                           prob = (1 - goCont_function(TSF = TSFcat)) * corr_seed_surv,
                                                                                                           size = 1)
          
          
          # Size of the germinated seedlings
          
          # print("SEEDS SIZE")
          
          # Seedling size parameters (mean, sd, and degrees of freedom)
          seedling_size_parameters = seedling_size_function(prevwinterT = seedlingSize_prevwinterT_value,
                                                            abLarge = seedlingSize_dens_value,
                                                            TSFcont = seedlingSize_TSF_value,
                                                            year = year_obs, 
                                                            population = population)
          
          # Get seedling size by sampling a truncated Student-t distribution
          # to match the model family, using the corresponding mean, sd, and df
          seed_produced$rownb = seq(1, nrow(seed_produced))
          
          seed_produced$size_unscaled[which(seed_produced$goCont == 1)] = apply(seed_produced[which(seed_produced$goCont == 1), ], 
                                                                                1, 
                                                                                function(x) rtt(1, 
                                                                                                location = seedling_size_parameters$mean, 
                                                                                                scale = seedling_size_parameters$sd, 
                                                                                                df = seedling_size_parameters$df, 
                                                                                                left = 0))
          
          # Assign zero to individuals with negative size
          if(any(seed_produced$size_unscaled[which(seed_produced$goCont == 1)] == Inf, na.rm = T)) seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled == Inf)] = max(seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled != Inf)], na.rm = T)
          
          if(any(is.na(seed_produced$size_unscaled[which(seed_produced$goCont == 1)]), na.rm = T)) seed_produced$size_unscaled[which(seed_produced$goCont == 1 & seed_produced$size_unscaled < 0)] = 0
          
          
          # Seedling ID
          
          # print("SEEDS ID")
          
          seed_produced$ID = paste("Seed", 
                                   seq(max_seed_ID + 1, max_seed_ID + 1 + nrow(seed_produced) - 1), 
                                   sep = "_")
          
          seed_produced = seed_produced[, colnames(seed_produced)[-which(colnames(seed_produced) == "rownb")]]
          
          # Update max seed ID
          max_seed_ID = max(as.numeric(unlist(lapply(strsplit(seed_produced$ID, split = "_"), function(x) x[2]))))
        }
      }
      
      ### SURVIVAL ###
      
      # print("SURVIVAL")
      
      data_droso$survival[which(data_droso$size < 4.5)] = rbinom(n = nrow(data_droso[which(data_droso$size < 4.5)]),
                                                                 prob = non_repro_survival_function(size = nrSurv_size_value,
                                                                                                    abLarge = nrSurv_dens_value,
                                                                                                    fallR = nrSurv_fallR_value,
                                                                                                    summerT = nrSurv_summerT_value,
                                                                                                    TSFcont = nrSurv_TSF_value,
                                                                                                    year = year_obs,
                                                                                                    population = population),
                                                                 size = 1)
      
      data_droso$survival[which(data_droso$size >= 4.5)] = rbinom(n = nrow(data_droso[which(data_droso$size >= 4.5)]),
                                                                  prob = repro_survival_function(size = rSurv_size_value,
                                                                                                 abLarge = rSurv_dens_value,
                                                                                                 fallR = rSurv_fallR_value,
                                                                                                 summerT = rSurv_summerT_value,
                                                                                                 TSFcont = rSurv_TSF_value,
                                                                                                 year = year_obs,
                                                                                                 population = population),
                                                                  size = 1)
      
      ### GROWTH ###
      
      # print("GROWTH")
      
      data_droso$sizeNext = NA
      
      # Growth parameters (mean, sd, and degrees of freedom)
      growth_parameters = growth_function(size = growth_size_value, 
                                          fallR = growth_fallR_value, 
                                          abLarge = growth_dens_value,
                                          TSFcont = growth_TSF_value,
                                          year = year_obs, 
                                          population = population)
      
      # Get size at next timestep by sampling a truncated Student-t distribution
      # to match the model family, using the corresponding mean, sd, and df
      data_droso$rownb = seq(1, nrow(data_droso))
      
      data_droso$sizeNext[which(data_droso$survival == 1)] = apply(data_droso[which(data_droso$survival == 1), ],
                                                                   1,
                                                                   function(x) rtt(1,
                                                                                   location = growth_parameters$mean[as.numeric(x[grep("rownb", names(x))])],
                                                                                   scale = growth_parameters$sd,
                                                                                   df = growth_parameters$df,
                                                                                   left = 0, 
                                                                                   right = max(data_droso$size_unscaled)))
      
      
      # Assign max size in current dataset to individual with infinite size
      if(any(data_droso$sizeNext[which(data_droso$survival == 1)] == Inf, na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext == Inf)] = max(data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext != Inf)], na.rm = T)
      
      if(any(is.na(data_droso$sizeNext[which(data_droso$survival == 1)]), na.rm = T)) data_droso$sizeNext[which(data_droso$survival == 1 & data_droso$sizeNext < 0)] = 0
      
      # Format dataset
      data_droso = data_droso[which(data_droso$survival == 1), colnames(data_droso)[-which(colnames(data_droso) %in% c("rownb"))]]
    }
    
    
    if(nrow(data_SB) > 0){
      
      ### SEEDBANK ###
      
      # print("SEEDBANK")
      
      # Seeds germinating from the seedbank 
      
      # print("SEEDS GERMINATING")
      
      data_SB$outSB = rbinom(n = nrow(data_SB),
                             prob = outSB_function(TSF = TSFcat,
                                                   seedSurv = corr_seed_surv),
                             size = 1)
      
      # Assign a size to the germinated seeds, add them to the dataset 
      # of aboveground individuals and remove them from the seedbank data
      
      # print("SEEDS SIZE")
      
      # Seedling size parameters (mean, sd, and degrees of freedom)
      seedling_size_parameters = seedling_size_function(prevwinterT = seedlingSize_prevwinterT_value,
                                                        abLarge = seedlingSize_dens_value,
                                                        TSFcont = seedlingSize_TSF_value,
                                                        year = year_obs, 
                                                        population = population)
      
      # Get seedling size by sampling a truncated Student-t distribution
      # to match the model family, using the corresponding mean, sd, and df
      data_SB$rownb = seq(1, nrow(data_SB))
      
      if(length(which(data_SB$outSB == 1)) > 0){
        
        data_SB$size_unscaled[which(data_SB$outSB == 1)] = apply(data_SB[which(data_SB$outSB == 1), ], 
                                                                 1, 
                                                                 function(x) rtt(1, 
                                                                                 location = seedling_size_parameters$mean, 
                                                                                 scale = seedling_size_parameters$sd, 
                                                                                 df = seedling_size_parameters$df, 
                                                                                 left = 0))
        
        # Assign zero in current dataset to individual with negative size
        if(any(is.na(data_SB$size_unscaled[which(data_SB$outSB == 1)]), na.rm = T)) data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled < 0)] = 0
        
        if(any(data_SB$size_unscaled[which(data_SB$outSB == 1)] == Inf, na.rm = T)) data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled == Inf)] = max(data_SB$size_unscaled[which(data_SB$outSB == 1 & data_SB$size_unscaled != Inf)], na.rm = T)
        
      }
      
      
      # Preparing the seedbank data to merge with the continuous germination data
      
      data_SB$goCont = data_SB$outSB
      data_SB$goSB = NA
      
      # Merge datasets 
      seed_produced = rbind(seed_produced, data_SB[which(data_SB$outSB == 1), colnames(data_SB)[which(colnames(data_SB) %in% colnames(seed_produced))]])
      
      # Remove germinated seeds from seedbank data
      data_SB = data_SB[which(data_SB$outSB == 0), ]
      
      data_SB = data_SB[, colnames(data_SB)[-which(colnames(data_SB) %in% c("rownb", "outSB"))]]
      
      
      # Seeds staying in and going to the seedbank
      
      # print("SEEDS STAYING")
      
      data_SB$staySB = rbinom(n = nrow(data_SB),
                              prob = staySB_function(TSF = TSFcat),
                              size = 1)
      
      # Keep only seeds staying in the seedbank
      data_SB = data_SB[which(data_SB$staySB == 1), ]
    }
    
    # Merge seeds going to SB to seedbank data
    seed_produced$staySB = seed_produced$goSB
    
    if(nrow(seed_produced[which(seed_produced$staySB == 1), ]) > 0){
      
      data_SB = rbind(data_SB, seed_produced[which(seed_produced$staySB == 1), ])
    }
    
    if(nrow(data_SB) > 0){
      
      # Format seedbank data to match other datasets
      data_SB = data_SB[, colnames(data_SB)[-which(colnames(data_SB) %in% c("staySB", "outSB", "rownb", "goCont", "goSB"))]]
    }
    
    seedbank_size_vec[i] = nrow(data_SB)
    
    if(nrow(seed_produced) > 0){
      
      # Keep only seeds germinating 
      seed_produced = seed_produced[which(seed_produced$goCont == 1), colnames(seed_produced)[-which(colnames(seed_produced) %in% c("goCont", "goSB", "staySB", "outSB", "rownb"))]]
      
      # Cap the number of recruits if needed
      if(!is.null(recruitCap) & nrow(seed_produced) > 0){
        
        # Get number of seedlings per quadrat
        quadratsAboveMaxSeedlings = aggregate(ID ~ quadratID, 
                                              data = seed_produced,
                                              FUN = function(x) length(x))
        
        # Keep quadrats where the number of seedlings is above the threshold
        quadratsAboveMaxSeedlings = quadratsAboveMaxSeedlings[which(quadratsAboveMaxSeedlings$ID > recruitCap$ID[which(recruitCap$TSFcont == TSFcat)]), ]
        
        # If quadrats are above the threshold, sample the seedlings that
        # will be kept
        if(nrow(quadratsAboveMaxSeedlings) > 0){
          
          for(quadrat in quadratsAboveMaxSeedlings$quadratID){
            
            recruitsKept_ID = sample(seq(1, nrow(seed_produced[which(seed_produced$quadratID == quadrat), ])), 
                                     size = recruitCap$ID[which(recruitCap$TSFcont == TSFcat)], replace = F)
            
            recruitsKept = seed_produced[which(seed_produced$quadratID == quadrat)[recruitsKept_ID], ]
            
            seed_produced = seed_produced[-which(seed_produced$quadratID == quadrat), ]
            seed_produced = rbind(seed_produced, recruitsKept)
          }
        }
      }
      
      
      # Adding new seedlings to the population
      
      # print("ADDING NEW SEEDLINGS TO POPULATION")
      
      data_droso = rbind(data_droso[-1, ], seed_produced)
    }
    
    log_lambda[i] = log((nrow(data_droso) + nrow(data_SB))/(nrow(sim_data[which(sim_data$time_sim == i-1), ]) + seedbank_size_vec[i-1])) # Calculate log lambda (with seedbank)
    
    
    ## SAVE DATA ###
    
    # print("SAVE DATA")
    
    time_sim = i # Update timestep
    
    # Merge yearly individual data with full individual data
    if(nrow(data_droso) > 0) data_droso = cbind(data_droso, time_sim)
    
    sim_data = rbind(sim_data, data_droso[, colnames(data_droso)[which(colnames(data_droso) %in% colnames(sim_data))]])
    
    # Format yearly data
    if(nrow(data_droso) > 0) data_droso = data_droso[, -ncol(data_droso)]
    
    # Assess population quasi-extinction
    if(nrow(data_droso) < 5 & nrow(data_SB) < 50){
      
      extinction = 1
      # break
      
    }
    
    
    ##### NEW  DATA FOR T + 1 #####
    
    # print("NEW DATA FOR T+1")
    
    data_droso$size_unscaled[which(!is.na(data_droso$sizeNext))] = data_droso$sizeNext[which(!is.na(data_droso$sizeNext))] # Assign new size to individuals
    
  }
  
  # Create summary data
  # data_agg = aggregate(size ~ time_sim, data = sim_data, function(x) length(x))
  # data_agg_repro = aggregate(size ~ time_sim, data = sim_data[which(sim_data$size > 4.5), ], function(x) length(x))
  # 
  # data_agg$run = sim
  # data_agg_repro$run = sim
  
  # ibm_data = rbind(ibm_data, data_agg)
  # ibm_data_repro = rbind(ibm_data_repro, data_agg_repro)
  
  
  return(list(# pop_data = sim_data,
    # pop_size = ibm_data,
    # pop_size_repro = ibm_data_repro,
    log_lambda = log_lambda,
    extinction = extinction))
}


## 5.2. Simulation function ----
# -------------------------

ibm = function(n_sim = 1000,               # Number of simulations
                       n_years = 50,               # Number of years per simulation
                       first_year = sample(seq(2016, 2021)),          # Year of start of simulation
                       years_RE = seq(2016, 2021), # Observed years available for random year effect
                       data_initial,               # Initial dataset
                       covariates_values,          # Values of covariates for various scenarios
                       focal_vr, 
                       focal_cov,                  # Focal covariate
                       focal_min_max,              # Focal covariate fixed value (min or max)
                       seedbank_size = 10000,      # Initial seedbank size
                       recruitCap,                 # Maximum number of recruits per quadrat
                       max_nbFlowers,              # Maximum number of flowers per individual
                       population){          # Fire frequency
  
  # Initialization - Prepare storing objects
  
  extinction_vector = rep(0, n_sim)
  log_lambda_array = array(NA, dim = c(n_sim, n_years))
  # pop_size_list = vector(mode = "list", length = n_sim)
  # pop_size_repro_list = vector(mode = "list", length = n_sim)
  # pop_data_list = vector(mode = "list", length = n_sim)
  
  # ibm_data = NULL
  # ibm_data_repro = NULL
  
  # Calculate number of flowers
  data_initial$nbFlowers = data_initial$fs * data_initial$fps
  
  # Format initial data
  data_initial = data_initial[, c("site", "quadratID", "ID", "TSFcont", "size_unscaled", "sizeNext", 
                                  "fl", "nbFlowers", 
                                  "surv", "time", "abLarge")]
  
  colnames(data_initial) = c("site", "quadratID", "ID", "TSFcont", "size_unscaled", "sizeNext", 
                             "flowering", "nbFlowers", 
                             "survival", "time", "abLarge")
  
  # Prepare initial seedbank data
  seedbank_initial = data.frame(site = population, quadratID = NA, ID = paste("Seed", seq(1, seedbank_size), sep = "_"),
                                TSFcont = NA, size_unscaled = NA, sizeNext = NA, 
                                flowering = NA, nbFlowers = NA,
                                survival = NA, time = NA, abLarge = NA)
  
  
  # Projection
  
  # For each simulation, project the population
  # If there is an error (exploding population), restart the simulation
  # to have another sequence of years
  
  for(sim in 1:n_sim){
    
    print(paste("Iteration", sim))
    
    while(TRUE){
      
      ibm_sim_result =
        try(ibm_sim(n_sim = n_sim,
                    n_years = n_years,
                    sim = sim,
                    first_year = first_year,
                    years_RE = years_RE,
                    covariates_values = covariates_values,
                    focal_vr = focal_vr,
                    focal_cov = focal_cov,
                    focal_min_max = focal_min_max,
                    seedbank_size = seedbank_size,
                    recruitCap = recruitCap,
                    max_nbFlowers,
                    population = population,
                    # ibm_data = ibm_data,
                    # ibm_data_repro = ibm_data_repro,
                    data_initial = data_initial,
                    seedbank_initial = seedbank_initial),
            silent = TRUE)
      
      if(!is(ibm_sim_result, 'try-error')) break
    }
    
    # Fill in result objects
    log_lambda_array[sim, ] = ibm_sim_result$log_lambda
    extinction_vector[sim] = ibm_sim_result$extinction
    # pop_size_list[[sim]] = ibm_sim_result$pop_size
    # pop_size_repro_list[[sim]] = ibm_sim_result$pop_size_repro
    # pop_data_list[[sim]] = ibm_sim_result$pop_data
  }
  
  return(list(# pop_data = pop_data_list,
    # pop_size = pop_size_list,
    # pop_size_repro = pop_size_repro_list,
    log_lambda = log_lambda_array,
    extinction = extinction_vector))
}




###########################################################################
#
# 6. Projections ----
#
###########################################################################

n_sim = 100
n_years = 50

## 6.1. Vertedero ----
# ---------------

population = "Vertedero"

recruitCap = aggregate(ID ~ quadratID + time + TSFcont_unscaled,
                       data = droso[which(droso$site == population &
                                                    droso$stage == "SD"), ], FUN = function(x) length(x))
colnames(recruitCap)[which(colnames(recruitCap) == "TSFcont_unscaled")] = "TSFcont"

recruitCap$TSFcont[which(recruitCap$TSFcont > 4)] = 4

recruitCap = aggregate(ID ~ TSFcont,
                       data = recruitCap, FUN = max)

recruitCap$ID = round(recruitCap$ID)

recruitCap = rbind(recruitCap, 
                   data.frame(TSFcont = c(0, 1, 2),
                              ID = c(round(max(recruitCap$ID) * 1.5), max(recruitCap$ID), round(mean(recruitCap$ID[which(recruitCap$TSFcont != 1)])))))

recruitCap = recruitCap[order(recruitCap$TSFcont), ]

max_nbFlowers = aggregate(fs * fps ~ TSFcont_unscaled,
                          data = droso[which(droso$site == population), ], FUN = max)
colnames(max_nbFlowers) = c("TSFcont", "nbFlowers")

max_nbFlowers$TSFcont[which(max_nbFlowers$TSFcont > 4)] = 4

max_nbFlowers = aggregate(nbFlowers ~ TSFcont,
                          data = max_nbFlowers, FUN = max)


ncpus <- 4

startIBM <- function(sim){
  
  if (sim > 1) {      #slave, give it other random seed,
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
  }
  
  ibm(n_sim = n_sim/ncpus, 
              n_years = n_years, 
              data_initial = droso[which(droso$site == population), ],
              covariates_values = covariates_values,
              focal_vr = focal_vr,
              focal_cov = focal_cov,
              focal_min_max = focal_min_max,
              recruitCap = recruitCap,
              max_nbFlowers = max_nbFlowers,
              population = population)
  
}


## 6.1.1. Focal vital rate = non-reproductive survival focal covariate = summerT, focal value = min ----
# -------------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "summerT"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_summerT_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_summerT_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.2. Focal vital rate = non-reproductive survival, focal covariate = summerT, focal value = max ----
# --------------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "summerT"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_summerT_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_summerT_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.3. Focal vital rate = non-reproductive survival, focal covariate = fallR, focal value = min ----
# ------------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "fallR"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_fallR_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_fallR_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.4. Focal vital rate = non-reproductive survival, focal covariate = fallR, focal value = max ----
# ------------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "fallR"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_fallR_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_fallR_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.5. Focal vital rate = non-reproductive survival, focal covariate = dens, focal value = min ----
# -----------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "dens"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_dens_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_dens_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.6. Focal vital rate = non-reproductive survival ,focal covariate = dens, focal value = max ----
# -----------------------------------------------------------------------------------------------

focal_vr = "non_repro_survival"
focal_cov = "dens"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_nrSurv_dens_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_nrSurv_dens_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.7. Focal vital rate = reproductive survival focal covariate = summerT, focal value = min ----
# ---------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "summerT"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_summerT_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_summerT_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.8. Focal vital rate = reproductive survival, focal covariate = summerT, focal value = max ----
# ----------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "summerT"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_summerT_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_summerT_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.9. Focal vital rate = reproductive survival, focal covariate = fallR, focal value = min ----
# --------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "fallR"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_fallR_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_fallR_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.10. Focal vital rate = reproductive survival, focal covariate = fallR, focal value = max ----
# ---------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "fallR"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_fallR_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_fallR_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.11. Focal vital rate = reproductive survival, focal covariate = dens, focal value = min ----
# --------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "dens"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_dens_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_dens_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.12. Focal vital rate = reproductive survival ,focal covariate = dens, focal value = max ----
# --------------------------------------------------------------------------------------------

focal_vr = "repro_survival"
focal_cov = "dens"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_rSurv_dens_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_rSurv_dens_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.13. Focal vital rate = flowering probability, focal covariate = prevwinterT, focal value = min ----
# ---------------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "prevwinterT"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_prevwinterT_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_prevwinterT_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.14. Focal vital rate = flowering probability, focal covariate = prevwinterT, focal value = max ----
# ---------------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "prevwinterT"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_prevwinterT_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_prevwinterT_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.15. Focal vital rate = flowering probability, focal covariate = prevfallR, focal value = min ----
# -------------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "prevfallR"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_prevfallR_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_prevfallR_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.16. Focal vital rate = flowering probability, focal covariate = prevfallR, focal value = max ----
# -------------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "prevfallR"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_prevfallR_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_prevfallR_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.17. Focal vital rate = flowering probability, focal covariate = dens, focal value = min ----
# --------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "dens"
focal_min_max = "min"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_dens_min <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_dens_min, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()


## 6.1.18. Focal vital rate = flowering probability, focal covariate = dens, focal value = max ----
# --------------------------------------------------------------------------------------------

focal_vr = "flowering"
focal_cov = "dens"
focal_min_max = "max"


## PROCESS ## 

# Set initial values on all cpus
# Make a list for every parameter (or vector if parameters are just single numbers).
# All lists have as many entries as there are chains.
# If this is the first run, you provide initial parameter values, if this is not the first run,
# read in as initial values the last values of the previous run

# Set up parallel environment
sfInit(parallel = TRUE, cpus = ncpus, slaveOutfile = "Output/IBM_Progress.txt")
sfExport(list = c(ls(), ".Random.seed"))

sfLibrary("mgcv", character.only=TRUE)
sfLibrary("crch", character.only=TRUE)
sfLibrary("snowfall", character.only=TRUE)
# Every iteration is over "restart" number of steps at the end of each iteration,
# Results of each set of iterations are saved to a file
ibm_vertedero_flowering_dens_max <- sfClusterApplyLB(1:ncpus, startIBM)

save(ibm_vertedero_flowering_dens_max, file = paste0("Output/IBM_", population, "_", focal_vr, "_", focal_cov, "_", focal_min_max, ".RData"))
sfStop()