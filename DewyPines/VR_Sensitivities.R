####################################

# This script calculates the scaled sensitivities of all three populations of dewy pine in southern Spain according to Morris et al. 2020
# The lambdas were provided by Eva Conquet

# Date: 07.06.2024
# Author: Esin Ickin

######################################


#  POP: VERTEDERO ######################################################

# 0) Prepare session ##################
rm(list=ls())

library(dplyr)

# setwd



# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)
str(df)

levels(df$focal_cov) # covariates
levels(df$focal_vr) # vital rates
levels(df$focal_min_max) # max or min covariates
levels(df$population) # now we have 3 populations
# maybe calculate sens of three pop separately then take the mean

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_WithCovariates.csv")
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")


# Calculate mean density per square to standardize covariates
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")
nbSquares = aggregate(quadratID ~ time + site, data = droso, FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, data = droso, function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, data = density_per_square, function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Climate variables time series to standardize covariates
summerT_timeseries     = aggregate(summerT_unscaled ~ site + time, data = droso, mean)
prevwinterT_timeseries = aggregate(prevwinterT_unscaled ~ site + time, data = droso, mean)
fallR_timeseries       = aggregate(fallR_unscaled ~ site + time, data = droso, mean)
prevfallR_timeseries   = aggregate(prevfallR_unscaled ~ site + time, data = droso, mean)


# Subset data to the three focal natural populations
droso_natural = droso[which(droso$site %in% c("Vertedero", "SierraCarboneraY5", "SierraRetinY5")), ]


# Seedbank data
droso_seedbank = read.csv("Data/droso_SeedBank_NoDormancyLoss.csv")


# Number of flowers
seeds_per_flower = 9.8



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate

population = "Vertedero"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)


# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)



# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)


## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
V.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))

## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
V.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# filter df for prevfallR and population and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & focal_vr == "flowering" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevfallR
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevfallR.flowering <- data.frame(variable = "prevfallR",
                                  lambda_max = max.prevfallR$lambda_max,
                                  lambda_min = min.prevfallR$lambda_min,
                                  vr="flowering"
)

# Sensitivities
V.prevfallR.flowering <- abs((prevfallR.flowering$lambda_max-prevfallR.flowering$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# filter df for prevwinterT and population and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & focal_vr == "flowering" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevwinterT
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterT.flowering <- data.frame(variable = "prevwinterT",
                                    lambda_max = max.prevwinterT$lambda_max,
                                    lambda_min = min.prevwinterT$lambda_min,
                                    vr="flowering"
)

# Sensitivities of population V to fall rain
V.prevwinterT.flowering <- abs((prevwinterT.flowering$lambda_max-prevwinterT.flowering$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))


## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
V.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities
V.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
V.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
V.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "Vertedero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities
V.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))




# POP: SIERRA CARBONERA Y5 ##############################################

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)
str(df)

levels(df$focal_cov) # covariates
levels(df$focal_min_max) # max or min covariates
levels(df$population) # now we have 3 populations
# maybe calculate sens of three pop separately then take the mean

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_WithCovariates.csv")
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")


# Calculate mean density per square to standardize covariates
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")
nbSquares = aggregate(quadratID ~ time + site, data = droso, FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, data = droso, function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, data = density_per_square, function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Climate variables time series to standardize covariates
summerT_timeseries     = aggregate(summerT_unscaled ~ site + time, data = droso, mean)
prevwinterT_timeseries = aggregate(prevwinterT_unscaled ~ site + time, data = droso, mean)
fallR_timeseries       = aggregate(fallR_unscaled ~ site + time, data = droso, mean)
prevfallR_timeseries   = aggregate(prevfallR_unscaled ~ site + time, data = droso, mean)


# Subset data to the three focal natural populations
droso_natural = droso[which(droso$site %in% c("Vertedero", "SierraCarboneraY5", "SierraRetinY5")), ]


# Seedbank data
droso_seedbank = read.csv("Data/droso_SeedBank_NoDormancyLoss.csv")


# Number of flowers
seeds_per_flower = 9.8



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "SierraCarboneraY5"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)


# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# 3) Sensitivity to fallR of population "SierraCarboneraY5" (SC) ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)


## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities of population V to fall rain
SC.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities of population V to fall rain
SC.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR of population "SierraCarboneraY5" (SC) ############

## 4.1 flowering -----------------------------
# filter df for prevfallR and population and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & focal_vr == "flowering" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevfallR
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevfallR.flowering <- data.frame(variable = "prevfallR",
                                  lambda_max = max.prevfallR$lambda_max,
                                  lambda_min = min.prevfallR$lambda_min,
                                  vr="flowering"
)

# Sensitivities of population V to fall rain
SC.prevfallR.flowering <- abs((prevfallR.flowering$lambda_max-prevfallR.flowering$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# filter df for prevwinterT and population and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & focal_vr == "flowering" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevwinterT
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterT.flowering <- data.frame(variable = "prevwinterT",
                                    lambda_max = max.prevwinterT$lambda_max,
                                    lambda_min = min.prevwinterT$lambda_min,
                                    vr="flowering"
)

# Sensitivities of population V to fall rain
SC.prevwinterT.flowering <- abs((prevwinterT.flowering$lambda_max-prevwinterT.flowering$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))


## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
SC.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))




## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
SC.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
SC.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
SC.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities of population V to fall rain
SC.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))




# POP: SIERRA RETIN Y5 ##########################################


# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_WithCovariates.csv")
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")


# Calculate mean density per square to standardize covariates
droso$quadratID = paste(droso$transect, droso$subQuadrat, sep = "_")
nbSquares = aggregate(quadratID ~ time + site, data = droso, FUN = function(x) length(unique(x)))
density_per_square = aggregate(abLarge_unscaled ~ quadratID + time + site, data = droso, function(x) unique(x))
yearly_density_per_square = aggregate(abLarge_unscaled ~ time + site, data = density_per_square, function(x) sum(x))
yearly_density_per_square$abLarge_unscaled = yearly_density_per_square$abLarge_unscaled/nbSquares$quadratID


# Climate variables time series to standardize covariates
summerT_timeseries     = aggregate(summerT_unscaled ~ site + time, data = droso, mean)
prevwinterT_timeseries = aggregate(prevwinterT_unscaled ~ site + time, data = droso, mean)
fallR_timeseries       = aggregate(fallR_unscaled ~ site + time, data = droso, mean)
prevfallR_timeseries   = aggregate(prevfallR_unscaled ~ site + time, data = droso, mean)


# Subset data to the three focal natural populations
droso_natural = droso[which(droso$site %in% c("Vertedero", "SierraCarboneraY5", "SierraRetinY5")), ]


# Seedbank data
droso_seedbank = read.csv("Data/droso_SeedBank_NoDormancyLoss.csv")


# Number of flowers
seeds_per_flower = 9.8



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "SierraRetinY5"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)


# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)



# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities of population V to fall rain
SR.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities of population V to fall rain
SR.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR of population "SierraRetinY5" (SR) ############

## 4.1 flowering -----------------------------
# filter df for prevfallR and population and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & focal_vr == "flowering" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevfallR
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevfallR.flowering <- data.frame(variable = "prevfallR",
                                  lambda_max = max.prevfallR$lambda_max,
                                  lambda_min = min.prevfallR$lambda_min,
                                  vr="flowering"
)

# Sensitivities of population V to fall rain
SR.prevfallR.flowering <- abs((prevfallR.flowering$lambda_max-prevfallR.flowering$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# filter df for prevwinterT and population and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & focal_vr == "flowering" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max prevwinterT
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterT.flowering <- data.frame(variable = "prevwinterT",
                                    lambda_max = max.prevwinterT$lambda_max,
                                    lambda_min = min.prevwinterT$lambda_min,
                                    vr="flowering"
)

# Sensitivities of population V to fall rain
SR.prevwinterT.flowering <- abs((prevwinterT.flowering$lambda_max-prevwinterT.flowering$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))

## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
SR.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))



## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
SR.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
SR.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
SR.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities
SR.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))




# POP: Buejo ##########################################

# no effect of TSF
# flowering is affected by prevwinterR but in natural populations it's affected by prevfallR and prevwinterT

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_disturbed.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_DormancyLoss.csv")

seeds_per_flower = 9.8

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
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso, mean)


# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "Bujeo"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)

# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# prevwinterR
min_prevwinterR=min(droso_pop$prevwinterR_unscaled, na.rm = T)
max_prevwinterR=max(droso_pop$prevwinterR_unscaled, na.rm = T)

sd_prevwinterR=sd(droso_pop$prevwinterR_unscaled, na.rm = T)
mean_prevwinterR=mean(droso_pop$prevwinterR_unscaled, na.rm = T)

# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities
B.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
B.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# no prevfallR

## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# no prevwinterT
## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
B.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
B.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
B.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities
B.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities of population V to fall rain
B.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))



# 8) Sensitivity to prevwinterR ##########################

## 8.1 flowering -----------------------------
# filter df for prevwinter and population and only select cols that we need for now
prevwinterR <- filter(df, focal_cov == "prevwinterR" & focal_vr == "flowering"  & population == "Bujeo") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max
max.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterR)[colnames(max.prevwinterR)=="lambda"] <- "lambda_max"

# filter for min prevwinterR
min.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterR)[colnames(min.prevwinterR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterR.flowering <- data.frame(variable = "prevwinterR",
                             lambda_max = max.prevwinterR$lambda_max,
                             lambda_min = min.prevwinterR$lambda_min,
                             vr="flowering"
)

# Sensitivities
B.prevwinterR.flowering <- abs((prevwinterR.flowering$lambda_max-prevwinterR.flowering$lambda_min)/((max_prevwinterR-min_prevwinterR)/sd_prevwinterR))




# POP: MonteraTorero ##########################################

# no effect of TSF
# flowering is affected by prevwinterR but in natural populations it's affected by prevfallR and prevwinterT

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_disturbed.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_DormancyLoss.csv")

seeds_per_flower = 9.8

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
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso, mean)



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "MonteraTorero"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)

# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# prevwinterR
min_prevwinterR=min(droso_pop$prevwinterR_unscaled, na.rm = T)
max_prevwinterR=max(droso_pop$prevwinterR_unscaled, na.rm = T)

sd_prevwinterR=sd(droso_pop$prevwinterR_unscaled, na.rm = T)
mean_prevwinterR=mean(droso_pop$prevwinterR_unscaled, na.rm = T)

# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities
MT.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
MT.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# no prevfallR

## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# no prevwinterT
## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
MT.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
MT.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
MT.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities
MT.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))

## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities of population V to fall rain
MT.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))



# 8) Sensitivity to prevwinterR ##########################

## 8.1 flowering -----------------------------
# filter df for prevwinter and population and only select cols that we need for now
prevwinterR <- filter(df, focal_cov == "prevwinterR" & focal_vr == "flowering"  & population == "MonteraTorero") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max
max.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterR)[colnames(max.prevwinterR)=="lambda"] <- "lambda_max"

# filter for min prevwinterR
min.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterR)[colnames(min.prevwinterR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterR.flowering <- data.frame(variable = "prevwinterR",
                                    lambda_max = max.prevwinterR$lambda_max,
                                    lambda_min = min.prevwinterR$lambda_min,
                                    vr="flowering"
)

# Sensitivities
MT.prevwinterR.flowering <- abs((prevwinterR.flowering$lambda_max-prevwinterR.flowering$lambda_min)/((max_prevwinterR-min_prevwinterR)/sd_prevwinterR))



# POP: Prisoneros ##########################################

# no effect of TSF
# flowering is affected by prevwinterR but in natural populations it's affected by prevfallR and prevwinterT

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_disturbed.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_DormancyLoss.csv")

seeds_per_flower = 9.8

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
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso, mean)



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "Prisoneros"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)

# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# prevwinterR
min_prevwinterR=min(droso_pop$prevwinterR_unscaled, na.rm = T)
max_prevwinterR=max(droso_pop$prevwinterR_unscaled, na.rm = T)

sd_prevwinterR=sd(droso_pop$prevwinterR_unscaled, na.rm = T)
mean_prevwinterR=mean(droso_pop$prevwinterR_unscaled, na.rm = T)

# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities
P.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
P.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# no prevfallR

## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# no prevwinterT
## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
P.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
P.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
P.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities
P.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))

## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities of population V to fall rain
P.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))




# 8) Sensitivity to prevwinterR ##########################

## 8.1 flowering -----------------------------
# filter df for prevwinter and population and only select cols that we need for now
prevwinterR <- filter(df, focal_cov == "prevwinterR" & focal_vr == "flowering"  & population == "Prisoneros") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max
max.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterR)[colnames(max.prevwinterR)=="lambda"] <- "lambda_max"

# filter for min prevwinterR
min.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterR)[colnames(min.prevwinterR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterR.flowering <- data.frame(variable = "prevwinterR",
                                    lambda_max = max.prevwinterR$lambda_max,
                                    lambda_min = min.prevwinterR$lambda_min,
                                    vr="flowering"
)

# Sensitivities
P.prevwinterR.flowering <- abs((prevwinterR.flowering$lambda_max-prevwinterR.flowering$lambda_min)/((max_prevwinterR-min_prevwinterR)/sd_prevwinterR))



# POP: Retin ##########################################

# no effect of TSF
# flowering is affected by prevwinterR but in natural populations it's affected by prevfallR and prevwinterT

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_disturbed.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_DormancyLoss.csv")

seeds_per_flower = 9.8

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
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso, mean)



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "Retin"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)

# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# prevwinterR
min_prevwinterR=min(droso_pop$prevwinterR_unscaled, na.rm = T)
max_prevwinterR=max(droso_pop$prevwinterR_unscaled, na.rm = T)

sd_prevwinterR=sd(droso_pop$prevwinterR_unscaled, na.rm = T)
mean_prevwinterR=mean(droso_pop$prevwinterR_unscaled, na.rm = T)

# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities
R.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
R.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# no prevfallR

## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# no prevwinterT
## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
R.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities of population V to fall rain
R.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
R.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities
R.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))

## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities
R.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))



# 8) Sensitivity to prevwinterR ##########################

## 8.1 flowering -----------------------------
# filter df for prevwinter and population and only select cols that we need for now
prevwinterR <- filter(df, focal_cov == "prevwinterR" & focal_vr == "flowering"  & population == "Retin") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max
max.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterR)[colnames(max.prevwinterR)=="lambda"] <- "lambda_max"

# filter for min prevwinterR
min.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterR)[colnames(min.prevwinterR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterR.flowering <- data.frame(variable = "prevwinterR",
                                    lambda_max = max.prevwinterR$lambda_max,
                                    lambda_min = min.prevwinterR$lambda_min,
                                    vr="flowering"
)

# Sensitivities
R.prevwinterR.flowering <- abs((prevwinterR.flowering$lambda_max-prevwinterR.flowering$lambda_min)/((max_prevwinterR-min_prevwinterR)/sd_prevwinterR))


# POP: SCarbDist ##########################################

# no effect of TSF
# flowering is affected by prevwinterR but in natural populations it's affected by prevfallR and prevwinterT

# 1) Load data ###########################
df=read.csv("VR_results_df.csv",stringsAsFactors = T)

# 2) Covariates ###########################################
# code provided by Eva Conquet

# Dewy-pine data
droso = read.csv("Data/droso_disturbed.csv")
droso_full = read.csv("Data/dataDroso2022.csv")
droso_full$quadratID = paste(droso_full$transect, droso_full$subQuadrat, sep = "_")

droso_seedbank = read.csv("Data/droso_SeedBank_DormancyLoss.csv")

seeds_per_flower = 9.8

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
prevwinterR_timeSeries = aggregate(prevwinterR_unscaled ~ time + site, 
                                   data = droso, mean)



# Get temperature (next summer and previous winter),
# rainfall (next and previous fall), density, size, and TSF values 
# for each projection scenario and each vital rate
levels(factor(droso$site))

population = "SCarbDist"
droso_pop = droso[which(droso$site == population & droso$time != 2022), ]

#  Min, max, standard deviation, and mean of each covariate

# Next summer mean max. daily temperature
min_summerT = min(droso_pop$summerT_unscaled, na.rm = T)
max_summerT = max(droso_pop$summerT_unscaled, na.rm = T)

mean_summerT  = mean(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)
sd_summerT    = sd(summerT_timeseries$summerT_unscaled[which(summerT_timeseries$site == population)], na.rm = T)


# Previous winter mean max. daily temperature
min_prevwinterT = min(droso_pop$prevwinterT_unscaled, na.rm = T)
max_prevwinterT = max(droso_pop$prevwinterT_unscaled, na.rm = T)

mean_prevwinterT  = mean(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)
sd_prevwinterT    = sd(prevwinterT_timeseries$prevwinterT_unscaled[which(prevwinterT_timeseries$site == population)], na.rm = T)


# Next fall cumulative rainfall
min_fallR = min(droso_pop$fallR_unscaled, na.rm = T)
max_fallR = max(droso_pop$fallR_unscaled, na.rm = T)

mean_fallR  = mean(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)
sd_fallR    = sd(fallR_timeseries$fallR_unscaled[which(fallR_timeseries$site == population)], na.rm = T)


# Previous fall cumulative rainfall
min_prevfallR = min(droso_pop$prevfallR_unscaled, na.rm = T)
max_prevfallR = max(droso_pop$prevfallR_unscaled, na.rm = T)

mean_prevfallR  = mean(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)
sd_prevfallR    = sd(prevfallR_timeseries$prevfallR_unscaled[which(prevfallR_timeseries$site == population)], na.rm = T)


# Density
min_dens = min(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)         
max_dens = max(yearly_density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)     

sd_dens   = sd(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)
mean_dens = mean(density_per_square$abLarge_unscaled[which(yearly_density_per_square$site == population)], na.rm = T)


# Size
# We first calculate the mean size for each year and take the min and max
# to avoid having multiple years where the min and max occur.

yearly_mean_size = aggregate(size_unscaled ~ time, data = droso_pop, mean)

min_size = min(yearly_mean_size$size_unscaled, na.rm = T)
max_size = max(yearly_mean_size$size_unscaled, na.rm = T)

sd_size   = sd(droso_pop$size_unscaled, na.rm = T)
mean_size = mean(droso_pop$size_unscaled, na.rm = T)

# TSF
min_TSF = min(droso_pop$TSFcont_unscaled, na.rm = T)
max_TSF = max(droso_pop$TSFcont_unscaled, na.rm = T)

sd_TSF   = sd(droso_pop$TSFcont_unscaled, na.rm = T)
mean_TSF = mean(droso_pop$TSFcont_unscaled, na.rm = T)

# prevwinterR
min_prevwinterR=min(droso_pop$prevwinterR_unscaled, na.rm = T)
max_prevwinterR=max(droso_pop$prevwinterR_unscaled, na.rm = T)

sd_prevwinterR=sd(droso_pop$prevwinterR_unscaled, na.rm = T)
mean_prevwinterR=mean(droso_pop$prevwinterR_unscaled, na.rm = T)

# 3) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_vr)
levels(df$population)
levels(df$focal_min_max)

## 3.1 flowering -----------------------------
# no fallR in flowering

## 3.2 non_repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "non_repro_survival" & population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.non_repro_surv <- data.frame(variable = "fallR",
                                   lambda_max = max.fallR$lambda_max,
                                   lambda_min = min.fallR$lambda_min,
                                   vr="non_repro_survival"
)

# Filter out rows with infinite values in lambda columns
#fallR.non_repro_surv <- fallR.non_repro_surv[is.finite(fallR.non_repro_surv$lambda_max) & is.finite(fallR.non_repro_surv$lambda_min), ]

# Sensitivities
SCD.fallR.non_repro_surv <- abs((fallR.non_repro_surv$lambda_max-fallR.non_repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 3.3 repro_survival -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & focal_vr == "repro_survival" & population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max fallR
max.fallR <- fallR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# filter for min fallR 
min.fallR <- fallR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
fallR.repro_surv <- data.frame(variable = "fallR",
                               lambda_max = max.fallR$lambda_max,
                               lambda_min = min.fallR$lambda_min,
                               vr="repro_survival"
)

# Sensitivities
SCD.fallR.repro_surv <- abs((fallR.repro_surv$lambda_max-fallR.repro_surv$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 4) Sensitivity to prevfallR ############

## 4.1 flowering -----------------------------
# no prevfallR

## 4.2 non_repro_survival -----------------------
# no prevfall in this vital rate

## 4.3 repro_survival -----------------------
# no prevfall in this vital rate


# 5) Sensitivity to prevwinterT ##########################

## 5.1 flowering -----------------------------
# no prevwinterT
## 5.2 non_repro_survival--------------------
# no prevwinterT in this vital rate

## 5.3 repro_survival--------------------------
# no prevwinterT in this vital rate


# 6) Sensitivity to summerT ##########################

## 6.1 flowering -----------------------------
# no summer T in this vital rate

## 6.2 non_repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "non_repro_survival" &  population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.non_repro_surv <- data.frame(variable = "summerT",
                                     lambda_max = max.summerT$lambda_max,
                                     lambda_min = min.summerT$lambda_min,
                                     vr="non_repro_survival"
)

# Sensitivities of population V to fall rain
SCD.summerT.non_repro_surv <- abs((summerT.non_repro_surv$lambda_max-summerT.non_repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


## 6.3 repro_survival------------------------
# filter df for summerT and population and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & focal_vr == "repro_survival" &  population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max summerT
max.summerT <- summerT %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# filter for min summerT 
min.summerT <- summerT %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
summerT.repro_surv <- data.frame(variable = "summerT",
                                 lambda_max = max.summerT$lambda_max,
                                 lambda_min = min.summerT$lambda_min,
                                 vr="repro_survival"
)

# Sensitivities
SCD.summerT.repro_surv <- abs((summerT.repro_surv$lambda_max-summerT.repro_surv$lambda_min)/((max_summerT-min_summerT)/sd_summerT))


# 7) Sensitivity to density ##########################

## 7.1 flowering -----------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "flowering"  & population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.flowering <- data.frame(variable = "dens",
                             lambda_max = max.dens$lambda_max,
                             lambda_min = min.dens$lambda_min,
                             vr="flowering"
)

# Sensitivities of population V to fall rain
SCD.dens.flowering <- abs((dens.flowering$lambda_max-dens.flowering$lambda_min)/((max_dens-min_dens)/sd_dens))


## 7.2 non_repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "non_repro_survival"  & population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min dens
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.non_repro_surv <- data.frame(variable = "dens",
                                  lambda_max = max.dens$lambda_max,
                                  lambda_min = min.dens$lambda_min,
                                  vr="non_repro_survival"
)

# Sensitivities
SCD.dens.non_repro_surv <- abs((dens.non_repro_surv$lambda_max-dens.non_repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))

## 7.3 repro_survival------------------------
# filter df for dens and population and only select cols that we need for now
dens <- filter(df, focal_cov == "dens" & focal_vr == "repro_survival" &  population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max dens
max.dens <- dens %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# filter for min 
min.dens <- dens %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
dens.repro_surv <- data.frame(variable = "dens",
                              lambda_max = max.dens$lambda_max,
                              lambda_min = min.dens$lambda_min,
                              vr="repro_survival"
)

# Sensitivities
SCD.dens.repro_surv <- abs((dens.repro_surv$lambda_max-dens.repro_surv$lambda_min)/((max_dens-min_dens)/sd_dens))




# 8) Sensitivity to prevwinterR ##########################

## 8.1 flowering -----------------------------
# filter df for prevwinter and population and only select cols that we need for now
prevwinterR <- filter(df, focal_cov == "prevwinterR" & focal_vr == "flowering"  & population == "SCarbDist") %>%
  select(focal_cov, focal_min_max, lambda, population,focal_vr)

# filter for max
max.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "max")

# rename col from lambda to lambda_max
colnames(max.prevwinterR)[colnames(max.prevwinterR)=="lambda"] <- "lambda_max"

# filter for min prevwinterR
min.prevwinterR <- prevwinterR %>%
  filter(focal_min_max == "min")

# rename col from lambda to lambda_min
colnames(min.prevwinterR)[colnames(min.prevwinterR)=="lambda"] <- "lambda_min"

# merge lambda_max and lambda_min
prevwinterR.flowering <- data.frame(variable = "prevwinterR",
                                    lambda_max = max.prevwinterR$lambda_max,
                                    lambda_min = min.prevwinterR$lambda_min,
                                    vr="flowering"
)

# Sensitivities
SCD.prevwinterR.flowering <- abs((prevwinterR.flowering$lambda_max-prevwinterR.flowering$lambda_min)/((max_prevwinterR-min_prevwinterR)/sd_prevwinterR))



################################################################################################################################################################################


# SAVE ALL SENSITIVITIES ##############################

# group sensitivities according to the drivers and vital rates

# and take the means across the three sites

# fallR non_repro_surv
df_fallR_non_repro_surv=data.frame(V=V.fallR.non_repro_surv,
                                   SC=SC.fallR.non_repro_surv,
                                   SR=SR.fallR.non_repro_surv,
                                   B=B.fallR.non_repro_surv,
                                   MT=MT.fallR.non_repro_surv,
                                   P=P.fallR.non_repro_surv,
                                   R=R.fallR.non_repro_surv,
                                   SCD=SCD.fallR.non_repro_surv)

df_fallR_non_repro_surv$sens=rowMeans(df_fallR_non_repro_surv[,c("V","SC","SR",
                                                                 "B","MT","P","R","SCD")],)

# fallR repro_surv
df_fallR_repro_surv=data.frame(V=V.fallR.repro_surv,
                               SC=SC.fallR.repro_surv,
                               SR=SR.fallR.repro_surv,
                               B=B.fallR.repro_surv,
                               MT=MT.fallR.repro_surv,
                               P=P.fallR.repro_surv,
                               R=R.fallR.repro_surv,
                               SCD=SCD.fallR.repro_surv)

df_fallR_repro_surv$sens=rowMeans(df_fallR_repro_surv[,c("V","SC","SR",
                                                         "B","MT","P","R","SCD")],)


# prevfallR flowering
df_prevfallR_flowering=data.frame(V=V.prevfallR.flowering,
                                  SC=SC.prevfallR.flowering,
                                  SR=SR.prevfallR.flowering)

df_prevfallR_flowering$sens=rowMeans(df_prevfallR_flowering[,c("V","SC","SR")],)

# prevwinterT flowering
df_prevwinterT_flowering=data.frame(V=V.prevwinterT.flowering,
                                    SC=SC.prevwinterT.flowering,
                                    SR=SR.prevwinterT.flowering)

df_prevwinterT_flowering$sens=rowMeans(df_prevwinterT_flowering[,c("V","SC","SR")],)

# prevwinterR flowering
df_prevwinterR_flowering=data.frame(B=B.prevwinterR.flowering,
                                    MT=MT.prevwinterR.flowering,
                                    P=P.prevwinterR.flowering,
                                    R=R.prevwinterR.flowering,
                                    SCD=SCD.prevwinterR.flowering)
df_prevwinterR_flowering$sens=rowMeans(df_prevwinterR_flowering[,c("B","MT","P","R","SCD")])

# summerT non_repro_surv
df_summerT_non_repro_surv=data.frame(V=V.summerT.non_repro_surv,
                                     SC=SC.summerT.non_repro_surv,
                                     SR=SR.summerT.non_repro_surv,
                                     
                                     B=B.summerT.non_repro_surv,
                                     MT=MT.summerT.non_repro_surv,
                                     P=P.summerT.non_repro_surv,
                                     R=R.summerT.non_repro_surv,
                                     SCD=SCD.summerT.non_repro_surv)

df_summerT_non_repro_surv$sens=rowMeans(df_summerT_non_repro_surv[,c("V","SC","SR","B","MT","P","R","SCD")],)

# summerT repro_surv
df_summerT_repro_surv=data.frame(V=V.summerT.repro_surv,
                                 SC=SC.summerT.repro_surv,
                                 SR=SR.summerT.repro_surv,
                                 
                                 B=B.summerT.repro_surv,
                                 MT=MT.summerT.repro_surv,
                                 P=P.summerT.repro_surv,
                                 R=R.summerT.repro_surv,
                                 SCD=SCD.summerT.repro_surv)

df_summerT_repro_surv$sens=rowMeans(df_summerT_repro_surv[,c("V","SC","SR","B","MT","P","R","SCD")],)

# dens flowering
df_dens_flowering=data.frame(V=V.dens.flowering,
                             SC=SC.dens.flowering,
                             SR=SR.dens.flowering,
                             
                             B=B.dens.flowering,
                             MT=MT.dens.flowering,
                             P=P.dens.flowering,
                             R=R.dens.flowering,
                             SCD=SCD.dens.flowering)

df_dens_flowering$sens=rowMeans(df_dens_flowering[,c("V","SC","SR","B","MT","P","R","SCD")],)

# dens non_repro_surv
df_dens_non_repro_surv=data.frame(V=V.dens.non_repro_surv,
                                  SC=SC.dens.non_repro_surv,
                                  SR=SR.dens.non_repro_surv,
                                  
                                  B=B.dens.non_repro_surv,
                                  MT=MT.dens.non_repro_surv,
                                  P=P.dens.non_repro_surv,
                                  R=R.dens.non_repro_surv,
                                  SCD=SCD.dens.non_repro_surv)

df_dens_non_repro_surv$sens=rowMeans(df_dens_non_repro_surv[,c("V","SC","SR","B","MT","P","R","SCD")],)

# dens repro_surv
df_dens_repro_surv=data.frame(V=V.dens.repro_surv,
                              SC=SC.dens.repro_surv,
                              SR=SR.dens.repro_surv,
                              
                              B=B.dens.repro_surv,
                              MT=MT.dens.repro_surv,
                              P=P.dens.repro_surv,
                              R=R.dens.repro_surv,
                              SCD=SCD.dens.repro_surv)

df_dens_repro_surv$sens=rowMeans(df_dens_repro_surv[,c("V","SC","SR","B","MT","P","R","SCD")],)



VR.sens=data.frame(study.doi="Conquet et al. in prep",
                   year.of.publication=2024,
                   group="Plants",
                   species="Drosophyllum lusitanicum",
                   continent="Europe",
                   driver=rep(c("rain","rain",
                                "rain","rain",
                                "temperature","temperature",
                                "density","density","density","rain"),each=100),
                   driver.type=rep(c("C","C",
                                     "C","C",
                                     "C","C",
                                     "D","D","D","C"),each=100),
                   stage.age=rep(c("non-reproductive","reproductive",
                                   "reproductive","reproductive",
                                   "non-reproductive","reproductive",
                                   "non-reproductive","non-reproductive","reproductive","reproductive"),each=100),
                   vital.rates=rep(c("survival","survival",
                                     "flowering","flowering",
                                     "survival","survival",
                                     "flowering","survival","survival",
                                     "flowering"),each=100),
                   sens=c(df_fallR_non_repro_surv$sens,
                          df_fallR_repro_surv$sens,
                          df_prevfallR_flowering$sens,
                          df_prevwinterT_flowering$sens,
                          df_summerT_non_repro_surv$sens,
                          df_summerT_repro_surv$sens,
                          df_dens_flowering$sens,
                          df_dens_non_repro_surv$sens,
                          df_dens_repro_surv$sens,
                          df_prevwinterR_flowering$sens),
                   mat=2, # age at sexual maturity [Source: Paniw et al.]
                   n.vr=5, # number of vital rates with covariates
                   n.pam=156, # number of parameters
                   dens=1, # density dependence in it?
                   biotic_interactions=0, # any other biotic interactions?
                   lambda.sim=1, # was lambda calculated analytically (0) or using simulation (1)?
                   study.length=11)


write.csv(VR.sens,"SensVR_DewyPines.csv",row.names = F)


df=read.csv("/Users/esinickin/Desktop/DewyPines/SensVR_DewyPines.csv")


