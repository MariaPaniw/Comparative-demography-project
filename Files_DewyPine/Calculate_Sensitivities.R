####################################

# This script calculates the scaled sensitivities of all three populations of dewy pine in southern Spain according to Morris et al. 2020
# The lambdas were provided by Eva Conquet

# Date: 15.03.2024
# Author: Esin Ickin

######################################


#  POP: VERTEDERO ######################################################

# 0) Prepare session ##################
rm(list=ls())

library(dplyr)

# setwd
setwd("~/OneDrive - Universität Zürich UZH/Master Thesis/pert_analyses/DewyPines")


# 1) Load data ###########################
df=read.csv("results_df.csv",stringsAsFactors = T)
str(df)

levels(df$focal_cov) # covariates
levels(df$focal_full) # covariates but with full meaning
levels(df$focal_min_max) # max or min covariates
levels(df$other_mean_obs) # covariation included or not (mean)

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


# Correction factors
pop_year_corr = read.csv("Data/PopYearCorrection_Natural.csv")

corr = data.frame(sigmaS = c(0.45, 0.45, 0.35, 0.35, 0.33, 
                             0.84, 0.84, 0.82, 0.8, 0.82),
                  TSF = rep(c("0", "1", "2", "3", "4"), 2),
                  LS = c(rep("HLS", 5), rep("LLS", 5)))
corr = corr[which(corr$LS == "LLS"), ]



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





# 3) Sensitivity to fallR of population "Vertedero" (V) ######################
levels(df$focal_cov)
levels(df$focal_full)

levels(df$population)

## 3.1 no covariation -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.no.cov <- data.frame(variable = "fallR",
                           cov = 0,
                           lambda_max = max.fallR$lambda_max,
                           lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.no.cov <- fallR.no.cov[is.finite(fallR.no.cov$lambda_max) & is.finite(fallR.no.cov$lambda_min), ]

# Sensitivities of population SC to fall rain
V.fallR <- abs((fallR.no.cov$lambda_max-fallR.no.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))
hist(V.fallR)
sum(is.na(V.fallR))
sum(is.infinite(V.fallR))


## 3.2 covariation -----------------------------
# filter df for fallR and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.cov <- data.frame(variable = "fallR",
                        cov = 1,
                        lambda_max = max.fallR$lambda_max,
                        lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.cov <- fallR.cov[is.finite(fallR.cov$lambda_max) & is.finite(fallR.cov$lambda_min), ]

# save sensitivities with covariation of population SC to fall rain
V.fallR.cov <- abs((fallR.cov$lambda_max-fallR.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))
hist(V.fallR)







# 4) Sensitivity to prevfallR of population "Vertedero" (V) ############
## 4.1 no covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.no.cov <- data.frame(variable = "prevfallR",
                               cov = 0,
                               lambda_max = max.prevfallR$lambda_max,
                               lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.no.cov <- prevfallR.no.cov[is.finite(prevfallR.no.cov$lambda_max) & is.finite(prevfallR.no.cov$lambda_min), ]

V.prevfallR <- abs((prevfallR.no.cov$lambda_max-prevfallR.no.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR"  & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.cov <- data.frame(variable = "prevfallR",
                            cov = 1,
                            lambda_max = max.prevfallR$lambda_max,
                            lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.cov <- prevfallR.cov[is.finite(prevfallR.cov$lambda_max) & is.finite(prevfallR.cov$lambda_min), ]

V.prevfallR.cov <- abs((prevfallR.cov$lambda_max-prevfallR.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))




# 5) Sensitivity to prevwinterT of population "Vertedero" (V) ##########################
levels(df$focal_cov)


## 5.1 no covariation ---------------
# filter df for prevwinterT and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT no cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT no cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.no.cov <- data.frame(variable = "prevwinterT",
                                 cov = 0,
                                 lambda_max = max.prevwinterT$lambda_max,
                                 lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.no.cov <- prevwinterT.no.cov[is.finite(prevwinterT.no.cov$lambda_max) & is.finite(prevwinterT.no.cov$lambda_min), ]

V.prevwinterT <- abs((prevwinterT.no.cov$lambda_max-prevwinterT.no.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))

## 5.2 covariation ---------------
# filter df for prevwinterT and pop and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.cov <- data.frame(variable = "prevwinterT",
                              cov = 1,
                              lambda_max = max.prevwinterT$lambda_max,
                              lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.cov <- prevwinterT.cov[is.finite(prevwinterT.cov$lambda_max) & is.finite(prevwinterT.cov$lambda_min), ]

V.prevwinterT.cov <- abs((prevwinterT.cov$lambda_max-prevwinterT.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))



# 6) Sensitivity to summerT of population "Vertedero" (V) ##########################
levels(df$focal_cov)

## 6.1 no covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population =="Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT no cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.no.cov <- data.frame(variable = "summerT",
                             cov = 0,
                             lambda_max = max.summerT$lambda_max,
                             lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.no.cov <- summerT.no.cov[is.finite(summerT.no.cov$lambda_max) & is.finite(summerT.no.cov$lambda_min), ]

V.summerT <- abs((summerT.no.cov$lambda_max-summerT.no.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))

## 6.2 covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population == "Vertedero") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.cov <- data.frame(variable = "summerT",
                          cov = 1,
                          lambda_max = max.summerT$lambda_max,
                          lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.cov <- summerT.cov[is.finite(summerT.cov$lambda_max) & is.finite(summerT.cov$lambda_min), ]

V.summerT.cov <- abs((summerT.cov$lambda_max-summerT.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))



# POP: Sierra Carbonera Y5 ##############################################


# 1) Load data ###########################
df=read.csv("results_df.csv",stringsAsFactors = T)
str(df)

levels(df$focal_cov) # covariates
levels(df$focal_full) # covariates but with full meaning
levels(df$focal_min_max) # max or min covariates
levels(df$other_mean_obs) # covariation included or not (mean)

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


# Correction factors
pop_year_corr = read.csv("Data/PopYearCorrection_Natural.csv")

corr = data.frame(sigmaS = c(0.45, 0.45, 0.35, 0.35, 0.33, 
                             0.84, 0.84, 0.82, 0.8, 0.82),
                  TSF = rep(c("0", "1", "2", "3", "4"), 2),
                  LS = c(rep("HLS", 5), rep("LLS", 5)))
corr = corr[which(corr$LS == "LLS"), ]



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



# 3) Sensitivity to fallR of population "SierracarboneraY5" (SC) ######################
levels(df$focal_cov)
levels(df$focal_full)

levels(df$population)

## 3.1 no covariation -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.no.cov <- data.frame(variable = "fallR",
                           cov = 0,
                           lambda_max = max.fallR$lambda_max,
                           lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.no.cov <- fallR.no.cov[is.finite(fallR.no.cov$lambda_max) & is.finite(fallR.no.cov$lambda_min), ]

# Sensitivities of population SC to fall rain
SC.fallR <- abs((fallR.no.cov$lambda_max-fallR.no.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))
hist(SC.fallR)
sum(is.na(SC.fallR))
sum(is.infinite(SC.fallR))


## 3.2 covariation -----------------------------
# filter df for fallR and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.cov <- data.frame(variable = "fallR",
                        cov = 1,
                        lambda_max = max.fallR$lambda_max,
                        lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.cov <- fallR.cov[is.finite(fallR.cov$lambda_max) & is.finite(fallR.cov$lambda_min), ]

# save sensitivities with covariation of population SC to fall rain
SC.fallR.cov <- abs((fallR.cov$lambda_max-fallR.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR)) # if we filter, which we have to because otherwise there are NAs and infinite values, then there's only 3
sum(is.na(SC.fallR.cov))
sum(is.infinite(SC.fallR.cov))
hist(SC.fallR.cov)


# 4) Sensitivity to prevfallR of population "SierracarboneraY5" (SC) ############
## 4.1 no covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.no.cov <- data.frame(variable = "prevfallR",
                               cov = 0,
                               lambda_max = max.prevfallR$lambda_max,
                               lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.no.cov <- prevfallR.no.cov[is.finite(prevfallR.no.cov$lambda_max) & is.finite(prevfallR.no.cov$lambda_min), ]

SC.prevfallR <- abs((prevfallR.no.cov$lambda_max-prevfallR.no.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR"  & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.cov <- data.frame(variable = "prevfallR",
                            cov = 1,
                            lambda_max = max.prevfallR$lambda_max,
                            lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.cov <- prevfallR.cov[is.finite(prevfallR.cov$lambda_max) & is.finite(prevfallR.cov$lambda_min), ]

SC.prevfallR.cov <- abs((prevfallR.cov$lambda_max-prevfallR.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))



# 5) Sensitivity to prevwinterT of population "SierraCarboneraY5" (SC) ##########################
levels(df$focal_cov)


## 5.1 no covariation ---------------
# filter df for prevwinterT and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT no cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT no cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.no.cov <- data.frame(variable = "prevwinterT",
                                 cov = 0,
                                 lambda_max = max.prevwinterT$lambda_max,
                                 lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.no.cov <- prevwinterT.no.cov[is.finite(prevwinterT.no.cov$lambda_max) & is.finite(prevwinterT.no.cov$lambda_min), ]

SC.prevwinterT <- abs((prevwinterT.no.cov$lambda_max-prevwinterT.no.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))

## 5.2 covariation ---------------
# filter df for prevwinterT and pop and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.cov <- data.frame(variable = "prevwinterT",
                              cov = 1,
                              lambda_max = max.prevwinterT$lambda_max,
                              lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.cov <- prevwinterT.cov[is.finite(prevwinterT.cov$lambda_max) & is.finite(prevwinterT.cov$lambda_min), ]

SC.prevwinterT.cov <- abs((prevwinterT.cov$lambda_max-prevwinterT.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))




# 6) Sensitivity to summerT of pop SC ##########################
levels(df$focal_cov)

## 6.1 no covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population =="SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT no cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.no.cov <- data.frame(variable = "summerT",
                             cov = 0,
                             lambda_max = max.summerT$lambda_max,
                             lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.no.cov <- summerT.no.cov[is.finite(summerT.no.cov$lambda_max) & is.finite(summerT.no.cov$lambda_min), ]

SC.summerT <- abs((summerT.no.cov$lambda_max-summerT.no.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))

## 6.2 covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population == "SierraCarboneraY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.cov <- data.frame(variable = "summerT",
                          cov = 1,
                          lambda_max = max.summerT$lambda_max,
                          lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.cov <- summerT.cov[is.finite(summerT.cov$lambda_max) & is.finite(summerT.cov$lambda_min), ]

SC.summerT.cov <- abs((summerT.cov$lambda_max-summerT.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))




# POP: Sierra Retin Y5 ##########################################


# 1) Load data ###########################
df=read.csv("results_df.csv",stringsAsFactors = T)
str(df)

levels(df$focal_cov) # covariates
levels(df$focal_full) # covariates but with full meaning
levels(df$focal_min_max) # max or min covariates
levels(df$other_mean_obs) # covariation included or not (mean)

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


# Correction factors
pop_year_corr = read.csv("Data/PopYearCorrection_Natural.csv")

corr = data.frame(sigmaS = c(0.45, 0.45, 0.35, 0.35, 0.33, 
                             0.84, 0.84, 0.82, 0.8, 0.82),
                  TSF = rep(c("0", "1", "2", "3", "4"), 2),
                  LS = c(rep("HLS", 5), rep("LLS", 5)))
corr = corr[which(corr$LS == "LLS"), ]



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



# 3) Sensitivity to fallR of population "SierraRetinY5" (SR) ######################
levels(df$focal_cov)
levels(df$focal_full)

levels(df$population)

## 3.1 no covariation -----------------------------
# filter df for fallR and population and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.no.cov <- data.frame(variable = "fallR",
                           cov = 0,
                           lambda_max = max.fallR$lambda_max,
                           lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.no.cov <- fallR.no.cov[is.finite(fallR.no.cov$lambda_max) & is.finite(fallR.no.cov$lambda_min), ]

# Sensitivities of population SC to fall rain
SR.fallR <- abs((fallR.no.cov$lambda_max-fallR.no.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))
hist(SC.fallR)
sum(is.na(SC.fallR))
sum(is.infinite(SC.fallR))


## 3.2 covariation -----------------------------
# filter df for fallR and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda, population)

# filter for max fallR no cov
max.fallR <- fallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.fallR)[colnames(max.fallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.fallR <- select(max.fallR, - focal_min_max)

# filter for min fallR no cov
min.fallR <- fallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.fallR)[colnames(min.fallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.fallR <- select(min.fallR, - focal_min_max)

# merge lambda_max and lambda_min
fallR.cov <- data.frame(variable = "fallR",
                        cov = 1,
                        lambda_max = max.fallR$lambda_max,
                        lambda_min = min.fallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
fallR.cov <- fallR.cov[is.finite(fallR.cov$lambda_max) & is.finite(fallR.cov$lambda_min), ]

# save sensitivities with covariation of population SC to fall rain
SR.fallR.cov <- abs((fallR.cov$lambda_max-fallR.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR)) # if we filter, which we have to because otherwise there are NAs and infinite values, then there's only 3
sum(is.na(SC.fallR.cov))
sum(is.infinite(SC.fallR.cov))
hist(SC.fallR.cov)


# 4) Sensitivity to prevfallR of population "SierraRetinY5" (SR) ############
## 4.1 no covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.no.cov <- data.frame(variable = "prevfallR",
                               cov = 0,
                               lambda_max = max.prevfallR$lambda_max,
                               lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.no.cov <- prevfallR.no.cov[is.finite(prevfallR.no.cov$lambda_max) & is.finite(prevfallR.no.cov$lambda_min), ]

SR.prevfallR <- abs((prevfallR.no.cov$lambda_max-prevfallR.no.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 4.2 covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR"  & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevfallR no cov
max.prevfallR <- prevfallR %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevfallR)[colnames(max.prevfallR)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevfallR <- select(max.prevfallR, - focal_min_max)

# filter for min prevfallR no cov
min.prevfallR <- prevfallR %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevfallR)[colnames(min.prevfallR)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevfallR <- select(min.prevfallR, - focal_min_max)

# merge lambda_max and lambda_min
prevfallR.cov <- data.frame(variable = "prevfallR",
                            cov = 1,
                            lambda_max = max.prevfallR$lambda_max,
                            lambda_min = min.prevfallR$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevfallR.cov <- prevfallR.cov[is.finite(prevfallR.cov$lambda_max) & is.finite(prevfallR.cov$lambda_min), ]

SR.prevfallR.cov <- abs((prevfallR.cov$lambda_max-prevfallR.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))



# 5) Sensitivity to prevwinterT of population "SierraRetinY5" (SR) ##########################
levels(df$focal_cov)


## 5.1 no covariation ---------------
# filter df for prevwinterT and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT no cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT no cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.no.cov <- data.frame(variable = "prevwinterT",
                                 cov = 0,
                                 lambda_max = max.prevwinterT$lambda_max,
                                 lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.no.cov <- prevwinterT.no.cov[is.finite(prevwinterT.no.cov$lambda_max) & is.finite(prevwinterT.no.cov$lambda_min), ]

SR.prevwinterT <- abs((prevwinterT.no.cov$lambda_max-prevwinterT.no.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))

## 5.2 covariation ---------------
# filter df for prevwinterT and pop and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max prevwinterT cov
max.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.prevwinterT)[colnames(max.prevwinterT)=="lambda"] <- "lambda_max"

# remove max_min column
max.prevwinterT <- select(max.prevwinterT, - focal_min_max)

# filter for min prevwinterT cov
min.prevwinterT <- prevwinterT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.prevwinterT)[colnames(min.prevwinterT)=="lambda"] <- "lambda_min"

# remove max_min column
min.prevwinterT <- select(min.prevwinterT, - focal_min_max)

# merge lambda_max and lambda_min
prevwinterT.cov <- data.frame(variable = "prevwinterT",
                              cov = 1,
                              lambda_max = max.prevwinterT$lambda_max,
                              lambda_min = min.prevwinterT$lambda_min
)
# Filter out rows with infinite values in lambda columns
prevwinterT.cov <- prevwinterT.cov[is.finite(prevwinterT.cov$lambda_max) & is.finite(prevwinterT.cov$lambda_min), ]

SR.prevwinterT.cov <- abs((prevwinterT.cov$lambda_max-prevwinterT.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))




# 6) Sensitivity to summerT of pop SC ##########################
levels(df$focal_cov)

## 6.1 no covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population =="SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT no cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.no.cov <- data.frame(variable = "summerT",
                             cov = 0,
                             lambda_max = max.summerT$lambda_max,
                             lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.no.cov <- summerT.no.cov[is.finite(summerT.no.cov$lambda_max) & is.finite(summerT.no.cov$lambda_min), ]

SR.summerT <- abs((summerT.no.cov$lambda_max-summerT.no.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))

## 6.2 covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT" & population == "SierraRetinY5") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max summerT cov
max.summerT <- summerT %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.summerT)[colnames(max.summerT)=="lambda"] <- "lambda_max"

# remove max_min column
max.summerT <- select(max.summerT, - focal_min_max)

# filter for min summerT no cov
min.summerT <- summerT %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.summerT)[colnames(min.summerT)=="lambda"] <- "lambda_min"

# remove max_min column
min.summerT <- select(min.summerT, - focal_min_max)

# merge lambda_max and lambda_min
summerT.cov <- data.frame(variable = "summerT",
                          cov = 1,
                          lambda_max = max.summerT$lambda_max,
                          lambda_min = min.summerT$lambda_min
)
# Filter out rows with infinite values in lambda columns
summerT.cov <- summerT.cov[is.finite(summerT.cov$lambda_max) & is.finite(summerT.cov$lambda_min), ]

SR.summerT.cov <- abs((summerT.cov$lambda_max-summerT.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))





# SAVE ALL SENSITIVITIES ##############################

# group sensitivities according to the drivers and whether covariation is included or not

# fall (autumn) rain and no covariation
df_fallR=data.frame(sens=c(SC.fallR,
                           V.fallR
                           #SR.fallR,
                           ),
                    driver="fallR",
                    cov=0)

# fall rain with covariation
df_fallR.cov=data.frame(sens=c(SC.fallR.cov,
                           V.fallR.cov
                           #SR.fallR.cov,
),
driver="fallR",
cov=1)

# previous fall rain and no covariation
df_prevfallR=data.frame(sens=c(#SC.prevfallR,
                                SR.prevfallR,
                                V.prevfallR),
                      driver="prevfallR",
                      cov=0)

# previous fall rain with covariation
df_prevfallR.cov=data.frame(sens=c(#SC.prevfallR.cov,
                            SR.prevfallR.cov,
                            V.prevfallR.cov),
                            driver="prevfallR",
                            cov=1)

# previous winter temperature and no covariation
df_prevwinterT=data.frame(sens=c(SC.prevwinterT,
                                 SR.prevwinterT,
                                 V.prevwinterT),
                          driver="prevwinterT",
                          cov=0)

# previous winter temperature with covariation
df_prevwinterT.cov=data.frame(sens=c(SC.prevwinterT.cov,
                                 SR.prevwinterT.cov,
                                 V.prevwinterT.cov),
                          driver="prevwinterT",
                          cov=1)

# summer temperature and no covariation
df_summerT=data.frame(sens=c(SC.summerT,
                             #SR.summerT,
                             V.summerT),
                      driver="summerT",
                      cov=0)

# summer temperature with covariation
df_summerT.cov=data.frame(sens=c(SC.summerT.cov,
                             #SR.summerT.cov,
                             V.summerT.cov),
                      driver="summerT",
                      cov=1)

# SR.fallR.cov is empty, so omit SR.fallR too
# SC.prevfallR.cov is empty, so omit SC.prevfallR too
# SR.summerT.cov is empty, so omit SR.summerT too

# only population V has at least some sensitivities to all the drivers
# after talking with author Conquet, we can proceed with just this one population

# but I will omit V.prevfallR too because it only has 2 resampled sensitivities
# and it's similar to V.fallR
# same with summerT but it's similar to prevwinterT anyways

Sens_DewyPines=data.frame(study.doi="Conquet et al. in prep",
                          year.of.publication=2024,
                          group="Plants",
                          species="Drosophyllum lusitanicum",
                          continent="Europe",
                          driver=rep(c("rain","temperature"),each=20),
                          driver.type="C",
                          stage.age="all",
                          vital.rates="all",
                          sens=c(V.fallR[1:10],V.fallR.cov[1:10],
                                 V.prevwinterT[1:10],V.prevwinterT.cov[1:10]),
                          cov=rep(c(0,1),each=10),
                          gen.time=1, # generation time
                          n.vr=5, # number of vital rates with covariates
                          n.pam=156, # number of parameters
                          dens=1, # density dependence in it?
                          biotic_interactions=0) # any other biotic interactions?



write.csv(Sens_DewyPines,"Sens_DewyPines.csv",row.names = F)









