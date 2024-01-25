####################################

# This script calculates the scaled sensitivities of the dewy pine according to Morris et al. 2020
# The lambdas were provided by Eva Conquet

# Date: 9.11.2023
# Author: Esin Ickin

######################################

# 0) Prepare session ##################
rm(list=ls())

library(dplyr)

setwd("~/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/DewyPines")

# 1) Load data ###########################
df=read.csv("results_df.csv",stringsAsFactors = T)
str(df)
names(df$focal_cov)

levels(df$focal_cov) # covariates
levels(df$focal_full) # covariates but with full meaning
levels(df$focal_min_max) # max or min covariates
levels(df$other_mean_obs) # covariation included or not (mean)

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






# 3) Sensitivity to density ################
# aka "above ground density of nlarge individuals"

## 3.1 no covariation -----------------------------
# filter df for dens and only select cols that we need for now
dens <- filter(df, focal_cov == "dens") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max dens no cov
max.dens <- dens %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# remove max_min column
max.dens <- select(max.dens, - focal_min_max)


# filter for min dens no cov
min.dens <- dens %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# remove max_min column
min.dens <- select(min.dens, - focal_min_max)

# merge lambda_max and lambda_min
dens.no.cov <- data.frame(variable = "dens",
                   cov = 0,
                   lambda_max = max.dens$lambda_max,
                   lambda_min = min.dens$lambda_min
                   )
# Filter out rows with infinite values in lambda columns
dens.no.cov <- dens.no.cov[is.finite(dens.no.cov$lambda_max) & is.finite(dens.no.cov$lambda_min), ]

dens.no.cov$sens <- abs((dens.no.cov$lambda_max-dens.no.cov$lambda_min)/((max_dens-min_dens)/sd_dens))






## 3.2 covariation -----------------------------
# filter df for dens and only select cols that we need for now
dens <- filter(df, focal_cov == "dens") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max dens cov
max.dens <- dens %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.dens)[colnames(max.dens)=="lambda"] <- "lambda_max"

# remove max_min column
max.dens <- select(max.dens, - focal_min_max)


# filter for min dens cov
min.dens <- dens %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.dens)[colnames(min.dens)=="lambda"] <- "lambda_min"

# remove max_min column
min.dens <- select(min.dens, - focal_min_max)

# merge lambda_max and lambda_min
dens.cov <- data.frame(variable = "dens",
                   cov = 1,
                   lambda_max = max.dens$lambda_max,
                   lambda_min = min.dens$lambda_min
)
# Filter out rows with infinite values in lambda columns
dens.cov <- dens.cov[is.finite(dens.cov$lambda_max) & is.finite(dens.cov$lambda_min), ]

dens.cov$sens <- abs((dens.cov$lambda_max-dens.cov$lambda_min)/((max_dens-min_dens)/sd_dens))

# 4) Sensitivity to fallR ######################
levels(df$focal_cov)
levels(df$focal_full)

## 4.1 no covariation -----------------------------
# filter df for fallR and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

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

fallR.no.cov$sens <- abs((fallR.no.cov$lambda_max-fallR.no.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


## 4.2 covariation -----------------------------
# filter df for fallR and only select cols that we need for now
fallR <- filter(df, focal_cov == "fallR") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

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

fallR.cov$sens <- abs((fallR.cov$lambda_max-fallR.cov$lambda_min)/((max_fallR-min_fallR)/sd_fallR))


# 5) Sensitivity to prevfallR ############
## 5.1 no covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR") %>%
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

prevfallR.no.cov$sens <- abs((prevfallR.no.cov$lambda_max-prevfallR.no.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))


## 5.2 covariation ----------------
levels(df$focal_cov)
levels(df$focal_full)

# filter df for prevfallR and only select cols that we need for now
prevfallR <- filter(df, focal_cov == "prevfallR") %>%
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

prevfallR.cov$sens <- abs((prevfallR.cov$lambda_max-prevfallR.cov$lambda_min)/((max_prevfallR-min_prevfallR)/sd_prevfallR))

# 6) Sensitivity to prevwinterT ##########################
levels(df$focal_cov)

## 6.1 no covariation ---------------
# filter df for prevwinterT and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT") %>%
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

prevwinterT.no.cov$sens <- abs((prevwinterT.no.cov$lambda_max-prevwinterT.no.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))

## 6.2 covariation ---------------
# filter df for prevwinterT and only select cols that we need for now
prevwinterT <- filter(df, focal_cov == "prevwinterT") %>%
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

prevwinterT.cov$sens <- abs((prevwinterT.cov$lambda_max-prevwinterT.cov$lambda_min)/((max_prevwinterT-min_prevwinterT)/sd_prevwinterT))


# 7) Sensitivity to summerT ##########################
levels(df$focal_cov)

## 7.1 no covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT") %>%
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

summerT.no.cov$sens <- abs((summerT.no.cov$lambda_max-summerT.no.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))

## 7.2 covariation ---------------
# filter df for summerT and only select cols that we need for now
summerT <- filter(df, focal_cov == "summerT") %>%
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

summerT.cov$sens <- abs((summerT.cov$lambda_max-summerT.cov$lambda_min)/((max_summerT-min_summerT)/sd_summerT))

# 8) Sensitivity to TSF ##########################
levels(df$focal_cov)

## 8.1 no covariation ---------------
# filter df for TSF and only select cols that we need for now
TSF <- filter(df, focal_cov == "TSF") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max TSF no cov
max.TSF <- TSF %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.TSF)[colnames(max.TSF)=="lambda"] <- "lambda_max"

# remove max_min column
max.TSF <- select(max.TSF, - focal_min_max)

# filter for min TSF no cov
min.TSF <- TSF %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.TSF)[colnames(min.TSF)=="lambda"] <- "lambda_min"

# remove max_min column
min.TSF <- select(min.TSF, - focal_min_max)

# merge lambda_max and lambda_min
TSF.no.cov <- data.frame(variable = "TSF",
                             cov = 0,
                             lambda_max = max.TSF$lambda_max,
                             lambda_min = min.TSF$lambda_min
)
# Filter out rows with infinite values in lambda columns
TSF.no.cov <- TSF.no.cov[is.finite(TSF.no.cov$lambda_max) & is.finite(TSF.no.cov$lambda_min), ]

TSF.no.cov$sens <- abs((TSF.no.cov$lambda_max-TSF.no.cov$lambda_min)/((max_TSF-min_TSF)/sd_TSF))


## 8.2 covariation ---------------
# filter df for TSF and only select cols that we need for now
TSF <- filter(df, focal_cov == "TSF") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max TSF cov
max.TSF <- TSF %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_max
colnames(max.TSF)[colnames(max.TSF)=="lambda"] <- "lambda_max"

# remove max_min column
max.TSF <- select(max.TSF, - focal_min_max)

# filter for min TSF cov
min.TSF <- TSF %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "obs")

# rename col from lambda to lambda_min
colnames(min.TSF)[colnames(min.TSF)=="lambda"] <- "lambda_min"

# remove max_min column
min.TSF <- select(min.TSF, - focal_min_max)

# merge lambda_max and lambda_min
TSF.cov <- data.frame(variable = "TSF",
                         cov = 1,
                         lambda_max = max.TSF$lambda_max,
                         lambda_min = min.TSF$lambda_min
)
# Filter out rows with infinite values in lambda columns
TSF.cov <- TSF.no.cov[is.finite(TSF.cov$lambda_max) & is.finite(TSF.cov$lambda_min), ]

TSF.cov$sens <- abs((TSF.cov$lambda_max-TSF.no.cov$lambda_min)/((max_TSF-min_TSF)/sd_TSF))

## 9.1 no covariation ---------------
# filter df for size and only select cols that we need for now
size <- filter(df, focal_cov == "size") %>%
  select(focal_cov, focal_min_max, other_mean_obs, lambda)

# filter for max TSF no cov
max.TSF <- TSF %>%
  filter(focal_min_max == "max") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_max
colnames(max.TSF)[colnames(max.TSF)=="lambda"] <- "lambda_max"

# remove max_min column
max.TSF <- select(max.TSF, - focal_min_max)

# filter for min TSF no cov
min.TSF <- TSF %>%
  filter(focal_min_max == "min") %>%
  filter(other_mean_obs == "mean")

# rename col from lambda to lambda_min
colnames(min.TSF)[colnames(min.TSF)=="lambda"] <- "lambda_min"

# remove max_min column
min.TSF <- select(min.TSF, - focal_min_max)

# merge lambda_max and lambda_min
TSF.no.cov <- data.frame(variable = "TSF",
                         cov = 0,
                         lambda_max = max.TSF$lambda_max,
                         lambda_min = min.TSF$lambda_min
)
# Filter out rows with infinite values in lambda columns
TSF.no.cov <- TSF.no.cov[is.finite(TSF.no.cov$lambda_max) & is.finite(TSF.no.cov$lambda_min), ]

TSF.no.cov$sens <- abs((TSF.no.cov$lambda_max-TSF.no.cov$lambda_min)/((max_TSF-min_TSF)/sd_TSF))












# 10) Merge dataframes together, from each pick 10 rows ####################################
rm(final.sens)
# had to exlude sensitivities to TSF because too many NAs
# same with sensitivities to prevfallR
final.sens <- rbind(dens.no.cov[1:10,],dens.cov[1:10,],
                    fallR.no.cov[1:10,],fallR.cov[1:10,],
                    prevwinterT.no.cov[1:10,],prevwinterT.cov[1:10,],
                    summerT.no.cov[1:10,],summerT.cov[1:10,])

# remove lambda columns
final.sens <- final.sens %>%
  select(-lambda_max) %>%
  select(-lambda_min)

# add columns
final.sens$study.doi <- "Conquet et al. in prep."
final.sens$year.of.publication <- "2024"
final.sens$group <- "Plants"
final.sens$species <- "Drosophyllum lusitanicum"
final.sens$continent <- "Europe"
final.sens$stage.age <- "all"
final.sens$vital.rates <- "all"

# rename column name from variable to driver
colnames(final.sens)[colnames(final.sens)=="variable"] <- "driver"
# density is biotic driver type
# the rest is climatic aka abiotic
final.sens$driver.type <- ifelse(final.sens$driver == "dens", "B","A")

# rearrange columns so that it fits to the rest of the species
final.sens <- final.sens %>%
  select(study.doi,year.of.publication,group,species,continent,driver,driver.type,stage.age,vital.rates,sens,cov)

# 11) Save output #####################
write.csv(final.sens, "DewyPines_Sens.csv", row.names = F)






