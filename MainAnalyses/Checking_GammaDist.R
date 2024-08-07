###
# Here we check whether the Gamma distribution is the best option for our main analysis

# Author of this script: Esin Ickin
# Date: 07.08.2024
###


# 0. Prepare session ############################################

rm(list=ls())

# set wd
setwd("/Users/esinickin/Desktop/MainAnalyses")

# load libraries
library(lme4)


# 1. Load data #########################################################
df=read.csv("AllSens.csv")

### 2. EDA & edit data ##################

str(df)
summary(df)

#hist(df$sens)
#hist(df$n.vr)
#hist(df$mat)
#hist(df$n.pam)


sum(is.na(df$sens))
sum(is.infinite(df$sens))

# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform age at sexual maturity
df$mat=log(df$mat)

# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# 3. GLMM: Sensitivities to all climatic variables ##################################################################

# filter df so that we only include sensitivities to climatic variables
climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","Winterlength")) # remove SAM and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see species Marmota flaviventris by Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
#climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any

climate_df$driver=factor(climate_df$driver)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

# GLMM 1: Gamma distribution
m1.gamma <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

#summary(m1.gamma)

# GLMM 2: Log-normal distribution
m2.lognormal <- lmer(log(sens) ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|group/species) , data = climate_df)

#summary(m2.lognormal)

AIC(m1.gamma,m2.lognormal)


# DHARMa ####################
library(DHARMa)

# GLMM 1:
# Simulate residuals
simulationOutput.gamma <- simulateResiduals(fittedModel = m1.gamma, n = 250)

# Residual vs. Fitted Values Plot
plotResiduals(simulationOutput.gamma)

# Q-Q Plot
plotQQunif(simulationOutput.gamma)

# Residuals vs. Key Predictors
plotResiduals(simulationOutput.gamma, form = climate_df$cov)

# Overall Residual Plot
plot(simulationOutput.gamma)


# GLMM 2:

# GLMM 1:
# Simulate residuals
simulationOutput.lognormal <- simulateResiduals(fittedModel = m2.lognormal, n = 250)

# Residual vs. Fitted Values Plot
plotResiduals(simulationOutput.lognormal)

# Q-Q Plot
plotQQunif(simulationOutput.lognormal)

# Residuals vs. Key Predictors
plotResiduals(simulationOutput.lognormal, form = climate_df$cov)

# Overall Residual Plot
plot(simulationOutput.lognormal)

