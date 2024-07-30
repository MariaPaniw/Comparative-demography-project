###
# This is a script containing the main analysis (GLMM 1 of GlobalAnalysis) but without the African striped mouse (Rhabdomys pumilio)
# to test whether the exclusion of this species (that has a monthly matrix population model and not an annual one like other species) makes a difference

# Author of this script: Esin Ickin
# Date: 30.07.2024
###


# 0. Prepare session ############################################

rm(list=ls())

# set wd
setwd("~/Documents/Master Thesis/pert_analyses")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)

## 0.1 Load data #########################################################

# load df
df=read.csv("AllSens.csv")

## 0.2 EDA & edit data ##################
str(df)

# remove Rhabdomys pumilio
# because the MPM is monthly and not annually
df=df[df$species!="Rhabdomys pumilio",]


summary(df)
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


# 1. GLMM: Sensitivities to all climatic variables ##################################################################

climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","PET","Winterlength")) # remove SAM and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

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

m1 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m1)
# as you can see, the results don't change dramatically compared to the original main analysis 
# which can be found in the script GlobalAnalysis.R



# 2. GLMM: Sens to temperature vs rain ##################################################################
# split climate further to temperature and rain


levels(factor(climate_df$driver))

# group the drivers into rain and temperature
climate_df$driver[climate_df$driver %in% c("Pbr", "Pwn", "rain","fallR","Precipitation","precipitation","Rain","Rainfall","ROS")]="rain"

climate_df$driver[climate_df$driver %in% c("lagged temperature","prevwinterT","sea ice","SST","SSTA_b","SSTA_b_l","SSTA_m","SSTA_m_l","summerT","Tat","Tbr","temperature","Temperature","Twn","SST1","SST2","SST3")]="temperature"

climate_df=filter(climate_df, !driver %in% c("SAM","Q","PET","Winterlength"))
# remove SAM (a driver that influences climate like rain and temperature, so not really a climatic driver)
# remove Q because now I can't sort it into rain or temperature
# remove winterlength again for the same reasons as above

# remove sens == 0 if there is 
min(climate_df$sens)
#climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution

# from chr to factors
climate_df$driver=factor(climate_df$driver)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species) # 36 species

length(levels(climate_df$species))
levels(climate_df$sp)


m2 <- glmer(sens ~ cov*dens + cov*driver + dens*driver + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m2)


# 3. GLMM: Vital-rate specific sensitivities ########
rm(list=ls())
dfVR=read.csv("AllSensVR.csv")


# remove Rhabdomys pumilio
# because the MPM is monthly and not annually
dfVR=dfVR[dfVR$species!="Rhabdomys pumilio",]

### 3.1 EDA & edit data ###################
str(dfVR)


sum(is.na(dfVR$sens))
sum(is.infinite(dfVR$sens))

hist(dfVR$sens)
# very skewed, mostly between 0 and 1, positive and continuous
# continuous
# always positive
# very right skewed
# probably best dist is gamma dist with log link in the GLMM

hist(dfVR$mat)
max(dfVR$mat)
min(dfVR$mat)

hist(dfVR$n.vr)
hist(dfVR$n.pam)

# add new column with mean number of parameters per vital rate
dfVR$par.per.vr=dfVR$n.pam/dfVR$n.vr

# log-transform age at sexual maturity
dfVR$mat=log(dfVR$mat)
# log-transform also number of vital rates and parameters per vital rate
dfVR$n.vr=log(dfVR$n.vr)
dfVR$par.per.vr=log(dfVR$par.per.vr)


### 3.2 group vital rates #################

# group vital rates into reproduction, survival, or trait change
levels(factor(dfVR$vital.rates))
dfVR[dfVR=="breeding probability"]="reproduction"
dfVR[dfVR=="litter probability"]="reproduction"
dfVR[dfVR=="litter size"]="reproduction"
dfVR[dfVR=="denning survival"]="survival"
dfVR[dfVR=="maturation probability"]="trait change"
dfVR[dfVR=="fecundity"]="reproduction"
dfVR[dfVR=="growth"]="trait change"
dfVR[dfVR=="recruitment"]="reproduction"
dfVR[dfVR=="offspring mass"]="trait change"
dfVR[dfVR=="transition"]="trait change"
dfVR[dfVR=="flowering"]="reproduction"
levels(factor(dfVR$vital.rates))


# group (st)ages into non-reproductive and reproductive individuals
levels(factor(dfVR$stage.age))
# non-reproductive
dfVR[dfVR=="babies"]="non-reproductive"
dfVR[dfVR=="immature"]="non-reproductive"
dfVR[dfVR=="juvenile"]="non-reproductive"
dfVR[dfVR=="non-reproductive adult"]="non-reproductive"
dfVR[dfVR=="sapling"]="non-reproductive"
dfVR[dfVR=="juveniles"]="non-reproductive"
dfVR[dfVR=="calves"]="non-reproductive"
dfVR[dfVR=="yearling"]="non-reproductive"
dfVR[dfVR=="female juvenile"]="non-reproductive"
dfVR[dfVR=="male juvenile"]="non-reproductive"
dfVR[dfVR=="pups"]="non-reproductive"
dfVR[dfVR=="subadult"]="non-reproductive"

# reproductive
dfVR[dfVR=="adult"]="reproductive"
dfVR[dfVR=="previous breeder"]="reproductive"
dfVR[dfVR=="previous non-breeder"]="reproductive"
dfVR[dfVR=="reproductive adult"]="reproductive"
dfVR[dfVR=="all"]="reproductive" 
dfVR[dfVR=="adults"]="reproductive"
dfVR[dfVR=="female adult"]="reproductive"
dfVR[dfVR=="male adult"]="reproductive"
dfVR[dfVR=="dominant adult"]="reproductive"
dfVR[dfVR=="helper adult"]="reproductive"

# reindeers and dippers have age classes from 1 to 6 and 1 to 4 respectively
# reindeer stages 2-5 are reproductive
dfVR$stage.age[dfVR$stage.age=="age class 2" & dfVR$species=="Rangifer tarandus"]="reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 3" & dfVR$species=="Rangifer tarandus"]="reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 4" & dfVR$species=="Rangifer tarandus"]="reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 5" & dfVR$species=="Rangifer tarandus"]="reproductive"
# reindeer stage 1 & 6 is non-reproductive
dfVR$stage.age[dfVR$stage.age=="age class 1" & dfVR$species=="Rangifer tarandus"]="non-reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 6" & dfVR$species=="Rangifer tarandus"]="non-reproductive"

# dipper 
dfVR$stage.age[dfVR$stage.age=="age class 1" & dfVR$species=="Cinclus cinclus"]="non-reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 2" & dfVR$species=="Cinclus cinclus"]="reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 3" & dfVR$species=="Cinclus cinclus"]="reproductive"
dfVR$stage.age[dfVR$stage.age=="age class 4" & dfVR$species=="Cinclus cinclus"]="reproductive"

# check
levels(factor(dfVR$stage.age))
levels(factor(dfVR$vital.rates))


# now put stage.age and vital.rates into one column
dfVR$new_vitalrates=paste(dfVR$stage.age, dfVR$vital.rates)
dfVR$new_vitalrates=factor(dfVR$new_vitalrates)

# pool all reproduction together
levels(dfVR$new_vitalrates)[c(1,4)] =c("reproduction","reproduction")
levels(dfVR$new_vitalrates)


## 3.3 make standardized label for drivers ##########
# for the GLMM
# rain-related drivers
dfVR$driver[dfVR$driver %in% c("Pbr", "Pwn", "rain","fallR","Precipitation","precipitation","Rain","Rainfall","ROS")]="rain"

# temperature stuff
dfVR$driver[dfVR$driver %in% c("lagged temperature","prevwinterT","sea ice","SST","SSTA_b","SSTA_b_l","SSTA_m","SSTA_m_l","summerT","Tat","Tbr","temperature","Temperature","Twn","SST1","SST2","SST3","Q")]="temperature"

# density stuff
dfVR$driver[dfVR$driver %in% c("density","intraD","IntraDens","lagged density")]="density"

# biotic drivers
dfVR$driver[dfVR$driver %in% c("Chla","food","goose abundance","interD","InterDens","lagged food","reindeer carcass availability")]="biotic"


# check
levels(factor(dfVR$driver))


# remove SAM and winterlength
dfVR=dfVR[!dfVR$driver %in% c("SAM","winterlength"),]


levels(factor(dfVR$driver))

# filter df so that we only include climatic covariates 
climate_df=filter(dfVR,driver %in% c("rain","temperature"))
levels(factor(climate_df$driver))
levels(factor(climate_df$species))
levels(factor(climate_df$new_vitalrates))

growth=climate_df[climate_df$new_vitalrates %in% c("non-reproductive trait change","reproductive trait change"),]
levels(factor(growth$species))

# remove trait change because this is what gave huge uncertainties
climate_df=droplevels(climate_df[!climate_df$new_vitalrates%in%c("reproductive trait change","non-reproductive trait change"),])

levels(climate_df$new_vitalrates)

## 3.4 GLMM: Sens to all Climate Variables ####################################

# remove sens == 0 if there is 
min(climate_df$sens)
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution

# re-define factors
climate_df$species=factor(climate_df$species) # now we have 23 species
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group)
climate_df$driver=factor(climate_df$driver)
climate_df$driver.type=factor(climate_df$driver.type)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)

m3 <- glmer(sens ~ dens*new_vitalrates + mat*new_vitalrates + n.vr + par.per.vr + (1+new_vitalrates|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m3)


