###########################################################

# This script contains the main analysis but 
# with species only at one site or with only one population

#############################################################
rm(list=ls())

# set wd
 setwd("/Users/esinickin/Documents/Master Thesis/MainAnalyses/SciAdv")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)
library(MuMIn)

df=read.csv("AllSens.csv") # 41

head(df)
unique(df$species)

### EDA & edit data ##################

# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform generation time
df$mat=log(df$mat)
# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# remove species that have spatially-explicit models (multiple sites and/or multiple populations)

multi_sites_spp <- c("Pyrrhula pyrrhula", "Lophophanes cristatus", "Certhia familiaris", "Sitta europaea", "Prunella modularis", "Linaria cannabina", "Turdus torquatus", "Prunella collaris", "Perisoreus infaustus", # birds
                   "Dracocephalum austriacum",
                   "Drosophyllum lusitanicum",
                   "Fagus sylvatica", "Pinus halepensis", "Pinus nigra", "Pinus pinaster", "Pinus pinea", "Pinus sylvestris", "Pinus uncinata", "Quercus faginea", "Quercus ilex", "Quercus robur/petraea", "Quercus suber" # plants
                   )
length(multi_sites_spp) # 22 species to remove?!

# or remove only species with multiple populations

multi_pop_spp <- c("Perisoreus infaustus","Dracocephalum austriacum",
                    "Drosophyllum lusitanicum")

# now filter both
no_multi_sites_spp_df <- df %>%
  filter(!species %in% multi_sites_spp)

no_multi_pop_spp_df <- df %>%
  filter(!species %in% multi_pop_spp)



# GLMM: Sens to all Climate Variables (no species with multiple sites) ##################################################################
# climate and not splitting by rain/temp
# cov*dens + mat + par + random effects

# filter no_multi_sites_spp_df so that we only include climatic covariates (what do you think?)
climate_df=filter(no_multi_sites_spp_df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","PET","Winterlength")) # remove SAM, PET, and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any


climate_df$driver=factor(climate_df$driver)

climate_df$densV2=climate_df$dens
climate_df$densV2[climate_df$biotic_interactions==1]=1

climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

climate_df$lambda.sim=factor(climate_df$lambda.sim)


# aggregate sensitivities modeled per lambda type and whether density was included in models to check distribution

sum=aggregate(sens~lambda.sim+dens,function(x) length(unique(x)), data=climate_df)

ggplot(sum, aes(fill=lambda.sim, y=sens, x=dens)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Density effects")+
  ylab("# of |S| values")


# Change order of facto levels

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))


m1 <- glmer(sens ~ cov*dens + mat+ n.vr +par.per.vr+ (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)


summary(m1) 

length(unique(climate_df$species))

2.088/(2.088+0.256+0.108+0.001+0.718)
0.256/(2.088+0.256+0.108+0.001+0.718)
0.108/(2.088+0.256+0.108+0.001+0.718)
0.001/(2.088+0.256+0.108+0.001+0.718)
0.718/(2.088+0.256+0.108+0.001+0.718)


r.squaredGLMM(m1)


# GLMM: Sens to all Climate Variables (no species with multiple populations) ##################################################################
# climate and not splitting by rain/temp
# cov*dens + mat + par + random effects

# filter no_multi_pop_spp_df so that we only include climatic covariates (what do you think?)
climate_df=filter(no_multi_pop_spp_df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","PET","Winterlength")) # remove SAM, PET, and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any


climate_df$driver=factor(climate_df$driver)

climate_df$densV2=climate_df$dens
climate_df$densV2[climate_df$biotic_interactions==1]=1

climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

climate_df$lambda.sim=factor(climate_df$lambda.sim)


# aggregate sensitivities modeled per lambda type and whether density was included in models to check distribution

sum=aggregate(sens~lambda.sim+dens,function(x) length(unique(x)), data=climate_df)

ggplot(sum, aes(fill=lambda.sim, y=sens, x=dens)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Density effects")+
  ylab("# of |S| values")


# Change order of facto levels

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))


m1 <- glmer(sens ~ cov*dens + mat+ n.vr +par.per.vr+ (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)


summary(m1) 

length(unique(climate_df$species))

r.squaredGLMM(m1)


1.843/(1.843+0.189+0.001+0.001+0.818)
0.189/(1.843+0.189+0.001+0.001+0.818)
0.001/(1.843+0.189+0.001+0.001+0.818)
0.001/(1.843+0.189+0.001+0.001+0.818)
0.818/(1.843+0.189+0.001+0.001+0.818)


