# This script calculates the population growth rates and sensitivities to temperature and rainfall with and without accounting for covariation
# For 13 tree species from Garcia-Callejas et al. 2016

# Species 1: Conifers (group)
# Species 2: Deciduous (group)
# Species 3: Fagus sylvatica
# Species 4: Juniperus thurifera
# Species 5: Pinus halpensis
# Species 6: Pinus nigra
# Species 7: Pinus pinaster
# Species 8: Pinus pinea
# Species 9: Pinus sylvestris
# Species 10: Pinus uncinata
# Species 11: Quercus faginea
# Species 12: Quercus ilex
# Species 13: Quercus pyrenaica
# Species 14: Quercus robur/petraea
# Species 15: Quercus suber
# Species 16: Sclerophyllous (group)

# remove species 4 and 13 because they don't work

# Prepare session #####################
# set wd
setwd("~/Files_SpanishTrees")

rm(list=ls())

library(dplyr)

# climate covariates
clim=read.csv("ClimateCovariates.csv")

# SIMULATION 1: ###########################################
## 1) TEMPERATURE #####################
df.temp=NULL # empty df to store sensitivities
df.mu.temp=NULL # empty df to store sensitivities averaging across space


for (i in c(1:3,5:12,14:16)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
  keep=which(l.mu.minTcov<3&l.mu.minTcov>-3)
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<3&l.mu.maxTcov>-3)
    
  keep=keep[keep%in%keep2]

  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<3&l.mu.minTnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T------------
  df.t=read.csv(paste0("MaxTempNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<3&l.mu.maxTnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp=rbind(df.temp,data.frame(sens.cov,sens.no.cov,
                         species=paste0("sp", i)))
  ## averaging across space --------
  df.mu.temp=rbind(df.mu.temp,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                               sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                            species=paste0("sp", i)))
}
## add variables -------
df.temp$type="Temperature"
df.mu.temp$type="Temperature"
df.mu.temp$sim=1 # this was for the first run


## 2) RAINFALL ############################
# empty df to store sens
df.rain=NULL
df.mu.rain=NULL # sens averaged across space

rm(sens.cov,sens.no.cov) # remove items above

for (i in c(1:3,5:12,14:16)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<3&l.mu.minPcov>-3)
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<3&l.mu.maxPcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<3&l.mu.minPnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<3&l.mu.maxPnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain=rbind(df.rain,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
 
   ## averaging across space -------------
  df.mu.rain=rbind(df.mu.rain,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
}

#### add variables -----------
df.rain$type="Precipitation"
df.mu.rain$type="Precipitation"
df.mu.rain$sim=1


# SIMULATION 2:  ###########################################
## 1) TEMPERATURE #####################
# define what we want to keep and remove the rest:
df_keep=c("df.mu.temp","df.mu.rain")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

df.temp2=NULL # empty df to store sensitivities
df.mu.temp2=NULL # empty df to store sensitivities averaging across space

for (i in c(1:3,5:12,14:16)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run2/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<3&l.mu.minTcov>-3)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<3&l.mu.maxTcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<3&l.mu.minTnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<3&l.mu.maxTnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp2=rbind(df.temp2,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
  ## averaging across space --------
  df.mu.temp2=rbind(df.mu.temp2,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
}
## add variables -------
df.temp2$type="Temperature"
df.mu.temp2$type="Temperature"
df.mu.temp2$sim=2





## 2) RAINFALL ############################
# empty df to store sens
df.rain2=NULL
df.mu.rain2=NULL # sens averaged across space

rm(sens.cov,sens.no.cov) # remove items above

for (i in c(1:3,5:12,14:16)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<3&l.mu.minPcov>-3)
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<3&l.mu.maxPcov>-3)
    
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<3&l.mu.minPnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<3&l.mu.maxPnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain2=rbind(df.rain2,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain2=rbind(df.mu.rain2,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
}

#### add variables -----------
df.rain2$type="Precipitation"
df.mu.rain2$type="Precipitation"
df.mu.rain2$sim=2

# SIMULATION 3: ###########################################
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp3=NULL # empty df to store sensitivities
df.mu.temp3=NULL # empty df to store sensitivities averaging across space

for (i in c(1:3,5:12,14:16)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run3/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<3&l.mu.minTcov>-3)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<3&l.mu.maxTcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<3&l.mu.minTnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<3&l.mu.maxTnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp3=rbind(df.temp3,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  ## averaging across space --------
  df.mu.temp3=rbind(df.mu.temp3,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}
## add variables -------
df.temp3$type="Temperature"
df.mu.temp3$type="Temperature"
df.mu.temp3$sim=3

## 2) RAINFALL ############################
# empty df to store sens
df.rain3=NULL
df.mu.rain3=NULL # sens averaged across space

rm(sens.cov,sens.no.cov) # remove items above

for (i in c(1:3,5:12,14:16)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<3&l.mu.minPcov>-3)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<3&l.mu.maxPcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<3&l.mu.minPnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<3&l.mu.maxPnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain3=rbind(df.rain3,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain3=rbind(df.mu.rain3,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}

#### add variables -----------
df.rain3$type="Precipitation"
df.mu.rain3$type="Precipitation"
df.mu.rain3$sim=3


# SIMULATION 4:#####################################
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2",
          "df.mu.temp3","df.mu.rain3")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp4=NULL # empty df to store sensitivities
df.mu.temp4=NULL # empty df to store sensitivities averaging across space

for (i in c(1:3,5:12,14:16)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run4/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<3&l.mu.minTcov>-3)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<3&l.mu.maxTcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<3&l.mu.minTnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<3&l.mu.maxTnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp4=rbind(df.temp4,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  ## averaging across space --------
  df.mu.temp4=rbind(df.mu.temp4,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}
## add variables -------
df.temp4$type="Temperature"
df.mu.temp4$type="Temperature"
df.mu.temp4$sim=4

## 2) RAINFALL ############################
# empty df to store sens
df.rain4=NULL
df.mu.rain4=NULL # sens averaged across space

rm(sens.cov,sens.no.cov) # remove items above

for (i in c(1:3,5:12,14:16)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<3&l.mu.minPcov>-3)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<3&l.mu.maxPcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<3&l.mu.minPnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<3&l.mu.maxPnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain4=rbind(df.rain4,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain4=rbind(df.mu.rain4,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}

#### add variables -----------
df.rain4$type="Precipitation"
df.mu.rain4$type="Precipitation"
df.mu.rain4$sim=4


# SIMULATION 5:#####################################
# doesn't work yet because one folder is empty
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2",
          "df.mu.temp3","df.mu.rain3",
          "df.mu.temp4","df.mu.rain4")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp5=NULL # empty df to store sensitivities
df.mu.temp5=NULL # empty df to store sensitivities averaging across space

for (i in c(1:3,5:12,14:16)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run5/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<3&l.mu.minTcov>-3)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<3&l.mu.maxTcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<3&l.mu.minTnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<3&l.mu.maxTnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp5=rbind(df.temp5,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  ## averaging across space --------
  df.mu.temp5=rbind(df.mu.temp5,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}
## add variables -------
df.temp5$type="Temperature"
df.mu.temp5$type="Temperature"
df.mu.temp5$sim=5

## 2) RAINFALL ############################
# empty df to store sens
df.rain5=NULL
df.mu.rain5=NULL # sens averaged across space

rm(sens.cov,sens.no.cov) # remove items above

for (i in c(1:3,5:12,14:16)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<3&l.mu.minPcov>-3)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<3&l.mu.maxPcov>-3)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<3&l.mu.minPnoCov>-3)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<3&l.mu.maxPnoCov>-3)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain5=rbind(df.rain5,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain5=rbind(df.mu.rain5,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
}

#### add variables -----------
df.rain5$type="Precipitation"
df.mu.rain5$type="Precipitation"
df.mu.rain5$sim=5




# Rearrange dataframes ###########################
df <- rbind(df.mu.temp,df.mu.temp2,df.mu.temp3,df.mu.temp4,df.mu.temp5,
                 df.mu.rain,df.mu.rain2,df.mu.rain3,df.mu.rain4,df.mu.rain5)

df_cov <- df %>% 
  select(-sens.no.cov)

colnames(df_cov)[colnames(df_cov)=="sens.cov"] <- "sens"
df_cov$cov <- 1

df_no_cov <- df %>%
  select(-sens.cov)

colnames(df_no_cov)[colnames(df_no_cov)=="sens.no.cov"] <- "sens"
df_no_cov$cov <- 0

df2 <- rbind(df_no_cov,df_cov)

df2 <- select(df2, species, type, sens, cov,sim)
df2$full.species <- c("Conifers","Deciduous","Fagus sylvatica","Pinus halepensis","Pinus nigra","Pinus pinaster","Pinus pinea","Pinus sylvestris","Pinus uncinata",
                            "Quercus faginea","Quercus ilex",
                            "Quercus robur/petraea","Quercus suber","Sclerophyllous")

df2$study.doi <- "10.1093/jpe/rtw081"
df2$group <- "Plants"
df2$year.of.publication <- "2016"
df2$driver.type <- "A"
df2$continent <- "Europe"
df2$stage.age <- "all"
df2$vital.rates <- "all"
colnames(df2)[colnames(df2)=="type"] <- "driver"

df2 <- df2 %>% select(study.doi,year.of.publication,group,full.species,continent,driver,driver.type,stage.age,vital.rates,sens,cov)
colnames(df2)[colnames(df2)=="full.species"] <- "species"

# SAVE OUTPUT ####################
write.csv(df2, "SensSpainTrees.csv",row.names = F)
