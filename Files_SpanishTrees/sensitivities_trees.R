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

# discard "species" 1, 2 and 16
# remove species 4 because no lambdas
# discard species 13 too same reason

# Prepare session #####################
# set wd
setwd("~/Documents/Master Thesis/pert_analyses/SpainTrees")

rm(list=ls())

library(dplyr)

# climate covariates
clim=read.csv("ClimateCovariates.csv")

# SIMULATION 1: ###########################################
## 1) TEMPERATURE #####################
df.temp=NULL # empty df to store sensitivities
df.mu.temp=NULL # empty df to store sensitivities averaging across space

df.l_ratio.temp=NULL # empty df to store log ratios
df.mu.l_ratio.temp=NULL # empty df to store log ratios averaging across species


# species 4
# run 1
df.t1=read.csv("MinTempCov/results.sp4.csv")
# calculate lambdas
l1=df.t1$adults_2070/df.t1$adults_2060
l2=df.t1$adults_2080/df.t1$adults_2070
l3=df.t1$adults_2090/df.t1$adults_2080

sp4_1=na.omit(cbind(l1,l2,l3))
sp4_1.l.mu.minTcov=apply(sp4_1,1,mean)

# run 2
df.t2=read.csv("MinTempCov/run2/results.sp4.csv")
# calculate lambdas
l1=df.t2$adults_2070/df.t2$adults_2060
l2=df.t2$adults_2080/df.t2$adults_2070
l3=df.t2$adults_2090/df.t2$adults_2080

sp4_2=na.omit(cbind(l1,l2,l3))
sp4_2.l.mu.minTcov=apply(sp4_2,1,mean)


# run 3
df.t3=read.csv("MinTempCov/run3/results.sp4.csv")
# calculate lambdas
l1=df.t3$adults_2070/df.t3$adults_2060
l2=df.t3$adults_2080/df.t3$adults_2070
l3=df.t3$adults_2090/df.t3$adults_2080

sp4_3=na.omit(cbind(l1,l2,l3))
sp4_3.l.mu.minTcov=apply(sp4_3,1,mean)


# run 4
df.t4=read.csv("MinTempCov/run4/results.sp4.csv")
# calculate lambdas
l1=df.t4$adults_2070/df.t4$adults_2060
l2=df.t4$adults_2080/df.t4$adults_2070
l3=df.t4$adults_2090/df.t4$adults_2080

sp4_4=na.omit(cbind(l1,l2,l3))
sp4_4.l.mu.minTcov=apply(sp4_4,1,mean)

# run 5
df.t5=read.csv("MinTempCov/run5/results.sp4.csv")
# calculate lambdas
l1=df.t5$adults_2070/df.t5$adults_2060
l2=df.t5$adults_2080/df.t5$adults_2070
l3=df.t5$adults_2090/df.t5$adults_2080

sp4_5=na.omit(cbind(l1,l2,l3))
sp4_5.l.mu.minTcov=apply(sp4_5,1,mean)

# all zeros :(

# try now with species 13
# run 1
df.t1=read.csv("MinTempCov/results.sp13.csv")
# calculate lambdas
l1=df.t1$adults_2070/df.t1$adults_2060
l2=df.t1$adults_2080/df.t1$adults_2070
l3=df.t1$adults_2090/df.t1$adults_2080

sp13_1=na.omit(cbind(l1,l2,l3))
sp13_1.l.mu.minTcov=apply(sp13_1,1,mean)

# run 2
df.t2=read.csv("MinTempCov/run2/results.sp13.csv")
# calculate lambdas
l1=df.t2$adults_2070/df.t2$adults_2060
l2=df.t2$adults_2080/df.t2$adults_2070
l3=df.t2$adults_2090/df.t2$adults_2080

sp13_2=na.omit(cbind(l1,l2,l3))
sp13_2.l.mu.minTcov=apply(sp13_2,1,mean)


# run 3
df.t3=read.csv("MinTempCov/run3/results.sp13.csv")
# calculate lambdas
l1=df.t3$adults_2070/df.t3$adults_2060
l2=df.t3$adults_2080/df.t3$adults_2070
l3=df.t3$adults_2090/df.t3$adults_2080

sp13_3=na.omit(cbind(l1,l2,l3))
sp13_3.l.mu.minTcov=apply(sp13_3,1,mean)


# run 4
df.t4=read.csv("MinTempCov/run4/results.sp13.csv")
# calculate lambdas
l1=df.t4$adults_2070/df.t4$adults_2060
l2=df.t4$adults_2080/df.t4$adults_2070
l3=df.t4$adults_2090/df.t4$adults_2080

sp13_4=na.omit(cbind(l1,l2,l3))
sp13_4.l.mu.minTcov=apply(sp13_4,1,mean)

# run 5
df.t5=read.csv("MinTempCov/run5/results.sp13.csv")
# calculate lambdas
l1=df.t5$adults_2070/df.t5$adults_2060
l2=df.t5$adults_2080/df.t5$adults_2070
l3=df.t5$adults_2090/df.t5$adults_2080

sp13_5=na.omit(cbind(l1,l2,l3))
sp13_5.l.mu.minTcov=apply(sp13_5,1,mean)

###start##############################

for (i in c(3,5:12,14:15)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
  keep=which(l.mu.minTcov<1.005&l.mu.minTcov>0.75)
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<1.005&l.mu.maxTcov>0.75)
    
  keep=keep[keep%in%keep2]

  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<1.005&l.mu.minTnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T------------
  df.t=read.csv(paste0("MaxTempNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<1.005&l.mu.maxTnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp=rbind(df.temp,data.frame(sens.cov,sens.no.cov,
                         species=paste0("sp", i)))
  
  # log ratios
  l_ratio.cov=abs(log(l.mu.maxTcov[keep]/l.mu.minTcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxTnoCov[keep]/l.mu.minTnoCov[keep]))
    
  df.l_ratio.temp=rbind(df.l_ratio.temp,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                   species=paste0("sp",i)))
  
  
  ## averaging across space --------
  df.mu.temp=rbind(df.mu.temp,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                               sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                            species=paste0("sp", i)))
  
  df.mu.l_ratio.temp=rbind(df.mu.l_ratio.temp,
                           data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                      l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                 species=paste0("sp", i)))
  
  
}

## add variables -------
df.temp$type="Temperature"
df.mu.temp$type="Temperature"
df.mu.temp$sim=1 # this was for the first run

df.l_ratio.temp$type="Temperature"
df.mu.l_ratio.temp$type="Temperature"
df.mu.l_ratio.temp$sim=1 # this was for the first run



sens.all.sites=df.temp
sens.all.sites$sim=1
levels(factor(sens.all.sites$species))

sens.all.sites <- sens.all.sites %>%
  group_by(species) %>%
  mutate(sites = row_number())


## 2) RAINFALL ############################
# empty df to store sens
df.rain=NULL
df.mu.rain=NULL # sens averaged across space

df.l_ratio.rain=NULL # empty df to store log ratios
df.mu.l_ratio.rain=NULL # empty df to store log ratios averaging across species

rm(sens.cov,sens.no.cov,
   l_ratio.cov,l_ratio.no.cov) # remove items above

for (i in c(3,5:12,14:15)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<1.005&l.mu.minPcov>0.75)
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<1.005&l.mu.maxPcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<1.005&l.mu.minPnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<1.005&l.mu.maxPnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain=rbind(df.rain,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
  
  # log ratios
  l_ratio.cov=abs(log(l.mu.maxPcov[keep]/l.mu.minPcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxPnoCov[keep]/l.mu.minPnoCov[keep]))

 df.l_ratio.rain=rbind(df.l_ratio.rain,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                  species=paste0("sp",i)))
  
   ## averaging across space -------------
  df.mu.rain=rbind(df.mu.rain,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
 
 df.mu.l_ratio.rain=rbind(df.mu.l_ratio.rain,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                        l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                        species=paste0("sp", i)))
}

#### add variables -----------
df.rain$type="Precipitation"
df.mu.rain$type="Precipitation"
df.mu.rain$sim=1

df.rain$sim=1

df.l_ratio.rain$type="Precipitation"
df.mu.l_ratio.rain$type="Precipitation"
df.mu.l_ratio.rain$sim=1
df.l_ratio.rain$sim=1


df.rain <- df.rain %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.rain)

# SIMULATION 2:  ###########################################
## 1) TEMPERATURE #####################
# define what we want to keep and remove the rest:
df_keep=c("df.mu.temp","df.mu.rain","sens.all.sites",
          "df.mu.l_ratio.temp","df.mu.l_ratio.rain")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

df.temp2=NULL # empty df to store sensitivities
df.mu.temp2=NULL # empty df to store sensitivities averaging across space

df.l_ratio.temp2=NULL
df.mu.l_ratio.temp2=NULL

for (i in c(3,5:12,14:15)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run2/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<1.005&l.mu.minTcov>0.75)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<1.005&l.mu.maxTcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<1.005&l.mu.minTnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<1.005&l.mu.maxTnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp2=rbind(df.temp2,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
  
  # log ratios
  l_ratio.cov=abs(log(l.mu.maxTcov[keep]/l.mu.minTcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxTnoCov[keep]/l.mu.minTnoCov[keep]))
  
  ## averaging across space --------
  df.mu.temp2=rbind(df.mu.temp2,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
  
  df.mu.l_ratio.temp2=rbind(df.mu.l_ratio.temp2,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                           l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                           species=paste0("sp", i)))
  
}
## add variables -------
df.temp2$type="Temperature"
df.temp2$sim=2
df.mu.temp2$type="Temperature"
df.mu.temp2$sim=2

# log ratios
df.l_ratio.temp2$type="Temperature"
df.l_ratio.temp2$sim=2
df.mu.l_ratio.temp2$type="Temperature"
df.mu.l_ratio.temp2$sim=2


df.temp2 <- df.temp2 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.temp2)




## 2) RAINFALL ############################
# empty df to store sens
df.rain2=NULL
df.mu.rain2=NULL # sens averaged across space

df.l_ratio.rain2=NULL
df.mu.l_ratio.rain2=NULL

rm(sens.cov,sens.no.cov,
   l_ratio.cov,l_ratio.no.cov) # remove items above

for (i in c(3,5:12,14:15)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<1.005&l.mu.minPcov>0.75)
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<1.005&l.mu.maxPcov>0.75)
    
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<1.005&l.mu.minPnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run2/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<1.005&l.mu.maxPnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain2=rbind(df.rain2,data.frame(sens.cov,sens.no.cov,
                                   species=paste0("sp", i)))
  
  l_ratio.cov=abs(log(l.mu.maxPcov[keep]/l.mu.minPcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxPnoCov[keep]/l.mu.minPnoCov[keep]))
  
  df.l_ratio.rain2=rbind(df.l_ratio.rain2,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain2=rbind(df.mu.rain2,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                         sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                         species=paste0("sp", i)))
  
  df.mu.l_ratio.rain2=rbind(df.mu.l_ratio.rain2,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                           l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                           species=paste0("sp", i)))
}

#### add variables -----------
df.rain2$type="Precipitation"
df.rain2$sim=2
df.mu.rain2$type="Precipitation"
df.mu.rain2$sim=2

df.l_ratio.rain2$type="Precipitation"
df.l_ratio.rain2$sim=2
df.mu.l_ratio.rain2$type="Precipitation"
df.mu.l_ratio.rain2$sim=2

df.rain2 <- df.rain2 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.rain2)




# SIMULATION 3: ###########################################
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2","sens.all.sites",
          "df.mu.l_ratio.temp","df.mu.l_ratio.rain",
          "df.mu.l_ratio.temp2","df.mu.l_ratio.rain2")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp3=NULL # empty df to store sensitivities
df.mu.temp3=NULL # empty df to store sensitivities averaging across space

df.l_ratio.temp3=NULL
df.mu.l_ratio.temp3=NULL

for (i in c(3,5:12,14:15)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run3/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<1.005&l.mu.minTcov>0.75)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<1.005&l.mu.maxTcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<1.005&l.mu.minTnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<1.005&l.mu.maxTnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp3=rbind(df.temp3,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  # log ratios
  l_ratio.cov=abs(log(l.mu.maxTcov[keep]/l.mu.minTcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxTnoCov[keep]/l.mu.minTnoCov[keep]))
  
  df.l_ratio.temp3=rbind(df.l_ratio.temp3,data.frame(l_ratio.cov,l_ratio.no.cov,
                                     species=paste0("sp", i)))
  
  

  
  ## averaging across space --------
  df.mu.temp3=rbind(df.mu.temp3,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
  
  df.mu.l_ratio.temp3=rbind(df.mu.l_ratio.temp3,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                          l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                          species=paste0("sp", i)))
  
}

## add variables -------
df.temp3$type="Temperature"
df.temp3$sim=3
df.mu.temp3$type="Temperature"
df.mu.temp3$sim=3

df.l_ratio.temp3$type="Temperature"
df.l_ratio.temp3$sim=3
df.mu.l_ratio.temp3$type="Temperature"
df.mu.l_ratio.temp3$sim=3


df.temp3 <- df.temp3 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.temp3)


## 2) RAINFALL ############################
# empty df to store sens
df.rain3=NULL
df.mu.rain3=NULL # sens averaged across space

df.l_ratio.rain3=NULL
df.mu.l_ratio.rain3=NULL

rm(sens.cov,sens.no.cov,
   l_ratio.cov,l_ratio.no.cov) # remove items above

for (i in c(3,5:12,14:15)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<1.005&l.mu.minPcov>0.75)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<1.005&l.mu.maxPcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<1.005&l.mu.minPnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run3/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<1.005&l.mu.maxPnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain3=rbind(df.rain3,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  l_ratio.cov=abs(log(l.mu.maxPcov[keep]/l.mu.minPcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxPnoCov[keep]/l.mu.minPnoCov[keep]))
  
  df.l_ratio.rain3=rbind(df.l_ratio.rain3,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain3=rbind(df.mu.rain3,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))

  df.mu.l_ratio.rain3=rbind(df.mu.l_ratio.rain3,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                           l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                           species=paste0("sp", i)))
  }

#### add variables -----------
df.rain3$type="Precipitation"
df.rain3$sim=3
df.mu.rain3$type="Precipitation"
df.mu.rain3$sim=3

df.l_ratio.rain3$type="Precipitation"
df.l_ratio.rain3$sim=3
df.mu.l_ratio.rain3$type="Precipitation"
df.mu.l_ratio.rain3$sim=3


df.rain3 <- df.rain3 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.rain3)


# SIMULATION 4:#####################################
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2","sens.all.sites",
          "df.mu.temp3","df.mu.rain3",
          "df.mu.l_ratio.temp","df.mu.l_ratio.rain",
          "df.mu.l_ratio.temp2","df.mu.l_ratio.rain2",
          "df.mu.l_ratio.temp3","df.mu.l_ratio.rain3")
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp4=NULL # empty df to store sensitivities
df.mu.temp4=NULL # empty df to store sensitivities averaging across space

df.l_ratio.temp4=NULL
df.mu.l_ratio.temp4=NULL

for (i in c(3,5:12,14:15)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run4/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<1.005&l.mu.minTcov>0.75)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<1.005&l.mu.maxTcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<1.005&l.mu.minTnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<1.005&l.mu.maxTnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp4=rbind(df.temp4,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  # log ratios
  l_ratio.cov=abs(log(l.mu.maxTcov[keep]/l.mu.minTcov[keep]))
  
  l_ratio.no.cov=abs(log(l.mu.maxTnoCov[keep]/l.mu.minTnoCov[keep]))
  
  df.l_ratio.temp4=rbind(df.l_ratio.temp4,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                     species=paste0("sp", i)))
  
  ## averaging across space --------
  df.mu.temp4=rbind(df.mu.temp4,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
  
  df.mu.l_ratio.temp4=rbind(df.mu.l_ratio.temp4,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                           l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                           species=paste0("sp", i)))
}



## add variables -------
df.temp4$type="Temperature"
df.temp4$sim=4
df.mu.temp4$type="Temperature"
df.mu.temp4$sim=4

df.l_ratio.temp4$type="Temperature"
df.l_ratio.temp4$sim=4
df.mu.l_ratio.temp4$type="Temperature"
df.mu.l_ratio.temp4$sim=4



df.temp4 <- df.temp4 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.temp4)


## 2) RAINFALL ############################
# empty df to store sens
df.rain4=NULL
df.mu.rain4=NULL # sens averaged across space

df.l_ratio.rain4=NULL
df.mu.l_ratio.rain4=NULL

rm(sens.cov,sens.no.cov,l_ratio.cov,l_ratio.no.cov) # remove items above

for (i in c(3,5:12,14:15)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<1.005&l.mu.minPcov>0.75)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<1.005&l.mu.maxPcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<1.005&l.mu.minPnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run4/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<1.005&l.mu.maxPnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain4=rbind(df.rain4,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  l_ratio.cov=abs(log(l.mu.maxPcov[keep]/l.mu.minPcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxPnoCov[keep]/l.mu.minPnoCov[keep]))
  
  df.l_ratio.rain4=rbind(df.l_ratio.rain4,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                     species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain4=rbind(df.mu.rain4,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
  
  df.mu.l_ratio.rain4=rbind(df.mu.l_ratio.rain4,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                           l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                           species=paste0("sp", i)))
}

#### add variables -----------
df.rain4$type="Precipitation"
df.rain4$sim=4
df.mu.rain4$type="Precipitation"
df.mu.rain4$sim=4

df.l_ratio.rain4$type="Precipitation"
df.l_ratio.rain4$sim=4
df.mu.l_ratio.rain4$type="Precipitation"
df.mu.l_ratio.rain4$sim=4


df.rain4 <- df.rain4 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.rain4)


# SIMULATION 5:#####################################
# doesn't work yet because one folder is empty
df_keep=c("df.mu.temp","df.mu.rain",
          "df.mu.temp2","df.mu.rain2","sens.all.sites",
          "df.mu.temp3","df.mu.rain3",
          "df.mu.temp4","df.mu.rain4",
          "df.mu.l_ratio.temp","df.mu.l_ratio.rain",
          "df.mu.l_ratio.temp2","df.mu.l_ratio.rain2",
          "df.mu.l_ratio.temp3","df.mu.l_ratio.rain3",
          "df.mu.l_ratio.temp4","df.mu.l_ratio.rain4"
          )
all_objects=ls()
objects_to_remove=setdiff(all_objects, df_keep)
rm(list=objects_to_remove)

## 1) TEMPERATURE #####################
df.temp5=NULL # empty df to store sensitivities
df.mu.temp5=NULL # empty df to store sensitivities averaging across space

df.l_ratio.temp5=NULL
df.mu.l_ratio.temp5=NULL

for (i in c(3,5:12,14:15)) {
  
  ### min T cov ----------
  df.t=read.csv(paste0("MinTempCov/run5/results.sp", i, ".csv"))
  
  # calculate lambdas
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTcov=apply(l,1,mean)
  
 keep=which(l.mu.minTcov<1.005&l.mu.minTcov>0.75)
  
  
  minT=df.t$temperature
  
  ## max T cov ------------
  df.t=read.csv(paste0("MaxTempCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxTcov<1.005&l.mu.maxTcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxT=df.t$temperature
  
  
  ## min T ------------
  df.t=read.csv(paste0("MinTempNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minTnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minTnoCov<1.005&l.mu.minTnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max T ------------
  df.t=read.csv(paste0("MaxTempNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxTnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxTnoCov<1.005&l.mu.maxTnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens -------
  sens.cov=abs((l.mu.maxTcov[keep]-l.mu.minTcov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxTnoCov[keep]-l.mu.minTnoCov[keep])/((maxT[keep]-minT[keep])/0.5))
  
  df.temp5=rbind(df.temp5,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
  l_ratio.cov=abs(log(l.mu.maxTcov[keep]/l.mu.minTcov[keep]))
  l_ratio.no.cov=abs(log(l.mu.maxTnoCov[keep]/l.mu.minTnoCov[keep]))
  
  df.l_ratio.temp5=rbind(df.l_ratio.temp5,data.frame(l_ratio.cov,l_ratio.no.cov,species=paste0("sp",i)))
  
  ## averaging across space --------
  df.mu.temp5=rbind(df.mu.temp5,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
  
  df.mu.l_ratio.temp5=rbind(df.mu.l_ratio.temp5,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                          l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                          species=paste0("sp", i)))
}
## add variables -------
df.temp5$type="Temperature"
df.temp5$sim=5
df.mu.temp5$type="Temperature"
df.mu.temp5$sim=5

df.l_ratio.temp5$type="Temperature"
df.l_ratio.temp5$sim=5
df.mu.l_ratio.temp5$type="Temperature"
df.mu.l_ratio.temp5$sim=5



df.temp5 <- df.temp5 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.temp5)

## 2) RAINFALL ############################
# empty df to store sens
df.rain5=NULL
df.mu.rain5=NULL # sens averaged across space

df.l_ratio.rain5=NULL
df.mu.l_ratio.rain5=NULL

rm(sens.cov,sens.no.cov,l_ratio.cov,l_ratio.no.cov) # remove items above

for (i in c(3,5:12,14:15)) {
  
  ### min P cov -----------------
  df.t=read.csv(paste0("MinPrecipCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPcov=apply(l,1,mean)
  
  keep=which(l.mu.minPcov<1.005&l.mu.minPcov>0.75)
  
  
  minP=df.t$precipitation
  
  ## max P cov -------------------
  df.t=read.csv(paste0("MaxPrecipCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPcov=apply(l,1,mean)
  
  keep2=which(l.mu.maxPcov<1.005&l.mu.maxPcov>0.75)
  
  keep=keep[keep%in%keep2]
  
  maxP=df.t$precipitation
  
  
  ## min P ------------------
  df.t=read.csv(paste0("MinPrecipNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.minPnoCov=apply(l,1,mean)
  
  keep3=which(l.mu.minPnoCov<1.005&l.mu.minPnoCov>0.75)
  
  keep=keep[keep%in%keep3]
  
  
  ## max P  ----------------
  df.t=read.csv(paste0("MaxPrecipNoCov/run5/results.sp", i, ".csv"))
  
  l1=df.t$adults_2070/df.t$adults_2060
  l2=df.t$adults_2080/df.t$adults_2070
  l3=df.t$adults_2090/df.t$adults_2080
  
  l=cbind(l1,l2,l3)
  l.mu.maxPnoCov=apply(l,1,mean)
  
  keep4=which(l.mu.maxPnoCov<1.005&l.mu.maxPnoCov>0.75)
  
  keep=keep[keep%in%keep4]
  
  ## calculate sens ------------
  sens.cov=abs((l.mu.maxPcov[keep]-l.mu.minPcov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  sens.no.cov=abs((l.mu.maxPnoCov[keep]-l.mu.minPnoCov[keep])/((maxP[keep]-minP[keep])/0.5))
  
  df.rain5=rbind(df.rain5,data.frame(sens.cov,sens.no.cov,
                                     species=paste0("sp", i)))
  
 l_ratio.cov=abs(log(l.mu.maxPcov[keep]/l.mu.minPcov[keep]))
 l_ratio.no.cov=abs(log(l.mu.maxPnoCov[keep]/l.mu.minPnoCov[keep]))
 
 df.l_ratio.rain5=rbind(df.l_ratio.rain5,data.frame(l_ratio.cov,l_ratio.no.cov,
                                                    species=paste0("sp", i)))
  
  ## averaging across space -------------
  df.mu.rain5=rbind(df.mu.rain5,data.frame(sens.cov=mean(sens.cov[!is.nan(sens.cov)&is.finite(sens.cov)]),
                                           sens.no.cov=mean(sens.no.cov[!is.nan(sens.no.cov)&is.finite(sens.no.cov)]),
                                           species=paste0("sp", i)))
 
 df.mu.l_ratio.rain5=rbind(df.mu.l_ratio.rain5,data.frame(l_ratio.cov=mean(l_ratio.cov[!is.nan(l_ratio.cov)&is.finite(l_ratio.cov)]),
                                                          l_ratio.no.cov=mean(l_ratio.no.cov[!is.nan(l_ratio.no.cov)&is.finite(l_ratio.no.cov)]),
                                                          species=paste0("sp", i)))
}

#### add variables -----------
df.rain5$type="Precipitation"
df.rain5$sim=5
df.mu.rain5$type="Precipitation"
df.mu.rain5$sim=5

df.l_ratio.rain5$type="Precipitation"
df.l_ratio.rain5$sim=5
df.mu.l_ratio.rain5$type="Precipitation"
df.mu.l_ratio.rain5$sim=5

df.rain5 <- df.rain5 %>%
  group_by(species) %>%
  mutate(sites = row_number())

sens.all.sites=rbind(sens.all.sites,df.rain5)

# get mean max and min covariates 
max.P=mean(maxP)
min.P=mean(minP)

max.T=mean(maxT)
min.T=mean(minT)

# Rearrange dataframes ###########################
df <- rbind(df.mu.temp,df.mu.temp2,df.mu.temp3,df.mu.temp4,df.mu.temp5,
                 df.mu.rain,df.mu.rain2,df.mu.rain3,df.mu.rain4,df.mu.rain5)

df.l_ratio <- rbind(df.mu.l_ratio.temp,
                    df.mu.l_ratio.temp2,
                    df.mu.l_ratio.temp3,
                    df.mu.l_ratio.temp4,
                    df.mu.l_ratio.temp5,
                    df.mu.l_ratio.rain,
                    df.mu.l_ratio.rain2,
                    df.mu.l_ratio.rain3,
                    df.mu.l_ratio.rain4,
                    df.mu.l_ratio.rain5)

df.l_ratio_cov <- df.l_ratio %>% 
  select(-l_ratio.no.cov)

colnames(df.l_ratio_cov)[colnames(df.l_ratio_cov)=="l_ratio.cov"] <- "l_ratio"
df.l_ratio_cov$cov <- 1

df.l_ratio_no_cov <- df.l_ratio %>%
  select(-l_ratio.cov)

colnames(df.l_ratio_no_cov)[colnames(df.l_ratio_no_cov)=="l_ratio.no.cov"] <- "l_ratio"
df.l_ratio_no_cov$cov <- 0

df2_l_ratio <- rbind(df.l_ratio_no_cov,df.l_ratio_cov)

# lkasfmaed
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
df2$species

df2_l_ratio <- select(df2_l_ratio, species, type, l_ratio, cov, sim)

df2
df2$full.species <- c("Fagus sylvatica","Pinus halepensis","Pinus nigra","Pinus pinaster","Pinus pinea","Pinus sylvestris","Pinus uncinata",
                            "Quercus faginea","Quercus ilex",
                            "Quercus robur/petraea","Quercus suber")

df2$full.species

df2$study.doi <- "10.1093/jpe/rtw081"
df2$group <- "Plants"
df2$year.of.publication <- "2016"
df2$driver.type <- "C"
df2$continent <- "Europe"
df2$stage.age <- "all"
df2$vital.rates <- "all"
df2$l_ratio <- df2_l_ratio$l_ratio

colnames(df2)[colnames(df2)=="type"] <- "driver"

df2 <- df2 %>% select(study.doi,year.of.publication,group,full.species,continent,driver,driver.type,stage.age,vital.rates,sens,cov,l_ratio)

colnames(df2)[colnames(df2)=="full.species"] <- "species"

levels(factor(df2$species))

# add species-specific age at maturity
df.mat=data.frame(species=c("Fagus sylvatica",
                            "Pinus halepensis",
                            "Pinus nigra",
                            "Pinus pinaster",
                            "Pinus pinea",
                            "Pinus sylvestris",
                            "Pinus uncinata",
                            "Quercus faginea",
                            "Quercus ilex",
                            "Quercus robur/petraea",
                            "Quercus suber"),
                      mat=c(40,
                            5,
                            27.5,
                            12.5,
                            22.5,
                            27.5,
                            15,
                            30,
                            30,
                            30,
                            30)) # see supporting information for sources
df2=merge(df2, df.mat, by="species",all.x=T)

df2$n.vr=4  # number of vital rates with covariates
df2$n.pam=21 # number of parameters
df2$dens=0 # are there any density effects?
df2$biotic_interactions=0 # or other biotic interactions? 
df2$lambda.sim=1
df2$study.length=22

df2 <- df2 %>% relocate(l_ratio, .after = last_col())

# SAVE OUTPUT ####################
#write.csv(df2, "SensSpainTrees.csv",row.names = F)

# rearrange new df #############

# old df
#old.sens <- read.csv("SensSpainTrees.csv")

# new
sens.no.cov <- select(sens.all.sites, sens.no.cov, species, sim, type, sites)

names(sens.no.cov)[names(sens.no.cov) == "sens.no.cov"] <- "sens"  
names(sens.no.cov)[names(sens.no.cov) == "type"] <- "driver" 
sens.no.cov$cov=0

sens.cov <- select(sens.all.sites, sens.cov, species, sim, type, sites)
names(sens.cov)[names(sens.cov) == "sens.cov"] <- "sens"  
names(sens.cov)[names(sens.cov) == "type"] <- "driver" 
sens.cov$cov=1

new.sens <- rbind(sens.no.cov,sens.cov)
levels(factor(new.sens$species))

species_lookup <- data.frame(
  species_code = c("sp3", "sp4", "sp5", "sp6", "sp7", "sp8", "sp9", 
                   "sp10", "sp11", "sp12", "sp14", "sp15"),
  latin_name = c("Fagus sylvatica", "Juniperus thurifera", "Pinus halpensis", 
                 "Pinus nigra", "Pinus pinaster", "Pinus pinea", "Pinus sylvestris", 
                 "Pinus uncinata", "Quercus faginea", "Quercus ilex", 
                 "Quercus robur/petraea", "Quercus suber")
)

new.sens <- new.sens %>%
  left_join(species_lookup, by = c("species" = "species_code")) %>%
  mutate(species = latin_name) %>%
  select(-latin_name)

# add species-specific age at maturity
df.mat=data.frame(species=c("Fagus sylvatica",
                            "Pinus halepensis",
                            "Pinus nigra",
                            "Pinus pinaster",
                            "Pinus pinea",
                            "Pinus sylvestris",
                            "Pinus uncinata",
                            "Quercus faginea",
                            "Quercus ilex",
                            "Quercus robur/petraea",
                            "Quercus suber"),
                  mat=c(40,
                        5,
                        27.5,
                        12.5,
                        22.5,
                        27.5,
                        15,
                        30,
                        30,
                        30,
                        30)) # see supporting information for sources
new.sens=merge(new.sens, df.mat, by="species",all.x=T)

new.sens$study.doi <- "10.1093/jpe/rtw081"
new.sens$group <- "Plants"
new.sens$year.of.publication <- "2016"
new.sens$driver.type <- "C"
new.sens$continent <- "Europe"
new.sens$stage.age <- "all"
new.sens$vital.rates <- "all"

new.sens$n.vr=4  # number of vital rates with covariates
new.sens$n.pam=21 # number of parameters
new.sens$dens=0 # are there any density effects?
new.sens$biotic_interactions=0 # or other biotic interactions? 
new.sens$lambda.sim=1
new.sens$study.length=22

colnames(new.sens)[colnames(new.sens)=="sites"] <- "population"

# rearrange cols
Sens_Trees_sites <- select(new.sens, species, study.doi, year.of.publication, group, continent, driver, driver.type, stage.age, vital.rates, sens, cov, mat, n.vr, n.pam, dens, biotic_interactions, lambda.sim, study.length, population)

#write.csv(Sens_Trees_sites, "Sens_Trees_sites.csv",row.names = F)
