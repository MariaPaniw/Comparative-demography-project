# This script performs analyses as part of Paniw et al. XXXX
# Please read the README file on GitHub for a workflow of the analyses

# Created by: Maria Paniw
# Date created: June 20, 2022

rm(list = ls())

library(plyr)
library(ggplot2)
library(ggnewscale)
library(boot)
library(patchwork)
library(popbio)

setwd("~/Desktop/Comparative-demography-project-main/Files_Shrubs")

### Load MCMC parameters

load("allShrubs.Rdata")

### Parameters 

# !!!! NOTE: although we have rain for all years, we using only values for years where we also have density observations

rain=c(0.1392206, -0.1514142,  NA,  1.5511706, NA, NA, -1.1548692,  NA, NA,  NA,
       0.6807192, NA, -0.5835423, -0.6776162)


#### INITIAL DENSITIES 

### INTERSPECIFIC
load("Cistus libanotis.Rdata")

neigh1=C

load("Halimium commutatum.Rdata")

neigh2=C

load("Lavandula stoechas.Rdata")

neigh3=C


load("Rosmarinus officinalis.Rdata")

neigh4=C


load("Halimium halimifolium.Rdata")

neigh5=C

neighHal1=log(neigh1+neigh2+neigh3+neigh4)

neighCis=log(neigh2+neigh3+neigh4+neigh5)

#### INTRASPECIFIC

num=read.csv("shrub_number.csv")

ad.h1=log(num[num$species%in%"Halimium halimifolium",c("X.adults")]+0.001)

ad.c=log(num[num$species%in%"Cistus libanotis",c("X.adults")]+0.001)

##### Environmental min, mean, max

min.rain=(-1.1548692)
max.rain=1.5511706

mean.rain=mean(rain,na.rm=T)
sd.rain=sd(rain,na.rm=T)

min.inter.h1=min(neighHal1)
max.inter.h1=max(neighHal1)

mean.inter.h1=mean(neighHal1)
sd.inter.h1=sd(neighHal1)

min.inter.c=min(neighCis)
max.inter.c=max(neighCis)

mean.inter.c=mean(neighCis)
sd.inter.c=sd(neighCis)


ad.h1=matrix(ad.h1,18,8,byrow = F)

min.intra.h1=min(ad.h1)
max.intra.h1=max(ad.h1)

mean.intra.h1=mean(ad.h1,na.rm=T)
sd.intra.h1=sd(ad.h1,na.rm=T)


ad.c=matrix(ad.c,18,8,byrow = F)

min.intra.c=min(ad.c,na.rm=T)
max.intra.c=max(ad.c,na.rm=T)

mean.intra.c=mean(ad.c,1,mean,na.rm=T)
sd.intra.c=sd(ad.c,na.rm=T)

# Halilimium halimifolium ###########

survM.h1 <- function(rain,neigh,pu){
  
  mu=inv.logit(out1$sims.list$a0.h1[pu] +out1$sims.list$a1.h1[pu]*rain
               + out1$sims.list$a3.h1[pu]*neigh)
  
  return(mu)
}


rec.h1=1.471


survSap.h1 <- function(rain,ad,pu){
  
  mu=inv.logit(out1$sims.list$b0.h1[pu] + out1$sims.list$b2.h1[pu]*ad )
  
  return(mu)
}

survS.h1 = 0.517/5


# Cistus libanotis ##########################

survM.c <- function(rain,ad,neigh,pu){
  
  mu=inv.logit(out1$sims.list$a0.c[pu] +out1$sims.list$a1.c[pu]*rain
               + out1$sims.list$a2.c[pu]*ad + out1$sims.list$a3.c[pu]*neigh)
  
  return(mu)
}



rec.c=1.18


survSap.c <- out1$sims.list$phiS.c

survS.c = 1.677/5

# Sensitivity Analysis ################################################
n.stage=3

#create a big loop
# I sample 100 from the 2100 MCMC posterior samples
par.samp=sample(1:2100,100,replace=F)


min.ad.h1=mean(log(num[num$species%in%"Halimium halimifolium"&num$year%in%2013,c("X.adults")]+0.001))
max.ad.h1=mean(log(num[num$species%in%"Halimium halimifolium"&num$year%in%2010,c("X.adults")]+0.001))

min.bi.h1=mean(neighHal1[,7])
max.bi.h1=mean(neighHal1[,4])

min.ad.c=mean(log(num[num$species%in%"Cistus libanotis"&num$year%in%2013,c("X.adults")]+0.001))
max.ad.c=mean(log(num[num$species%in%"Cistus libanotis"&num$year%in%2010,c("X.adults")]+0.001))


min.bi.c=mean(neighCis[,7])
max.bi.c=mean(neighCis[,4])

# rainfall
deltaR.h1=NULL
deltaR.c=NULL
deltaR.h1.V2=NULL
deltaR.c.V2=NULL

# density
deltaInter.h1=NULL
deltaInter.c=NULL
deltaInter.h1.V2=NULL
deltaInter.c.V2=NULL
deltaIntra.h1=NULL
deltaIntra.c=NULL
deltaIntra.h1.V2=NULL
deltaIntra.c.V2=NULL

# log ratios
# rainfall
l_ratio.R.h1=NULL
l_ratio.R.c=NULL
l_ratio.R.h1.V2=NULL
l_ratio.R.c.V2=NULL

# density
l_ratio.Inter.h1=NULL
l_ratio.Inter.c=NULL
l_ratio.Inter.h1.V2=NULL
l_ratio.Inter.c.V2=NULL
l_ratio.Intra.h1=NULL
l_ratio.Intra.c=NULL
l_ratio.Intra.h1.V2=NULL
l_ratio.Intra.c.V2=NULL


for(i in 1:length(par.samp)){
  
  ## RAINFALL ####################
  
  ### 1. Sensitivity to rain assuming mean dens
  
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=max.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=max.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=max.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=min.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=min.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=min.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaR.h1[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  l_ratio.R.h1[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  rm(mpm.min)
  rm(mpm.max)
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=max.rain,ad=mean.intra.c,neigh=mean.inter.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=min.rain,ad=mean.intra.c,neigh=mean.inter.c,par.samp[i])),n.stage,n.stage)
  
  deltaR.c[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  l_ratio.R.c[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  ### 2. Sensitivity to rain assuming covarying dens

  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=max.rain,ad=max.ad.h1,par.samp[i]),0.21*survSap.h1(rain=max.rain,ad=max.ad.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=max.rain,neigh=max.bi.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min =  matrix(c(0,survS.h1,0,
                      0,(1-0.21)*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),0.21*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),
                      rec.h1,0,survM.h1(rain=min.rain,neigh=min.bi.h1,par.samp[i])),n.stage,n.stage)
  
  deltaR.h1.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  l_ratio.R.h1.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  
  ### get deltas for each of 18 sites:
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=max.rain,ad=max.ad.c,neigh=max.bi.c,par.samp[i])),n.stage,n.stage)
  
  mpm.min =  matrix(c(0,survS.c,0,
                      0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                      rec.c,0,survM.c(rain=min.rain,ad=min.ad.c,neigh=min.bi.c,par.samp[i])),n.stage,n.stage)
  
  deltaR.c.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  l_ratio.R.c.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  
# INTERSPECIFIC DENSITY ##############

  ### 1. Sensitivity to inter assuming mean rain/intra
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=max.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=min.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaInter.h1[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.h1-min.inter.h1)/sd.inter.h1)
  l_ratio.Inter.h1[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean.rain,ad=mean.intra.c,neigh=max.inter.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean.rain,ad=mean.intra.c,neigh=min.inter.c,par.samp[i])),n.stage,n.stage)
  
  deltaInter.c[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.c-min.inter.c)/sd.inter.c)
  l_ratio.Inter.c[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  ### 2. Sensitivity to inter assuming covarying rain/intra
  
  # Halimium
  year.max.n=c(4,4,1,4,4,4,4,4,1,7,4,4,4,7,4,4,4,4)
  year.min.n=c(1,14,11,7,7,11,2,11,2,11,2,2,2,14,2,2,14,14)
  
  year.max=c(2010,2010,2007,2010,2010,2010,2010,2010,2007,2013,2010,2010,2010,2013,2010,2010,2010,2010)
  year.min=c(2007,2020,2017,2013,2013,2017,2008,2017,2008,2017,2008,2008,2008,2020,2008,2008,2020,2020)
  
  
  max.ad.h1=mean(c(ad.h1[1,3],ad.h1[2,3],ad.h1[3,1],ad.h1[4,3],ad.h1[5,3],ad.h1[6,3],ad.h1[7,3],ad.h1[8,3],ad.h1[9,1],
                   ad.h1[10,4],ad.h1[11,3],ad.h1[12,3],ad.h1[13,3],ad.h1[14,4],ad.h1[15,3],ad.h1[16,3],ad.h1[17,3],ad.h1[18,3]))
  
  min.ad.h1=mean(c(ad.h1[1,1],ad.h1[2,7],ad.h1[3,5],ad.h1[4,4],ad.h1[5,4],ad.h1[6,5],ad.h1[7,2],ad.h1[8,5],ad.h1[9,2],
                   ad.h1[10,5],ad.h1[11,2],ad.h1[12,2],ad.h1[13,2],ad.h1[14,7],ad.h1[15,2],ad.h1[16,2],ad.h1[17,7],ad.h1[18,7]))
  
  
  ### get deltas for each of 18 sites:
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.max.n]),ad=max.ad.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.max.n]),ad=max.ad.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.max.n]),neigh=max.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.min.n]),ad=min.ad.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.min.n]),ad=min.ad.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.min.n]),neigh=min.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaInter.h1.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.h1-min.inter.h1)/sd.inter.h1)
  l_ratio.Inter.h1.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  # Cistus
  year.max.n=c(4,1,1,1,1,4,4,2,1,1,4,4,1,7,4,4,4,7)
  year.min.n=c(2,14,11,7,7,11,2,11,14,14,2,2,2,14,13,2,14,14)
  
  year.max=c(2010,2007,2007,2007,2007,2010,2010,2008,2007,2007,2010,2010,2007,2013,2010,2010,2010,2013)
  year.min=c(2008,2020,2017,2013,2013,2017,2008,2017,2020,2020,2008,2008,2008,2020,2019,2008,2020,2020)
  
  
  max.ad.c=mean(c(ad.c[1,3],ad.c[2,1],ad.c[3,1],ad.c[4,1],ad.c[5,1],ad.c[6,3],ad.c[7,3],ad.c[8,2],ad.c[9,1],
                  ad.c[10,1],ad.c[11,3],ad.c[12,3],ad.c[13,1],ad.c[14,4],ad.c[15,3],ad.c[16,3],ad.c[17,3],ad.c[18,4]))
  
  min.ad.c=mean(c(ad.c[1,2],ad.c[2,7],ad.c[3,5],ad.c[4,4],ad.c[5,4],ad.c[6,5],ad.c[7,2],ad.c[8,5],ad.c[9,7],
                  ad.c[10,7],ad.c[11,2],ad.c[12,2],ad.c[13,2],ad.c[14,7],ad.c[15,6],ad.c[16,2],ad.c[17,7],ad.c[18,7]))
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.max.n]),ad=max.ad.c,neigh=max.inter.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.min.n]),ad=min.ad.c,neigh=min.inter.c,par.samp[i])),n.stage,n.stage)
  
  deltaInter.c.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.c-min.inter.c)/sd.inter.c)
  l_ratio.Inter.c.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  # INTRASPECIFIC DENSITY #############################
  
  ### 1. Sensitivity to intra assuming mean rain/inter
  
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=max.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=max.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=min.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=min.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.h1[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.h1-min.intra.h1)/sd.intra.h1)
  l_ratio.Intra.h1[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean.rain,ad=max.intra.c,neigh=mean.inter.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean.rain,ad=min.intra.c,neigh=mean.inter.c,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.c[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.c-min.intra.c)/sd.intra.c)
  l_ratio.Intra.c[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
 
  ### 2. Sensitivity to intra assuming covarying rain/inter
  
  # Halimium
  year.max.n=c(4,7,1,1,1,4,1,2,2,1,4,4,1,1,4,1,1,2)
  year.min.n=c(2,14,11,14,14,14,14,14,14,14,14,11,14,14,14,14,14,14)
  
  
  max.bi.h1=mean(c(neighHal1[1,4],neighHal1[2,7],neighHal1[3,1],neighHal1[4,1],neighHal1[5,1],neighHal1[6,4],neighHal1[7,1],neighHal1[8,2],neighHal1[9,2],
                   neighHal1[10,1],neighHal1[11,4],neighHal1[12,4],neighHal1[13,1],neighHal1[14,1],neighHal1[15,4],neighHal1[16,1],neighHal1[17,1],neighHal1[18,2]))
  
  min.bi.h1=mean(c(neighHal1[1,2],neighHal1[2,14],neighHal1[3,11],neighHal1[4,14],neighHal1[5,14],neighHal1[6,14],neighHal1[7,14],neighHal1[8,14],neighHal1[9,14],
                   neighHal1[10,14],neighHal1[11,14],neighHal1[12,11],neighHal1[13,14],neighHal1[14,14],neighHal1[15,14],neighHal1[16,14],neighHal1[17,14],neighHal1[18,14]))
  
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.max.n]),ad=max.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.max.n]),ad=max.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.max.n]),neigh=max.bi.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.min.n]),ad=min.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.min.n]),ad=min.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.min.n]),neigh=min.bi.h1,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.h1.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.h1-min.intra.h1)/sd.intra.h1)
  l_ratio.Intra.h1.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
  
  # Cistus
  year.max.n=c(4,4,2,1,4,2,4,4,1,13,4,4,7,7,4,4,4,4)
  year.min.n=c(11,11,11,14,14,11,2,14,2,14,2,1,1,14,14,1,14,2)
  
  
  max.bi.c=mean(c(neighCis[1,4],neighCis[2,4],neighCis[3,2],neighCis[4,1],neighCis[5,4],neighCis[6,2],neighCis[7,4],neighCis[8,4],neighCis[9,1],
                  neighCis[10,13],neighCis[11,4],neighCis[12,4],neighCis[13,7],neighCis[14,7],neighCis[15,4],neighCis[16,4],neighCis[17,4],neighCis[18,4]))
  
  min.bi.c=mean(c(neighCis[1,11],neighCis[2,11],neighCis[3,11],neighCis[4,14],neighCis[5,14],neighCis[6,11],neighCis[7,2],neighCis[8,14],neighCis[9,2],
                  neighCis[10,14],neighCis[11,2],neighCis[12,1],neighCis[13,1],neighCis[14,14],neighCis[15,14],neighCis[16,1],neighCis[17,14],neighCis[18,2]))
  
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.max.n]),ad=max.intra.c,neigh=max.bi.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.min.n]),ad=min.intra.c,neigh=min.bi.c,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.c.V2[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.c-min.intra.c)/sd.intra.c)
  l_ratio.Intra.c.V2[i]=abs(log(lambda(mpm.max)/lambda(mpm.min)))
}

  ### save output ########################
# sensitivities of all vital rates

# Halimium
df.h1=data.frame(study.doi="10.1098/rspb.2022.1494",
                 year.of.publication="2023",
                 group="Plants",
                 species="Halimium halimifolium",
                 continent="Europe",
                 driver=rep(c("rain","intraD","interD"),each=200),
                 driver.type=rep(c("C","D","B"),each=200), # intradensity = density and interdensity = biotic
                 stage.age="all",
                 vital.rates="all",
                 sens=c(deltaR.h1,deltaR.h1.V2,
                        deltaIntra.h1,deltaIntra.h1.V2,
                        deltaInter.h1,deltaInter.h1.V2),
                 cov=rep(c(0,1),each=100),
                 mat=5.8, # age at maturity
                 n.vr=2, # number of vital rates with covariates
                 n.pam=5, # number of parameters
                 dens=1, # are there intraspecific density effects?
                 biotic_interactions=1, # are there other biotic interactions?
                 lambda.sim=0,
                 study.length=8,
                 l_ratio=c(l_ratio.R.h1,l_ratio.R.h1.V2,
                           l_ratio.Intra.h1,l_ratio.Intra.h1.V2,
                           l_ratio.Inter.h1,l_ratio.Inter.h1.V2))

write.csv(df.h1, "Sens_Halimium.csv", row.names = F)


# Cistus
df.c=data.frame(study.doi="10.1098/rspb.2022.1494",
                 year.of.publication="2023",
                 group="Plants",
                 species="Cistus libanotis",
                 continent="Europe",
                 driver=rep(c("rain","intraD","interD"),each=200),
                 driver.type=rep(c("C","D","B"),each=200), # intradensity = density and interdensity = biotic
                 stage.age="all",
                 vital.rates="all",
                sens=c(deltaR.c,deltaR.c.V2,
                         deltaIntra.c,deltaIntra.c.V2,
                         deltaInter.c,deltaInter.c.V2),
                cov=rep(c(0,1),each=100),
                mat=4.6,
                n.vr=1, # number of vital rates with covariates
                n.pam=4, # number of parameters of these vital rates
                dens=1,
                biotic_interactions=1,
                lambda.sim=0,
                study.length=8,
                l_ratio=c(l_ratio.R.c,l_ratio.R.c.V2,
                          l_ratio.Intra.c,l_ratio.Intra.c.V2,
                          l_ratio.Inter.c,l_ratio.Inter.c.V2))

write.csv(df.c, "Sens_Cistus.csv", row.names = F)



# Sensitivities per vital rate ###########################
n.stage=3

#create a big loop
# I sample 100 from the 2100 MCMC posterior samples
par.samp=sample(1:2100,100,replace=F)

min.ad.h1=mean(log(num[num$species%in%"Halimium halimifolium"&num$year%in%2013,c("X.adults")]+0.001))
max.ad.h1=mean(log(num[num$species%in%"Halimium halimifolium"&num$year%in%2010,c("X.adults")]+0.001))

min.bi.h1=mean(neighHal1[,7])
max.bi.h1=mean(neighHal1[,4])

min.ad.c=mean(log(num[num$species%in%"Cistus libanotis"&num$year%in%2013,c("X.adults")]+0.001))
max.ad.c=mean(log(num[num$species%in%"Cistus libanotis"&num$year%in%2010,c("X.adults")]+0.001))

min.bi.c=mean(neighCis[,7])
max.bi.c=mean(neighCis[,4])

deltaR.h1.survSap=NULL
deltaR.h1.survM=NULL
deltaR.c.survM=NULL
deltaInter.h1.survM=NULL
deltaInter.c.survM=NULL
deltaIntra.h1.survSap=NULL
deltaIntra.c.survM=NULL

for(i in 1:length(par.samp)){
  
## RAINFALL ####################
  
  # 1. Sapling survival (survSap)
  # Halimium
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=max.rain,ad=max.ad.h1,par.samp[i]),0.21*survSap.h1(rain=max.rain,ad=max.ad.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min =  matrix(c(0,survS.h1,0,
                      0,(1-0.21)*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),0.21*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),
                      rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaR.h1.survSap[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  
  
  # Cistus
  # doesn't have a function for survSap
  
  # 2. Adult survival (survM)
  # Halimium
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=max.rain,neigh=max.bi.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min =  matrix(c(0,survS.h1,0,
                      0,(1-0.21)*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),0.21*survSap.h1(rain=min.rain,ad=min.ad.h1,par.samp[i]),
                      rec.h1,0,survM.h1(rain=min.rain,neigh=min.bi.h1,par.samp[i])),n.stage,n.stage)
  
  deltaR.h1.survM[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  
  # Cistus
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=max.rain,ad=max.ad.c,neigh=max.bi.c,par.samp[i])),n.stage,n.stage)
  
  mpm.min =  matrix(c(0,survS.c,0,
                      0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                      rec.c,0,survM.c(rain=min.rain,ad=min.ad.c,neigh=min.bi.c,par.samp[i])),n.stage,n.stage)
  
  deltaR.c.survM[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/sd.rain)
  
  
# INTERSPECIFIC DENSITY ###############################
  
  # 1. Adult survival (survM)
  
  # Halimium
  year.max.n=c(4,4,1,4,4,4,4,4,1,7,4,4,4,7,4,4,4,4)
  year.min.n=c(1,14,11,7,7,11,2,11,2,11,2,2,2,14,2,2,14,14)
  
  year.max=c(2010,2010,2007,2010,2010,2010,2010,2010,2007,2013,2010,2010,2010,2013,2010,2010,2010,2010)
  year.min=c(2007,2020,2017,2013,2013,2017,2008,2017,2008,2017,2008,2008,2008,2020,2008,2008,2020,2020)
  
  
  max.ad.h1=mean(c(ad.h1[1,3],ad.h1[2,3],ad.h1[3,1],ad.h1[4,3],ad.h1[5,3],ad.h1[6,3],ad.h1[7,3],ad.h1[8,3],ad.h1[9,1],
                   ad.h1[10,4],ad.h1[11,3],ad.h1[12,3],ad.h1[13,3],ad.h1[14,4],ad.h1[15,3],ad.h1[16,3],ad.h1[17,3],ad.h1[18,3]))
  
  min.ad.h1=mean(c(ad.h1[1,1],ad.h1[2,7],ad.h1[3,5],ad.h1[4,4],ad.h1[5,4],ad.h1[6,5],ad.h1[7,2],ad.h1[8,5],ad.h1[9,2],
                   ad.h1[10,5],ad.h1[11,2],ad.h1[12,2],ad.h1[13,2],ad.h1[14,7],ad.h1[15,2],ad.h1[16,2],ad.h1[17,7],ad.h1[18,7]))
  
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.max.n]),neigh=max.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean.rain,ad=mean.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean(rain[year.min.n]),neigh=min.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaInter.h1.survM[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.h1-min.inter.h1)/sd.inter.h1)
  
  # Cistus
  year.max.n=c(4,1,1,1,1,4,4,2,1,1,4,4,1,7,4,4,4,7)
  year.min.n=c(2,14,11,7,7,11,2,11,14,14,2,2,2,14,13,2,14,14)
  
  year.max=c(2010,2007,2007,2007,2007,2010,2010,2008,2007,2007,2010,2010,2007,2013,2010,2010,2010,2013)
  year.min=c(2008,2020,2017,2013,2013,2017,2008,2017,2020,2020,2008,2008,2008,2020,2019,2008,2020,2020)
  
  
  max.ad.c=mean(c(ad.c[1,3],ad.c[2,1],ad.c[3,1],ad.c[4,1],ad.c[5,1],ad.c[6,3],ad.c[7,3],ad.c[8,2],ad.c[9,1],
                  ad.c[10,1],ad.c[11,3],ad.c[12,3],ad.c[13,1],ad.c[14,4],ad.c[15,3],ad.c[16,3],ad.c[17,3],ad.c[18,4]))
  
  min.ad.c=mean(c(ad.c[1,2],ad.c[2,7],ad.c[3,5],ad.c[4,4],ad.c[5,4],ad.c[6,5],ad.c[7,2],ad.c[8,5],ad.c[9,7],
                  ad.c[10,7],ad.c[11,2],ad.c[12,2],ad.c[13,2],ad.c[14,7],ad.c[15,6],ad.c[16,2],ad.c[17,7],ad.c[18,7]))
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.max.n]),ad=max.ad.c,neigh=max.inter.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.min.n]),ad=min.ad.c,neigh=min.inter.c,par.samp[i])),n.stage,n.stage)
  
  deltaInter.c.survM[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.inter.c-min.inter.c)/sd.inter.c)
  
  
# INTRASPECIFIC DENSITY #################################
  # 1. Sapling survival (survSap)
  # Halimium
  year.max.n=c(4,7,1,1,1,4,1,2,2,1,4,4,1,1,4,1,1,2)
  year.min.n=c(2,14,11,14,14,14,14,14,14,14,14,11,14,14,14,14,14,14)
  
  max.bi.h1=mean(c(neighHal1[1,4],neighHal1[2,7],neighHal1[3,1],neighHal1[4,1],neighHal1[5,1],neighHal1[6,4],neighHal1[7,1],neighHal1[8,2],neighHal1[9,2],
                   neighHal1[10,1],neighHal1[11,4],neighHal1[12,4],neighHal1[13,1],neighHal1[14,1],neighHal1[15,4],neighHal1[16,1],neighHal1[17,1],neighHal1[18,2]))
  
  min.bi.h1=mean(c(neighHal1[1,2],neighHal1[2,14],neighHal1[3,11],neighHal1[4,14],neighHal1[5,14],neighHal1[6,14],neighHal1[7,14],neighHal1[8,14],neighHal1[9,14],
                   neighHal1[10,14],neighHal1[11,14],neighHal1[12,11],neighHal1[13,14],neighHal1[14,14],neighHal1[15,14],neighHal1[16,14],neighHal1[17,14],neighHal1[18,14]))
  
  mpm.max = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.max.n]),ad=max.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.max.n]),ad=max.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  mpm.min = matrix(c(0,survS.h1,0,
                     0,(1-0.21)*survSap.h1(rain=mean(rain[year.min.n]),ad=min.intra.h1,par.samp[i]),0.21*survSap.h1(rain=mean(rain[year.min.n]),ad=min.intra.h1,par.samp[i]),
                     rec.h1,0,survM.h1(rain=mean.rain,neigh=mean.inter.h1,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.h1.survSap[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.h1-min.intra.h1)/sd.intra.h1)
  
  # 1. Adult survival (survM)
  # Cistus
  year.max.n=c(4,4,2,1,4,2,4,4,1,13,4,4,7,7,4,4,4,4)
  year.min.n=c(11,11,11,14,14,11,2,14,2,14,2,1,1,14,14,1,14,2)
  
  max.bi.c=mean(c(neighCis[1,4],neighCis[2,4],neighCis[3,2],neighCis[4,1],neighCis[5,4],neighCis[6,2],neighCis[7,4],neighCis[8,4],neighCis[9,1],
                  neighCis[10,13],neighCis[11,4],neighCis[12,4],neighCis[13,7],neighCis[14,7],neighCis[15,4],neighCis[16,4],neighCis[17,4],neighCis[18,4]))
  
  min.bi.c=mean(c(neighCis[1,11],neighCis[2,11],neighCis[3,11],neighCis[4,14],neighCis[5,14],neighCis[6,11],neighCis[7,2],neighCis[8,14],neighCis[9,2],
                  neighCis[10,14],neighCis[11,2],neighCis[12,1],neighCis[13,1],neighCis[14,14],neighCis[15,14],neighCis[16,1],neighCis[17,14],neighCis[18,2]))
  
  mpm.max = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.max.n]),ad=max.intra.c,neigh=max.bi.c,par.samp[i])),n.stage,n.stage)
  
  
  mpm.min = matrix(c(0,survS.c,0,
                     0,(1-0.336)*survSap.c[par.samp[i]],0.336*survSap.c[par.samp[i]],
                     rec.c,0,survM.c(rain=mean(rain[year.min.n]),ad=min.intra.c,neigh=min.bi.c,par.samp[i])),n.stage,n.stage)
  
  deltaIntra.c.survM[i]=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.intra.c-min.intra.c)/sd.intra.c)
  
}

## save output ####################

# sens per vital rate

# Halimium
df.h2=data.frame(study.doi="10.1098/rspb.2022.1494",
                 year.of.publication="2023",
                 group="Plants",
                 species="Halimium halimifolium",
                 continent="Europe",
                 driver=rep(c("rain","rain","intraD","interD"),each=100),
                 driver.type=rep(c("C","C","D","B"),each=100), # intradensity = density and interdensity = biotic
                 stage.age=rep(c("sapling","adult"),each=100),
                 vital.rates="survival",
                 sens=c(deltaR.h1.survSap,
                        deltaR.h1.survM,
                        deltaIntra.h1.survSap,
                        deltaInter.h1.survM),
                 mat=5.8, # generation time
                 n.vr=2, # number of vital rates with covariates
                 n.pam=5, # number of parameters
                 dens=1,
                 biotic_interactions=1,
                 lambda.sim=0,
                 study.length=8)

write.csv(df.h2, "Sens_VR_Halimium.csv",row.names = F)


# Cistus
df.c2=data.frame(study.doi="10.1098/rspb.2022.1494",
                 year.of.publication="2023",
                 group="Plants",
                 species="Cistus libanotis",
                 continent="Europe",
                 driver=rep(c("rain","intraD","interD"),each=100),
                 driver.type=rep(c("C","D","B"),each=100), # intradensity = density and interdensity = biotic
                 stage.age="adult",
                 vital.rates="survival",
                 sens=c(deltaR.c.survM,
                        deltaIntra.c.survM,
                        deltaInter.c.survM),
                 mat=4.6,
                 n.vr=1, # number of vital rates with covariates
                 n.pam=4, # number of parameters of these vital rates
                 dens=1,
                 biotic_interactions=1,
                 lambda.sim=0,
                 study.length=8)

write.csv(df.c2, "Sens_VR_Cistus.csv",row.names = F)
