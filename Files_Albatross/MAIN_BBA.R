################################################################

# Original author: Stephanie Jenouvrier
# Study DOI: 10.1111/1365-2656.12827
# This code was converted from MATLAB to R via ChatGPT-3.5
# Edited by: Esin Ickin
# Date: 5.03.2024

#################################################################


# 0) Prepare session ###########################################

rm(list=ls())

# load packages
library(R.matlab)
library(popbio)

# set wd
setwd("~/Documents/Master Thesis/pert_analyses/Albatross")

# 1) Load data etc. #################################################

#load SST data
data <- readMat("SST.mat")

yearSSTobs <- 1982:2014

# Functions
source("invlogit.R")
source("invlogitG.R")
source("parameter_ENVSTO.R")
source("popmat.R")

# DESCRIPTION OF INPUTS

# SST = Sea Surface Temperature. This is a matrix of dimension ny*3. This is SST* in Figure 1 in Jenouvrier et al. 2018
# It contains SST for each of the three sectors for each year, so that SST is of dimension ny*3, with ny = the number of years.
# Col 1 = Juvenile sector, Col 2 = Winter sector and Col 3= Breeding sector.
SST <- as.data.frame(data$SSTsd2)
colnames(SST)=c("SST1","SST2","SST3")

# TEMP = Temperature. This is a vector of dimension ny*1. This is the temperature recorded by GLS for each year. It comes from devices set on the legs of individuals. SSTG in Figure 1 in
# Jenouvrier et al. 2018
# the data is from Figure S1b
SSTG <- read.csv("SSTG.csv")


# 2) Covariates ######################################

# Plot SST observations
#plot(yearSSTobs, SST[,1], type = "l", xlab = "Year", ylab = "SST", main = "SST Observations",col="darkgreen")
#lines(yearSSTobs, SST[,2], col = "darkblue")
#lines(yearSSTobs, SST[,3], col = "purple")
#legend("topleft", legend=c("Juvenile", "Wintering", "Breeding"), col=c("darkgreen", "darkblue", "purple"), lty=1)
# 1: SST* in the juvenile sector during the wintering season (May to Aug)
# 2: SST* in the wintering sector of adults (July to Sept)
# 3: SST* in the breeding sector (October of year t to March of year t+1)

# also: TEMP == SST[,2] because in the original code it stated that "temperature recorded by GLS correlates with SSTobs[, 2]"
# let's check:
#plot(yearSSTobs,SST$SST1,type="l",col="green")
#lines(SSTG$year,SSTG$SSTG,type="l",col="blue")
# problem: only years 2006 - 2013 of SSTG
# now get only SSTs for this period because otherwise we can't do covariations

SST$year=yearSSTobs
SST=SST[SST$year>= 2006,]
SST$SSTG=SSTG$SSTG
# remove NAs
SST=na.omit(SST)

# SST1
max.SST1=max(SST$SST1)
min.SST1=min(SST$SST1)
mean.SST1=mean(SST$SST1)
sd.SST1=sd(SST$SST1)

# SST2
max.SST2=max(SST$SST2)
min.SST2=min(SST$SST2)
mean.SST2=mean(SST$SST2)
sd.SST2=sd(SST$SST2)

# SST3
max.SST3=max(SST$SST3)
min.SST3=min(SST$SST3)
mean.SST3=mean(SST$SST3)
sd.SST3=sd(SST$SST3)

# SSTG
max.SSTG=max(SST$SSTG)
min.SSTG=min(SST$SSTG)
mean.SSTG=mean(SST$SSTG)
sd.SSTG=sd(SST$SSTG)

# covariation
# SST1
SST2_when_SST1_max=SST$SST2[which(SST$SST1==max(SST$SST1))][1]
SST2_when_SST1_min=SST$SST2[which(SST$SST1==min(SST$SST1))][1]

SST3_when_SST1_max=SST$SST3[which(SST$SST1==max(SST$SST1))][1]
SST3_when_SST1_min=SST$SST3[which(SST$SST1==min(SST$SST1))][1]

SSTG_when_SST1_max=SST$SSTG[which(SST$SST1==max(SST$SST1))][1]
SSTG_when_SST1_min=SST$SSTG[which(SST$SST1==min(SST$SST1))][1]

# SST2
SST1_when_SST2_max=SST$SST1[which(SST$SST2==max(SST$SST2))][1]
SST1_when_SST2_min=SST$SST1[which(SST$SST2==min(SST$SST2))][1]

SST3_when_SST2_max=SST$SST3[which(SST$SST2==max(SST$SST2))][1]
SST3_when_SST2_min=SST$SST3[which(SST$SST2==min(SST$SST2))][1]

SSTG_when_SST2_max=SST$SSTG[which(SST$SST2==max(SST$SST2))][1]
SSTG_when_SST2_min=SST$SSTG[which(SST$SST2==min(SST$SST2))][1]


# SST3
SST1_when_SST3_max=SST$SST1[which(SST$SST3==max(SST$SST3))][1]
SST1_when_SST3_min=SST$SST1[which(SST$SST3==min(SST$SST3))][1]

SST2_when_SST3_max=SST$SST2[which(SST$SST3==max(SST$SST3))][1]
SST2_when_SST3_min=SST$SST2[which(SST$SST3==min(SST$SST3))][1]

SSTG_when_SST3_max=SST$SSTG[which(SST$SST3==max(SST$SST3))][1]
SSTG_when_SST3_min=SST$SSTG[which(SST$SST3==min(SST$SST3))][1]


# SSTG
SST1_when_SSTG_max=SST$SST1[which(SST$SSTG==max(SST$SSTG))][1]
SST1_when_SSTG_min=SST$SST1[which(SST$SSTG==min(SST$SSTG))][1]

SST2_when_SSTG_max=SST$SST2[which(SST$SSTG==max(SST$SSTG))][1]
SST2_when_SSTG_min=SST$SST2[which(SST$SSTG==min(SST$SSTG))][1]

SST3_when_SSTG_max=SST$SST3[which(SST$SSTG==max(SST$SSTG))][1]
SST3_when_SSTG_min=SST$SST3[which(SST$SSTG==min(SST$SSTG))][1]


mean.theta=parameter_ENV(SST1 = mean.SST1, SST2 = mean.SST2, SST3 = mean.SST3, TEMP = mean.SSTG, ny = 7)
mean.theta=na.omit(mean.theta)
mean.A=popmat(mean.theta)
mean.lambda=lambda(mean.A)
mean.lambda

# 3) Sensitivity analysis ##############
# according to Morris et al. 2020

## 3.1 SST1 #########
delta.SST1=NULL
delta.SST1.cov=NULL
delta.SST2=NULL
delta.SST2.cov=NULL
delta.SST3=NULL
delta.SST3.cov=NULL
delta.SSTG=NULL
delta.SSTG.cov=NULL



# Start loop ##########
ny=7 # number of years
for(i in 1:100){
  ### max SST1 no cov ###################
  max.theta=parameter_ENV(SST1=max.SST1,SST2=mean.SST2,SST3=mean.SST3,TEMP=mean.SSTG,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SST1 no cov #########################
  min.theta=parameter_ENV(SST1=min.SST1,SST2=mean.SST2,SST3=mean.SST3,TEMP=mean.SSTG,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  #### Sens to SST1 no cov ------
  delta.SST1[i]=abs((max.lam-min.lam)/((max.SST1-min.SST1)/sd.SST1))
  
  ### max SST1 cov ###########################
  max.theta=parameter_ENV(SST1=max.SST1,SST2=SST2_when_SST1_max,SST3=SST3_when_SST1_max,TEMP=SSTG_when_SST1_max,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam.cov=lambda(max.A)
  
  ### min SST1 cov ##############################
  min.theta=parameter_ENV(SST1=min.SST1,SST2=SST2_when_SST1_min,SST3=SST3_when_SST1_min,TEMP=SSTG_when_SST1_min,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam.cov=lambda(min.A)
  
  #### Sens to SST1 cov ------
  delta.SST1.cov[i]=abs((max.lam.cov-min.lam.cov)/((max.SST1-min.SST1)/sd.SST1))
  
  
  ## 3.2 SST2 #########
  ### max SST2 no cov ###################
  max.theta=parameter_ENV(SST1=mean.SST1,SST2=max.SST2,SST3=mean.SST3,TEMP=mean.SSTG,ny=ny)
  #max.A=popmat(mean(t(max.theta)))
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SST2 no cov ###################
  min.theta=parameter_ENV(SST1=mean.SST1,SST2=min.SST2,SST3=mean.SST3,TEMP=mean.SSTG,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  #### Sens to SST2 no cov ------
  delta.SST2[i]=abs((max.lam-min.lam)/((max.SST2-min.SST2)/sd.SST2))
  
  
  ### max SST2 cov ###################
  max.theta=parameter_ENV(SST1=SST1_when_SST2_max,SST2=max.SST2,SST3=SST3_when_SST2_max,TEMP=SSTG_when_SST2_max,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  
  ### min SST2 cov ###################
  min.theta=parameter_ENV(SST1=SST1_when_SST2_min,SST2=min.SST2,SST3=SST3_when_SST2_min,TEMP=SSTG_when_SST2_min,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  #### Sens to SST2 cov ------
  delta.SST2.cov[i]=abs((max.lam-min.lam)/((max.SST2-min.SST2)/sd.SST2))
  
  ## 3.3 SSTG #########
  ### max SST3 no cov ###################
  max.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=max.SST3,TEMP=mean.SSTG,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SST3 no cov ###################
  min.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=min.SST3,TEMP=mean.SSTG,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  #### Sens To SST3 no cov --------
  delta.SST3[i]=abs((max.lam-min.lam)/((max.SST3-min.SST3)/sd.SST3))
  
  
  ### max SST3 cov ###################
  max.theta=parameter_ENV(SST1=SST1_when_SST3_max,SST2=SST2_when_SST3_max,SST3=max.SST3,TEMP=SSTG_when_SST3_max,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SST3 cov ###################
  min.theta=parameter_ENV(SST1=SST1_when_SST3_min,SST2=SST2_when_SST3_min,SST3=min.SST3,TEMP=SSTG_when_SST3_min,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  
  #### Sens To SST3 cov --------
  delta.SST3.cov[i]=abs((max.lam-min.lam)/((max.SST3-min.SST3)/sd.SST3))
  
  
  ## 3.4 SSTG #########
  ### max SSTG no cov ###################
  max.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=mean.SST3,TEMP=max.SSTG,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SSTG no cov ###################
  min.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=mean.SST3,TEMP=min.SSTG,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  #### Sens To SSTG no cov --------
  delta.SSTG[i]=abs((max.lam-min.lam)/((max.SSTG-min.SSTG)/sd.SSTG))
  
  
  ### max SSTG cov ###################
  max.theta=parameter_ENV(SST1=SST1_when_SSTG_max,SST2=SST2_when_SSTG_max,SST3=SST3_when_SSTG_max,TEMP=max.SSTG,ny=ny)
  max.A=popmat(na.omit(max.theta))
  max.lam=lambda(max.A)
  
  ### min SSTG cov ###################
  min.theta=parameter_ENV(SST1=SST1_when_SSTG_min,SST2=SST2_when_SSTG_min,SST3=SST3_when_SSTG_min,TEMP=min.SSTG,ny=ny)
  min.A=popmat(na.omit(min.theta))
  min.lam=lambda(min.A)
  
  
  #### Sens To SSTG cov --------
  delta.SSTG.cov[i]=abs((max.lam-min.lam)/((max.SSTG-min.SSTG)/sd.SSTG))
  

}

# 4) Output ##########
hist(delta.SST1)
hist(delta.SST1.cov)
hist(delta.SST2)
hist(delta.SST2.cov)
hist(delta.SST3)
hist(delta.SST3.cov)
hist(delta.SSTG)
hist(delta.SSTG.cov)

# save output
Results=data.frame(study.doi="110.1111/1365-2656.12827",
           year.of.publication="2018",
           group="Birds",
           species="Thalassarche melanophris",
           continent="Sub-Antarctic",
           driver=c(rep(c("SST1","SST2","SST3","SSTG"),each=200)),
           driver.type="C",
           stage.age="all",
           vital.rates="all",
           sens=c(delta.SST1,delta.SST1.cov,
                  delta.SST2,delta.SST2.cov,
                  delta.SST3,delta.SST3.cov,
                  delta.SSTG,delta.SSTG.cov),
           cov=rep(c(0,1),each=100),
           mat=10.6, # age at first reproduction, source: Myhrvold et al. 2015
           n.vr=10, # number of vital rates with covariates
           n.pam=47, # number of total parameters in these vital rates
           dens=0, # density dependence in it?
           biotic_interactions=0, # any biotic interactions?
           lambda.sim=0, # was lambda calculated analytically or using simulation?
           study.length=35) # study length in years

write.csv(Results,"Sens_Albatross.csv",row.names = F)
