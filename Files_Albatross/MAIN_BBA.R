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
setwd("~/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/Albatross")

# 1) Load data etc. #################################################

#load SST data
data <- readMat("SST.mat")
SST <- as.data.frame(data$SSTsd2)
colnames(SST)=c("SST1","SST2","SST3")

yearSSTobs <- 1982:2014

# Functions
source("invlogit.R")
source("invlogitG.R")
source("parameter_ENV.R")
source("popmat.R")

# 2) Covariates ######################################

# Plot SST observations
plot(yearSSTobs, SST[,1], type = "l", xlab = "Year", ylab = "SST", main = "SST Observations",col="darkgreen")
lines(yearSSTobs, SST[,2], col = "darkblue")
lines(yearSSTobs, SST[,3], col = "purple")
legend("topleft", legend=c("Juvenile", "Wintering", "Breeding"), col=c("darkgreen", "darkblue", "purple"), lty=1)
# 1: SST* in the juvenile sector during the wintering season (May to Aug)
# 2: SST* in the wintering sector of adults (July to Sept)
# 3: SST* in the breeding sector (October of year t to March of year t+1)

# also: TEMP == SST[,2] because in the original code it stated that "temperature recorded by GLS correlates with SSTobs[, 2]"

# SST1
max.SST1=max(SST[,1])
min.SST1=min(SST[,1])
mean.SST1=mean(SST[,1])
sd.SST1=sd(SST[,1])

# SST2
max.SST2=max(SST[,2])
min.SST2=min(SST[,2])
mean.SST2=mean(SST[,2])
sd.SST2=sd(SST[,2])

# SST3
max.SST3=max(SST[,3])
min.SST3=min(SST[,3])
mean.SST3=mean(SST[,3])
sd.SST3=sd(SST[,3])

# covariation
# SST1
SST2_when_SST1_max=SST$SST2[which(SST$SST1==max(SST$SST1))][1]
SST2_when_SST1_min=SST$SST2[which(SST$SST1==min(SST$SST1))][1]

SST3_when_SST1_max=SST$SST3[which(SST$SST1==max(SST$SST1))][1]
SST3_when_SST1_min=SST$SST3[which(SST$SST1==min(SST$SST1))][1]

# SST2
SST1_when_SST2_max=SST$SST1[which(SST$SST2==max(SST$SST2))][1]
SST1_when_SST2_min=SST$SST1[which(SST$SST2==min(SST$SST2))][1]

SST3_when_SST2_max=SST$SST3[which(SST$SST2==max(SST$SST2))][1]
SST3_when_SST2_min=SST$SST3[which(SST$SST2==min(SST$SST2))][1]

# SST3
SST1_when_SST3_max=SST$SST1[which(SST$SST3==max(SST$SST3))][1]
SST1_when_SST3_min=SST$SST1[which(SST$SST3==min(SST$SST3))][1]

SST2_when_SST3_max=SST$SST2[which(SST$SST3==max(SST$SST3))][1]
SST2_when_SST3_min=SST$SST2[which(SST$SST3==min(SST$SST3))][1]

# new function for parameter_ENV
# adapted from parameter_ENV.R
parameter_ENV <- function(SST1,SST2,SST3,ny=1,i){ 
  # in the original code ny = length(SST), but why?
  # but any other number (i.e. 1) gives the same values
  # i for resampling and getting the uncertainties
  
  # INITIALIZE theta matrix
  theta <- matrix(0, nrow = 68, ncol = ny)
  
  # SURVIVAL for the 25 classes
  # First year survival (sigma 1)
  theta[1, ] <- inv.logit(-0.23 - 0.24 * SST1 - 0.32 * SST1^2)
  
  # Survival for the other 24 classes (sigma 2 - sigma 25)
  theta[2:10, ] <- rep(0.93, ny)
  theta[11:16, ] <- rep(0.94, ny)
  theta[17:22, ] <- rep(0.92, ny)
  theta[23, ] <- rep(0.94, ny)
  theta[24:25, ] <- rep(0.92, ny)
  
  # BREEDING PROBABILITIES AND SUCCESS
  
  # PRE-BREEDERS
  # Recruitment (6 parameters)
  # !!!!! theta[26:31, ] <- sapply(4:9, function(i) plogis(-4.891594061 + (-1.382420815 + 0.563784374 * i) * SST3))
  
  # Function of age [different intercepts for each age] 
  theta[26, ] = inv.logit(-4.891594061); # age 4; beta 5
  theta[27, ] = inv.logit(-3.556583792); # age 5; beta 6
  theta[28, ] = inv.logit(-2.247152916); # age 6; beta 7
  theta[29, ] = inv.logit(-1.511223273 ); # age 7; beta 8
  theta[30, ] = inv.logit(-1.108869401); # age 8; beta 9
  theta[31, ] = inv.logit(-0.772677029 ); # age 9 +; beta 10
  # Breeding success (6 parameters)
  theta[32, ] <- inv.logit(-0.300619322 - 0.154409282 * SST2 - 0.049391794 * SST2^2 +
                             0.135000818 * SST3 + 0.135000818 * SST3^2)
  theta[33,] <- theta[32]
  theta[34,] <- theta[32]
  theta[35,] <- theta[32]
  theta[36,] <- theta[32]
  theta[37,] <- theta[32]
  
  # FIRST-TIME BREEDERS
  # Breeding probability
  bs <- 1.284297702 - 0.154409282 * SST2 - 0.049391794 * SST2^2 +
    0.135000818 * SST3 + 0.124121331 * SST3^2
  cs <- 0.65367776 + 0.030113582 * SST2 - 0.057345669 * SST2^2 +
    0.2029 * SST3 + 0.0783 * SST3^2
  theta[38:43, ] <- invlogitG(bs,cs) + invlogitG(cs,bs)
  
  bf <- 0.81440189 - 0.154409282 * SST2 - 0.049391794 * SST2^2 +
    0.135000818 * SST3 + 0.124121331 * SST3^2
  cf <- 0.390548647 + 0.030113582 * SST2 - 0.057345669 * SST2^2 +
    0.2029 * SST3 + 0.0783 * SST3^2
  theta[44:49, ] <- invlogitG(bf,cf) + invlogitG(cf,bf)
  
  # Breeding success
  RETs <- (exp(rnorm(100,mean=5.494,sd=0.01)[i] + rnorm(100,mean=0.05,sd=0.01)[i] - rnorm(100,mean=0.015,sd=0.006)[i] * SST2) - 248.8625954) / 10.59881224
  theta[53:58, ] <- plogis(0.9785 + 0.4366 * RETs + 0.3161)
  # maybe this is return date
  # if so, then some of the variables might be from Table S2.4b of the study Desprez et al. 2018
  # then the intercept: 5.494 (SE=0.01)
  # breeding status: 0.05 (SE=0.01)
  # sex: doesn't match (Table: -0.024 SE=0.01)
  # Temp2: doesn't match perfectly (Table: -0.012, SE=0.006)
  
  
  RETf <- (exp(rnorm(100,mean=5.494,sd=0.01)[i] - rnorm(100,mean=0.015,sd=0.006)[i] * SST2) - 248.8625954) / 10.59881224
  theta[59:64, ] <- plogis(0.9785 + 0.4366 * RETf + 0.3161)
  # intercept same as above and temp similar
  
  
  # EXPERIENCED BREEDERS
  theta[50, ] <- theta[38, ]
  theta[51, ] <- theta[44, ]
  
  bnb <- 0.013026711 - 0.154409282 * SST2 - 0.049391794 * SST2^2 +
    0.135000818 * SST3 + 0.124121331 * SST3^2
  cnb <- -0.2112 + 0.030113582 * SST2 - 0.057345669 * SST2^2 +
    0.2029 * SST3 + 0.0783 * SST3^2
  theta[52, ] <- invlogitG(bnb,cnb) + invlogitG(cnb,bnb)
  
  theta[65:67, ] <- theta[53:55, ]
  
  # Fecundity per females = sex ratio
  theta[68, ] <- rep(0.50, ny)
  
  return(theta)
}

# 3) Sensitivity analysis ##############
# according to Morris et al. 2020 (DOI:...)

## 3.1 SST1 #########
delta.SST1=NULL
delta.SST1.cov=NULL
delta.SST2=NULL
delta.SST2.cov=NULL
delta.SST3=NULL
delta.SST3.cov=NULL

# Start loop ##########
for(i in 1:100){
  ### max SST1 no cov ###################
  max.theta=parameter_ENV(SST1=max.SST1,SST2=mean.SST2,SST3=mean.SST3,i=i)
  max.A=popmat(max.theta[,1])
  max.lam=lambda(max.A)
  
  ### min SST1 no cov #########################
  min.theta=parameter_ENV(SST1=min.SST1,SST2=mean.SST2,SST3=mean.SST3,i=i)
  min.A=popmat(min.theta)
  min.lam=lambda(min.A)
  
  #### Sens to SST1 no cov ------
  delta.SST1[i]=abs((max.lam-min.lam)/((max.SST1-min.SST1)/sd.SST1))
  
  ### max SST1 cov ###########################
  max.theta=parameter_ENV(SST1=max.SST1,SST2=SST2_when_SST1_max,SST3=SST3_when_SST1_max,i=i)
  max.A=popmat(max.theta)
  max.lam.cov=lambda(max.A)
  
  ### min SST1 cov ##############################
  min.theta=parameter_ENV(SST1=min.SST1,SST2=SST2_when_SST1_min,SST3=SST3_when_SST1_min,i=i)
  min.A=popmat(min.theta)
  min.lam.cov=lambda(min.A)
  
  #### Sens to SST1 cov ------
  delta.SST1.cov[i]=abs((max.lam.cov-min.lam.cov)/((max.SST1-min.SST1)/sd.SST1))
  
  
  ## 3.2 SST2 #########
  ### max SST2 no cov ###################
  max.theta=parameter_ENV(SST1=mean.SST1,SST2=max.SST2,SST3=mean.SST3,i=i)
  #max.A=popmat(mean(t(max.theta)))
  max.A=popmat(max.theta)
  max.lam=lambda(max.A)
  
  ### min SST2 no cov ###################
  min.theta=parameter_ENV(SST1=mean.SST1,SST2=min.SST2,SST3=mean.SST3,i=i)
  min.A=popmat(min.theta)
  min.lam=lambda(min.A)
  
  #### Sens to SST2 no cov ------
  delta.SST2[i]=abs((max.lam-min.lam)/((max.SST2-min.SST2)/sd.SST2))
  
  
  ### max SST2 cov ###################
  max.theta=parameter_ENV(SST1=SST1_when_SST2_max,SST2=max.SST2,SST3=SST3_when_SST2_max,i=i)
  #max.A=popmat(mean(t(max.theta)))
  #max.lam=lambda(max.A)
  
  #max.A=popmat(t(max.theta))
  #max.lam=lambda(max.A)
  #max.lam # way too small, doesn't make sense
  
  max.A=popmat(max.theta)
  max.lam=lambda(max.A)
  
  
  ### min SST2 cov ###################
  min.theta=parameter_ENV(SST1=SST1_when_SST2_min,SST2=min.SST2,SST3=SST3_when_SST2_min,i=i)
  min.A=popmat(min.theta)
  min.lam=lambda(min.A)
  
  #### Sens to SST2 cov ------
  delta.SST2.cov[i]=abs((max.lam-min.lam)/((max.SST2-min.SST2)/sd.SST2))
  
  
  
  ## 3.3 SST3 #########
  ### max SST3 no cov ###################
  max.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=max.SST3,i=i)
  max.A=popmat(max.theta)
  max.lam=lambda(max.A)
  
  ### min SST3 no cov ###################
  min.theta=parameter_ENV(SST1=mean.SST1,SST2=mean.SST2,SST3=min.SST3,i=i)
  min.A=popmat(min.theta)
  min.lam=lambda(min.A)
  
  #### Sens To SST3 no cov --------
  delta.SST3[i]=abs((max.lam-min.lam)/((max.SST3-min.SST3)/sd.SST3))
  
  
  ### max SST3 cov ###################
  max.theta=parameter_ENV(SST1=SST1_when_SST3_max,SST2=SST2_when_SST3_max,SST3=max.SST3,i=i)
  max.A=popmat(max.theta)
  max.lam=lambda(max.A)
  
  ### min SST3 cov ###################
  min.theta=parameter_ENV(SST1=SST1_when_SST3_min,SST2=SST2_when_SST3_min,SST3=min.SST3,i=i)
  min.A=popmat(min.theta)
  min.lam=lambda(min.A)
  
  
  #### Sens To SST3 cov --------
  delta.SST3.cov[i]=abs((max.lam-min.lam)/((max.SST3-min.SST3)/sd.SST3))
  
  
  
}

# 4) Output ##########
hist(delta.SST1)
hist(delta.SST1.cov)
hist(delta.SST2)
hist(delta.SST2.cov)
hist(delta.SST3)
hist(delta.SST3.cov)

# save output
Results=data.frame(study.doi="110.1111/1365-2656.12827",
           year.of.publication="2018",
           group="Birds",
           species="Thalassarche melanophris",
           continent="Sub-Antarctic",
           driver=c(rep(c("SST1","SST2","SST3"),each=200)),
           driver.type="C",
           stage.age="all",
           vital.rates="all",
           sens=c(delta.SST1,delta.SST1.cov,
                  delta.SST2,delta.SST2.cov,
                  delta.SST3,delta.SST3.cov),
           cov=rep(c(0,1),each=100),
           gen.time=9.8, # generation time [average age at first breeding is 9·8 years, Weimerskirch, Clobert & Jouventin 1987]
           n.vr=10, # number of vital rates with covariates
           n.pam=47, # number of total parameters in these vital rates
           dens=0, # density dependence in it?
           biotic_interactions=0) # any biotic interactions?

write.csv(Results,"Sens_Albatross.csv",row.names = F)
