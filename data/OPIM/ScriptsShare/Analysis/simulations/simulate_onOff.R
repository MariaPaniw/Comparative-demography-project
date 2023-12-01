#Stochastic simulations of populationi projection models whose vital rate
#correlations are either all 'on' or all 'off', whereby 'off' means set to zero. 
#This is done for:
#1. Simulations using mean joint posterior parameters.
#2. Simulations using 100 joint posterior draws.

#Load packages and the random sequences of integers to seed stochastic simulations
library(lme4) ; library(msm) ; library(Matrix) ; library(mvtnorm)
rseed.vec     <- read.csv("Data/rseed.vec100K.csv")[,1] #seed for stochastic simulations

##1. Simulations using "mean joint posterior parameters"#########################
meanParameters=NULL
source("Analysis/populationModels/helianthella_MPM.R")
onH=projection_onOff("on","helianthella")
offH=projection_onOff("off","helianthella")

source("Analysis/populationModels/opuntia_IPM.R")
onC=projection_onOff("on","opuntia")
offC=projection_onOff("off","opuntia")

source("Analysis/populationModels/orchis_IPM.R")
onO=projection_onOff("on","orchis")
offO=projection_onOff("off","orchis")

#Write out stochastic lambdas and lambda_t variances  
lambdaS      <- c(exp(mean(onH)),exp(mean(offH)),
                  exp(mean(onC)),exp(mean(offC)),
                  exp(mean(onO)),exp(mean(offO)))
varLambda_t  <- c(var(exp(onH)),var(exp(offH)),
                  var(exp(onC)),var(exp(offC)),
                  var(exp(onO)),var(exp(offO)))
meanParametersLambda=data.frame(correlation=rep(c("on","off"),3),
                                species=c("helianthella","opuntia","orchis"),
                                lambdaS,varLambda_t)
write.csv(meanParametersLambda,"Results/simulations/meanParamLam.csv",row.names=F)


##2. Simulations using 100 joint posterior draws#################################

##Helianthella---------------------------------------------------------------------
on_sim=off_sim=matrix(NA,45000,100) #series of 45000 realized lambdas (lambda_t in the text)
source("Analysis/populationModels/helianthella_MPM.R")
for(samID in 1:100){
  
  parametersHeli(samID)
  on_sim[,samID]=projection_onOff("on","helianthella") #series of log(lambdas_t)
  off_sim[,samID]=projection_onOff("off","helianthella") #series of log(lambdas_t)

}
#Write results out
write.csv(on_sim,"Results/simulations/heli_lambdaS_on.csv",row.names=F)  
write.csv(off_sim,"Results/simulations/heli_lambdaS_off.csv",row.names=F) 
rm(list=c("on_sim","off_sim"))


#Opuntia----------------------------------------------------------------------------
on_sim=off_sim=matrix(NA,45000,100) #series of 45000 realized lambdas (lambda_t in the text)
source("Analysis/populationModels/opuntia_IPM.R")
for(samID in 1:100){
  
  parametersOpuntia(samID)
  on_sim[,samID]=projection_onOff("on","opuntia")#series of log(lambdas_t)
  off_sim[,samID]=projection_onOff("off","opuntia")#series of log(lambdas_t)

}
#Write results out
write.csv(on_sim,"Results/simulations/opuntia_lambdaS_on.csv",row.names=F)  
write.csv(off_sim,"Results/simulations/opuntia_lambdaS_off.csv",row.names=F) 
rm(list=c("on_sim","off_sim"))


##Orchis-----------------------------------------------------------------------------
on_sim=off_sim=matrix(NA,45000,100) #series of 45000 realized lambdas (lambda_t in the text)
source("Analysis/populationModels/orchis_IPM.R")
for(samID in 1:100){
  
  parametersOrchis(samID)
  on_sim[,samID]=projection_onOff("on","orchis") #series of log(lambdas_t)
  off_sim[,samID]=projection_onOff("off","orchis") #series of log(lambdas_t)

}
#Write results out
write.csv(on_sim,"Results/simulations/orchis_lambdaS_on.csv",row.names=F)  
write.csv(off_sim,"Results/simulations/orchis_lambdaS_off.csv",row.names=F) 
rm(list=c("on_sim","off_sim"))

