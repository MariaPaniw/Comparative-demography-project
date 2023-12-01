#Stochastic simulations of populationi projection models with 
#only one pairwise vital rate correlation turned off at a time.
#the function projection_pairwiseOnOff() calculates var(lambda_t) for
#each 45000 time step simulation run. 

#Load packages and the random sequences of integers to seed stochastic simulations
library(lme4) ; library(msm) ; library(Matrix) ; library(mvtnorm)
rseed.vec     <- read.csv("Data/rseed.vec100K.csv")[,1] #seed for stochastic simulations

#HELIANTHELLA-------------------------------------------------------------------------
heli_CorrOnOff=as.data.frame(matrix(NA,100,10)) #Placeholders for variance values 
source("Analysis/populationModels/helianthella_MPM.R")
#Loop for Orchis
for(samID in 1:100){
  
  parametersHeli(samID) #Read specific posterior sample values
  heli_CorrOnOff[simI,]  <- projection_pairwiseOnOff('helianthella')

}

#FORMAT DATA: Column names are the "vital rate combinations" (vrCombination)
vrCombination=which(upper.tri(matrix(0,5,5)),arr.ind=T)
vrCombination[vrCombination==1]="g" ; vrCombination[vrCombination==2]="s" ; vrCombination[vrCombination==3]="f" 
vrCombination[vrCombination==4]="fe" ; vrCombination[vrCombination==5]="fr" 

#Names of combinations
names(heli_CorrOnOff)=paste0(vrCombination[,1],vrCombination[,2])
#Introduce sample number
heli_CorrOnOff$sample=c(1:100)

#make files "long"
heli_CorrOnOff=melt(heli_CorrOnOff,id.vars="sample",measure.vars = names(heli_CorrOnOff)[1:10])
names(heli_CorrOnOff)[3]="lambdaVar" 
names(heli_CorrOnOff)[2]="comb"

#write output
write.csv(heli_CorrOnOff,"Results/simulations/heli_pairwise_onOff.csv",row.names=F)


#OPUNTIA-------------------------------------------------------------------------
opuntia_CorrOnOff=as.data.frame(matrix(NA,100,6)) #Placeholders for variance
source("Analysis/populationModels/opuntia_IPM.R")
#Loop for Orchis
for(samID in 1:100){
  
  parametersOpuntia(samID)
  opuntia_CorrOnOff[simI,]  <- projection_pairwiseOnOff('opuntia')

}

#FORMAT DATA: Column names are the "vital rate combinations" (vrCombination)
vrCombination=which(upper.tri(matrix(0,4,4)),arr.ind=T)
vrCombination[vrCombination==1]="g" ; vrCombination[vrCombination==2]="s" ; vrCombination[vrCombination==3]="f" ; vrCombination[vrCombination==4]="fe"

#Names of combinations
names(opuntia_CorrOnOff)=paste0(vrCombination[,1],vrCombination[,2])
#Introduce sample number
opuntia_CorrOnOff$sample=c(1:100)

#make files "long"
opuntia_CorrOnOff=melt(opuntia_CorrOnOff,id.vars="sample",measure.vars = names(opuntia_CorrOnOff)[1:6])
names(opuntia_CorrOnOff)[3]="lambdaVar" 
names(opuntia_CorrOnOff)[2]="comb"

#write output
write.csv(opuntia_CorrOnOff,"Results/simulations/opuntia_pairwise_onOff.csv",row.names=F)


#ORCHIS-------------------------------------------------------------------------
orchis_CorrOnOff=as.data.frame(matrix(NA,100,15)) #Placeholders for variance
source("Analysis/populationModels/orchis_IPM.R")
#Loop for Orchis
for(samID in 1:100){
  
  parametersOrchis(samID)
  orchis_CorrOnOff[simI,] <- projection_pairwiseOnOff('orchis')

}

#FORMAT DATA: Column names are the "vital rate combinations" (vrCombination)
vrCombination=which(upper.tri(matrix(0,6,6)),arr.ind=T)
vrCombination[vrCombination==1]="g" ; vrCombination[vrCombination==2]="s" ; vrCombination[vrCombination==3]="f" ; vrCombination[vrCombination==4]="fe"
vrCombination[vrCombination==4]="fe" ; vrCombination[vrCombination==5]="fr" ; vrCombination[vrCombination==6]="d" 

#Names of combinations
names(orchis_CorrOnOff)=paste0(vrCombination[,1],vrCombination[,2])
#Introduce sample number
orchis_CorrOnOff$sample=c(1:100)

#make files "long"
orchis_CorrOnOff=melt(orchis_CorrOnOff,id.vars="sample",measure.vars = names(orchis_CorrOnOff)[1:15])
names(orchis_CorrOnOff)[3]="lambdaVar"
names(orchis_CorrOnOff)[2]="comb"

#write file out
write.csv(orchis_CorrOnOff,"Results/simulations/orchis_pairwise_onOff.csv",row.names=F)
