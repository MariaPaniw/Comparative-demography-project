setwd("C:/Users/ac79/Downloads/Dropbox/Bayes group project/ms/Revision/ScriptsShare")

#Vital rate model fitting (these might take a few hours to fit)########################
source("Analysis/vitalRates/heli_vr_estimation.R")
source("Analysis/vitalRates/opuntia_vr_estimation.R")
source("Analysis/vitalRates/orchis_vr_estimation.R")


#Load functions needed to run simulations##############################################

#Functions that load one of the 100 joint posterior draws------------------------------
source("Analysis/simulations/parametersHeli.R")     #Helianthella posterior draws
source("Analysis/simulations/parametersOpuntia.R")  #Opuntia posterior draws
source("Analysis/simulations/parametersOrchis.R")   #Orchis posterior draws

#Functions to run stochastic simulations-----------------------------------------------
#Function runs stochastic simulations where 
#correlations are either all "on" or all "off"
source("Analysis/simulations/projection_onOff.R")
#Function runs stochastic simulations where only one 
#pairwise correlations at a time is turned "off" 
source("Analysis/simulations/projection_pairwiseOnOff.R")


#Simulation runs#######################################################################

#Each of these two scripts will take up to several days to run.
#If possible, we recommend running as many separate simulations as possible on a supercomputer.
source("Analysis/simulations/simulate_onOff.R") #simulations for Figure 4 in the main text 
source("Analysis/simulations/simulate_pairwiseOnOff.R")  #simulations for Figure 3 in the main text
