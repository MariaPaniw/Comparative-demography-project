# Script to reorganize sensitivities of Dracocephalum #############

# 0. Prepare session ##########################

# empty R's brain
rm(list=ls())

# load packages
library(tidyverse)

# set wd
setwd("/Users/esinickin/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/Dracocephalum")


# 1. Load data ####################
sens_no_cov=read.csv("SensNoCov.csv")
head(sens_no_cov)

sens_cov=read.csv("SensCov.csv")
head(sens_cov)

hist(sens_no_cov$sens.pet.no.cov-sens_cov$sens.pet.cov)
hist(sens_no_cov$sens.precip.no.cov-sens_cov$sens.precip.cov)

pet_no_cov=sens_no_cov$sens.pet.no.cov
pr_no_cov=sens_no_cov$sens.precip.no.cov
pet_cov=sens_cov$sens.pet.cov
pr_cov=sens_cov$sens.precip.cov


# 2. Create new dataframe #################

Sens_Dracocephalum=data.frame(species="Dracocephalum austriacum",
                               study.doi="Evers et al. in Prep",
                               year.of.publication="2024",
                               group="Plants",
                               continent="Europe",
                               driver=rep(c("PET","precipitation"),each=200),
                               driver.type="C",
                               stage.age="all",
                               vital.rates="all",
                               sens=c(pet_no_cov,pet_cov,pr_no_cov,pr_cov),
                               cov=rep(c(0,1),each=100),
                               gen.time=2, # Source: Dostálek & Münzbergová 2013
                               n.vr=5, # number of total vital rates
                               n.pam=45, #length of params df (number of parameters)
                               dens=0,
                               biotic_interactions=0)

write.csv(Sens_Dracocephalum, "Sensitivities_Dracocephalum.csv",row.names = F)
