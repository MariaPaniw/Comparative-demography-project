
# Sensitivity analysis for Magellanic Penguins ##

# Original script by: T.J. Clark-Wolf, P. Dee Boersma, Ginger A. Rebstock, and Briana Abrahms
# Study: Clark-Wolf et al. 2023
# Study DOI: https://www.pnas.org/doi/10.1073/pnas.2209821120

# This script contains an IPM
# and sensitivity analysis done by Esin Ickin

# Date: 31.01.2024


# 0) Prepare session ###########################

# FIRST: run the script "CRData.R"

# set directory
setwd("~/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/MagellanicPenguin/Data")

# load packages
library(tidyverse)
library(nimble)
library(ggplot2)
library(MCMCvis)
library(lubridate)


# 1) Parameters ########################################3

#### Load population data ###############

# REPRODUCTION: NUMNESTS = number of nests, numfledged = number fledged!
repro <- read.csv("ReproSuccess_9.21.21.csv", header = F)
repro$V1[1] <- "1984" # fixed!
colnames(repro) <- c("bookyear", "numnests", "numeggs", "clutchsize", "numhatch", "meanhatch", "numfledged", "rs")
repro$bookyear <- as.Date(as.character(repro$bookyear), "%Y")
repro <- repro %>% arrange(bookyear)
# missing 1982, 2011

# number of surveyed nests
R <- as.numeric(c(NA, repro$numnests[1:28], NA, repro$numnests[29:36]))

# number of fledges
J <- as.numeric(c(NA, repro$numfledged[1:28], NA, repro$numfledged[29:36]))

# POPULATION DATA (using Ginger's normal adjustment, not sliding scale...)
# NOTE: I have 2020 pop data, going to exclude for now...
pop <- read.csv("StakeSurveyTOMACT22_9.21.21.csv", header = F)
pop[1,1] <- 1987 # fix weird stuff
colnames(pop) <- c("bookyear","00N01E", "00N02E", "00N03E", "00N04E", "00N05E", "00N06E", "00N07E", "00N08E",
                   "00N09E", "00N12E", "00N14E", "00N15E", "00N16E", "00N17E", "00N18E", "01N14E", "01S14E",
                   "02N14E", "02S14E", "03S14E", "04S14E", "05S14E", "total")
pop$density <- pop$total/22 # number of 100 meter survey areas = 22

# convert by 1.36
pop$blah <- c(pop$density[1:28]/1.36, pop$density[29:33])
pop$popest2 <- pop$blah*35240 # total habitat reported in the paper...(Table 1, 2012 survey)

y <- c(rep(NA, 5), pop$popest2[1:24], NA, pop$popest2[25:32])

# SURVIVAL DATA - see "CRData.R" to load data
# drop banded adults never sexed (we assume all juveniles are not sexed)
# drop banded juveniles that were recaptured and never sexed
combo$pengid <- as.character(combo$pengid)
combo$age <- as.character(combo$age)
combo_b <- combo %>%
  mutate(sum = rowSums(across(where(is.numeric)))) # sum across capture histories

combo2 <- combo_b %>%
  filter(age == 2, sex == "NULL") %>%
  mutate(sex = sample(c("M","F"), n(), prob = c(0.5,0.5), replace = T))

combo2.5 <- combo_b %>%
  filter(age == 2, sex != "NULL")

combo3 <- combo_b %>%
  filter(age == 0) %>%
  bind_rows(combo2) %>%
  bind_rows(combo2.5)

# randomly assign sex to recaptured, unsexed individuals (n = 784)
combo4 <- combo3 %>%
  filter(age == 0, sum > 1, sex == "NULL") %>%
  mutate(sex = sample(c("M","F"), n(), prob = c(0.5,0.5), replace = T))

# retain unrecaptured, unsexed individuals
combo4.5 <- combo3 %>%
  filter(age == 0, sum <= 1, sex == "NULL")

# retain unsexed individuals
combo5 <- combo3 %>%
  filter(age == 0, sex != "NULL")

combo6 <- combo3 %>%
  filter(age == 2) %>%
  bind_rows(combo4) %>%
  bind_rows(combo4.5) %>%
  bind_rows(combo5)

# individuals marked as juveniles
CH.J <- combo6 %>% filter(age == 0) %>%
  dplyr::select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)

# individuals marked as adults, sexed as females
CH.A <- combo6 %>% filter(age == 2) %>%
  filter(sex == "F") %>%
  dplyr::select(-pengid, -sex, -age, -sum) %>%
  as.matrix(.)

# reformat into pure numerics
CH.J2 <- matrix(as.numeric(CH.J), nrow(CH.J), ncol(CH.J))
CH.A2 <- matrix(as.numeric(CH.A), nrow(CH.A), ncol(CH.A))

# set-up
cap <- apply(CH.J2, 1, sum) # total # of juvenile captures
ind <- which(cap >=2) # which inds were recaptured; only 3321
CH.J.R <- CH.J2[ind,] # juveniles recaptured at least once
CH.J.N <- CH.J2[-ind,] # juveniles not recaptured

# remove first capture of re-captured individuals (b/c they are adults with re-capture)
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A2, CH.J.R1)
CH.A.marray <- R2ucare::marray(CH.A.m, freq = c(rowSums(CH.J2),rowSums(CH.A2))) #before the function was marray() but didn't work, then I did R2ucare::marray() but frequency is missing
# frequency: is a vector with the number of individuals having the corresponding encounter history
# I think this one is right, maybe?
# R = number of released individuals, m = m-array with upper triangly filled only and never = number of individuals never recaptured


# create CH matrix for juveniles, ignoring subsequent recaptures (those who were recaptured)
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

# create m-array for these
CH.J.R.marray <- R2ucare::marray(CH.J.R2,freq=c(rowSums(CH.J.R1),rowSums(CH.J.R2))) # marray same as above, added freq, but not sure if correct

# the last column ought to show # of juveniles not recaptured again,
# and should all be zeros, since all were released as adults
#CH.J.R.marray[,dim(CH.J2)[2]] <- 0 
# but this doesn't work
# so,
# look at dim
dim(CH.J.R.marray$m)
# here it should be 38 but max is 37: CH.J.R.marray$m[,dim(CH.J2)[2]] 
# so we do 37
CH.J.R.marray$m[,37,] <- 0
# now it works but not sure if correct

# create m-array for juvs never recaptured and add to previous array
CH.J.N.marray <- R2ucare::marray(CH.J.N, freq=c(rowSums(CH.J.N)))
CH.J.marray <- CH.J.R.marray$m + CH.J.N.marray$m # can add them only if $m included which wasn't there in the originial code

#### Load environmental data ##############

# PRECIPITATION

weather <- read.csv("Site Weather.csv")

# convert NAs to 0s in precipitation
weather1 <- weather %>% 
  mutate(Precipitation = replace(Precipitation, which(is.na(Precipitation)), 0.0))

# convert character dates to actual dates
weather1$WeathDate <- as.Date(weather1$WeathDate)

# calculate total precipitation between Oct 15 and Dec 15
rain_60 <-  weather1 %>%
  filter(WeathDate >= as.Date(paste(year(WeathDate), 10, 15, sep = "-")),
         WeathDate <= as.Date(paste(year(WeathDate), 12, 15, sep = "-"))) %>%
  group_by(BookYear) %>%
  filter(!is.na(Precipitation)) %>%
  summarize(amt_rain = sum(Precipitation))


# TEMPERATURE

# number of days with max temp > 25C per breeding season
num_days_25 <- weather1 %>% group_by(BookYear) %>%
  filter(!is.na(MaxTemp)) %>%
  summarize(n_days = n(),
            n_gt25 = sum(MaxTemp > 25),
            p_gt25 = n_gt25/n_days)

# OCEAN - MIGRATION - SSTA
ssta_migration <- read.csv("ssta.PC1.migration.11.12.21.csv")

# OCEAN - BREEDING - SSTA
ssta_breeding <- read.csv("ssta_summarized.12.2.21.csv")


## Create covariates dataframe #############
# of scaled covariates (adopted from the original code)
# for sensitivity analyses later

cov=data.frame(rain = as.vector(scale(c(36.6, rain_60$amt_rain[1:28], 36.6, rain_60$amt_rain[29:36]))),
               temp = as.vector(scale(c(0.383, num_days_25$p_gt25[1:28], 0.383, num_days_25$p_gt25[29:36]))),
               ssta_b = as.vector(ssta_breeding$av_s[1:38]),
               ssta_b_l = as.vector(c(ssta_breeding$av_s[2:38],0)),
               ssta_m = as.vector(ssta_migration$mean.2[1:38]),
               ssta_m_l = as.vector(c(ssta_migration$mean.2[2:39])))

# write.csv(cov, "cov.csv", row.names = F)

# load cov if already saved
cov=read.csv("cov.csv")

### COVARIATION #######################
# to calculate sensitivities with covariation (see below under section 3)

# rain
temp_when_rain_max=cov$temp[which(cov$rain==max(cov$rain))][1]
temp_when_rain_min=cov$temp[which(cov$rain==min(cov$rain))][1]

ssta_b_when_rain_max=cov$ssta_b[which(cov$rain==max(cov$rain))][1]
ssta_b_when_rain_min=cov$ssta_b[which(cov$rain==min(cov$rain))][1]

ssta_b_l_when_rain_max=cov$ssta_b_l[which(cov$rain==max(cov$rain))][1]
ssta_b_l_when_rain_min=cov$ssta_b_l[which(cov$rain==min(cov$rain))][1]

ssta_m_when_rain_max=cov$ssta_m[which(cov$rain==max(cov$rain))][1]
ssta_m_when_rain_min=cov$ssta_m[which(cov$rain==min(cov$rain))][1]

ssta_m_l_when_rain_max=cov$ssta_m_l[which(cov$rain==max(cov$rain))][1]
ssta_m_l_when_rain_min=cov$ssta_m_l[which(cov$rain==min(cov$rain))][1]

# temp
rain_when_temp_max=cov$rain[which(cov$temp==max(cov$temp))][1]
rain_when_temp_min=cov$rain[which(cov$temp==min(cov$temp))][1]

ssta_b_when_temp_max=cov$ssta_b[which(cov$temp==max(cov$temp))][1]
ssta_b_when_temp_min=cov$ssta_b[which(cov$temp==min(cov$temp))][1]

ssta_b_l_when_temp_max=cov$ssta_b_l[which(cov$temp==max(cov$temp))][1]
ssta_b_l_when_temp_min=cov$ssta_b_l[which(cov$temp==min(cov$temp))][1]

ssta_m_when_temp_max=cov$ssta_m[which(cov$temp==max(cov$temp))][1]
ssta_m_when_temp_min=cov$ssta_m[which(cov$temp==min(cov$temp))][1]

ssta_m_l_when_temp_max=cov$ssta_m_l[which(cov$temp==max(cov$temp))][1]
ssta_m_l_when_temp_min=cov$ssta_m_l[which(cov$temp==min(cov$temp))][1]

# ssta_b
rain_when_sstab_max=cov$rain[which(cov$ssta_b==max(cov$ssta_b))][1]
rain_when_sstab_min=cov$rain[which(cov$ssta_b==min(cov$ssta_b))][1]

temp_when_sstab_max=cov$temp[which(cov$ssta_b==max(cov$ssta_b))][1]
temp_when_sstab_min=cov$temp[which(cov$ssta_b==min(cov$ssta_b))][1]

ssta_m_when_sstab_max=cov$ssta_m[which(cov$ssta_b==max(cov$ssta_b))][1]
ssta_m_when_sstab_min=cov$ssta_m[which(cov$ssta_b==min(cov$ssta_b))][1]

ssta_m_l_when_sstab_max=cov$ssta_m_l[which(cov$ssta_b==max(cov$ssta_b))][1]
ssta_m_l_when_sstab_min=cov$ssta_m_l[which(cov$ssta_b==min(cov$ssta_b))][1]



# 1) Sensitivity Analysis ###############################

# We want to calculate the sensitivities of lambda to climatic drivers rain and temperature
# We use the method adopted from Morris et al. 2020 where they calculated the "scaled" sensitivities
# which is calculated like that: Sens_to_precip = abs((lambda_max - lambda_min) / ((max.precip - min.precip)/sd.precip))
# for that we need to get lambda_max and lambda_min, where the climatic driver in question (in this example precip) is first set to its maximum value and then lambda is calculated, and then set to its minimum value and then get lambda
# But we want 2 types of sensitivities, once without covariation and once with

# WITHOUT COVARIATION:
# means that the other climatic driver (in this example temperature) is set to its mean value

# WITH COVARIATION:
# means that the other climatic driver is set to its observed value when precip was at its maximum (then get lambda_max) and when it was at its minimum (to get lambda_min)

# for this analysis see script "SensAnalysis_MagellanicPenguins.R"


# 2) Sensitivities to Rain ##################

# Now we do the same thing but this time only perturbing the covariates in one vital rate at a time.
# The vital rates are: fecundity, juvenile survival, and adult survival.

# Perturb rain in fecundity
# Rain is only a covariate in the fecundity model and not in the others

#### 2.1 Fecundity #########

pert.fec <- nimbleCode({ # function of IPM but fecundity is perturbed
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,) 
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) 
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) 
  
  # mean demographic parameters 
  l.mphij ~ dnorm(0, 0.34)
  l.mphia ~ dnorm(0, 0.34)
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) 
  
  # covariates
  for (i in 1:11){alpha[i] ~ dnorm(0, 0.34)}
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    # fecundity
    log(fec[t]) <-  l.mfec + (alpha[1] * pert.rain[t]) + (alpha[2] * pert.temp[t]) + (alpha[3] * pert.ssta_b[t]) + (alpha[4] * pert.ssta_m[t])
    # juvenile survival
    logit(phi.juv[t]) <- l.mphij + (alpha[5] * ssta_b[t]) + (alpha[6] * ssta_m_l[t])
    # adult survival
    logit(phi.ad[t]) <-  l.mphia + (alpha[7] * temp[t]) + (alpha[8] * ssta_b[t]) + (alpha[9] * ssta_m_l[t])
    
    log(omega[t]) <-  l.mim + (alpha[10] * ssta_b[t]) + (alpha[11] * ssta_m_l[t]) #mixed
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  mphij <- exp(l.mphij) / (1+exp(l.mphij))
  mphia <- exp(l.mphia) / (1+exp(l.mphia))
  mfec <- exp(l.mfec)
  mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears,] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears,] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

# my data
CH.A.marray <- CH.A.marray$m
my.data <- list(marr.j = CH.J.marray,
                marr.a = CH.A.marray,
                rel.j = rowSums(CH.J.marray),
                rel.a = rowSums(CH.A.marray),
                y = y,
                J = J,#jz,
                R = R)#R # before it was rz but that didn't work

# inits
initial.values <- function(){list(l.mphij = rnorm(1, 0.2, 0.5),
                                  l.mphia = rnorm(1, 0.2, 0.5),
                                  l.mfec = rnorm(1, 0.2, 0.5),
                                  l.mim = rnorm(1, 0.2, 0.5),
                                  l.pj = rnorm(1, 0.2, 1),
                                  l.pa = rnorm(1, 0.2, 1),
                                  alpha = rnorm(11, 0.2, 1),
                                  sigma.y = runif(1,0,10),
                                  N1 = rpois(dim(CH.J.marray)[2], 25000), 
                                  Nimm = rpois(dim(CH.J.marray)[2], 20000), 
                                  Nad = rpois(dim(CH.J.marray)[2], 240000))} 

###### Max rain #########

# First run IPM when max rain and covariation among covariates in fecundity

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_rain_max,38),
                     pert.rain = rep(max(cov$rain),38), # perturb rain
                     pert.ssta_b = rep(ssta_b_when_rain_max,38),
                     pert.ssta_m = rep(ssta_m_when_rain_max,38))

# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                            constants = my.constants,
                            data = my.data,
                            inits = initial.values,
                            monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                         "Ntot", "sigma.y", "alpha",
                                         "geomean.lambda"),
                            niter = 300000,
                            nburnin = 100000,
                            nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxRainFec_lambda = mean_lambda # these are the lambdas

###### Min rain #######

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_rain_min,38),
                     pert.rain = rep(min(cov$rain),38), # perturb rain
                     pert.ssta_b = rep(ssta_b_when_rain_min,38),
                     pert.ssta_m = rep(ssta_m_when_rain_min,38))

# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                       constants = my.constants,
                       data = my.data,
                       inits = initial.values,
                       monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                    "Ntot", "sigma.y", "alpha",
                                    "geomean.lambda"),
                       niter = 300000,
                       nburnin = 100000,
                       nchains = 1) # one chain is enough, right? before it was 3


# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinRainFec_lambda = mean_lambda # these are the lambdas

## Sensitivities ###########

# now calculate the sensitivities of fecundity to rain

SensFec.Rain=NULL

for(i in 1:100){
  SensFec.Rain[i] = abs((MaxRainFec_lambda[i] - MinRainFec_lambda[i])/(max(cov$rain)-min(cov$rain)/1))
}

hist(SensFec.Rain)




# 3) Sensitivities to Temperature ##################

# Now we do the same thing but this time only perturbing the covariates in one vital rate at a time.
# The vital rates are: fecundity, juvenile survival, and adult survival.



#### 3.1 Fecundity #################

# Perturb temperature in fecundity

###### Max temp #########

# First run IPM when max temp and covariation among covariates in fecundity (and then later when min temp with covariation... you get the idea)

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(max(cov$temp),38), # perturb temp
                     pert.rain = rep(rain_when_temp_max,38),
                     pert.ssta_b = rep(ssta_b_when_temp_max,38),
                     pert.ssta_m = rep(ssta_m_when_temp_max,38))

# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                        "Ntot", "sigma.y", "alpha",
                                        "geomean.lambda"),
                           niter = 300000,
                           nburnin = 100000,
                           nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxTempFec_lambda = mean_lambda # these are the lambdas

##### Min temp ############

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(min(cov$temp),38), # perturb temp
                     pert.rain = rep(rain_when_temp_min,38),
                     pert.ssta_b = rep(ssta_b_when_temp_min,38),
                     pert.ssta_m = rep(ssta_m_when_temp_min,38))

# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                        "Ntot", "sigma.y", "alpha",
                                        "geomean.lambda"),
                           niter = 300000,
                           nburnin = 100000,
                           nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinTempFec_lambda = mean_lambda # these are the lambdas

## Sensitivities ###########

# now calculate the sensitivities of fecundity to temperature

SensFec.Temp=NULL

for(i in 1:100){
  SensFec.Temp[i] = abs((MaxTempFec_lambda[i] - MinTempFec_lambda[i])/(max(cov$temp)-min(cov$temp)/1))
}

hist(SensFec.Temp)



#### 3.2 Adult Survival #########


pert.survA <- nimbleCode({ # function of IPM but adult survival is perturbed
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,) 
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) 
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) 
  
  # mean demographic parameters 
  l.mphij ~ dnorm(0, 0.34)
  l.mphia ~ dnorm(0, 0.34)
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) 
  
  # covariates
  for (i in 1:11){alpha[i] ~ dnorm(0, 0.34)}
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    # fecundity
    log(fec[t]) <-  l.mfec + (alpha[1] * rain[t]) + (alpha[2] * temp[t]) + (alpha[3] * ssta_b[t]) + (alpha[4] * ssta_m[t])
    # juvenile survival
    logit(phi.juv[t]) <- l.mphij + (alpha[5] * ssta_b[t]) + (alpha[6] * ssta_m_l[t])
    # adult survival
    logit(phi.ad[t]) <-  l.mphia + (alpha[7] * pert.temp[t]) + (alpha[8] * pert.ssta_b[t]) + (alpha[9] * pert.ssta_m_l[t])
    
    log(omega[t]) <-  l.mim + (alpha[10] * ssta_b[t]) + (alpha[11] * ssta_m_l[t]) #mixed
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  mphij <- exp(l.mphij) / (1+exp(l.mphij))
  mphia <- exp(l.mphia) / (1+exp(l.mphia))
  mfec <- exp(l.mfec)
  mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears,] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears,] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

###### Max temp #########

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(max(cov$temp),38), # perturb temp
                     pert.ssta_b = rep(ssta_b_when_temp_max,38),
                     pert.ssta_m_l = rep(ssta_m_l_when_temp_max, 38)
                    )


# run code!
pert.survA.out <- nimbleMCMC(code = pert.survA,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                        "Ntot", "sigma.y", "alpha",
                                        "geomean.lambda"),
                           niter = 300000,
                           nburnin = 100000,
                           nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survA.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxTempSurvA_lambda = mean_lambda # these are the lambdas

###### Min temp #######

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(min(cov$temp),38), # perturb temp
                     pert.ssta_b = rep(ssta_b_when_temp_min,38),
                     pert.ssta_m_l = rep(ssta_m_l_when_temp_min, 38)
)


# run code!
pert.survA.out <- nimbleMCMC(code = pert.survA,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                          "Ntot", "sigma.y", "alpha",
                                          "geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 1) # one chain is enough, right? before it was 3


# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survA.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinTempSurvA_lambda = mean_lambda # these are the lambdas

## Sensitivities ###########

# now calculate the sensitivities of adult survival to temp

SensA.Rain=NULL

for(i in 1:100){
  SensA.Rain[i] = abs((MaxTempSurvA_lambda[i] - MinTempSurvA_lambda[i])/(max(cov$temp)-min(cov$temp)/1))
}

hist(SensA.Rain)


# 4) Sensitivities to SSTA_b ##################

#### 4.1 Fecundity #################

# Perturb SSTA_b in fecundity

###### Max ssta_b #########

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_sstab_max,38),
                     pert.rain = rep(rain_when_sstab_max,38),
                     pert.ssta_b = rep(max(cov$ssta_b),38),
                     pert.ssta_m = rep(ssta_m_when_sstab_max,38))

# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                        "Ntot", "sigma.y", "alpha",
                                        "geomean.lambda"),
                           niter = 300000,
                           nburnin = 100000,
                           nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxSSTABFec_lambda = mean_lambda # these are the lambdas

##### Min ssta_b ############

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_sstab_min,38),
                     pert.rain = rep(rain_when_sstab_min,38),
                     pert.ssta_b = rep(min(cov$ssta_b),38),
                     pert.ssta_m = rep(ssta_m_when_sstab_min,38))


# run code!
pert.fec.out <- nimbleMCMC(code = pert.fec,
                           constants = my.constants,
                           data = my.data,
                           inits = initial.values,
                           monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                        "Ntot", "sigma.y", "alpha",
                                        "geomean.lambda"),
                           niter = 300000,
                           nburnin = 100000,
                           nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.fec.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinSSTABFec_lambda = mean_lambda # these are the lambdas

## Sensitivities ###########

# now calculate the sensitivities of fecundity to ssta_b

SensFec.SSTAB=NULL

for(i in 1:100){
  SensFec.SSTAB[i] = abs((MaxSSTABFec_lambda[i] - MinSSTABFec_lambda[i])/(max(cov$ssta_b)-min(cov$ssta_b)/1))
}

hist(SensFec.SSTAB)

#### 4.2 Juvenile Survival #########


pert.survJ <- nimbleCode({ # function of IPM but juvenile survival is perturbed
  
  ####
  # 1. Define priors for the parameters - very vague priors for everything...
  ####
  
  # initial population sizes
  N1[1] ~ T(dnorm(25000, 4e-08),0,) 
  Nad[1] ~ T(dnorm(240000, 5e-08),0,) 
  Nimm[1] ~ T(dnorm(20000, 4e-08),0,) 
  
  # mean demographic parameters 
  l.mphij ~ dnorm(0, 0.34)
  l.mphia ~ dnorm(0, 0.34)
  l.mfec ~ dnorm(0, 0.34)
  l.mim ~ dnorm(0, 0.34)
  l.pj ~ dnorm(0, 0.34) 
  l.pa ~ dnorm(0, 0.34)
  
  # observation error
  tauy <- pow(sigma.y, -2) # (precision)
  sigma.y ~ dunif(0, 50000) # (sd) 
  
  # covariates
  for (i in 1:11){alpha[i] ~ dnorm(0, 0.34)}
  
  ####
  # 2. Constrain parameters
  ####
  
  for (t in 1:(nyears-1)){
    # fecundity
    log(fec[t]) <-  l.mfec + (alpha[1] * rain[t]) + (alpha[2] * temp[t]) + (alpha[3] * ssta_b[t]) + (alpha[4] * ssta_m[t])
    # juvenile survival
    logit(phi.juv[t]) <- l.mphij + (alpha[5] * pert.ssta_b[t]) + (alpha[6] * pert.ssta_m_l[t])
    # adult survival
    logit(phi.ad[t]) <-  l.mphia + (alpha[7] * temp[t]) + (alpha[8] * ssta_b[t]) + (alpha[9] * ssta_m_l[t])
    
    log(omega[t]) <-  l.mim + (alpha[10] * ssta_b[t]) + (alpha[11] * ssta_m_l[t]) #mixed
  }
  
  # detection constraints - juveniles
  for (t in 1:29){logit(pj[t]) <- l.pj}
  pj[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pj[t]) <- l.pj}
  
  # detection constraints - adults
  for (t in 1:29){logit(pa[t]) <- l.pa}
  pa[30] <- 0.00001 
  for (t in 31:(nyears-1)){logit(pa[t]) <- l.pa}
  
  ####
  # 3. Derived parameters
  ####
  
  mphij <- exp(l.mphij) / (1+exp(l.mphij))
  mphia <- exp(l.mphia) / (1+exp(l.mphia))
  mfec <- exp(l.mfec)
  mim <- exp(l.mim)
  mpj <- exp(l.pj) / (1+exp(l.pj))
  mpa <- exp(l.pa) / (1+exp(l.pa))
  
  # population growth rate
  for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1]/Ntot[t]
    l.lambda[t] <- log(lambda[t])
  }
  
  # geometric mean
  geomean.lambda <- exp((1/(nyears-1))*sum(l.lambda[1:(nyears-1)]))
  
  ####
  # 4. Likelihoods
  ####
  
  # 4.1 Likelihoods for population count data
  for (t in 2:nyears){
    mean1[t] <- 0.5 * fec[t-1] * phi.juv[t-1] * Nad[t-1]
    N1[t] ~ dpois(mean1[t])
    Nad[t] ~ dbin(phi.ad[t-1], Ntot[t-1])
    mpo[t] <- Ntot[t-1] * omega[t-1]
    Nimm[t] ~ dpois(mpo[t])
  }
  for (t in 1:nyears){
    Ntot[t] <- N1[t] + Nad[t] + Nimm[t]
  }
  
  # observation process
  for (t in 1:nyears){
    y[t] ~ T(dnorm(Ntot[t], tauy),0,)
  }
  
  # 4.2 Likelihoods for capture-recapture data
  # define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr.j[t,1:nyears,] ~ dmulti(pr.j[t,1:nyears], rel.j[t])
    marr.a[t,1:nyears,] ~ dmulti(pr.a[t,1:nyears], rel.a[t])
  }
  
  # m-array cell probabilities for juveniles and adults
  for (t in 1:(nyears-1)){
    qa[t] <- 1 - pa[t] # probability of recapture, adults
    qj[t] <- 1 - pj[t] # probability of recapture, juveniles
    
    # main diagonal
    pr.j[t,t] <- phi.juv[t] * pj[t]
    pr.a[t,t] <- phi.ad[t] * pa[t]
    
    # above main diagonal
    for (j in (t+1):(nyears-1)){
      pr.j[t,j] <- phi.juv[t] * prod(phi.ad[(t+1):j]) * qj[t] * prod(qa[t:(j-1)]) * pa[j]/qa[t]
      pr.a[t,j] <- prod(phi.ad[t:j]) * prod(qa[t:(j-1)]) * pa[j]
    }
    
    # below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  
  for (t in 1:(nyears-1)){
    # last column: probability of non-recapture
    pr.j[t,nyears] <- 1 - sum(pr.j[t, 1:(nyears-1)])
    pr.a[t,nyears] <- 1 - sum(pr.a[t, 1:(nyears-1)])
  }
  
  # 4.3 Likelihoods for productivity data
  for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t]*fec[t]
  }
})

###### Max ssta_b #########

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.ssta_b = rep(max(cov$ssta_b),38),
                     pert.ssta_m_l = rep(ssta_m_l_when_sstab_max,38))


# run code!
pert.survJ.out <- nimbleMCMC(code = pert.survJ,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                          "Ntot", "sigma.y", "alpha",
                                          "geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survJ.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxSSTABSurvJ_lambda = mean_lambda # these are the lambdas

###### Min ssta_b #######

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.ssta_b = rep(min(cov$ssta_b),38),
                     pert.ssta_m_l = rep(ssta_m_l_when_sstab_min,38))


# run code!
pert.survJ.out <- nimbleMCMC(code = pert.survJ,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                          "Ntot", "sigma.y", "alpha",
                                          "geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 1) # one chain is enough, right? before it was 3


# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survJ.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinSSTABSurvJ_lambda = mean_lambda # these are the lambdas


## Sensitivities ###########

# now calculate the sensitivities of juvenile survival to ssta_b

SensJ.SSTAB=NULL

for(i in 1:100){
  SensJ.SSTAB[i] = abs((MaxSSTABSurvJ_lambda[i] - MinSSTABSurvJ_lambda[i])/(max(cov$ssta_b)-min(cov$ssta_b)/1))
}

hist(SensJ.SSTAB)



#### 4.3 Adult Survival ########

###### Max ssta_b #########

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_sstab_max,38),
                     pert.ssta_b = rep(max(cov$ssta_b),38),
                     pert.ssta_m_l = rep(ssta_m_l_when_sstab_max,38)
)


# run code!
pert.survA.out <- nimbleMCMC(code = pert.survA,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                          "Ntot", "sigma.y", "alpha",
                                          "geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 1) # one chain is enough, right? before it was 3

# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survA.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MaxSSTABSurvA_lambda = mean_lambda # these are the lambdas

###### Min ssta_b #######

# my.constants
my.constants <- list(nyears = dim(CH.J.marray)[2],
                     rain = rep(0,38), 
                     temp = rep(0,38),
                     ssta_b = rep(0,38),
                     ssta_b_l = rep(0,38),
                     ssta_m = rep(0,38),
                     ssta_m_l = rep(0,38),
                     pert.temp =  rep(temp_when_sstab_min,38),
                     pert.ssta_b = rep(min(cov$ssta_b),38),
                     pert.ssta_m_l = rep(ssta_m_l_when_sstab_min,38)
)


# run code!
pert.survA.out <- nimbleMCMC(code = pert.survA,
                             constants = my.constants,
                             data = my.data,
                             inits = initial.values,
                             monitors = c("mpj", "mpa", "mphij", "mphia", "mfec", "mim",
                                          "Ntot", "sigma.y", "alpha",
                                          "geomean.lambda"),
                             niter = 300000,
                             nburnin = 100000,
                             nchains = 1) # one chain is enough, right? before it was 3


# get lambdas

# get only 100 iterations, to get the uncertainties around our estimates
coeffs=pert.survA.out[201:300,] # now we have 100 rep.

lambda=NULL
mean_lambda=NULL

for(i in 1:100){ # 100 iterations for uncertainties
  
  for(t in 1:37){ # lambda = Ntot[t+1] / Ntot[t]
    
    lambda[t] <- coeffs[i,t+1]/coeffs[i,t]
  }
  mean_lambda[i] <- mean(lambda[17:37]) # get mean but only of the last 20, treat the first bunch as "burn-ins"
  
}

MinSSTABSurvA_lambda = mean_lambda # these are the lambdas

## Sensitivities ###########

# now calculate the sensitivities of adult survival to ssta_b

SensA.SSTAB=NULL

for(i in 1:100){
  SensA.SSTAB[i] = abs((MaxSSTABSurvA_lambda[i] - MinSSTABSurvA_lambda[i])/(max(cov$ssta_b)-min(cov$ssta_b)/1))
}

hist(SensA.SSTAB)

# Save output #############
New_df=data.frame(species="Spheniscus magellanicus",
                  study.doi="10.1073/pnas.2209821120",
                  year.of.publication="2022",
                  group="Birds",
                  continent="South America",
                  driver=rep(c("Rain","SSTA_b",
                               "Rain","SSTA_b","Temperature",
                               "SSTA_b"),each=100),
                  driver.type="C",
                  stage.age=rep(c("adult","adult",
                                  "adult","adult","adult",
                                  "juvenile"),each=100),
                  vital.rates=rep(c("survival","survival",
                                    "fecundity","fecundity","fecundity",
                                    "survival"),each=100),
                  sens=c(SensA.Rain,SensA.SSTAB,
                         SensFec.Rain,SensFec.SSTAB,SensFec.Temp,
                         SensJ.SSTAB),
                  mat=2.8, # Myhrvold et al. 2015
                  n.vr=3, # number of vital rates with coviariates
                  n.pam=12, # number of total parameters of these vital rates
                  dens=0,
                  biotic_interactions=0,
                  lambda.sim=1, # lambda calculated using simulations
                  study.length=37)

write.csv(New_df,"SensVR_MPenguins.csv",row.names = F)

