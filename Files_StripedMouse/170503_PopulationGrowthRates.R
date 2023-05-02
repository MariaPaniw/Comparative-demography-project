library(lme4)
library(ggplot2)
library(popbio)
library(plyr)
library(gridExtra)

# Load workspace containing fitted models
load('170503_allmodels_final.RData')

#################################
#### CMR PREDICTION FORMULAS ####
#################################

# SURVIVAL

# Provides the y value (linear predictor) of S, where logit(S) = y and y = beta1*x1 + beta2*x2
S.y = function(stage,temp,food,dens_F){
  
  X <- c(1, as.numeric(stage=='I'), as.numeric(stage=='P'), temp, food, dens_F,
         (stage=='I')*temp, (stage=='P')*temp,
         (stage=='I')*food, (stage=='P')*food,
         (stage=='I')*dens_F, (stage=='P')*dens_F,
         temp*food, food*dens_F,
         (stage=='I')*temp*food, (stage=='P')*temp*food,
         (stage=='I')*food*dens_F, (stage=='P')*food*dens_F)
  
  X %*% est$estimate[1:18]
}

# Transforms the y value with the correct link function
S.predict = function(stage,temp,food,dens_F){
  return(plogis(S.y(stage,temp,food,dens_F)))
}

# MATURATION
Psi.y = function(stage,temp,food,dens_F){
  X <- c(as.numeric(stage=='I'), as.numeric(stage=='P'),
         temp, food, dens_F,
         food*dens_F, temp*dens_F, temp*food)
  X %*% est$estimate[57:64]
}

Psi.predict = function(stage,temp,food,dens_F){
  return(1/(1+2*exp(-Psi.y(stage,temp,food,dens_F))))    
}


######################################
#### Preparing Environmental Data ####
######################################

# Already in the workspace as dataframe "env"

# Extract scaling parameters for population density
m_dens <- 19.22219
sd_dens <- 8.199321

# Set "baseline" area
area.fix <- 19343.53/10000 # 1.9434353

# Add "lagged" environmental covariates (t-1)
env$temp_t <- NA
env$food_t <- NA
env$dens_F_t <- NA

for(i in 2:nrow(env)){
  env$temp_t[i] <- env$temp[i-1]
  env$food_t[i] <- env$food[i-1]
  env$dens_F_t[i] <- env$dens_F[i-1]
}

############################################
#### Calculating 'Real' Population Size ####
############################################

# Reading in population counts
pop.count <- read.csv('170110_allfemales_stage.csv')
pop.count <- subset(pop.count, !(session%in%c(6,8,16,17,18,38,39)))

pop.count$total <- pop.count$breeders + pop.count$philopatrics + pop.count$immatures + pop.count$unknown

# Correcting observed numbers with recapture probabilities
pop.count <- merge(pop.count, recap, by = 'session', all.x = TRUE)

pop.count$breeders_corr <- pop.count$breeders / pop.count$pB
pop.count$philopatrics_corr <- pop.count$philopatrics / pop.count$pP
pop.count$immatures_corr <- pop.count$immatures / pop.count$pP

pop.count$total_corr <- pop.count$breeders_corr + pop.count$philopatrics_corr + pop.count$immatures_corr

# Reading in data on area size over time
data <- read.csv('161222_soc_fac_ad.csv')
area.t <- data$area_h_corr

##############################################
#### Building Matrix Model for Projection ####
##############################################

my.matrix = function(temp, temp_t, food, food_t, dens_F, dens_F_t){
  
  Si <- S.predict('I',temp,food,dens_F)
  Sp <- S.predict('P',temp,food,dens_F)
  Sb <- S.predict('B',temp,food,dens_F)
  
  PsiIB <- Psi.predict('I',temp,food,dens_F)
  PsiPB <- Psi.predict('P',temp,food,dens_F)
  
  newd <- data.frame(food_t = food_t, dens_F_t = dens_F_t, temp_t = temp_t) 
  
  Bp <- unname(predict(bpmodF, newdata = newd, re.form=NA, type='response'))
  
  Lp <- predict(lpmodF, type = "response", newdata = newd)
  LS <- predict(lsmod_Wlag, type = "response", newdata = newd)
  #LS <- 2.65
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}

timedens.matrix = function(t,dens){
  
  A <- my.matrix(env$temp[t],env$temp_t[t],
                 env$food[t],env$food_t[t],
                 dens[1],dens[2])
  return(A)
}

##########################################
#### CALCULATE DIFFERENT GROWTH RATES ####
##########################################

growth.rates = function(t){
  
  # Obtain the initial population vector
  N <- c(pop.count$immatures_corr[t],pop.count$philopatrics_corr[t],pop.count$breeders_corr[t])
  
  Nvec <- matrix(NA, 12000, 3)
  Nvec[1,] = c(pop.count$immatures_corr[t-1],pop.count$philopatrics_corr[t-1],pop.count$breeders_corr[t-1])/area.t[t-1]
  Nvec[2,] = c(pop.count$immatures_corr[t],pop.count$philopatrics_corr[t],pop.count$breeders_corr[t])/area.t[t]
  
  # Calculate the matrix
  A <- my.matrix(env$temp[t],env$temp_t[t],
               env$food[t],env$food_t[t],
               env$dens_F[t],env$dens_F_t[t])
  
  # Calculate 1-timestep transient growth rate
  T.lambda <- sum(A%*%N)/sum(N)
  
  # Calculate asymptotic growth rate (dominant right eigenvalue)
  A.lambda <- lambda(A)
  
  # Calculate pseudo-asymptotic growth rate (accounting for density feedback)
  #for(x in 2:11999){
  #  
  #  # Calculate corrected sum of philopatrics and breeders
  #  total.no <- Nvec[x,2] + Nvec[x,3]
  #  total.no_t <- Nvec[x-1,2] + Nvec[x-1,3]
  #  
  #  # Scale density
  #  scale.dens <- (total.no - m_dens) / sd_dens
  #  scale.dens_t <- (total.no_t - m_dens) / sd_dens
  #  
  #  dens <- c(scale.dens, scale.dens_t)
  #  
  #  # 4) Feed density into new matrix and project
  #  Nvec[x+1,] = timedens.matrix(t,dens) %*% Nvec[x,]
  #  
  #}
  
  #psA.lambda <- sum(Nvec[12000,])/sum(N)
  psA.lambda <- NA
  
  return(data.frame(time = t, A.lambda = A.lambda, T.lambda = T.lambda, psA.lambda = psA.lambda))
}

# Re-adding missing sessions in pop.count
sessions <- data.frame(session = c(1:117))
pop.count <- merge(pop.count, sessions, by = 'session', all.y = TRUE)

GR <- do.call("rbind", sapply(2:nrow(env), FUN = function(t) growth.rates(t), simplify = FALSE))


