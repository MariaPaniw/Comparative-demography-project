
rm(list=ls(all=TRUE))
library('popbio')
library(readr)
library(MASS)
library(boot)

setwd("/Users/maria/Dropbox/teaching/esin/Jays")

# Load covariates
env=read.csv("covariates_jays.csv")

# Load model outputs 
load("Survival_ModOuts_Jays.rdata")

surv.mod=results0[[1]]$results

load("Rep_ModOuts_Jays.rdata")

### Parameter definitions
# Srj  #Survival retained juvs
# Sdj  #Survival dispersed juvs
# Ssn  #Survival summer non- breeders
# Swn  #Survival winter non- breeders
# Ssb  #Survival summer breeders
# Swb  #Survival winter breeders
# Rsb  #Recruitment probability
# Csb  #Offspring number recruited 
# psi.rj #Transition probability from retained juv to summer breeder
# psi.dj #Transition probability from disp. juv to summer breeder
# psi.sn #Transition probability summer non-breeder to winter breeder
# psi.wn #Transition probability winter non-breeder to summer breeder

### Fitted models 
# Time vary-ing covariates = pop_d (population density), temp (temperature), snow (snow depth)
# Individual varying covariates = habitat - models predicted seperately for each habitat type (managed or unmanaged)

# Names of classes:
#' summer non-breeder (C) == sn
#' winter non-breeder(D) == wn
#' summer breeder (E) == sb
#' winter breeder (F) == wb
#' winter dispersed juvenile  (X) == dj
#' winter retained juvenile (Z) == rj 

## Survival - stage specific
# surv (S) ~ habitat + stage + snow + density + stage:density
betas_surv=surv.mod$beta[grep("S:",rownames(surv.mod$beta)),]
beta.est=betas_surv[,1]
names(beta.est)=row.names(betas_surv)
betas_surv= cbind(index=1:nrow(betas_surv),betas_surv)
coef.surv=beta.est


coef_Ssn = c(coef.surv[1],coef.surv[7],coef.surv[8],coef.surv[9]) # survival of summer non-breeder in managed forest, function of snow and pop_ad
coef_Swn = c(coef.surv[2],coef.surv[7],coef.surv[8],coef.surv[9],coef.surv[10]) # survival of winter non-breeder in managed forest, function of snow and pop_ad
coef_Ssb = c(coef.surv[3],coef.surv[7],coef.surv[8],coef.surv[9],coef.surv[11]) # survival of summer breeder in managed forest, function of snow and pop_ad
coef_Swb = c(coef.surv[4],coef.surv[7],coef.surv[8],coef.surv[9],coef.surv[12]) # survival of winter breeder  in managed forest, function of snow and pop_ad
coef_Sdj = c(coef.surv[5],coef.surv[7],coef.surv[8],coef.surv[9],coef.surv[13]) # survival of dispered juvenile  in managed forest, function of snow and pop_ad
coef_Srj = c(coef.surv[6],coef.surv[7],coef.surv[8],coef.surv[9],coef.surv[14]) # survival of retained juvenile  in managed forest, function of snow and pop_ad

vcv.Ssn=surv.mod$beta.vcv[c(1,7,8,9),c(1,7,8,9)]
vcv.Swn=surv.mod$beta.vcv[c(2,7,8,9,10),c(2,7,8,9,10)]
vcv.Ssb=surv.mod$beta.vcv[c(3,7,8,9,11),c(3,7,8,9,11)]
vcv.Swb=surv.mod$beta.vcv[c(4,7,8,9,12),c(4,7,8,9,12)]
vcv.Sdj=surv.mod$beta.vcv[c(5,7,8,9,13),c(5,7,8,9,13)]
vcv.Srj=surv.mod$beta.vcv[c(6,7,8,9,14),c(6,7,8,9,14)]

coefs.Ssn=mvrnorm(n = 100, mu = coef_Ssn, Sigma = vcv.Ssn)
coefs.Swn=mvrnorm(n = 100, mu = coef_Swn, Sigma = vcv.Swn)
coefs.Ssb=mvrnorm(n = 100, mu = coef_Ssb, Sigma = vcv.Ssb)
coefs.Swb=mvrnorm(n = 100, mu = coef_Swb, Sigma = vcv.Swb)
coefs.Sdj=mvrnorm(n = 100, mu = coef_Sdj, Sigma = vcv.Sdj)
coefs.Srj=mvrnorm(n = 100, mu = coef_Srj, Sigma = vcv.Srj)

colnames(coefs.Ssn) = c("stage_intercept","habitat","snow","pop_d")
colnames(coefs.Srj) = colnames(coefs.Sdj) = colnames(coefs.Swb) = colnames(coefs.Ssb) = colnames(coefs.Swn) = c("stage_intercept","habitat","snow","pop_d","stage_pop_d")

## Transition 
# psi = stagefromto + habitat + snow + habitat:snow
betas_tran=surv.mod$beta[grep("Psi:",rownames(surv.mod$beta)),]
beta.est=betas_tran[,1]
names(beta.est)=row.names(betas_tran)
(betas_tran=cbind(index=1:nrow(betas_tran),betas_tran))
coef_tran=beta.est

# coef_tran[5]*0 - habitat coef and not included for managed forest 
coef_psi.sn <- c(coef_tran[1], coef_tran[5], coef_tran[6], coef_tran[7]) #Transition probability summer non-breeder to winter breeder in managed forest
coef_psi.wn <- c(coef_tran[2], coef_tran[5], coef_tran[6], coef_tran[7]) #Transition probability winter non-breeder to summer breeder in managed forest
coef_psi.dj <- c(coef_tran[3], coef_tran[5], coef_tran[6], coef_tran[7]) #Transition probability from disp. juv to summer breeder in managed forest
coef_psi.rj <- c(coef_tran[4], coef_tran[5], coef_tran[6], coef_tran[7]) #Transition probability from retained juv to summer breeder in managed forest
colnames(surv.mod$beta.vcv)=rownames(surv.mod$beta.vcv) = rownames(surv.mod$beta)

vcv.psi.sn = surv.mod$beta.vcv[c(20,24:26),c(20,24:26)]
vcv.psi.wn = surv.mod$beta.vcv[c(21,24:26),c(21,24:26)]
vcv.psi.dj = surv.mod$beta.vcv[c(22,24:26),c(22,24:26)]
vcv.psi.rj = surv.mod$beta.vcv[c(23,24:26),c(23,24:26)]

coefs.psi.sn=mvrnorm(n = 100, mu = coef_psi.sn, Sigma = vcv.psi.sn)
coefs.psi.wn=mvrnorm(n = 100, mu = coef_psi.wn, Sigma = vcv.psi.wn)
coefs.psi.dj=mvrnorm(n = 100, mu = coef_psi.dj, Sigma = vcv.psi.dj)
coefs.psi.rj=mvrnorm(n = 100, mu = coef_psi.rj, Sigma = vcv.psi.rj)

## Reproduction 
# Rsb ~ habitat + pop_d + temp + habitat:temp + pop_d:temp 
# Csb ~ 1
coef.Rsb = c(fixef(mod_rep_nret)[1],fixef(mod_rep_nret)[2],fixef(mod_rep_nret)[3],fixef(mod_rep_nret)[4],fixef(mod_rep_nret)[5],fixef(mod_rep_nret)[6]) #Recruitment probability 
coefs.Rsb <- mvrnorm(n = 100, mu = coef.Rsb, Sigma = vcov(mod_rep_nret))
#
coefs.Csb <- mvrnorm(n = 100,mu = summary(mod_n_retained)$coefficients[1], Sigma = summary(mod_n_retained)$coefficients[2]) # (Intercept only model) #Offspring number recruited 

### Min, max of covariates:

min.snow=env$snow[env$snow==min(env$snow)][1]
max.snow=env$snow[env$snow==max(env$snow)][1]

sd.snow=sd(env$snow[seq(1,50,by=2)])
mean.snow=mean(env$snow[seq(1,50,by=2)])

temp.when.snow.min=env$temp[env$snow==min(env$snow)][1]
temp.when.snow.max=env$temp[env$snow==max(env$snow)][1]

# First value is managed, second is natural
pop.when.snow.min=env$pop_d[env$snow==min(env$snow)]
pop.when.snow.max=env$pop_d[env$snow==max(env$snow)]

#######

min.temp=env$temp[env$temp==min(env$temp)][1]
max.temp=env$temp[env$temp==max(env$temp)][1]

sd.temp=sd(env$temp[seq(1,50,by=2)])
mean.temp=mean(env$temp[seq(1,50,by=2)])

snow.when.temp.min=env$snow[env$temp==min(env$temp)][1]
snow.when.temp.max=env$snow[env$temp==max(env$temp)][1]

pop.when.temp.min=env$pop_d[env$temp==min(env$temp)]
pop.when.temp.max=env$pop_d[env$temp==max(env$temp)]

mean.popM=mean(env$pop_d[seq(1,50,by=2)])
mean.popN=mean(env$pop_d[seq(2,51,by=2)])

####################### For each coefficient sample, calculate lamba under different perturbations
sens.out.managed=NULL
sens.out.natural=NULL
for(i in 1:100){
  
  #----------------------------
  #MANAGED FOREST
  #----------------------------
  
  ### RETAINED JUVENILES
  
  ##min 

  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,pop.when.snow.min[1],pop.when.snow.min[1])))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,pop.when.snow.max[1],pop.when.snow.max[1])))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                     driver="Snow",
                                     driver.type="C",
                                     stage.age="retained juvenile",vital.rates="survival",
                                     sens=deltaSnow,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=1,
                                     sim=i))
  
  ### DISPERSED JUVENILES
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,pop.when.snow.min[1],pop.when.snow.min[1])))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,pop.when.snow.max[1],pop.when.snow.max[1])))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="dispersed juvenile",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### NON-BREEDING SUMMER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,pop.when.snow.min[1])))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,pop.when.snow.max[1])))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="non-breeder summer",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### NON-BREEDING WINTER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,pop.when.snow.min[1],pop.when.snow.min[1])))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,pop.when.snow.max[1],pop.when.snow.max[1])))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="non-breeder winter",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### BREEDING SUMMER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,pop.when.snow.min[1],pop.when.snow.min[1])))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,pop.when.snow.max[1],pop.when.snow.max[1])))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="breeder summer",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### BREEDING WINTER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,min.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,min.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,min.snow,pop.when.snow.min[1],pop.when.snow.min[1])))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,min.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,min.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,min.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,min.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,max.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,max.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,max.snow,pop.when.snow.max[1],pop.when.snow.max[1])))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,mean.popM,mean.temp,0,mean.popM*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,max.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,max.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,max.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,max.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="breeder winter",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### REPRODUCTION
 
  ##min 
  
  # habitat coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,mean.snow,mean.popM)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,pop.when.temp.min[1],min.temp,0,pop.when.temp.min[1]*min.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,mean.snow,0)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,mean.snow,0)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,mean.snow,0)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,mean.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,0,mean.snow,mean.popM)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,0,mean.snow,mean.popM,mean.popM)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,0,pop.when.temp.max[1],max.temp,0,pop.when.temp.max[1]*max.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,0,mean.snow,0)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,0,mean.snow,0)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,0,mean.snow,0)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,0,mean.snow,0)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaTemp=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp-min.temp)/sd.temp)
  
  sens.out.managed=rbind(sens.out.managed,data.frame(spec.driver="Breeding season temperature",
                                                     driver="Temperature",
                                                     driver.type="C",
                                                     stage.age="all",vital.rates="recruitment",
                                                     sens=deltaTemp,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  

  #----------------------------
  #NATURAL FOREST
  #----------------------------
  
  ### RETAINED JUVENILES
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,pop.when.snow.min[2],pop.when.snow.min[2])))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,pop.when.snow.max[2],pop.when.snow.max[2])))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="retained juvenile",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### DISPERSED JUVENILES
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,pop.when.snow.min[2],pop.when.snow.min[2])))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,pop.when.snow.max[2],pop.when.snow.max[2])))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="dispersed juvenile",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### NON-BREEDING SUMMER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,pop.when.snow.min[2])))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,pop.when.snow.max[2])))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="non-breeder summer",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### NON-BREEDING WINTER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,pop.when.snow.min[2],pop.when.snow.min[2])))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,pop.when.snow.max[2],pop.when.snow.max[2])))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="non-breeder winter",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### BREEDING SUMMER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,pop.when.snow.min[2],pop.when.snow.min[2])))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,pop.when.snow.max[2],pop.when.snow.max[2])))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="breeder summer",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### BREEDING WINTER
  
  ##min 
  
  # habita coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,min.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,min.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,min.snow,pop.when.snow.min[2],pop.when.snow.min[2])))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,min.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,min.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,min.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,min.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,max.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,max.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,max.snow,pop.when.snow.max[2],pop.when.snow.max[2])))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,mean.popN,mean.temp,1,mean.popN*mean.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,max.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,max.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,max.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,max.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaSnow=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.snow-min.snow)/sd.snow)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Snow depth",
                                                     driver="Snow",
                                                     driver.type="C",
                                                     stage.age="breeder winter",vital.rates="survival",
                                                     sens=deltaSnow,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  ### REPRODUCTION
  
  ##min 
  
  # habitat coef = 0 for managed for habitat
  Srj.min  = inv.logit(sum(coefs.Srj[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Sdj.min   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Ssn.min  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,mean.snow,mean.popN)))
  Swn.min   = inv.logit(sum(coefs.Swn[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Ssb.min   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Swb.min   = inv.logit(sum(coefs.Swb[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Rsb.min   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,pop.when.temp.min[2],min.temp,1,pop.when.temp.min[2]*min.temp)))
  Csb.min   = exp(sum(coefs.Csb[i,]))
  psi.rj.min = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,mean.snow,1)))
  psi.dj.min =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,mean.snow,1)))
  psi.sn.min  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,mean.snow,1)))
  psi.wn.min = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,mean.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.min*Csb.min 
  Mat[2,4] <- Rsb.min*Csb.min*0.8
  Mat[6,4] <- Ssb.min
  Mat[3,1] <- Srj.min*(1-psi.rj.min)
  Mat[4,1] <- Srj.min*psi.rj.min
  Mat[3,2] <- Sdj.min*(1-psi.dj.min)
  Mat[4,2] <- Sdj.min*psi.dj.min
  Mat[5,3] <- Ssn.min*(1-psi.sn.min)
  Mat[6,3] <- Ssn.min*psi.sn.min
  Mat[3,5] <- Swn.min*(1-psi.wn.min)
  Mat[4,5] <- Swn.min*psi.wn.min
  Mat[4,6] <- Swb.min
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.min=WinterMat%*% SummerMat
  
  ##max 
  
  # habita coef = 0 for managed for habitat
  Srj.max  = inv.logit(sum(coefs.Srj[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Sdj.max   = inv.logit(sum(coefs.Sdj[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Ssn.max  = inv.logit(sum(coefs.Ssn[i,]*c(1,1,mean.snow,mean.popN)))
  Swn.max   = inv.logit(sum(coefs.Swn[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Ssb.max   = inv.logit(sum(coefs.Ssb[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Swb.max   = inv.logit(sum(coefs.Swb[i,]*c(1,1,mean.snow,mean.popN,mean.popN)))
  Rsb.max   = inv.logit(sum(coefs.Rsb[i,]*c(1,1,pop.when.temp.max[2],max.temp,1,pop.when.temp.max[2]*max.temp)))
  Csb.max   = exp(sum(coefs.Csb[i,]))
  psi.rj.max = inv.logit(sum(coefs.psi.rj[i,]*c(1,1,mean.snow,1)))
  psi.dj.max =  inv.logit(sum(coefs.psi.dj[i,]*c(1,1,mean.snow,1)))
  psi.sn.max  = inv.logit(sum(coefs.psi.sn[i,]*c(1,1,mean.snow,1)))
  psi.wn.max = inv.logit(sum(coefs.psi.wn[i,]*c(1,1,mean.snow,1)))
  
  # ANNUAL PROJECTION MATRIX
  Mat <- matrix(0,6,6)
  Mat[1,4] <- 0.5*Rsb.max*Csb.max 
  Mat[2,4] <- Rsb.max*Csb.max*0.8
  Mat[6,4] <- Ssb.max
  Mat[3,1] <- Srj.max*(1-psi.rj.max)
  Mat[4,1] <- Srj.max*psi.rj.max
  Mat[3,2] <- Sdj.max*(1-psi.dj.max)
  Mat[4,2] <- Sdj.max*psi.dj.max
  Mat[5,3] <- Ssn.max*(1-psi.sn.max)
  Mat[6,3] <- Ssn.max*psi.sn.max
  Mat[3,5] <- Swn.max*(1-psi.wn.max)
  Mat[4,5] <- Swn.max*psi.wn.max
  Mat[4,6] <- Swb.max
  
  #ADD STAGE NAMES
  colnames(Mat) <- rownames(Mat) <- c("rj", "dj", "snb", "sb", "wnb", "wb")
  
  #SUMMER MATRIX
  SummerMat <- matrix(0,4,2)
  SummerMat[1,] <- Mat[1,3:4]
  SummerMat[2,] <- Mat[2,3:4]
  SummerMat[3,] <- Mat[5,3:4]
  SummerMat[4,] <- Mat[6,3:4]
  rownames(SummerMat) <- c("rj", "dj", "wnb", "wb")
  colnames(SummerMat) <- c("snb", "sb")
  
  #WINTER MATRIX
  WinterMat <- matrix(0,2,4)
  WinterMat[1:2,1:2] <- Mat[3:4,1:2]
  WinterMat[,3:4] <- Mat[3:4,5:6]
  colnames(WinterMat) <- c("rj", "dj", "wnb", "wb")
  rownames(WinterMat) <- c("snb", "sb")
  
  mpm.max=WinterMat%*% SummerMat
  
  deltaTemp=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp-min.temp)/sd.temp)
  
  sens.out.natural=rbind(sens.out.natural,data.frame(spec.driver="Breeding season temperature",
                                                     driver="Temperature",
                                                     driver.type="C",
                                                     stage.age="all",vital.rates="recruitment",
                                                     sens=deltaTemp,
                                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                                     cov=1,
                                                     sim=i))
  
  
  
}


sens.out.managed$population=1
sens.out.natural$population=2

### Different populations (managed/natural)

sens.out.pop=rbind(sens.out.managed,sens.out.natural)

ggplot(sens.out.pop, aes(x=sens, color=vital.rates)) +
  geom_density()+
  facet_grid(.~as.factor(population),scales = "free")

# Averaged over populations 

sens.out=aggregate(cbind(sens,l_ratio)~driver+stage.age+vital.rates+sim,mean,data=sens.out.pop)

ggplot(sens.out, aes(x=sens, color=vital.rates)) +
  geom_density()


sens_jays=data.frame(study.doi="10.1007/s00442-018-4100-z",year.of.publication=2018,
                      group="Birds",species="Perisoreus infaustus", continent="Europe",driver=sens.out$driver,driver.type="C",
                     stage.age=sens.out$stage.age,vital.rates=sens.out$vital.rates,sens=sens.out$sens,mat=1,n.vr=12,n.pam=28,dens=1,
                      biotic_interactions=0,lambda.sim=0,study.length=15)

write.csv(sens_jays,"sens_jays_vital_rates.csv",row.names = F)

sens_jays_pop=data.frame(study.doi="10.1007/s00442-018-4100-z",year.of.publication=2018,
                     group="Birds",species="Perisoreus infaustus",continent="Europe",driver=sens.out.pop$driver,driver.type="C",
                     stage.age=sens.out.pop$stage.age,vital.rates=sens.out.pop$vital.rates,sens=sens.out.pop$sens,mat=1,n.vr=12,n.pam=28,dens=1,
                     biotic_interactions=0,lambda.sim=0,study.length=15,population=sens.out.pop$population)

write.csv(sens_jays_pop,"sens_jays_pop_vital_rates.csv",row.names = F)
