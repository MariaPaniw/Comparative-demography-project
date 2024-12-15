
setwd("/Users/maria/Dropbox/teaching/esin/reindeer")

# SURVIVAL parameters
# load posterior samples
load("coeffic_lmer_withintercept.Rdata")
# "coeffic": all posterior samples for survival parameters (table S2)
# each row = one posterior sample (9,090 in total)
# each column = one parameter

#Column 1: intercept (including age 1)
age1.par=coeffic[,1]

#Column 2:  age 2
age2.par=coeffic[,2]

#Column 3: age 3
age3.par=coeffic[,3]

#Column 4: age 4
age4.par=coeffic[,4]

#Column 5: age 5
age5.par=coeffic[,5]

#Column 6: age 6
age6.par=coeffic[,6]

#Column 7: Nposthunt
Nposthunt.par=coeffic[,7]
#Column 8: year
#year.par=coeffic[,8]

#Column 9: ROS
ROS.par=coeffic[,9]
#Column 10: winter length
winterlength.par=coeffic[,10]

#Column 11: ROS effect age 2 (to be added to column 9)
ROS2.par=coeffic[,11]

#Column 12: ROS effect age 3 (to be added to column 9)
ROS3.par=coeffic[,12]

#Column 13: ROS effect age 4 (to be added to column 9)
ROS4.par=coeffic[,13]

#Column 14: ROS effect age 5 (to be added to column 9)
ROS5.par=coeffic[,14]

#Column 15: ROS effect age 6 (to be added to column 9)
ROS6.par=coeffic[,15]

#Column 16: k
k.par=coeffic[,16]

# SURVIVAL
# equation 1 in paper
# estimates are from supp tbl 2
# logit(S) = Age class 1 + N.posthunt + (ROS * exp(k*N.posthunt)) + winterlength  +  ROS'Age class 1
# Si = inv.logit(intercept+age + Nposthunt.par*dens + ROS.par*(ROS*exp(k.par*N.posthunt.par*dens)) + winterlength.par*winterlength + age*ROS.par*(ROS*exp(k.par*Nposthunt.par*dens)))

# age class 1
S1 <- function(dens,ROS,winterlength,age1.par,Nposthunt.par,ROS.par,k.par){
  s <- inv.logit(age1.par + Nposthunt.par*dens + (ROS.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age1.par*1))
  return(s)
}

#mean(S1(mean.dens,mean.ROS,mean.winterlength,age1.par,Nposthunt.par,ROS.par,k.par))

# age class 2
S2 <- function(dens,ROS,winterlength,age2.par,Nposthunt.par,ROS.par,ROS2.par,k.par){
  s <- inv.logit(age2.par*2 + Nposthunt.par*dens + (ROS.par+ROS2.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par+ROS2.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age2.par*2))
  return(s)
}

#mean(S2(mean.dens,mean.ROS,mean.winterlength,age2.par,Nposthunt.par,ROS.par,ROS2.par,k.par))



# age class 3
S3 <- function(dens,ROS,winterlength,age3.par,Nposthunt.par,ROS.par,ROS3.par,k.par){
  s <- inv.logit(age3.par*3 + Nposthunt.par*dens + (ROS.par+ROS3.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par+ROS3.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age3.par*3))
  return(s)
}

#mean(S3(mean.dens,mean.ROS,mean.winterlength,age3.par,Nposthunt.par,ROS.par,ROS3.par,k.par))

# age class 4
S4 <- function(dens,ROS,winterlength,age4.par,Nposthunt.par,ROS.par,ROS4.par,k.par){
  s <- inv.logit(age4.par*4 + Nposthunt.par*dens + (ROS.par+ROS4.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par+ROS4.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age4.par*4))
  return(s)
}

#mean(S4(mean.dens,mean.ROS,mean.winterlength,age4.par,Nposthunt.par,ROS.par,ROS4.par,k.par))

# age class 5
S5 <- function(dens,ROS,winterlength,age5.par,Nposthunt.par,ROS.par,ROS5.par,k.par){
  s <- inv.logit(age5.par*5 + Nposthunt.par*dens + (ROS.par+ROS5.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par+ROS5.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age5.par*5))
  return(s)
}

#mean(S5(mean.dens,mean.ROS,mean.winterlength,age5.par,Nposthunt.par,ROS.par,ROS5.par,k.par))

# age class 6
S6 <- function(dens,ROS,winterlength,age6.par,Nposthunt.par,ROS.par,ROS6.par,k.par){
  s <- inv.logit(age6.par*6 + Nposthunt.par*dens + (ROS.par+ROS6.par)*(ROS*exp(k.par*Nposthunt.par*dens)) + (ROS.par+ROS6.par)*(ROS*exp(k.par*Nposthunt.par*dens)*age6.par*6))
  return(s)
}

#mean(S6(mean.dens,mean.ROS,mean.winterlength,age6.par,Nposthunt.par,ROS.par,ROS6.par,k.par))



#############################################################

# FECUNDITY parameters
# load posterior samples
load("coefficF_lmer_withintercept.Rdata")

# intercept
#F.intercept=coefficF[,1]

# age class 2
F.age2.par=coefficF[,2]

# age class 3
F.age3.par=coefficF[,3]

# age class 4
F.age4.par=coefficF[,4]

# age class 5
F.age5.par=coefficF[,5]

# Nposthunt
F.Nposthunt=coefficF[,6]

# year
#F.year.par=coefficF[,7]

# ROS
F.ROS.par=coefficF[,8]

# winterlength
F.winterlength.par=coefficF[,9]

# ROS effect age 2 (to be addedto column 8)
F.ROS2.par=coefficF[,10]

# ROS effect age 3 (to be addedto column 8)
F.ROS3.par=coefficF[,11]

# ROS effect age 2 (to be addedto column 8)
F.ROS4.par=coefficF[,12]

# ROS effect age 2 (to be addedto column 8)
F.ROS5.par=coefficF[,13]

# k
F.k.par=coefficF[,14]


# FECUNDITY
# equation 2 in paper
# estimates are from supp tbl 2
# without year effect

# logit(F) = age class + N.posthunt + ROS * exp(k*N.posthunt)) + winter length + ROS'age class 

#Fi = inv.logit(F.agei.par*i + F.Nposthunt*dens + (F.ROS.par+F.ROSi.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)) + (F.ROS.par+F.ROSi.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)*F.agei.par*i))

# age class 2
F2 <- function(dens,ROS,winterlength,F.age2.par,F.Nposthunt,F.ROS.par,F.ROS2.par,F.k.par){
  f <- inv.logit(F.age2.par*2 + F.Nposthunt*dens + (F.ROS.par+F.ROS2.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)) + (F.ROS.par+F.ROS2.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)*F.age2.par*2))
  return(f)
  
}

#mean(F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par,F.Nposthunt,F.ROS.par,F.ROS2.par,F.k.par))

# age class 3
F3 <- function(dens,ROS,winterlength,F.age3.par,F.Nposthunt,F.ROS.par,F.ROS3.par,F.k.par){
  f <- inv.logit(F.age3.par*3 + F.Nposthunt*dens + (F.ROS.par+F.ROS3.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)) + (F.ROS.par+F.ROS3.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)*F.age3.par*3))
  return(f)
  
}
#mean(F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par,F.Nposthunt,F.ROS.par,F.ROS3.par,F.k.par))

# age class 4
F4 <- function(dens,ROS,winterlength,F.age4.par,F.Nposthunt,F.ROS.par,F.ROS4.par,F.k.par){
  f <- inv.logit(F.age4.par*4 + F.Nposthunt*dens + (F.ROS.par+F.ROS4.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)) + (F.ROS.par+F.ROS4.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)*F.age4.par*4))
  return(f)
  
}
#mean(F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par,F.Nposthunt,F.ROS.par,F.ROS4.par,F.k.par))

# age class 5
F5 <- function(dens,ROS,winterlength,F.age5.par,F.Nposthunt,F.ROS.par,F.ROS5.par,F.k.par){
  f <- inv.logit(F.age5.par*5 + F.Nposthunt*dens + (F.ROS.par+F.ROS5.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)) + (F.ROS.par+F.ROS5.par)*(ROS*exp(F.k.par*F.Nposthunt*dens)*F.age5.par*5))
  return(f)
  
}
#mean(F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par,F.Nposthunt,F.ROS.par,F.ROS5.par,F.k.par))

## 2) Covariates

cov <- read.csv("data_reindeer.txt",sep="")
# scaled pop density
cov$N <- (cov$N - mean(cov$N)) / sd(cov$N)

# ROS (rain-on-snow)
min.ROS=min(cov$ROS[cov$ROS>0])
max.ROS=max(cov$ROS)
sd.ROS=sd(cov$ROS)
mean.ROS=mean(cov$ROS)

# winter length
min.winterlength=min(cov$Winter.length)
max.winterlength=max(cov$Winter.length)
sd.winterlength=sd(cov$Winter.length)
mean.winterlength=mean(cov$Winter.length)

# population density (scaled)
min.dens=min(cov$N)
max.dens=max(cov$N)
sd.dens=1
mean.dens=0

# Covariation
## ROS
ROS_when_winterlength_max=cov$ROS[which(cov$Winter.length==max(cov$Winter.length))][1]
ROS_when_winterlength_min=cov$ROS[which(cov$Winter.length==min(cov$Winter.length))][1]
ROS_when_dens_max=cov$ROS[which(cov$N==max(cov$N))][1]
ROS_when_dens_min=cov$ROS[which(cov$N==min(cov$N))][1]
## winter length
winterlength_when_dens_max=cov$Winter.length[which(cov$N==max(cov$N))][1]
winterlength_when_dens_min=cov$Winter.length[which(cov$N==min(cov$N))][1]
winterlength_when_ROS_max=cov$Winter.length[which(cov$ROS==max(cov$ROS))][1]
winterlength_when_ROS_min=cov$Winter.length[which(cov$ROS==min(cov$ROS))][1]
## dens
dens_when_winterlength_max=cov$N[which(cov$Winter.length==max(cov$Winter.length))][1]
dens_when_winterlength_min=cov$N[which(cov$Winter.length==min(cov$Winter.length))][1]
dens_when_ROS_max=cov$N[which(cov$ROS==max(cov$ROS))][1]
dens_when_ROS_min=cov$N[which(cov$ROS==min(cov$ROS))][1]

lambda=NULL
for(i in 1:nrow(coeffic)){
  mpm <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                  F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                  F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                  F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                  F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                  0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  lambda[i] <- lambda(mpm)
}

Sens.ROS.NoCov=NULL
l_ratio.ROS.NoCov=NULL
which.coef=NULL
for(i in 1:nrow(coeffic)){
  # mpm when ROS max
  mpm.max <- matrix(c(0,S1(mean.dens,max.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,max.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,max.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,max.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,max.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,max.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,max.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,max.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,max.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,max.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  # mpm when ROS min
  mpm.min <- matrix(c(0,S1(mean.dens,min.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,min.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,min.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,min.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,min.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,min.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,min.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,min.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,min.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,min.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  # scaled sensitivities and log ratios
  
  if(lambda(mpm.min)>0.2){
    which.coef=c(which.coef,i)
    Sens.ROS.NoCov <- c(Sens.ROS.NoCov,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
    l_ratio.ROS.NoCov <- c(l_ratio.ROS.NoCov,abs(log(lambda(mpm.max)/lambda(mpm.min))))
  }
  

}

length(Sens.ROS.NoCov)
hist(Sens.ROS.NoCov)
mean(Sens.ROS.NoCov)
hist(l_ratio.ROS.NoCov)


Sens.ROS.Cov=NULL
l_ratio.ROS.Cov=NULL

for(i in which.coef){
  # mpm when ROS max
  mpm.max <- matrix(c(0,S1(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  # mpm when ROS min
  mpm.min <- matrix(c(0,S1(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  # scaled sensitivities
  Sens.ROS.Cov <- c(Sens.ROS.Cov,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
  l_ratio.ROS.Cov <- c(l_ratio.ROS.Cov,abs(log(lambda(mpm.max)/lambda(mpm.min))))
}

hist(Sens.ROS.Cov)
mean(Sens.ROS.Cov)
hist(l_ratio.ROS.Cov)

ScaledSens_Reindeer=data.frame(species="Rangifer tarandus",
                               study.doi="10.1038/s41467-019-09332-5",
                               year.of.publication="2019",
                               group="Mammals",
                               continent="Europe",
                               driver="ROS",
                               driver.type="C",
                               stage.age="all",
                               vital.rates="all",
                               sens=c(Sens.ROS.NoCov[1:100],Sens.ROS.Cov[1:100]),
                               cov=rep(c(0,1),each=100),
                               mat=1.9, # age at sexual maturity Source: Myhrvold et al. 2015
                               n.vr=10,
                               n.pam=49,
                               dens=1,
                               biotic_interactions=0,
                               lambda.sim=0,
                               study.length=20,
                               l_ratio=c(l_ratio.ROS.NoCov[1:100],l_ratio.ROS.Cov[1:100]))

# Change density plot line colors by groups
ggplot(ScaledSens_Reindeer, aes(x=sens, color=factor(cov))) +
  geom_density()

write.csv(ScaledSens_Reindeer,"SensReindeer.csv",row.names = F)

#### VITAL RATE:

#### 6.1 Survival for non-reproductive stages
## 6) ROS
# stages: 1-6
# reproductive stages: 2-5
# non-reproductive stages: 1 & 6

#### 6.1 Survival for non-reproductive stages

# stage 1
ROS.S1=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S1 <- c(ROS.S1,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
  
}

# stage 6
ROS.S6=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S6 <- c(ROS.S6,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
  
}



### 6.2 Survival for reproductive stages

# stage 2
ROS.S2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S2 <- c(ROS.S2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}

# stage 3
ROS.S3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S3 <- c(ROS.S3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}


# stage 4
ROS.S4=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S4 <- c(ROS.S4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}


# stage 5
ROS.S5=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.S5 <- c(ROS.S5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}



### 6.3 Fecundity

# stage 2
ROS.F2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.F2 <- c(ROS.F2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}

# stage 3
ROS.F3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.F3 <- c(ROS.F3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}


# stage 4
ROS.F4=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.F4 <- c(ROS.F4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}


# stage 5
ROS.F5=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_ROS_max,max.ROS,winterlength_when_ROS_max,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_ROS_min,min.ROS,winterlength_when_ROS_min,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  ROS.F5 <- c(ROS.F5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.ROS-min.ROS)/sd.ROS)))
}






## 7) Density

### 7.1 Survival for non-reproductive stages

# stage 1
DENS.S1=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S1 <- c(DENS.S1,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))
  
}

# stage 6
DENS.S6=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S6 <- c(DENS.S6,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))   
}



### 6.2 Survival for reproductive stages

# stage 2
DENS.S2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S2 <- c(DENS.S2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}

# stage 3
DENS.S3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S3 <- c(DENS.S3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))
}


# stage 4
DENS.S4=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S4 <- c(DENS.S4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}


# stage 5
DENS.S5=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(max.dens,ROS_when_dens_max,winterlength_when_dens_max,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(min.dens,ROS_when_dens_min,winterlength_when_dens_min,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.S5 <- c(DENS.S5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}



### 6.3 Fecundity

# stage 2
DENS.F2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(max.dens,ROS_when_dens_max,winterlength_when_dens_max,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(min.dens,ROS_when_dens_min,winterlength_when_dens_min,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.F2 <- c(DENS.F2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))
}

# stage 3
DENS.F3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(max.dens,ROS_when_dens_max,winterlength_when_dens_max,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(min.dens,ROS_when_dens_min,winterlength_when_dens_min,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.F3 <- c(DENS.F3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}


# stage 4
DENS.F4=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(max.dens,ROS_when_dens_max,winterlength_when_dens_max,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(min.dens,ROS_when_dens_min,winterlength_when_dens_min,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.F4 <- c(DENS.F4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}


# stage 5
DENS.F5=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(max.dens,ROS_when_dens_max,winterlength_when_dens_max,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(min.dens,ROS_when_dens_min,winterlength_when_dens_min,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  DENS.F5 <- c(DENS.F5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.dens-min.dens)/sd.dens)))}






## 8) Winterlength

### 7.1 Survival for non-reproductive stages

# stage 1
WL.S1=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S1 <- c(WL.S1,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
  
}

# stage 6
WL.S6=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S6 <- c(WL.S6,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
  
}



### 6.2 Survival for reproductive stages

# stage 2
WL.S2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S2 <- c(WL.S2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}

# stage 3
WL.S3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S3 <- c(WL.S3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}


# stage 4
WL.S4=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S4 <- c(WL.S4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}


# stage 5
WL.S5=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.S5 <- c(WL.S5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}



### 6.3 Fecundity

# stage 2
WL.F2=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.F2 <- c(WL.F2,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}

# stage 3
WL.F3=NULL
for(i in which.coef){
  
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.F3 <- c(WL.F3,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}


# stage 4
WL.F4=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(mean.dens,mean.ROS,mean.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.F4 <- c(WL.F4,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}


# stage 5
WL.F5=NULL
for(i in which.coef){
  mpm.max <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_winterlength_max,ROS_when_winterlength_max,max.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  
  mpm.min <- matrix(c(0,S1(mean.dens,mean.ROS,mean.winterlength,age1.par[i],Nposthunt.par[i],ROS.par[i],k.par[i]),0,0,0,0,
                      F2(mean.dens,mean.ROS,mean.winterlength,F.age2.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS2.par[i],F.k.par[i]),0,S2(mean.dens,mean.ROS,mean.winterlength,age2.par[i],Nposthunt.par[i],ROS.par[i],ROS2.par[i],k.par[i]),0,0,0,
                      F3(mean.dens,mean.ROS,mean.winterlength,F.age3.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS3.par[i],F.k.par[i]),0,0,S3(mean.dens,mean.ROS,mean.winterlength,age3.par[i],Nposthunt.par[i],ROS.par[i],ROS3.par[i],k.par[i]),0,0,
                      F4(mean.dens,mean.ROS,mean.winterlength,F.age4.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS4.par[i],F.k.par[i]),0,0,0,S4(mean.dens,mean.ROS,mean.winterlength,age4.par[i],Nposthunt.par[i],ROS.par[i],ROS4.par[i],k.par[i]),0,
                      F5(dens_when_winterlength_min,ROS_when_winterlength_min,min.winterlength,F.age5.par[i],F.Nposthunt[i],F.ROS.par[i],F.ROS5.par[i],F.k.par[i]),0,0,0,0,S5(mean.dens,mean.ROS,mean.winterlength,age5.par[i],Nposthunt.par[i],ROS.par[i],ROS5.par[i],k.par[i]),
                      0,0,0,0,0,S6(mean.dens,mean.ROS,mean.winterlength,age6.par[i],Nposthunt.par[i],ROS.par[i],ROS6.par[i],k.par[i])),6,6)
  
  WL.F5 <- c(WL.F5,abs((lambda(mpm.max)-lambda(mpm.min))/((max.winterlength-min.winterlength)/sd.winterlength)))
}



## 9) Save Output

Surv=data.frame(species="Rangifer tarandus",
                study.doi="10.1038/s41467-019-09332-5",
                year.of.publication="2019",
                group="Mammals",
                continent="Europe",
                driver=rep(c("ROS","density","winterlength"),each=600),
                driver.type=rep(c("C","D","C"),each=600),
                stage.age=rep(c(1,2,3,4,5,6),each=100),
                vital.rates="survival",
                sens=c(ROS.S1[1:100],ROS.S2[1:100],ROS.S3[1:100],ROS.S4[1:100],ROS.S5[1:100],ROS.S6[1:100],
                       DENS.S1[1:100],DENS.S2[1:100],DENS.S3[1:100],DENS.S4[1:100],DENS.S5[1:100],DENS.S6[1:100],
                       WL.S1[1:100],WL.S2[1:100],WL.S3[1:100],WL.S4[1:100],WL.S5[1:100],WL.S6[1:100]),
                mat=1.9, # age at sexual maturity Source: Myhrvold et al. 2015
                n.vr=10,
                n.pam=49,
                dens=1,
                biotic_interactions=0,
                lambda.sim=0,
                study.length=20)


Fec=data.frame(species="Rangifer tarandus",
               study.doi="10.1038/s41467-019-09332-5",
               year.of.publication="2019",
               group="Mammals",
               continent="Europe",
               driver=rep(c("ROS","density","winterlength"),each=400),
               driver.type=rep(c("C","D","C"),each=400),
               stage.age=rep(c(2,3,4,5),each=100),
               vital.rates="fecundity",
               sens=c(ROS.F2[1:100],ROS.F3[1:100],ROS.F4[1:100],ROS.F5[1:100],
                      DENS.F2[1:100],DENS.F3[1:100],DENS.F4[1:100],DENS.F5[1:100],
                      WL.F2[1:100],WL.F3[1:100],WL.F4[1:100],WL.F5[1:100]),
               mat=1.9, # age at sexual maturity Source: Myhrvold et al. 2015
               n.vr=10,
               n.pam=49,
               dens=1,
               biotic_interactions=0,
               lambda.sim=0,
               study.length=20)


SensVR=rbind(Surv,Fec)

write.csv(SensVR,"SensReindeer_VR.csv",row.names = F)



