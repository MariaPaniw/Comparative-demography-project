
# AIM: Build IPM for Lavandula stoechas from DNP


# Teresa Sánchez-Mejía 
# Last updated 11/11/2024


library(tidyverse) # gather / %>%
library(lme4) #MLM
library(mgcv)
library(gridExtra)
library(popbio)

setwd("~/Lavandula/")


# A. LOAD DATA AND GAM ####

load("Lavandula_data.RData")
load("Lavandula_GAM.RData")

# Lavandula data
# sajL  # survival continuous stage
# faL   # flowering continuous stage (reproductive adults)
# nfaL  # number of inflorescences continuous stage (reproductive adults)
# gajL  # growth continuous stage
# ssL   # survival discrete stage
# rsL   # recruitment (number seedlings/number reproductive adults)
# ss2L  # size seedling t+1

# B. DEFINE VITAL RATES FUNCTIONS #####

L_saj_cr = function(Size, Site,Year){ 
  
  beta <- coef(MFL_sajcr)
  Vb   <- vcov(MFL_sajcr)
  ilink <- family(MFL_sajcr)$linkinv
  Xp = predict(MFL_sajcr,                                               
               newdata=data.frame(Size=Size,Site=Site,Year=Year),type="lpmatrix")
  mrand <- mvrnorm(1, beta, Vb)
  
  p = as.numeric(ilink(Xp %*% mrand))
  return(p)
} # CONTINUOUS - Survival


L_fa_cr  = function(Size, P_t1){                                    
  
  beta <- coef(MFL_facr)
  Vb   <- vcov(MFL_facr)
  ilink <- family(MFL_facr)$linkinv
  
  Xp = predict(MFL_facr,                                               
               newdata=data.frame(Size=Size, P_t1= P_t1),type="lpmatrix")
  
  mrand <- mvrnorm(1, beta, Vb)
  
  TMe = as.numeric(ilink(Xp %*% mrand))
  
  return(TMe)
}   # CONTINUOUS - Flowering

L_gaj_cr = function(Size1, Size, Site, Year) {                                   
  
  beta <- coef(MFL_gajcr)
  Vb   <- vcov(MFL_gajcr)
  ilink <- family(MFL_gajcr)$linkinv
  
  Xp = predict(MFL_gajcr,                                               
               newdata=data.frame(Size=Size, Site=Site, Year=Year),type="lpmatrix")
  
  mrand <- mvrnorm(1, beta, Vb)
  
  mean = as.numeric(ilink(Xp %*% mrand))
  
  sd=sqrt(MFL_gajcr$sig2)
 
  return(dnorm(Size1,mean=mean,sd=sd))
}                   # CONTINUOUS - Growth    

L_rs_cr  = function(TMe,lintrad){                                        
  
  beta <- coef(MFL_rscr)
  Vb   <- vcov(MFL_rscr)
  ilink <- family(MFL_rscr)$linkinv
  
  Xp = predict(MFL_rscr,                                               
               newdata=data.frame(TMe=TMe, lintrad=lintrad),type="lpmatrix")
  
  mrand <- mvrnorm(1, beta, Vb)
  
  p = as.numeric(ilink(Xp %*% mrand))

  return(p)


  }                                   # DISCRETE - Recruitment

L_ss_cr  = function(TMaS_t1, lintrad) {                                           
  
  beta <- coef(MFL_sscr)
  Vb   <- vcov(MFL_sscr)
  ilink <- family(MFL_sscr)$linkinv
  
  Xp = predict(MFL_sscr,                                               
               newdata=data.frame(TMaS_t1=TMaS_t1, lintrad=lintrad),type="lpmatrix")
  
  mrand <- mvrnorm(1, beta, Vb)
  
  p = as.numeric(ilink(Xp %*% mrand))

  return(p)
}                                  # DISCRETE - Survival

L_ss2_cr = function(Size, TMaS_t1, Site){                                         
  
  beta <- coef(MFL_ss2cr)
  Vb   <- vcov(MFL_ss2cr)
  ilink <- family(MFL_ss2cr)$linkinv
  
  Xp = predict(MFL_ss2cr,                                               
               newdata=data.frame(TMaS_t1=TMaS_t1, Site=Site),type="lpmatrix")
  
  mrand <- mvrnorm(1, beta, Vb)
  
  mean = as.numeric(ilink(Xp %*% mrand))
  
  sd=sqrt(MFL_ss2cr$sig2)
  
  return(dnorm(Size,mean=mean,sd=sd))
}                          # DISCRETE - Seedling size 

# C. DEFINE KERNEL #####
# Define mid-point rule function and number of meshpoints

calc_midp = function(n,L,U){
  b = L+c(0:n)*(U-L)/n                                                          # Interval that each cell of the matrix covers 
  meshp = 0.5*(b[1:n]+b[2:(n+1)])                                               # Midpoint
  h=(U-L)/n
  return(list(meshp=meshp, h=h))
}

# Define Kernel

K.fnc_L_cr = function(midp, TMe, TMaS_t1,P_t1,lintrad, Site, Year) {
  
  K = array(0,c(n+1,n+1))                                                       # Define matrix
  
  ## Get functions (Survival, Growth and Fecundity) ##
  
  # Growth-Survival continuous stage
  G = outer(midp$meshp, midp$meshp, L_gaj_cr, Site=Site, Year=Year)    # Growth
  ## Control for eviction (equivalent to redistributing evicted sizes evenly among existing size classes)
  correct.G=1-apply(G,2,sum)                                                    # 1-sum of all columns (to check if the sum probabilities is less than 1)
  correct.G[correct.G<0]=0                                                      # just to prevent getting negative correction values (can happen when you play around and your IPM is too small)
  G[n,]=G[n,]+correct.G    
  
  S = diag(L_saj_cr(Size=midp$meshp,Site=Site, Year=Year))
  P = midp$h*as.matrix(G%*%S)   
  
  # Fecundity reproductive adults
  F_con =  L_fa_cr(Size=midp$meshp, P_t1=P_t1)
  F_dis = L_rs_cr(TMe=TMe, lintrad=lintrad)
  
  F = outer(F_dis,F_con,"*")
  
  # Survival discrete seedling stage
  Ss = matrix(as.vector(L_ss_cr(TMaS_t1=TMaS_t1, lintrad = lintrad)) *   # Size in which discrete stage is incorporated into continuous stage
                as.vector(L_ss2_cr(Size=midp$meshp, TMaS_t1=TMaS_t1,Site=Site)), 
              nrow = length(midp$meshp), ncol = 1)
  
  ## Place functions into matrix
  K[1,2:(n+1)] = F                                                              # Top row -> discrete seedling stage
  K[2:(n+1),1] = Ss                                                             # Continuous-Size stages (discrete seedling stage survival) 
  K[2:(n+1),2:(n+1)] = P                                                        # Continuous-Size stages (continuous-Size stage survival-growth )
  return(K)
}


# D. PERFORM IPM ####
## D.1. Set IPM parameters ####
# Set midpoint parameters
n=100
midp = calc_midp(
  n = 100,                                                                      
  L = 0.9*with(sajL, min(Size, na.rm =TRUE)),
  U = 1.1*with(sajL, max(Size, na.rm =TRUE)))

# Define a initial population vector (first year -2019- population)

n.site_L = unique(as.factor(rsL$Site)[!is.na(rsL$nadult)])                      

# Define covariates, RCD_t1 and lintrad (for all the years -2019:2022-)
# Define covariates, RCD_t1 and lintrad (for all the years -2019:2023-)
cov_L = rsL %>%                                                                 
  arrange(Site, Year) %>%
  select(Year, Site, lintrad, TMaS_t1,P_t1,TMe, lintrad)

cov_L=cov_L[!cov_L$Site%in%"MM4",]

#### Minimum and maximum of covariates 

#########
min.TMe=min(cov_L$TMe)
max.TMe=max(cov_L$TMe)

year.min.TMe=as.character(unique(cov_L$Year[cov_L$TMe==min.TMe]))
year.max.TMe=as.character(unique(cov_L$Year[cov_L$TMe==max.TMe]))

mean.TMe=mean(cov_L$TMe)
sd.TMe=sd(cov_L$TMe[!duplicated(cov_L$TMe)])

dens.when.TMe.min=mean(cov_L$lintrad[cov_L$TMe==min.TMe])
dens.when.TMe.max=mean(cov_L$lintrad[cov_L$TMe==max.TMe])

TMaS_t1.when.TMe.min=unique(cov_L$TMaS_t1[cov_L$TMe==min.TMe])
TMaS_t1.when.TMe.max=unique(cov_L$TMaS_t1[cov_L$TMe==max.TMe])

TMaS.when.TMe.min=unique(cov_L$TMaS[cov_L$TMe==min.TMe])
TMaS.when.TMe.max=unique(cov_L$TMaS[cov_L$TMe==max.TMe])

TMe_t1.when.TMe.min=unique(cov_L$TMe_t1[cov_L$TMe==min.TMe])
TMe_t1.when.TMe.max=unique(cov_L$TMe_t1[cov_L$TMe==max.TMe])

P_t1.when.TMe.min=unique(cov_L$P_t1[cov_L$TMe==min.TMe])
P_t1.when.TMe.max=unique(cov_L$P_t1[cov_L$TMe==max.TMe])

mean.dens=mean(cov_L$lintrad[cov_L$Site=="OM3"])
#########

min.TMaS_t1=min(cov_L$TMaS_t1)
max.TMaS_t1=max(cov_L$TMaS_t1)

year.min.TMaS_t1=as.character(unique(cov_L$Year[cov_L$TMaS_t1==min.TMaS_t1]))
year.max.TMaS_t1=as.character(unique(cov_L$Year[cov_L$TMaS_t1==max.TMaS_t1]))


mean.TMaS_t1=mean(cov_L$TMaS_t1[!duplicated(cov_L$TMaS_t1)])
sd.TMaS_t1=sd(cov_L$TMaS_t1[!duplicated(cov_L$TMaS_t1)])

dens.when.TMaS_t1.min=mean(cov_L$lintrad[cov_L$TMaS_t1==min.TMaS_t1])
dens.when.TMaS_t1.max=mean(cov_L$lintrad[cov_L$TMaS_t1==max.TMaS_t1])

TMe.when.TMaS_t1.min=unique(cov_L$TMe[cov_L$TMaS_t1==min.TMaS_t1])
TMe.when.TMaS_t1.max=unique(cov_L$TMe[cov_L$TMaS_t1==max.TMaS_t1])

TMaS.when.TMaS_t1.min=unique(cov_L$TMaS[cov_L$TMaS_t1==min.TMaS_t1])
TMaS.when.TMaS_t1.max=unique(cov_L$TMaS[cov_L$TMaS_t1==max.TMaS_t1])

TMe_t1.when.TMaS_t1.min=unique(cov_L$TMe_t1[cov_L$TMaS_t1==min.TMaS_t1])
TMe_t1.when.TMaS_t1.max=unique(cov_L$TMe_t1[cov_L$TMaS_t1==max.TMaS_t1])

P_t1.when.TMaS_t1.min=unique(cov_L$P_t1[cov_L$TMaS_t1==min.TMaS_t1])
P_t1.when.TMaS_t1.max=unique(cov_L$P_t1[cov_L$TMaS_t1==max.TMaS_t1])

#########
min.TMaS=min(cov_L$TMaS)
max.TMaS=max(cov_L$TMaS)

year.min.TMaS=as.character(unique(cov_L$Year[cov_L$TMaS==min.TMaS]))
year.max.TMaS=as.character(unique(cov_L$Year[cov_L$TMaS==max.TMaS]))


mean.TMaS=mean(cov_L$TMaS)
sd.TMaS=sd(cov_L$TMaS[!duplicated(cov_L$TMaS)])

dens.when.TMaS.min=mean(cov_L$lintrad[cov_L$TMaS==min.TMaS])
dens.when.TMaS.max=mean(cov_L$lintrad[cov_L$TMaS==max.TMaS])

TMaS_t1.when.TMaS.min=unique(cov_L$TMaS_t1[cov_L$TMaS==min.TMaS])
TMaS_t1.when.TMaS.max=unique(cov_L$TMaS_t1[cov_L$TMaS==max.TMaS])

TMe.when.TMaS.min=unique(cov_L$TMe[cov_L$TMaS==min.TMaS])
TMe.when.TMaS.max=unique(cov_L$TMe[cov_L$TMaS==max.TMaS])

TMe_t1.when.TMaS.min=unique(cov_L$TMe_t1[cov_L$TMaS==min.TMaS])
TMe_t1.when.TMaS.max=unique(cov_L$TMe_t1[cov_L$TMaS==max.TMaS])

P_t1.when.TMaS.min=unique(cov_L$P_t1[cov_L$TMaS==min.TMaS])
P_t1.when.TMaS.max=unique(cov_L$P_t1[cov_L$TMaS==max.TMaS])

#########
min.TMe_t1=min(cov_L$TMe_t1)
max.TMe_t1=max(cov_L$TMe_t1)

year.min.TMe_t1=as.character(unique(cov_L$Year[cov_L$TMe_t1==min.TMe_t1]))
year.max.TMe_t1=as.character(unique(cov_L$Year[cov_L$TMe_t1==max.TMe_t1]))


mean.TMe_t1=mean(cov_L$TMe_t1)
sd.TMe_t1=sd(cov_L$TMe_t1[!duplicated(cov_L$TMe_t1)])

dens.when.TMe_t1.min=mean(cov_L$lintrad[cov_L$TMe_t1==min.TMe_t1])
dens.when.TMe_t1.max=mean(cov_L$lintrad[cov_L$TMe_t1==max.TMe_t1])

TMaS_t1.when.TMe_t1.min=unique(cov_L$TMaS_t1[cov_L$TMe_t1==min.TMe_t1])
TMaS_t1.when.TMe_t1.max=unique(cov_L$TMaS_t1[cov_L$TMe_t1==max.TMe_t1])

TMaS.when.TMe_t1.min=unique(cov_L$TMaS[cov_L$TMe_t1==min.TMe_t1])
TMaS.when.TMe_t1.max=unique(cov_L$TMaS[cov_L$TMe_t1==max.TMe_t1])

TMe.when.TMe_t1.min=unique(cov_L$TMe[cov_L$TMe_t1==min.TMe_t1])
TMe.when.TMe_t1.max=unique(cov_L$TMe[cov_L$TMe_t1==max.TMe_t1])

P_t1.when.TMe_t1.min=unique(cov_L$P_t1[cov_L$TMe_t1==min.TMe_t1])
P_t1.when.TMe_t1.max=unique(cov_L$P_t1[cov_L$TMe_t1==max.TMe_t1])

#########
min.P_t1=min(cov_L$P_t1)
max.P_t1=max(cov_L$P_t1)

year.min.P_t1=as.character(unique(cov_L$Year[cov_L$P_t1==min.P_t1]))
year.max.P_t1=as.character(unique(cov_L$Year[cov_L$P_t1==max.P_t1]))


mean.P_t1=mean(cov_L$P_t1)
sd.P_t1=sd(cov_L$P_t1[!duplicated(cov_L$P_t1)])

dens.when.P_t1.min=mean(cov_L$lintrad[cov_L$P_t1==min.P_t1])
dens.when.P_t1.max=mean(cov_L$lintrad[cov_L$P_t1==max.P_t1])

TMaS_t1.when.P_t1.min=unique(cov_L$TMaS_t1[cov_L$P_t1==min.P_t1])
TMaS_t1.when.P_t1.max=unique(cov_L$TMaS_t1[cov_L$P_t1==max.P_t1])

TMaS.when.P_t1.min=unique(cov_L$TMaS[cov_L$P_t1==min.P_t1])
TMaS.when.P_t1.max=unique(cov_L$TMaS[cov_L$P_t1==max.P_t1])

TMe_t1.when.P_t1.min=unique(cov_L$TMe_t1[cov_L$P_t1==min.P_t1])
TMe_t1.when.P_t1.max=unique(cov_L$TMe_t1[cov_L$P_t1==max.P_t1])

TMe.when.P_t1.min=unique(cov_L$TMe[cov_L$P_t1==min.P_t1])
TMe.when.P_t1.max=unique(cov_L$TMe[cov_L$P_t1==max.P_t1])


## D.2. Do IPM  ####

library(MASS)

sens.out=NULL

for (i in 1:100){
  
  print(paste("Running iteration ",i))

     # TMe no cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = min.TMe,
        TMaS_t1 = 0,
        P_t1 = 0,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.TMe)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = max.TMe,
        TMaS_t1 = 0,
        P_t1 = 0,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.TMe)
      
      deltaTMe = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.TMe-min.TMe)/sd.TMe)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="TMe",
                                         driver="Temperature",
                                         driver.type="C",
                                         sens=deltaTMe,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=0,
                                         sim=i))


      # TMaS_t1 no cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = 0,
        TMaS_t1 = min.TMaS_t1,
        P_t1 = 0,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.TMaS_t1)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = 0,
        TMaS_t1 = max.TMaS_t1,
        P_t1 = 0,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.TMaS_t1)
      
      deltaTMaS_t1 = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.TMaS_t1-min.TMaS_t1)/sd.TMaS_t1)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="TMaS_t1",
                                         driver="Temperature",
                                         driver.type="C",
                                         sens=deltaTMaS_t1,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=0,
                                         sim=i))
      
      # P_t1 no cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = 0,
        TMaS_t1 = 0,
        P_t1 = min.P_t1,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.P_t1)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = 0,
        TMaS_t1 = 0,
        P_t1 = max.P_t1,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.P_t1)
      
      deltaP_t1 = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.P_t1-min.P_t1)/sd.P_t1)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="P_t1",
                                         driver="Rain",
                                         driver.type="C",
                                         sens=deltaP_t1,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=0,
                                         sim=i))
      
      ########################################################
      
      # TMe cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = min.TMe,
        TMaS_t1 = TMaS_t1.when.TMe.min,
        P_t1 = P_t1.when.TMe.min,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.TMe)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = max.TMe,
        TMaS_t1 = TMaS_t1.when.TMe.max,
        P_t1 = P_t1.when.TMe.max,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.TMe)
      
      
      deltaTMecov = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.TMe-min.TMe)/sd.TMe)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="TMe",
                                         driver="Temperature",
                                         driver.type="C",
                                         sens=deltaTMecov,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=1,
                                         sim=i))
      
      
      # TMaS_t1 cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = TMe.when.TMaS_t1.min,
        TMaS_t1 = min.TMaS_t1,
        P_t1 = P_t1.when.TMaS_t1.min,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.TMaS_t1)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = TMe.when.TMaS_t1.max,
        TMaS_t1 = max.TMaS_t1,
        P_t1 = P_t1.when.TMaS_t1.max,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.TMaS_t1)
      
      deltaTMaS_t1cov = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.TMaS_t1-min.TMaS_t1)/sd.TMaS_t1)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="TMaS_t1",
                                         driver="Temperature",
                                         driver.type="C",
                                         sens=deltaTMaS_t1cov,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=1,
                                         sim=i))
      
      # P_t1 cov
      # min
      mpm.min = K.fnc_L_cr (
        midp = midp, 
        TMe = TMe.when.P_t1.min,
        TMaS_t1 = TMaS_t1.when.P_t1.min,
        P_t1 = min.P_t1,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.min.P_t1)
      
      # max
      mpm.max = K.fnc_L_cr (
        midp = midp, 
        TMe = TMe.when.P_t1.max,
        TMaS_t1 = TMaS_t1.when.P_t1.max,
        P_t1 = max.P_t1,
        lintrad = mean.dens,
        Site = "MM3",
        Year=year.max.P_t1)
      
      deltaP_t1cov = abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.P_t1-min.P_t1)/sd.P_t1)
      
      
      sens.out=rbind(sens.out,data.frame(spec.driver="P_t1",
                                         driver="Rain",
                                         driver.type="C",
                                         sens=deltaP_t1cov,
                                         l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                         cov=1,
                                         sim=i))
    
}

ggplot(sens.out, aes(x=sens, color=factor(cov))) +
  geom_density()+
  facet_grid(driver~.,scales = "free")

ggplot(sens.out, aes(x=l_ratio, color=factor(cov))) +
  geom_density()+
  facet_grid(driver~.,scales = "free")



sens_lav=data.frame(species="Lavandula stoechas", study.doi="Sanchez Mejia in prep",year.of.publication=NA,
                           group="Plants",continent="Europe",driver=sens.out$driver,driver.type="C",
                           stage.age="all",vital.rates="all",sens=sens.out$sens,cov=sens.out$cov,mat=3,n.vr=7,n.pam=35,dens=0,
                           biotic_interactions=0,lambda.sim=0,study.length=5,l_ratio=sens.out$l_ratio)

write.csv(sens_lav,"sens_lavandula.csv",row.names = F)


