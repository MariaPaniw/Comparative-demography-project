## Code to build an Integral Projection Model (IPM) for Opuntia imbricata.

## IPM based on script from Aldo Compagnoni. Accompanying article https://doi.org/10.1002/ecm.1228
## 2 vital rate models are different. These are based on analysis from Sanne Evers. Accompanying article https://doi.org/10.1111/gcb.15519

# Edited by: Esin Ickin (14.09.2023) (step 7)
# Search for "ESIN" to find the edited lines
# In this code we're trying to estimate the uncertainties of the sensitivities. For that we resample 100 times from the posterior distributions.

library(tidyverse)
library(lme4)
setwd("~/OpuntiaImbricata")

s_mod <- readRDS("results/OPIM/OPIM_s_sliding2.rds")[[3]]$BestModel
fp_mod <- readRDS("results/gcb/OPIM_fp_month_result.rds")[[4]]$BestModel

# years climate dataframe
s_df <- readRDS("results/OPIM/OPIM_s_sliding2.rds")[[3]]$BestModelData %>%
  rename(tmin_n1_n9 = climate) %>% select(year, tavg_1_n5, tmin_n1_n9) %>% distinct
fp_df <- readRDS("results/gcb/OPIM_fp_month_result.rds")[[4]]$BestModelData %>%
  rename(tavg_17_12 = climate) %>% select(year, tavg_17_12) %>% distinct

clim_df <- left_join(s_df, fp_df)


#Useful functions
invlogit<-function(x){exp(x)/(1+exp(x))}
vol<-function(h,w,p){(1/3)*pi*h*(w/2)*(p/2)}


###1.Load vital rate model results#############################################################################

#"Variance-correlation matrix
vCor = read.csv("data/OPIM/ScriptsShare/Results/vitalRates/opuntia_correlation_mean.csv")
vCor = as.matrix(vCor)

#Coefficients
coefs1=read.csv("data/OPIM/ScriptsShare/Results/vitalRates/opuntia_coefficient_mean.csv")
#coefs=read.csv("Results/vitalRates/cholla_coef_ljk_(nbin2015).csv")
# ESIN:
coefs=read.csv("data/OPIM/ScriptsShare/Results/vitalRates/opuntia_coefficient_posterior.csv") # get posteriors and not the mean as above


#Posterior(s) of correlations
corPost=read.csv("data/OPIM/ScriptsShare/Results/vitalRates/opuntia_correlation_posterior.csv")

#Posterior of coefficients  
coefPost=read.csv("data/OPIM/ScriptsShare/Results/vitalRates/opuntia_coefficient_posterior.csv")


#Time-varying model parameters
fit.grow=coefs[,c("mu.int.G","beta.G","sigma.int.G","sigma")] # Growth (random year)
fit.surv=coefs[,c("mu.int.S","beta.S","sigma.int.S")]         # Survival (random year)
fit.flower=coefs[,c("mu.int.F","beta.F","sigma.int.F")]       # Probability of flowering
fit.buds=coefs[,c("mu.int.B","beta.B","sigma.int.B")]         # Size-dependent fertility


###2. Vital rates estimated 'from scratch'##########################################################

##A. Dispersal of seeds from maternal plant and survival until germination time
## (from fruit removal experiments)-----
fruit.dat<-read.csv("data/OPIM/opuntia_FruitSurvival.csv",T)
names(fruit.dat)

fruitsurv<-glm((Fr.on.grnd.not.chewed/Fr.on.plant~1),weights=Fr.on.plant,family="binomial",data=fruit.dat)
invlogit(coef(fruitsurv)[1])


##B. Germination----------------------
germ.dat<-read.csv("data/OPIM/opuntia_Germination.csv",T)
names(germ.dat)

## germination of 1-yo seeds
germ1yo<-glm(Seedlings04/Input~1,weights=Input,family="binomial",data=germ.dat)
invlogit(coef(germ1yo)[1])

## germination of 2-yo seeds
germ2yo<-glm(Seedlings05/(Input-Seedlings04)~1,weights=(Input-Seedlings04),family="binomial",data=germ.dat)
invlogit(coef(germ2yo)[1])


##C. Survival from germination to census------------------------
precensus.dat<-read.csv("data/OPIM/opuntia_PrecensusSurvival.csv",T)
names(precensus.dat)

precensus.surv<-glm(survive0405~1,family="binomial",data=precensus.dat)
invlogit(coef(precensus.surv)[1])


##D. Size distribution of recruits------------------------------
demog.dat<-read.csv("data/OPIM/opuntia_20042015.csv",T)

#MODIFICATIONS TO APPLY "vol" function#########################################################
#I change two "demog.dat$Perp_t" values that are == 0 to 0.2.
#Justified because these are actually small individuals (see below).
demog.dat[demog.dat$Perp_t==0,c("Height_t","Width_t","Perp_t")]
demog.dat$Perp_t[demog.dat$Perp_t==0]=0.2
#I also change "zeros" in year t1
demog.dat[c(394,3117,3258),c("Height_t1","Width_t1","Perp_t1")]
demog.dat$Perp_t1[c(3117,3258)]=0.2
#LIKELY MISTAKE: "Width_t1" is zero, but not in the previous years
#We therefore exclude this individual from analysis
demog.dat[394,]
demog.dat[demog.dat$TagID=="LJ12" & demog.dat$Plot=="T1",]
demog.dat=demog.dat[-394,]

#Data to calculate recruit size
recruits.dat <- subset(demog.dat,Recruit==1 & Year_t==2005)
names(recruits.dat)

kidsize.mean<-with(recruits.dat,
		mean(log(vol(h=Height_t,w=Width_t,p=Perp_t))))
kidsize.sd<-with(recruits.dat,
		sd(log(vol(h=Height_t,w=Width_t,p=Perp_t))))


#E. Calculate sizes for all individuals---------------------------
demog.dat$sizet<-vol(h=demog.dat$Height_t,w=demog.dat$Width_t,p=demog.dat$Perp_t)
demog.dat$sizet1<-vol(h=demog.dat$Height_t1,w=demog.dat$Width_t1,p=demog.dat$Perp_t1)


#F. Number of seeds per bud (from Johanna Ohm's data)---------------------------
seednum <- 139.4815 #Data comes from Ohm and Miller, published in Ecology in 2014 
  

###3.Collect IPM parameters#############################################################
cholla <- list()
cholla[[1]] <- fit.grow[,c("mu.int.G")]					## growth intercept
cholla[[2]] <- fit.grow[,c("beta.G")]					  ## growth slope
cholla[[3]] <- NA #fixef(fit.grow)[3]					  ## growth quadratic term
cholla[[4]] <- NA #resid.model$coef[1]					## intercept of size-dependent variance
cholla[[5]] <- NA #resid.model$coef[2]					## slope of size-dependent variance
cholla[[6]] <- fixef(s_mod)[1]					        ## survival intercept
cholla[[7]] <- fixef(s_mod)[2]					        ## survival slope
cholla[[8]]<- fixef(fp_mod)[1]				          ## pr_flower intercept
cholla[[9]]<- fixef(fp_mod)[2]				          ## pr_flower slope
cholla[[10]]<- fit.buds[,c("mu.int.B")]					## nfruits intercept
cholla[[11]]<- fit.buds[,c("beta.B")]					  ## nfruits slope
cholla[[12]]<- seednum						             ## seeds per fruit
cholla[[13]]<- invlogit(coef(fruitsurv)[1])^2 		## This is 6mo fruit survival, squared to estimate full year. But I think this way overestimates seed predation (I'm equating dispersal with predation).
											        ## Experiment with setting it to 1, i.e. all seeds get to the seed bank (note that germination rates account for seed predation)
cholla[[14]]<-invlogit(coef(germ1yo)[1])			## 1-year-old germination
cholla[[15]]<-invlogit(coef(germ2yo)[1])			## 2-year-old germination
cholla[[16]]<-invlogit(coef(precensus.surv)[1])	## pre-census recruit survival
cholla[[17]]<-kidsize.mean						## mean recruit size
cholla[[18]]<-kidsize.sd						## SD of recruit size
cholla[[19]]<-min(log(demog.dat$sizet),na.rm=T)	## minsize (floor)
cholla[[20]]<-max(log(demog.dat$sizet1),na.rm=T) ## maxsize (ceiling)
cholla[[21]]<-200								## matrix size (results seem v stable by 200)

## RANDOM EFFECTS
cholla[[22]] = fit.grow[,"sigma.int.G"]     ## sd of random year effects for growth 
cholla[[23]] = fit.surv[,"sigma.int.S"]     ## sd of random year effects for survival
cholla[[24]] = fit.flower[,"sigma.int.F"]   ## sd of random year effects for flowering
cholla[[25]] = fit.buds[,"sigma.int.B"]     ## sd of random year effects for nbuds

#Should I include "fixed effect SD" or the SD of residuals?
cholla[[26]] = fit.grow[,"sigma"]           #Only growth has residual variance


## CLIMATE EFFECTS
cholla[[27]] <- fixef(s_mod)[3]                 ## survival climate effect Tavg from 1 to -5
cholla[[28]] <- fixef(s_mod)[4]                 ## survival climate effect Tmin from -1 to -9
cholla[[29]] <- fixef(fp_mod)[3]                ## flower p climate effect Tavg from 17 to 12

#Run simulations, including "correlation on" sims. with "direct estimation of variance"
#diag(vCor)=cholla[22:25]^2
varCovar=function(corMat,varVec){
  out = corMat * (varVec %*% t(varVec))
  return(out)
}
#varVec=cholla[c(22:25)]
varVec=lapply(cholla[22:25], function(inner_list) inner_list[[1]])


###4. Define demographic functions###########################################################

#GROWTH FROM SIZE X TO Y
#Function produces a *probability density distribution* from each
#1. Size at time t, (the *x* argument) to each
#2. Size at time t+1 (the *y* argument)
gxy<-function(x,y,params,rfx,i){
	xb=pmin(pmax(x,params[[19]][[1]]),params[[20]][[1]]) #Transforms all values below/above limits in min/max size
	return(dnorm(y,mean=params[[1]][[i]]+params[[2]][[i]]*xb + rfx[1],sd=params[[26]][[i]]))
	}

#SURVIVAL AT SIZE X.
#This calculates the survival prob. of each size class 
sx<-function(x,params,rfx, s_tavg, s_tmin){
	xb=pmin(pmax(x,params[[19]][[1]]),params[[20]][[1]])
	return(invlogit(params[[6]][[1]] + params[[7]][[1]]*xb + params[[27]][[1]]*s_tavg + params[[28]][[1]]*s_tmin + rfx[2]))
	}

#SURVIVAL*GROWTH
pxy<-function(x,y,params,rfx, s_tavg, s_tmin){
	xb=pmin(pmax(x,params[[19]][[1]]),params[[20]][[1]])
	sx(xb,params,rfx, s_tavg, s_tmin)*gxy(xb,y,params,rfx,i)
	#gxy(xb,y,params,random) ## set survival equal to one and check colsums of T
	}

#PRODUCTION OF 1-YO SEEDS IN THE SEED BANK FROM X-SIZED MOMS
fx<-function(x,params,rfx, fp_tavg,i){
	xb=pmin(pmax(x,params[[19]][[1]]),params[[20]][[1]])
	p.flow<-invlogit(params[[8]][[1]] + params[[9]][[1]]*xb + params[[29]][[1]]*fp_tavg + rfx[3]) #Probability of flowering
	nfruits<-exp(params[[10]][[i]] + params[[11]][[i]]*xb + rfx[4])   #n. of flowers produced

  return(p.flow*nfruits*params[[12]][[1]]*params[[13]][[1]])  #Params[c(12,13)] are n. of seeds/fruit, 
                                                #plus seed survival to next year
}



#SIZE DISTRIBUTION OF RECRUITS
#If a recruit is produced, what is the likelihood of falling 
#being in a certain size (*y*) class?
recruits<-function(y,params){
	dnorm(x=y,mean=params[[17]][[1]],sd=params[[18]][[1]])
}


###5. Simulation parameters###############################################################
max.yrs       <- 50000
matsize       <- 200
extra.grid    <- 2        #How many row/col add to the part of the matrix referring to sizes?
random        <- T     
allparms      <- cholla
varCorr       <- vCor

#Eviction parameters
floor.extend   <- 0.5
ceiling.extend <- 2.25
lower          <- cholla[[19]][[1]] - floor.extend
upper          <- cholla[[20]][[1]] + ceiling.extend


##6. Construct IPM kernel################################################################################
bigmatrix<-function(params,random,sigma,lower,upper,rand.seed, s_tavg, s_tmin, fp_tavg,i){ 
  
	n<-params[[21]][[1]]
	L<-lower; U<-upper;
	#these are the upper and lower integration limits
	h<-(U-L)/n                   #Bin size
	b<-L+c(0:n)*h;               #Lower boundaries of bins 
	y<-0.5*(b[1:n]+b[2:(n+1)]);  #Bins' midpoints
	#these are the boundary points (b) and mesh points (y)
	
	# Random effects
	rfx <- matrix(0,1,4)  #Set random effect to 0 (in case random=F,corr=F) 
	if(random==T){        
	  set.seed(rand.seed)
	  rfx[1,1:4] = rmvnorm(n=1, mean=rep(0,4), sigma=sigma)
	}
 
	# Fertility matrix
	Fmat<-matrix(0,(n+2),(n+2))

	# Banked seeds go in top row
	Fmat[1,3:(n+2)]<-fx(y,params,rfx, fp_tavg,i)

	# Growth/survival transition matrix
	Tmat<-matrix(0,(n+2),(n+2))

	# Graduation to 2-yo seed bank = pr(not germinating as 1-yo)
	Tmat[2,1]<-1-params[[14]][[1]]

	# Graduation from 1-yo bank to cts size = germination * size distn * pre-census survival
	Tmat[3:(n+2),1]<-params[[14]][[1]]*recruits(y,params)*h   

	# Graduation from 2-yo bank to cts size = germination * size distn * pre-census survival
	Tmat[3:(n+2),2]<-params[[15]][[1]]*recruits(y,params)*h   

	# Growth/survival transitions among cts sizes
	Tmat[3:(n+2),3:(n+2)]<-t(outer(y,y,pxy,params=params,rfx=rfx, s_tavg = s_tavg, s_tmin = s_tmin))*h

	# Put it all together
	IPMmat<-Fmat+Tmat     #Full Kernel is simply a summation ot fertility
                        #and transition matrix
	return(list(IPMmat=IPMmat,Fmat=Fmat,Tmat=Tmat,meshpts=y,rfx=rfx))
}

## Calculate deterministic lambda
#i = 2                                   ## Set this i from 1:13 for the different years

#lambda<-Re(eigen(bigmatrix(params=cholla,random=F,
#                           lower=cholla[19],upper=cholla[20],rand.seed=42,
#                           s_tavg = clim_df$tavg_1_n5[i], s_tmin = clim_df$tmin_n1_n9[i], 
#                           fp_tavg = clim_df$tavg_17_12[i])$IPMmat)$values[1])
#lambda

# 7. Sensitivity analyses ################################################################################
# Here, we calculate scaled sensitivities, according to Morris et al. 2020 (DOI: <https://doi.org/10.1073/pnas.1918363117>)
# Note that this is a step that Esin will implement in her MS thesis. With the information given in 1-3, we should be able to run these analyses. 

#### Climate data ---------------------------------

# min, max, sd, mean of covariates
# tavg_1_n5
min_tavg_1_n5=min(clim_df$tavg_1_n5)
max_tavg_1_n5=max(clim_df$tavg_1_n5)
mean_tavg_1_n5=mean(clim_df$tavg_1_n5)
sd_tavg_1_n5=sd(clim_df$tavg_1_n5)

# tmin_n1_n9
min_tmin_n1_n9=min(clim_df$tmin_n1_n9)
max_tmin_n1_n9=max(clim_df$tmin_n1_n9)
mean_tmin_n1_n9=mean(clim_df$tmin_n1_n9)
sd_tmin_n1_n9=sd(clim_df$tmin_n1_n9)

# tavg_17_12
min_tavg_17_12=min(clim_df$tavg_17_12)
max_tavg_17_12=max(clim_df$tavg_17_12)
mean_tavg_17_12=mean(clim_df$tavg_17_12)
sd_tavg_17_12=sd(clim_df$tavg_17_12)

# covariation among covariates
# tavg_1_n5
tmin_n1_n9_when_tavg_1_n5_max=clim_df$tmin_n1_n9[which(clim_df$tavg_1_n5==max(clim_df$tavg_1_n5))][1]
tmin_n1_n9_when_tavg_1_n5_min=clim_df$tmin_n1_n9[which(clim_df$tavg_1_n5==min(clim_df$tavg_1_n5))][1]

tavg_17_12_when_tavg_1_n5_max=clim_df$tavg_17_12[which(clim_df$tavg_1_n5==max(clim_df$tavg_1_n5))][1]
tavg_17_12_when_tavg_1_n5_min=clim_df$tavg_17_12[which(clim_df$tavg_1_n5==min(clim_df$tavg_1_n5))][1]

# tmin_n1_n9
tavg_1_n5_when_tmin_n1_n9_max=clim_df$tavg_1_n5[which(clim_df$tmin_n1_n9==max(clim_df$tmin_n1_n9))][1]
tavg_1_n5_when_tmin_n1_n9_min=clim_df$tavg_1_n5[which(clim_df$tmin_n1_n9==min(clim_df$tmin_n1_n9))][1]

tavg_17_12_when_tmin_n1_n9_max=clim_df$tavg_17_12[which(clim_df$tmin_n1_n9==max(clim_df$tmin_n1_n9))][1]
tavg_17_12_when_tmin_n1_n9_min=clim_df$tavg_17_12[which(clim_df$tmin_n1_n9==min(clim_df$tmin_n1_n9))][1]

# tavg_17_12
tmin_n1_n9_when_tavg_17_12_max=clim_df$tmin_n1_n9[which(clim_df$tavg_17_12==max(clim_df$tavg_17_12))][1]
tmin_n1_n9_when_tavg_17_12_min=clim_df$tmin_n1_n9[which(clim_df$tavg_17_12==min(clim_df$tavg_17_12))][1]

tavg_1_n5_when_tavg_17_12_max=clim_df$tavg_1_n5[which(clim_df$tavg_17_12==max(clim_df$tavg_17_12))][1]
tavg_1_n5_when_tavg_17_12_min=clim_df$tavg_1_n5[which(clim_df$tavg_17_12==min(clim_df$tavg_17_12))][1]

#### IPM Kernel function ----------------------------------------------
kernel.fun<-function(tavg_1_n5,tmin_n1_n9,tavg_17_12,i){
 lam <-  Re(eigen(bigmatrix(params=cholla,random=F,
                     lower=cholla[[19]][[1]],upper=cholla[[20]][[1]],rand.seed=42,
                     s_tavg = tavg_1_n5, s_tmin = tmin_n1_n9, 
                     fp_tavg = tavg_17_12,i=i)$IPMmat)$values[1])
 return(lam)
  
}

#### Sensitivity to tavg_1_n5 -------------------------------------
# no covariation (mean covariates)
Delta.tavg_1_n5=NULL
for(i in 1:100){
  max_lam<-kernel.fun(max_tavg_1_n5,mean_tmin_n1_n9,mean_tavg_17_12,i=i)
  min_lam<-kernel.fun(min_tavg_1_n5,mean_tmin_n1_n9,mean_tavg_17_12,i=i)
  
  Delta.tavg_1_n5[i] <- abs((max_lam-min_lam)/((max_tavg_1_n5-min_tavg_1_n5)/sd_tavg_1_n5))
  
}

hist(Delta.tavg_1_n5)

# with covariation
Delta.tavg_1_n5.cov=NULL
for(i in 1:100){
  max_lam<-kernel.fun(max_tavg_1_n5,tmin_n1_n9_when_tavg_1_n5_max,tavg_17_12_when_tavg_1_n5_max,i=i)
  min_lam<-kernel.fun(min_tavg_1_n5,tmin_n1_n9_when_tavg_1_n5_min,tavg_17_12_when_tavg_1_n5_min,i=i)
  
  Delta.tavg_1_n5.cov[i] <- abs((max_lam-min_lam)/((max_tavg_1_n5-min_tavg_1_n5)/sd_tavg_1_n5))
  
}

hist(Delta.tavg_1_n5.cov)

# Differences between the two sensitivities (no cov-cov)
hist(Delta.tavg_1_n5-Delta.tavg_1_n5.cov)
mean(Delta.tavg_1_n5-Delta.tavg_1_n5.cov) #0.002312954

#### Sensitivity to tmin_n1_n9 -----------------------------------------------------------------
# no covariation
Delta.tmin_n1_n9=NULL
for(i in 1:100){
  max_lam<-kernel.fun(mean_tavg_1_n5,max_tmin_n1_n9,mean_tavg_17_12,i=i)
  min_lam<-kernel.fun(mean_tavg_1_n5,min_tmin_n1_n9,mean_tavg_17_12,i=i)
  
  Delta.tmin_n1_n9[i] <- abs((max_lam-min_lam)/((max_tmin_n1_n9-min_tmin_n1_n9)/sd_tmin_n1_n9))
}

hist(Delta.tmin_n1_n9)

# with covariation
Delta.tmin_n1_n9.cov=NULL
for(i in 1:100){
  max_lam<-kernel.fun(tavg_1_n5_when_tmin_n1_n9_max,max_tmin_n1_n9,tavg_17_12_when_tmin_n1_n9_max,i=i)
  min_lam<-kernel.fun(tavg_1_n5_when_tmin_n1_n9_min,min_tmin_n1_n9,tavg_17_12_when_tmin_n1_n9_min,i=i)
  
  Delta.tmin_n1_n9.cov[i] <- abs((max_lam-min_lam)/((max_tmin_n1_n9-min_tmin_n1_n9)/sd_tmin_n1_n9))
}

hist(Delta.tmin_n1_n9.cov)

# Differences (no cov - cov)
hist(Delta.tmin_n1_n9-Delta.tmin_n1_n9.cov)


#### Sensitivity to tavg_17_12 -----------------------------------------
# no covariation
Delta.tavg_17_12=NULL
for(i in 1:100){
  max_lam<-kernel.fun(mean_tavg_1_n5,mean_tmin_n1_n9,max_tavg_17_12,i=i)
  min_lam<-kernel.fun(mean_tavg_1_n5,mean_tmin_n1_n9,min_tavg_17_12,i=i)
  
  Delta.tavg_17_12[i] <- abs((max_lam-min_lam)/((max_tavg_17_12-min_tavg_17_12)/sd_tavg_17_12))
}

hist(Delta.tavg_17_12)

# with covariation
Delta.tavg_17_12.cov=NULL
for(i in 1:100){
  max_lam<-kernel.fun(tavg_1_n5_when_tavg_17_12_max,tmin_n1_n9_when_tavg_17_12_max,max_tavg_17_12,i=i)
  min_lam<-kernel.fun(tavg_1_n5_when_tavg_17_12_min,tmin_n1_n9_when_tavg_17_12_min,min_tavg_17_12,i=i)
  
  Delta.tavg_17_12.cov[i] <- abs((max_lam-min_lam)/((max_tavg_17_12-min_tavg_17_12)/sd_tavg_17_12))
}

hist(Delta.tavg_17_12.cov)

# Differences (no cov - cov)
hist(Delta.tavg_17_12 - Delta.tavg_17_12.cov)
mean(Delta.tavg_17_12 - Delta.tavg_17_12.cov)

# 8. Save Output -----------------------------------
delta_df <- data.frame(delta.tavg_1_n5 = Delta.tavg_1_n5,
                       delta.tavg_1_n5.cov = Delta.tavg_1_n5.cov,
                       delta.tavg_17_12 = Delta.tavg_17_12,
                       delta.tavg_17_12.cov = Delta.tavg_17_12.cov,
                       delta.tmin_n1_n9 = Delta.tmin_n1_n9,
                       delta.tmin_n1_n9.cov = Delta.tmin_n1_n9.cov)
#write.csv(delta_df, file = "Sens_w_resampling.csv", row.names = F)

