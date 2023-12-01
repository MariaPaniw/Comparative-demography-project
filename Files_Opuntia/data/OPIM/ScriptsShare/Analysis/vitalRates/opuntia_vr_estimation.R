#Script fits a Bayesian model for the four time-varying vital rates found in O. imbricata  

#Load packages
rm(list=ls())
library(rstan)

## define functions
invlogit<-function(x){exp(x)/(1+exp(x))}
vol<-function(h,w,p){(1/3)*pi*h*(w/2)*(p/2)} #update from BRAD - less biased VOL estimate 

## load demographic data
demog.dat<-read.csv("Data/opuntia_20042015.csv")


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


#Prepare data###############################################################################
demog.dat$sizet<-vol(h=demog.dat$Height_t,w=demog.dat$Width_t,p=demog.dat$Perp_t)
demog.dat$sizet1<-vol(h=demog.dat$Height_t1,w=demog.dat$Width_t1,p=demog.dat$Perp_t1)
demog.dat$Flower_t<-demog.dat$Goodbuds_t>0                  #Flower or not
#growth, survival, flowering, and fertility data
grow=na.omit(demog.dat[,c("Year_t1","sizet","sizet1","Survival_t1")])
surv=na.omit(demog.dat[,c("Year_t1","sizet","Survival_t1")])
flow=na.omit(demog.dat[,c("Year_t1","sizet","Flower_t")])
buds=na.omit(subset(demog.dat,Flower_t==T)[,c("Year_t1","sizet","Goodbuds_t")])


## Plot data
par(mfrow=c(2,2))
plot(log(grow$sizet), log(grow$sizet1))
plot(log(surv$sizet), surv$Survival_t1)
plot(log(flow$sizet), flow$Flower_t)
plot(log(buds$sizet), buds$Goodbuds_t)


#GROWTH - data for Stan 
years<-unique(grow$Year_t1)
NyearsG<-length(years)
nG<-nrow(grow)
xG<-log(grow$sizet)   #Log transform for normality
yG<-log(grow$sizet1)  #Log transform for normality
years<-grow$Year_t1-2004
newyearsG = as.numeric(factor(as.character(years)))


#SURVIVAL - data for Stan
years<-unique(surv$Year_t1)
NyearsS<-length(years)
nS<-nrow(surv)
xS<-log(surv$sizet)   #Log transform for normality
yS<-surv$Survival_t1  
years<-surv$Year_t1-2004
newyearsS = as.numeric(factor(as.character(years)))


#FLOWERING - data for Stan
years<-unique(flow$Year_t1)
NyearsF<-length(years)
nF<-nrow(flow)
xF<-log(flow$sizet)   #Log transform for normality
yF<- as.numeric(flow$Flower_t)  
years<-flow$Year_t1-2004
newyearsF = as.numeric(factor(as.character(years)))


#FERTILITY - data for Stan
years<-unique(buds$Year_t1)
NyearsB<-length(years)
nB<-nrow(buds)
xB<-log(buds$sizet)   #Log transform for normality
yB<-buds$Goodbuds_t  
years<-buds$Year_t1-2004
newyearsB = as.numeric(factor(as.character(years)))

#Test all data sets have the same number of years
NyearsS==NyearsG ; NyearsG==NyearsF ; NyearsF==NyearsB ; NyearsB==NyearsS
Nyears=NyearsB


### Data for Stan model
cholla_dat <- list(nG = nrow(grow), nS = nrow(surv), nF = nrow(flow), nB = nrow(buds),
                   nYears = NyearsG, #n. of years is same for all species
                   yearsG = newyearsG, yearsS = newyearsS, yearsF = newyearsF, yearsB = newyearsB,
                   yG = yG, yS = yS, yF = yF, yB = yB ,
                   xG = xG, xS = xS, xF = xF, xB = xB,
                   nVital=4,Omega = diag(4))


#Stan model
sink("Analysis/vitalRates/opuntia_vr_estimation.stan")
cat("
    data { 

      int<lower=1> nVital;         # N. of vital rates               
      int<lower=0> nYears;         # N. of years. Same for all vital rates.

      #Data for the growth model
      int<lower=0> nG;             # N. of data points for the growth model  
      int<lower=0> yearsG[nG];     # Index for years
      vector[nG] xG;               # log size at time t
      vector[nG] yG;               # log size at time t+1 
      
      #Data for the survival model
      int<lower=0> nS;             # N. of data points for the survival model  
      int<lower=0> yearsS[nS];     # Index for years
      vector[nS] xS;               # log size at time t
      int<lower=0,upper=1> yS[nS]; # Survival at time t+1. Values are either 0 or 1
      
      #Data for the flowering model
      int<lower=0> nF;             # N. of data points for the flowering model  
      int<lower=0> yearsF[nF];     # Index for years
      vector[nF] xF;               # log size at time t
      int<lower=0,upper=1> yF[nF]; # Flowering status at time t. Values are either 0 or 1
      
      #Data for the fertility model
      int<lower=0> nB;             # N. of data points for the fertility model  
      int<lower=0> yearsB[nB];     # Index for years
      vector[nB] xB;               # log size at time t
      int<lower=0> yB[nB];         # Number of flowers at time t

    }
    
    parameters {

      corr_matrix[nVital] Omega;    # prior correlation. This is a 4 by 4 matrix
      vector<lower=0>[nVital] tau;  # prior scale (st. dev.). A vector of length 4.
      
      vector[nVital] beta[nYears];  # intercept coeffs by year. A 10 by 4 matrix
      vector[nVital] u;             # intercept means. A vector of length 4.
      
      real<lower=0> sigma_y;        # Residual st. dev. for growth model
      real betaG;                   # Growth reg. slope
      real betaS;                   # Survival reg. slope
      real betaF;                   # Flowering reg. slope
      real betaB;                   # Fertility reg. slope
      real<lower=0> alphaB;         # Fertility dispersion parameter

    }
    
    model {
    
      int indGY;   #placeholders of indexes for the random year effects
      int indSY;
      int indFY;
      int indBY;
      real mG[nG]; #Placeholders for the growth (G), survival (S), flowering (F) and fertility (B) models
      real mS[nS];
      real mF[nF];
      real mB[nB];
      matrix[nVital,nVital] Sigma_beta; #Plance holder for var-cov matrix      

      #Hyperpriors
      tau ~ cauchy(0,2.5);   #Standard deviation of random effects
      Omega ~ lkj_corr(2);   #Correlation matrix
      u ~ normal(0,100);     #Means of random effects

      Sigma_beta <- quad_form_diag(Omega,tau);  #Generate var-cov matrix
      #Priors
      for (n in 1:nYears){ #Random year effects 
        beta[n] ~ multi_normal(u, Sigma_beta);
      }
      betaG ~ normal(0, 1000);   #Growth slope 
      betaS ~ normal(0, 1000);   #Survival slope
      betaF ~ normal(0, 1000);   #Flowering slope
      betaB ~ normal(0, 1000);   #Fertility slope 
      alphaB ~ inv_gamma(0.001, 0.001);   #Fertility dispersion parameter 
      sigma_y ~ inv_gamma(0.001, 0.001);  #Growth model residual st. dev.   

      #LIKELIHOODS###################################################
      #1. GROWTH MODEL
      for(ngrow in 1:nG){  # this (and subsequent 'for' loops) loops over random effects
        indGY <- yearsG[ngrow];
        mG[ngrow] <- beta[indGY,1] + betaG * xG[ngrow];
      }  
      yG ~ normal(mG, sigma_y);
          
      #2. SURVIVAL MODEL   
      for(nsurv in 1:nS){
        indSY <- yearsS[nsurv];
        mS[nsurv] <- beta[indSY,2] + betaS * xS[nsurv];
      }
      yS ~ bernoulli_logit(mS);
      
      #3. FLOWERING MODEL   
      for(nflow in 1:nF){
        indFY <- yearsF[nflow];
        mF[nflow] <- beta[indFY,3] + betaF * xF[nflow];
      }
      yF ~ bernoulli_logit(mF);
      
      #4. FERTILITY MODEL   
      for(nflow in 1:nB){
        indBY <- yearsB[nflow];
        mB[nflow] <- beta[indBY,4] + betaB * xB[nflow];
      }
      yB ~ neg_binomial_2_log(mB,alphaB);
    
    }
    ",fill=T)
sink()


#These two lines prompt stan to use four cores in parallel
rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores())
fit <- stan(file = 'Analysis/vitalRates/opuntia_vr_estimation.stan', data = cholla_dat, 
            iter = 5000, warmup = 1000, chains = 4)
out=extract(fit)



#MEAN PARAMETERS---------------------------------------------------------------------------------------------------------
#Correlation matrix
var_corrMat=apply(out$Omega,c(2,3),mean)
write.csv(var_corrMat,"Results/vitalRates/opuntia_correlation_mean.csv",row.names=F)

#Coefficients (NOTE: same means can be obtained using 'summary(fit)$summary')
betas=cbind(out$betaG,out$betaS,out$betaF,out$betaB,out$alphaB)
betas=apply(betas,2,mean)
muScaled=apply(out$u,2,mean)
sigma=mean(out$sigma_y)
sigScaled=apply(out$tau,2,mean)
outCoef=c(betas,NA,muScaled,sigma,sigScaled)
outCoef=as.data.frame(matrix(outCoef,1,15))
names(outCoef)=c("beta.G","beta.S","beta.F","beta.B","theta.B","deviance","mu.int.G","mu.int.S","mu.int.F",
                 "mu.int.B","sigma","sigma.int.G","sigma.int.S","sigma.int.F","sigma.int.B")
outCoef=outCoef[,c("beta.B","beta.F","beta.G","beta.S","deviance","mu.int.B","mu.int.F","mu.int.G","mu.int.S","sigma","sigma.int.B","sigma.int.F","sigma.int.G","sigma.int.S","theta.B")]

write.csv(outCoef,"Results/vitalRates/opuntia_coefficient_mean.csv",row.names=F)


#POSTERIOR OF PARAMETERS---------------------------------------------------------------------------------------------------------
#Correlation values
rhoSim=as.data.frame(cbind(out$Omega[,1,2],out$Omega[,1,3],out$Omega[,1,4],
                           out$Omega[,2,3],out$Omega[,2,4],out$Omega[,3,4]))
names(rhoSim)=c("growth-survival","growth-flowering","growth-fertility","survival-flowering","survival-fertility","flowering-fertility")
write.csv(rhoSim,"Results/vitalRates/opuntia_correlation_posterior.csv",row.names=F)


#ALL non-correlation simulations - Problem is aligning this to the old simulations file 
#Alpha names
alphaGrid=expand.grid(alpha="alpha",vr=c("B","F","G","S"),year=c(1:10))
alphaGrid=alphaGrid[order(alphaGrid[,1],alphaGrid[,2],alphaGrid[,3]),]
alphaGrid$nameChar=paste(alphaGrid$alpha,".",alphaGrid$vr,"[",alphaGrid$year,"]",sep="")
alphaSim=data.frame(cbind(out$beta[,,4],out$beta[,,3],out$beta[,,1],out$beta[,,2]))
names(alphaSim)=alphaGrid$nameChar
slopes=data.frame(out$betaB,out$betaF,out$betaG,out$betaS)
names(slopes)=c("beta.B","beta.F","beta.G","beta.S")
thetaB=data.frame(out$alphaB)
alphaMean=data.frame(out$u[,4],out$u[,3],out$u[,1],out$u[,2])
sigmaG=data.frame(out$sigma_y)
alphaSd=data.frame(out$tau[,4],out$tau[,3],out$tau[,1],out$tau[,2])
names(alphaMean)=c("mu.int.B","mu.int.F","mu.int.G","mu.int.S")
names(sigmaG)=c("sigma")
names(alphaSd)=c("sigma.int.B", "sigma.int.F", "sigma.int.G", "sigma.int.S")
names(thetaB)=c("theta.B")

simulations=cbind(alphaSim,slopes,NA,alphaMean,sigmaG,alphaSd,thetaB)

write.csv(simulations,"Results/vitalRates/opuntia_coefficient_posterior.csv",row.names=F)

