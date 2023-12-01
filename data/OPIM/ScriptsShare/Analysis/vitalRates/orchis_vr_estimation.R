#Script fits a Bayesian model for the six time-varying vital rates found in O. purpurea

## Load packages
rm(list=ls())
library(rstan) ; library(rstanmulticore)

## Define functions
invlogit<-function(x){exp(x)/(1+exp(x))}

## load demographic data
#Contact Hans Jacquemyn, hans.jacquemyn@kuleuven.be, for this data
orchid.demog<-read.csv("Data/orchis_2003-2014.csv")


#Prepare data###############################################################################
#growth, survival, flowering, fertility, flower-to-fruit, and dormancy data
grow=na.omit(orchid.demog[,c("begin.year","end.year","total.leaf.area","end.total.leaf.area","survival")])
surv=na.omit(orchid.demog[,c("begin.year","end.year","total.leaf.area","survival")])
flower=na.omit(orchid.demog[,c("begin.year","end.year","total.leaf.area","flowering")])
fert=na.omit(orchid.demog[,c("begin.year","end.year","total.leaf.area","number.flowers")])
fruit=na.omit(orchid.demog[,c("begin.year","end.year","number.fruits","number.flowers")])
dorm=na.omit(orchid.demog[,c("begin.year","end.year","total.leaf.area","dormant")])
#Take subset of data that only includes seedlings
seedling.size=orchid.demog[orchid.demog$number.leaves==-999 & orchid.demog$end.year<2013,]
seedling=na.omit(seedling.size[,c("begin.year","end.year","end.total.leaf.area")])


## Plot data
plot(log(grow$total.leaf.area),log(grow$end.total.leaf.area))
plot(log(surv$total.leaf.area),surv$survival)
plot(log(flower$total.leaf.area),flower$flowering)
plot(log(fert$total.leaf.area),fert$number.flowers)
plot(log(dorm$total.leaf.area),dorm$dormant)


#GROWTH - data for Stan 
years<-unique(grow$end.year)
NyearsG<-length(years)
nG<-nrow(grow)
xG<-log(grow$total.leaf.area)
yG<-log(grow$end.total.leaf.area)
newyearsG<-grow$end.year-2003


#SURVIVAL - data for Stan 
years<-unique(surv$end.year)
NyearsS<-length(years)
nS<-nrow(grow)
xS<-log(surv$total.leaf.area)
yS<-surv$survival
newyearsS<-surv$end.year-2003


#FLOWERING - data for Stan 
years<-unique(flower$end.year)
NyearsF<-length(years)
nF<-nrow(flower)
xF<-log(flower$total.leaf.area)
yF<-flower$flowering
newyearsF<-flower$end.year-2003


#FERTILITY - data for Stan 
years<-unique(fert$end.year)
NyearsR<-length(years)
nR<-nrow(fert)
xR<-log(fert$total.leaf.area)
yR<-fert$number.flowers
newyearsR<-fert$end.year-2003


#FLOWER-TO-FRUIT - data for Stan
#Note, no size 'x' values because this vital rates is size-independent
years<-unique(fruit$end.year)
NyearsP<-length(years)
nP<-nrow(fruit) 
yP<-fruit$number.fruits
Nfl<-fruit$number.flowers #Number of flowers for n in beta-binomial distribution
newyearsP<-fruit$end.year-2003


#DORMANCY - data for Stan
years<-unique(dorm$end.year)
NyearsD<-length(years)
nD<-nrow(dorm)
xD<-log(dorm$total.leaf.area)
yD<-dorm$dormant
newyearsD<-dorm$end.year-2003


### Data for Stan model
orchis_dat <- list(nG = nrow(grow), nS = nrow(surv), nF = nrow(flower), nR = nrow(fert), nP = nrow(fruit), nD = nrow(dorm),
                 nYears = NyearsG, #n. of years is same for all species
                 yearsG = newyearsG, yearsS = newyearsS, yearsF = newyearsF, yearsR = newyearsR, yearsP = newyearsP, yearsD = newyearsD, 
                 yG = yG, yS = yS, yF = yF, yR = yR, yD = yD, yP = yP,
                 xG = xG, xS = xS, xF = xF, xR = xR, xD = xD, Nfl= Nfl,
                 nVital=6,Omega = diag(6))


#Stan model
sink("Analysis/vitalRates/orchis_vr_estimation.stan")
cat("
    data {
      int<lower=1> nVital;         # N. of vital rates               
      int<lower=0> nYears;         # N. of years. Same for all vital rates

      int<lower=0> nG;             # N. of data points for the growth model
      int<lower=0> yearsG[nG];     # Index for years
      vector[nG] xG;               # log size at time t
      vector[nG] yG;               # log size at time t+1 
      
      int<lower=0> nS;             # N. of data points for the survival model 
      int<lower=0> yearsS[nS];     # Index for years
      vector[nS] xS;               # log size at time t
      int<lower=0,upper=1> yS[nS]; # Survival at time t. Values are either 0 or 1
      
      int<lower=0> nF;             # N. of data points for the flowering model      
      int<lower=0> yearsF[nF];     # Index for years
      vector[nF] xF;               # log size at time t
      int<lower=0,upper=1> yF[nF]; # Flowering at time t. Values are either 0 or 1
      
      int<lower=0> nR;             # N. of data points for the fertility model 
      int<lower=0> yearsR[nR];     # Index for years
      vector[nR] xR;               # log size at time t
      int<lower=0> yR[nR];         # Number of flowers at time t

      int<lower=0> nP;             # N. of data points for the flower-to-fruit model 
      int<lower=0> yearsP[nP];     # Index for years
      int<lower=0> Nfl[nP];        # Total number of flowers
      int<lower=0> yP[nP];         # N. of viable flowers

      int<lower=0> nD;             # N. of data points for the dormancy model
      int<lower=0> yearsD[nD];     # Index for years
      vector[nD] xD;               # log size at time t
      int<lower=0,upper=1> yD[nD]; # Dormancy status at time t+1. Values are either 0 or 1

    }
    
    parameters {

      corr_matrix[nVital] Omega;    # correlation matrix
      vector<lower=0>[nVital] tau;  # standard deviation of random effects
      
      vector[nVital] u;             # Mean of random year intercepts.
      vector[nVital] beta[nYears];  # Intercept coefficients by year
    
      real<lower=0> sigma_y;        # Residual for growth model
      real betaG;                   # Growth reg. slope
      real betaS;                   # Survival reg. slope
      real betaF;                   # Flowering reg. slope
      real betaR;                   # Fertility reg. slope
      real betaD;                   # Dormancy reg. slope
      real<lower=0> thetaR;         # Fertility dispersion parameter
      real<lower=0.1> lambda;       # Flower-to-fruit dispersion parameter

    }
    
    transformed parameters { #This code performs the re-parameterization of flower-to-fruit model

      #Reparameterization of beta binomial model
      real<lower=0> alpha[nYears];         
      real<lower=0> betaBin[nYears];
      real<lower=0,upper=1> phi[nYears];   #This is the inv_logit of the random year effect
  
      for (n in 1:nYears){                 #Loop through random year effects
        phi[n] <- inv_logit(beta[n,5]);
        alpha[n] <- lambda * phi[n];       #NOTE: lambda is the flower-to-fruit dispersion parameter
        betaBin[n] <- lambda * (1 - phi[n]);    
      }

    }
    
    model {
    
      #Variables to run models
      int indGY;   #placeholders of indexes for the random year effects
      int indSY;
      int indFY;
      int indRY;
      int indPY;
      int indDY;
      real mG[nG]; #Placeholders for the growth (G), survival (S), flowering (F),
      real mS[nS]; #fertility (B) models, and dormancy (D) models
      real mF[nF];
      real mR[nR];
      vector[nP] alphaP;   #place holder for Flower-to-fruit's alpha 
      vector[nP] betaP;    #place holder for Flower-to-fruit's beta  
      real mD[nD];
      matrix[nVital,nVital] Sigma_beta; #Plance holder for var-cov matrix

      #Hypepriors
      tau ~ cauchy(0,2.5);   #Standard deviation of random effects
      Omega ~ lkj_corr(2);   #Correlation matrix
      u ~ normal(0,100);     #Means of random effects
    
      #Priors
      Sigma_beta <- quad_form_diag(Omega,tau); #Generate Var-Cov matrix
      for (n in 1:nYears){ 
        beta[n] ~ multi_normal(u, Sigma_beta);
      }  
      betaG ~ normal(0, 1000);  #Growth slope
      betaS ~ normal(0, 1000);  #Survival slope
      betaF ~ normal(0, 1000);  #Flowering slope 
      betaR ~ normal(0, 1000);  #Fertility slope 
      betaD ~ normal(0, 1000);  #Dormancy slope   
      sigma_y ~ inv_gamma(0.001, 0.001);   #Growth model residual st. dev.   
      thetaR ~ inv_gamma(0.001, 0.001);    #Fertility dispersion parameter
      lambda ~ pareto(0.1,1.5);            #Flower-to-fruit dispersion parameter

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
      for(nflow in 1:nR){
        indRY <- yearsR[nflow];
        mR[nflow] <- beta[indRY,4] + betaR * xR[nflow];
      } 
      yR ~ neg_binomial_2_log(mR,thetaR);

      #5. FLOWER-TO-FRUIT MODEL   
      for(nfruit in 1:nP){
        indPY <- yearsP[nfruit];
        alphaP[nfruit] <- alpha[indPY];
        betaP[nfruit]  <- betaBin[indPY];
      } 
      yP ~ beta_binomial(Nfl, alphaP, betaP);

      #6. DORMANCY MODEL   
      for(d in 1:nD){
        indDY <- yearsD[d];
        mD[d] <- beta[indDY,6] + betaD * xD[d];
      }
      yD ~ bernoulli_logit(mD);

    }",fill=T)

sink()


rstan_options(auto_write = TRUE) #I am suggested to run these two lines 
options(mc.cores = parallel::detectCores())
fit <- stan(file = 'Analysis/vitalRates/orchis_vr_estimation.stan', data = orchis_dat, 
            iter = 5000, warmup = 1000, chains = 4,
            control=list(max_treedepth=15,adapt_delta=0.999))
out=extract(fit)


#MEAN PARAMETERS---------------------------------------------------------------------------------------------------------
#Mean Correlation matrix
var_corrMat=apply(out$Omega,c(2,3),mean)
write.csv(var_corrMat,"Results/vitalRates/orchid_correlation_mean.csv",row.names=F)

#Coefficients (NOTE: same means can be obtained using 'summary(fit)$summary')
betas=cbind(out$betaG,out$betaS,out$betaF,out$betaR,out$betaD)
betas=apply(betas,2,mean)
muScaled=apply(out$u,2,mean)
sigma=mean(out$sigma_y)
sigScaled=apply(out$tau,2,mean)
thetaR=mean(out$thetaR)
thetaP=mean(out$betaBin)
outCoef=c(betas,NA,muScaled,sigma,sigScaled,thetaR,thetaP)
outCoef=as.data.frame(matrix(outCoef,1,length(outCoef)))
names(outCoef)=c("beta.G","beta.S","beta.F","beta.R","beta.D","deviance","mu.intG",
                 "mu.intS","mu.intF","mu.intR","mu.intP","mu.intD","sigmaG","sigma.intG",
                 "sigma.intS","sigma.intF","sigma.intR","sigma.intP","sigma.intD",
                 "thetaR","thetaP")
outCoef=outCoef[,sort(names(outCoef))]

write.csv(outCoef,"Results/vitalRates/orchid_coefficient_mean.csv",row.names=F)


#POSTERIOR OF PARAMETERS---------------------------------------------------------------------------------------------------------
#Correlation values
rhoSim=as.data.frame(cbind(out$Omega[,1,2],out$Omega[,1,3],out$Omega[,1,4],out$Omega[,1,5],out$Omega[,1,6],
                           out$Omega[,2,3],out$Omega[,2,4],out$Omega[,2,5],out$Omega[,2,6],
                           out$Omega[,3,4],out$Omega[,3,5],out$Omega[,3,6],
                           out$Omega[,4,5],out$Omega[,4,6],
                           out$Omega[,5,6]))
names(rhoSim)=c("growth-survival","growth-flowering","growth-fertility","growth-fruiting","growth-dormancy",
                "survival-flowering","survival-fertility","survival-fruiting","survival-dormancy",
                "flowering-fertility","flowering-fruiting","flowering-dormancy",
                "fertility-fruiting","fertility-dormancy",
                "fruiting-dormancy")
write.csv(rhoSim,"Results/vitalRates/orchid_correlation_posterior.csv",row.names=F)


#NON-CORRELATION PARAMETERS
alphaGrid=expand.grid(alpha="alpha",vr=c("R","F","G","S","P","D"),year=c(1:11))
alphaGrid=alphaGrid[order(alphaGrid[,1],alphaGrid[,2],alphaGrid[,3]),]
alphaGrid$nameChar=paste(alphaGrid$alpha,".",alphaGrid$vr,"[",alphaGrid$year,"]",sep="")
alphaSim=data.frame(cbind(out$beta[,,4],out$beta[,,3],out$beta[,,1],out$beta[,,2],out$beta[,,5],out$beta[,,6]))
names(alphaSim)=alphaGrid$nameChar

slopes=data.frame(out$betaG,out$betaS,out$betaF,out$betaR,out$betaD)
names(slopes)=c("beta.G","beta.S","beta.F","beta.R","beta.D")
alphaMean=data.frame(out$u[,1],out$u[,2],out$u[,3],out$u[,4],out$u[,5],out$u[,6])
sigmaG=data.frame(out$sigma_y)
thetaR=data.frame(out$thetaR)
lambdaP=data.frame(out$lambda)
alphaSd=data.frame(out$tau[,1],out$tau[,2],out$tau[,3],out$tau[,4],out$tau[,5],out$tau[,6])
names(alphaMean)=c("mu.intG","mu.intS","mu.intF","mu.intR","mu.intP","mu.intD")
names(sigmaG)=c("sigmaG")
names(alphaSd)=c("sigma.intG", "sigma.intS","sigma.intF","sigma.intR","sigma.intP","sigma.intD")
names(thetaR)=c("thetaR")
names(lambdaP)=c("lambdaP")

simulations=cbind(alphaSim,slopes,NA,alphaMean,sigmaG,alphaSd,thetaR,lambdaP)

write.csv(simulations,"Results/vitalRates/orchid_coefficient_posterior.csv",row.names=F)
