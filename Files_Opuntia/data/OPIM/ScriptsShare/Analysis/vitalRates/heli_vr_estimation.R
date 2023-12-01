#Script fits a Bayesian model for the four time-varying vital rates found in O. imbricata  

#Load packages
rm(list=ls())
library(rstan)

## Demographic data ##
demog.dat<-read.csv("Data/helianthella_MAXFIELD_1998-2013.csv")

## Maniupulate data types ################################################################

#change class of clumps to numeric
demog.dat$clumpsT0<-as.numeric(as.character(demog.dat$clumpsT0))
demog.dat$clumpsT10<-as.numeric(as.character(demog.dat$clumpsT10))

#create column for year_t1 by adding 1 to year
demog.dat$year_t1<-demog.dat$year+1 # Shouldn't this be a factor?
demog.dat$totStalks <- demog.dat$OkstalksT + demog.dat$AbortStalks + demog.dat$BrokeStalks

#######################################################
## Subsetting data for different vital rate fxn's ##
#######################################################

#make four datasets- one for growth data, one for survival data, one for probability of flowering, one for fertility
#these will differ because different rows will have NA's depending on whether we need growth or survival data
#growth, survival, flowering, fertility, and flower-to-fruit data
grow=na.omit(demog.dat[,c("year","year_t1","clumpsT0","clumpsT10","surv")])
surv=na.omit(demog.dat[,c("year","year_t1","clumpsT0","surv")])

demog.dat.flow=subset(demog.dat,year!=2012)  #exclude last year of data for reproduction models
demog.dat.fert=subset(demog.dat.flow,pflowerT==1)  #only incluide individuals that flower
flow=na.omit(demog.dat.flow[,c("year","year_t1","clumpsT0","pflowerT")])
buds=na.omit(demog.dat.fert[,c("year","year_t1","clumpsT0","totStalks")])
abort=na.omit(demog.dat.fert[,c("year","year_t1","totStalks","OkstalksT")])

# Checking to see if all the years are present in each dataframe
lapply(list(grow, surv, flow, buds), function(x) unique(as.factor(x$year)))


#SET UP DATA#######################################################################################

#GROWTH - data for Stan 
nG<-nrow(grow) 
years<-unique(grow$year) 
NyearsG<-length(years) 
newyearsG<-grow$year-1997 #make new vector called 'newyearsG' that has '1' for 1998, '2' for 1999, '3' for 2000 etc. In the cholla dataset, we had to make a vector called years, and then use a for loop to take account of the missing year of data collection. The for loop populated a vector called newyearsG that has the year of the data collection, as here (starting with 1 for 1998). 
xG<-log(grow$clumpsT0)
yG<-grow$clumpsT10  

#SURVIVAL - data for Stan
nS<-nrow(surv) 
years<-unique(surv$year)  
NyearsS<-length(years)    
newyearsS<-surv$year-1997 
xS<-log(surv$clumpsT0)
yS<-surv$surv  
NyearsS==NyearsG #check to make sure that we have the same number of years of data for Survival and Growth data
Nyears=NyearsS   #For clarity, here I rename NyearsS 'Nyears'


#FLOWERING - data for Stan
nF<-nrow(flow)
years<-unique(flow$year)
NyearsF<-length(years)
newyearsF<-flow$year-1997 
xF<-log(flow$clumpsT0)
yF<-as.numeric(flow$pflowerT)
unique(flow$clumpsT0)
hist(log(flow$clumpsT0))


#FERTILITY - data for Stan
years<-unique(buds$year)
NyearsB<-length(years)
newyearsB<-buds$year-1997
nB<-nrow(buds)
xB<-log(buds$clumpsT0)
yB<-buds$totStalks 


#FLOWER-TO-FRUIT - data for Stan
nA<-nrow(abort)
years<-unique(abort$year)
NyearsA<-length(years)
newyearsA<-abort$year-1997
xA<-abort$totStalks 
yA<-abort$OkstalksT 


### Data for Stan model
heli_dat <- list(nG = nrow(grow), nS = nrow(surv), nF = nrow(flow), nB = nrow(buds), nA= nrow(abort),
                 nYears = NyearsG, #n. of years is same for all species
                 yearsG = newyearsG, yearsS = newyearsS, yearsF = newyearsF, yearsB = newyearsB, yearsA = newyearsA,
                 yG = yG, yS = yS, yF = yF, yB = yB, yA = yA,
                 xG = xG, xS = xS, xF = xF, xB = xB, xA = xA,
                 nVital=5,Omega = diag(5)) #, U=max(yG)+ceiling(0.1*max(yG))+1


#Stan model
sink("Analysis/vitalRates/heli_vr_estimation.stan")
cat("
    data {

      int<lower=1> nVital;         # N. of vital rates               
      int<lower=0> nYears;         # N. of years. Same for all vital rates.

      # Data for Growth model  
      int<lower=0> nG;             # N. of data points for the growth model  
      int<lower=0> yearsG[nG];     # Index for year random effect
      vector[nG] xG;               # log size (n. of clumps) at time t
      int yG[nG];                  # size (N. of clumps) size at time t+1 
      
      # Data for Survival model 
      int<lower=0> nS;               # N. of data points for the survival model
      int<lower=0> yearsS[nS];       # Index for year random effect
      vector[nS] xS;                 # log size (n. of clumps) at time t
      int<lower=0,upper=1> yS[nS];   # Survival at time t+1. Values are either 0 or 1
      
      # Data for Flowering model 
      int<lower=0> nF;              # N. of data points for the flowering model  
      int<lower=0> yearsF[nF];      # Index for year random effect
      vector[nF] xF;                # log size (n. of clumps) at time t
      int<lower=0,upper=1> yF[nF];  # Flowering status at time t. Values are either 0 or 1
      
      # Data for Fecundity model 
      int<lower=0> nB;            # N. of data points for the flowering model  
      int<lower=0> yearsB[nB];    # Index for year random effect
      vector[nB] xB;              # log size (n. of clumps) at time t
      int<lower=0> yB[nB];        # Total N. of flowering stalks at time t.

      # Data for Fruit-to-flower model 
      int<lower=0> nA;            # N. of data points for the flowering model    
      int<lower=0> yearsA[nA];    # Index for year random effect
      int<lower=0> xA[nA];        # Total number of stalks
      int<lower=0> yA[nA];        # N. of viable stalks

    }
    
    parameters {

      corr_matrix[nVital] Omega;    # correlation matrix
      vector<lower=0>[nVital] tau;  # standard deviation of random effects
      
      vector[nVital] u;             # Mean of random year intercepts.
      vector[nVital] beta[nYears];  # Intercept coefficients by year
      
      real betaG;                   # Growth reg. slope
      real betaS;                   # Survival reg. slope
      real betaF;                   # Flowering reg. slope
      real betaB;                   # Fertility reg. slope

      real<lower=0> nb_alphaG;       # Growth overdispersion paramter
      real<lower=0,upper=1> pMax;    # Asymptote for flowering model      
      real<lower=0> nb_alphaB;       # Fertility overdispersion paramter
      real<lower=0.1> lambdaA;       # Flower-to-fruit overdispersion parameter

    }

    transformed parameters {         #This code performs the re-parameterization of flower-to-fruit model
      
      #Reparameterization of beta binomial model
      real<lower=0> alpha[nYears];         
      real<lower=0> betaBin[nYears];
      real<lower=0,upper=1> phi[nYears];   #This is the inv_logit of the random year effect
  
      for (n in 1:nYears){                 #Loop for Beta Binomial year-specific parameters
        phi[n] <- inv_logit(beta[n,5]);
        alpha[n] <- lambdaA * phi[n];      #NOTE: lambdaA is the flower-to-fruit dispersion parameter
        betaBin[n] <- lambdaA * (1 - phi[n]);    
      }

    }

    model {

      int indGY; #placeholders of indexes for the random year effects
      int indSY;
      int indFY;
      int indBY;
      int indAY;
      vector[nG] mG; #Placeholders for the growth (G), survival (S), flowering (F) and fertility (B) models
      vector[nS] mS; 
      vector[nF] mF; 
      vector[nB] mB; 
      vector[nA] alphaA;   #place holder for Flower-to-fruit's alpha 
      vector[nA] betaA;    #place holder for Flower-to-fruit's beta  
      matrix[nVital,nVital] Sigma_beta; #Plance holder for var-cov matrix

      # Hyperpriors
      tau ~ cauchy(0,2.5);   #Standard deviation of random effects
      Omega ~ lkj_corr(2);   #Correlation matrix
      u ~ normal(0,100);     #Means of random effects
    
      # Priors
      Sigma_beta <- quad_form_diag(Omega,tau); #Generate Var-Cov matrix
      for (n in 1:nYears){
        beta[n] ~ multi_normal(u, Sigma_beta);
      }
      betaG ~ normal(0, 1000);   #Growth slope
      betaS ~ normal(0, 1000);   #Survival slope
      betaF ~ normal(0, 1000);   #Flowering slope
      betaB ~ normal(0, 1000);   #Fertility slope
      nb_alphaG ~ uniform(0, 100);   #Growth dispersion parameter
      nb_alphaB ~ uniform(0, 100);   #Fertility dispersion parameter
      lambdaA ~ pareto(0.1,1.5);   #Flower-to-fruit dispersion parameter
      pMax ~ uniform(0,1);         #Asymptote of flowering model

      #LIKELIHOODS###################################################
      #1. GROWTH MODEL
      for(ngrow in 1:nG){  # this (and subsequent 'for' loops) loops over random effects
        indGY<-yearsG[ngrow];
        mG[ngrow] <- beta[indGY,1] + betaG*xG[ngrow];
      }
      yG ~ neg_binomial_2_log(mG,nb_alphaG);
  
      #2. SURVIVAL MODEL  
      for(nsurv in 1:nS){
        indSY <- yearsS[nsurv];
        mS[nsurv] <- beta[indSY,2] + betaS*xS[nsurv];
      }
      yS ~ bernoulli_logit(mS);

      #3. FLOWERING MODEL   
      for(nflow in 1:nF){
        indFY <- yearsF[nflow];
        mF[nflow] <- inv_logit(beta[indFY,3] + betaF*xF[nflow])*pMax;
      }
      yF ~ bernoulli(mF);

      #4. FERTILITY MODEL   
      for(nflow in 1:nB){
        indBY <- yearsB[nflow];
        mB[nflow] <- beta[indBY,4] + betaB*xB[nflow];
      }
      yB ~ neg_binomial_2_log(mB, nb_alphaB);
    
      #5. FLOWER-TO-FRUIT MODEL
      for(nabor in 1:nA){
        indAY <- yearsA[nabor];
        alphaA[nabor] <- alpha[indAY];
        betaA[nabor]  <- betaBin[indAY];
      }
      yA ~ beta_binomial(xA, alphaA, betaA);

    }",fill=T)

sink()

rstan_options(auto_write = TRUE)  
options(mc.cores = parallel::detectCores())
fit <- stan(file = 'Analysis/vitalRates/heli_vr_estimation.stan', data = heli_dat, 
            iter = 5000, warmup = 1000, chains = 4)
out=extract(fit)


#MEAN PARAMETERS---------------------------------------------------------------------------------------------------------
#Variance-correlation matrix
var_corrMat=apply(out$Omega,c(2,3),mean)
write.csv(var_corrMat,"Results/vitalRates/heli_correlation_mean.csv",row.names=F)

#Coefficients (NOTE: same means can be obtained using 'summary(fit)$summary')
betas=cbind(out$betaG,out$betaG2,out$betaG3,out$betaS,out$betaS2,out$betaS3,
            out$betaF,out$betaF2,out$betaF3,out$betaB,out$betaB2,out$betaB3)
betas=apply(betas,2,mean)
muScaled=apply(out$u,2,mean)
sigScaled=apply(out$tau,2,mean)
thetaG=mean(out$nb_alphaG)
thetaB=mean(out$nb_alphaB)
thetaA=mean(out$lambdaA)
pMaxF=mean(out$pMax)
outCoef=c(betas,muScaled,sigScaled,thetaG,thetaB,thetaA,pMaxF)
outCoef=as.data.frame(matrix(outCoef,1,18))
names(outCoef)=c("beta.G","beta.S","beta.F","beta.B",
                 "mu.int.G","mu.int.S","mu.int.F","mu.int.B","mu.int.A",
                 "sigma.int.G","sigma.int.S","sigma.int.F","sigma.int.B","sigma.int.A",
                 "theta.G","theta.B","theta.A","pMaxF")
outCoef=outCoef[,c("beta.B","beta.F","beta.G","beta.S",
                   "mu.int.A","mu.int.B","mu.int.F","mu.int.G","mu.int.S",
                   "sigma.int.A","sigma.int.B","sigma.int.F","sigma.int.G","sigma.int.S",
                   "theta.A","theta.B","theta.G",
                   "pMaxF")]
write.csv(outCoef,"Results/vitalRates/heli_coefficient_mean.csv",row.names=F)


#POSTERIOR OF PARAMETERS---------------------------------------------------------------------------------------------------------
#Correlation values
rhoSim=as.data.frame(cbind(out$Omega[,1,2],out$Omega[,1,3],out$Omega[,1,4],out$Omega[,1,5],
                           out$Omega[,2,3],out$Omega[,2,4],out$Omega[,2,5],
                           out$Omega[,3,4],out$Omega[,3,5],out$Omega[,4,5]))
names(rhoSim)=c("growth-survival","growth-flowering","growth-fertility","growth-abortion",
                "survival-flowering","survival-fertility","survival-abortion",
                "flowering-fertility","flowering-abortion","fertility-abortion")
write.csv(rhoSim,"Results/vitalRates/heli_correlation_posterior.csv",row.names=F)


#ALL non-correlation posteriors - Problem is aligning this to the old simulations file 
#Alpha names
alphaGrid=expand.grid(alpha="alpha",vr=c("A","B","F","G","S"),year=c(1:14))
alphaGrid=alphaGrid[order(alphaGrid[,1],alphaGrid[,2],alphaGrid[,3]),]
alphaGrid$nameChar=paste(alphaGrid$alpha,".",alphaGrid$vr,"[",alphaGrid$year,"]",sep="")
alphaSim=data.frame(cbind(out$beta[,,5],out$beta[,,4],out$beta[,,3],out$beta[,,1],out$beta[,,2]))
names(alphaSim)=alphaGrid$nameChar
slopes=data.frame(out$betaB,out$betaF,out$betaG,out$betaS)
names(slopes)=c("beta.B","beta.F","beta.G","beta.S")
thetaA=data.frame(out$lambdaA)
thetaB=data.frame(out$nb_alphaB)
thetaG=data.frame(out$nb_alphaG)
pMaxF=data.frame(out$pMax)
alphaMean=data.frame(out$u[,5],out$u[,4],out$u[,3],out$u[,1],out$u[,2])
alphaSd=data.frame(out$tau[,5],out$tau[,4],out$tau[,3],out$tau[,1],out$tau[,2])
names(alphaMean)=c("mu.int.A","mu.int.B","mu.int.F","mu.int.G","mu.int.S")
names(alphaSd)=c("sigma.int.A","sigma.int.B", "sigma.int.F", "sigma.int.G", "sigma.int.S")
names(thetaA)=c("theta.A")
names(thetaB)=c("theta.B")
names(thetaG)=c("theta.G")
names(pMaxF)=c("pMaxF")
simulations=cbind(alphaSim,slopes,alphaMean,alphaSd,thetaA,thetaB,thetaG,pMaxF)#sigmaG,
write.csv(simulations,"Results/vitalRates/heli_coefficient_posterior.csv",row.names=F)
