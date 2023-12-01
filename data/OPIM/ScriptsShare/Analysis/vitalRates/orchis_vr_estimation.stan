
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
      for(ngrow in 1:nG){
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

    }
