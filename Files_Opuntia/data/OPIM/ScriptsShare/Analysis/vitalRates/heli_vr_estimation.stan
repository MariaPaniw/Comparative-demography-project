
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
      for(ngrow in 1:nG){
        indGY<-yearsG[ngrow];
        mG[ngrow] <- beta[indGY,1] + betaG*xG[ngrow]; # looping over R.E. to estimate params.
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

    }
