
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
      for(nflow in 1:nB){
        indBY <- yearsB[nflow];
        mB[nflow] <- beta[indBY,4] + betaB * xB[nflow];
      }
      yB ~ neg_binomial_2_log(mB,alphaB);
    
    }
    
