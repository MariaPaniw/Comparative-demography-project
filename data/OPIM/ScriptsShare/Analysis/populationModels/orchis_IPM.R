##Code to build an Integral Projection Model (IPM) for Orchis purpurea


###1.Load vital rate model results#############################################################################

#Coefficients means
coefs=read.csv("Results/vitalRates/orchid_coefficient_mean.csv",header=T)

#Correlation matrix
vCor = read.csv("Results/vitalRates/orchid_correlation_mean.csv")
vCor = as.matrix(forceSymmetric(as.matrix(vCor)))

#Posterior of coefficients  
coefPost=read.csv("Results/vitalRates/orchid_coefficient_posterior.csv")

#Posterior of correlations
corPost=read.csv("Results/vitalRates/orchid_correlation_posterior.csv")


##2.Collect IPM parameters################################################################################

##Orchid demographic data
#Contact Hans Jacquemyn, hans.jacquemyn@kuleuven.be, for this data
orchid.demog<-read.csv("Data/orchis_2003-2014.csv")

seeds   <- 6000                 # seeds per fruit; see text
eta     <- mean(c(0.007,0.023)) # seed to protocorm transition; see text
# protocorm survival probability
sigmapL <- 0.01485249           # light
# tuber survival probability
sigmatL <- 0.05940997           # light
#Size mean and SD
Dsize <- mean(log(orchid.demog$end.total.leaf.area[orchid.demog$number.leaves==0]),na.rm=T)
Dsd   <-   sd(log(orchid.demog$end.total.leaf.area[orchid.demog$number.leaves==0]),na.rm=T)
#Recruit size mean and SD
seedling.size=orchid.demog[orchid.demog$number.leaves==-999 & orchid.demog$end.year<2013,]
seedling=na.omit(seedling.size[,c("begin.year","end.year","end.total.leaf.area")])
seedlM  <- mean(log(seedling$end.total.leaf.area))
seedlSD <- sd(log(seedling$end.total.leaf.area))
matsize <- 200                  # Default size of the approximating matrix

##Parameters from models
fit.grow<-coefs[,c("mu.intG","beta.G","sigma.intG","sigmaG")]  # 5. Growth (random year)
fit.surv<-coefs[,c("mu.intS","beta.S","sigma.intS")]           # 6. Survival (random year)
fit.flower<-coefs[,c("mu.intF","beta.F","sigma.intF")]         # 7. Probability of flowering
fit.fert<-coefs[,c("mu.intR","beta.R","sigma.intR")]           # 8. Size-dependent fertility
fit.fruit<-coefs[,c("mu.intP","sigma.intP")]                   # 9. Probability that flower becomes fruit
fit.dorm<-coefs[,c("mu.intD","beta.D","sigma.intD")]           # 10. Probability of going dormant

##Calculate min and max size
minsize <- min(log(na.omit(orchid.demog$total.leaf.area)))
maxsize <- max(log(na.omit(orchid.demog$total.leaf.area)))

orchid<-c(0)

# [00] IPM
orchid[1]  <- minsize            # lower size limit
orchid[2]  <- maxsize            # upper size limit
orchid[3]  <- matsize            # matrix size
# [10] Full growth model ###Vegetative Growth
orchid[4] <- fit.grow["mu.intG"]               # Growth intercept  (baseline)
orchid[5] <- fit.grow["beta.G"]                       # Growth slope wrt size
orchid[6] <- fit.grow["sigmaG"]                       # Residual SD for growth
# [30] Survival
orchid[7] <- fit.surv["mu.intS"]                     # Survival intercept (baseline)
orchid[8] <- fit.surv["beta.S"]                     # Survival slope (size)
orchid[9] <- sigmapL                         # Protocorm survival
orchid[10] <- sigmatL                        # Tuber survival
# [40] Seeds
orchid[11] <- fit.fruit["mu.intP"]					  # Mean flower-to-fruits intercept
orchid[12] <- seeds                           # Seeds per fruit
# [50] Dormancy
orchid[13] <- fit.dorm["mu.intD"]                     # Dormancy intercept (baseline)
orchid[14] <- fit.dorm["beta.D"]                     # Dormancy slope (size)
orchid[15] <- Dsize                           # Mean size of plants emerging from dormancy
orchid[16] <- Dsd                             # SD   size of plants emerging from dormancy
# [60] Probability of Flowering
orchid[17] <- fit.flower["mu.intF"]                       # Flowering intercept (baseline)
orchid[18] <- fit.flower["beta.F"]                       # Dormancy slope (size)
# [70] Number of Flowers
orchid[19] <- fit.fert["mu.intR"]                     # Fertility intercept (baseline)
orchid[20] <- fit.fert["beta.R"]                     # Fertility slope (size)
# [80] Recruits
orchid[21] <- seedlM                     # Mean seedling size (constant)
orchid[22] <- seedlSD                      # Residual SD seedling size
orchid[23] <- eta                              # probability seed becomes protocorm

# Random effects
orchid[24] <- fit.grow["sigma.intG"]
orchid[25] <- fit.surv["sigma.intS"]
orchid[26] <- fit.flower["sigma.intF"]
orchid[27] <- fit.fert["sigma.intR"]
orchid[28] <- fit.fruit["sigma.intP"]
orchid[29] <- fit.dorm["sigma.intD"]

orchid<-as.numeric(orchid)

#FUNCTION: variance-covariance matrix ########################
#This function is used in the stochastic simulation script
varCovar=function(corMat,varVec){
  out = corMat * (varVec %*% t(varVec))
  return(out)
}
varVec=orchid[c(24:29)]

##3.Define demographic functions ################################################################################

# Inverse logic function
lgt <- function(x){ 
  y <- exp(x)/(1+exp(x))
  return(y)
}

# PRODUCTION OF PROTOCORMS BY Y-SIZED INDIVIDUALS
fx <- function(x,params,rfx){
  xb=pmin(pmax(x,params[1]),params[2])
  pflow <- lgt(params[17] + params[18]*xb + rfx[3])
  nflow <- exp(params[19] + params[20]*xb + rfx[4])
  nseed <- pflow*nflow*lgt(params[11]+rfx[5])*params[12] #params[12] = seed/fruit
  return(nseed)
}

## Growth
gxy <- function(x,y,params,rfx){
  xb=pmin(pmax(x,params[1]),params[2])
  return(dnorm(x=y,mean=c(params[4] + params[5]*xb + rfx[1]),sd=params[6]))
}

# SURVIVAL AT SIZE X
sx <- function(x,params,rfx){
  xb=pmin(pmax(x,params[1]),params[2])
  return(lgt(params[7] + params[8]*xb + rfx[2]))
  #return(1) # test 100% survival
}

# PROBABILITY OF DORMANCY AT SIZE X
dx <- function(x,params,rfx){
  xb=pmin(pmax(x,params[1]),params[2])
  return(lgt(params[13] + params[14]*xb + rfx[6]))
  #return(0) # test no plants going dormant
}

# (SURVIVAL) * (GROWTH) * (NOT GOING DORMANT)
pxy <- function(x,y,params,rfx){
  xb=pmin(pmax(x,params[1]),params[2])
  return(sx(xb,params,rfx)*(1-dx(xb,params,rfx))*gxy(xb,y,params,rfx))
}

# SIZE DISTRIBUTION OF PLANTS THAT EMERGE FROM TUBERS
recruits <- function(y,params){
  return(dnorm(x=y, mean=params[21], sd=params[22]))
}

# SIZE DISTRIBUTION OF PLANTS THAT EMERGE FROM DORMANCY
wakeup <- function(y,params){
  return(dnorm(x=y,mean=params[15],sd=params[16]))
}


###4.Simulation parameters######################################################
max.yrs       <- 50000
matsize       <- 200
extra.grid    <- 3        #How many row/col add to the part of the matrix referring to sizes?
random        <- T     
varCorr       <- vCor
allparms      <- orchid

#Eviction parameters
floor.extend   <- 0.14
ceiling.extend <- 0.71
lower          <- orchid[1] - floor.extend
upper          <- orchid[2] + ceiling.extend


###5. Construct Integral Projection Model (IPM)##########################################################
bigmatrix <- function(params,random,sigma,lower,upper,rand.seed){

  # params: Vital rate parameters
  # random: Use random effects? TRUE or FALSE
  # sigma: variance-covariance matrix for simulations
  # lower,upper: lower and upper bounds of size distribution (L and U in Eq. 8 of the main text)
  # rand.seed: sequence of values to use for random number generation during simulations
  #
  # The matrix looks like this:
  # [P= protocorms, T = tubers, D = dormant, G = growth]
  # [# = number, p = probability, f() = function]
  #
  #       1    2     3    4  ... N+4
  #
  #     _ P    T     D   G1  ... GN  _  <- FROM
  # P   | 0    0     0   #1S ... #NS |     1
  # T   | f(S) 0     0   0   ... 0   |     2
  # D   | 0    0     0   p1D ... pND |     4
  # G1  | 0    pR1  pD1  p11 ... pN1 |     5
  # ... | ..   ...  ...  ... ... ... |    ...
  # GN  | 0    pRN  pDN  p1N ... pNN |    N+4
  #     -                                 -
  # ^- TO
  
  ### DEFINE VARIABLES USED TO CREATE THE IPM ----------------------------------
  L          <- lower                    # Lower integration limit
  U          <- upper                    # Upper integration limit
  n          <- params[3]                # Matrix size
  h          <- (U-L)/n                  # Bin size (for size classes)
  b          <- L+c(0:n)*h               # Vector of lower boundaries to bins
  y          <- 0.5*(b[1:n]+b[2:(n+1)])  # Vector of midpoints to bins
  
  # Setup random effects 
  rfx    <- matrix(0,1,6)
  if(random==T){
    set.seed(rand.seed)
    rfx[1,1:6] = rmvnorm(n=1, mean=rep(0,6), sigma=sigma)  
  }
  
  ### Construct the matrix of reproductive transitions -------------------------
  
  # Initiate reproduction matrix
  Fmat <- matrix(0,(n+3),(n+3))
  
  # Production of new seeds goes in the top row
  Fmat[1,4:(n+3)] <- fx(y,params,rfx)
  
  ### Construct the matrix of non-reproductive transitions (growth, survival) --
  
  # Initiate transition matrix
  Tmat <- matrix(0,(n+3),(n+3))
  
  # Probability of dormancy (conditioned on survival)
  Tmat[3,4:(n+3)] <- sx(y,params,rfx) * dx(y,params,rfx)
  
  # Distribution of plants emerging from dormancy
  Tmat[4:(n+3),3] <- wakeup(y,params) * h
  
  # Growth function applied to the inner matrix
  Tmat[4:(n+3),4:(n+3)] <- t(outer(y,y,pxy,params=params,rfx=rfx)) * h

  #Tuber production: 'prob seed->protocorm' * 'prot. surv' * 'tuber surv' 
  Tmat[2,1] <- params[23] * params[9] * params[10]
 
  # Size distribution of plants (recruits) that emerge from tubers
  # survival has already been accounted for.
  Tmat[4:(n+3),2] <- recruits(y,params) * h
  
  ### Construct projection matrix ----------------------------------------------
  # The matrix is just the sum of the transition matrix and the reproduction 
  # matrix.  This is the discrete approximation to the continuous IPM kernel.
  IPM <- Tmat + Fmat
  
  ### Output everything --------------------------------------------------------
  return(list(IPMmat=IPM, Tmat=Tmat, Fmat=Fmat, rfx=rfx))
}

#Calculate deterministic lambda
#lambda<-Re(eigen(bigmatrix(params=orchid,random=F,
#lower=lower,upper=upper,sigma=vCor,rand.seed=42)$IPM)$values[1]);lambda