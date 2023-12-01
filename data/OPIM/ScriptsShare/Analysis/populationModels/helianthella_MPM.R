#Code to build a Matrix Projection Model (MPM) for Helianthella quinquenervis
options(stringsAsFactors=T)
#Useful functions
invlogit<-function(x){exp(x)/(1+exp(x))}


###1.Load vital rate model results#############################################################################

####Read in output from Bayesian estimation
#Coefficients
coefs <- read.csv("Results/vitalRates/heli_coefficient_mean.csv")
  
#"Variance-correlation matrix
vCor <- read.csv("Results/vitalRates/heli_correlation_mean.csv")
vCor <- as.matrix(vCor) 

#Posterior(s) of correlations
#a."Regular"
corPost=read.csv("Results/vitalRates/heli_correlation_posterior.csv") ##NEEDS UPDATING -TEX9/9/2015   ##read.csv("Results/vitalRates/heli_corPost_ljk.csv")

#Posterior of coefficients  
coefPost <- read.csv("Results/vitalRates/heli_coefficient_posterior.csv") ##NEEDS UPDATING -TEX9/9/2015   ##read.csv("Results/vitalRates/heli_simulations_ljk.csv")

#Time-varying model parameters
fit.grow=coefs[,c("mu.int.G","beta.G","sigma.int.G","theta.G")] #Growth 
fit.surv=coefs[,c("mu.int.S","beta.S","sigma.int.S")]           #Survival 
fit.flower=coefs[,c("mu.int.F","beta.F","sigma.int.F","pMaxF")] #Probability of Flowering
fit.buds=coefs[,c("mu.int.B","beta.B","sigma.int.B")]           #Fertility
fit.abort=coefs[,c("mu.int.A","sigma.int.A")]                   # Flower-to-fruit

###2. Vital rates estimated 'from scratch'################################################################

### Read in demography data
demog.dat<-read.csv("Data/helianthella_MAXFIELD_1998-2013.csv")

#Seedlings data (removes all adults)
seedling.demog.dat<-subset(demog.dat,demog.dat$seedling==1)  

#demog.dat$year_t1<-demog.dat$year+1 #create column for year_t1 by adding 1 to year
demog.dat$totStalks <- demog.dat$OkstalksT + demog.dat$AbortStalks + demog.dat$BrokeStalks

##Seedling production, survival and recruit size distribution############################################
#use only seedling data
heads.per.stalk<-2.17  			## flowerheads per stalk  ##hard-coded by Amy Iler
seeds.per.head<-64.1				## seeds per flowerhead   ##hard-coded by Amy Iler

#aggregate total seedlings observed by year
seedl=aggregate(seedling ~ year, data=demog.dat, FUN=sum)
##aggregate total OKstalks by year
OkStalk=aggregate(OkstalksT ~ year, data=demog.dat,FUN=sum,na.rm=T)
## combine seed and seedling data
OkStalkVsSeedl=merge(OkStalk,seedl)
##convert total OKstalks to total seeds
OkStalkVsSeedl$total.seeds=floor(seeds.per.head*heads.per.stalk*OkStalkVsSeedl$OkstalksT)

#Total seeds in year t-1 - we compare this to the number of seedlings at time t.
OkStalkVsSeedl$total.seeds.last.year<-NA #
OkStalkVsSeedl$total.seeds.last.year[2:15]<-OkStalkVsSeedl$total.seeds[1:14]

#7th year excluded because no stalks (stalks/seeds = infinite)
establish.dat<-subset(OkStalkVsSeedl,year!=1998)
establish.dat$P.recruit<-establish.dat$seedling / establish.dat$total.seeds.last.year
establish=glm(P.recruit ~ 1,weights=total.seeds.last.year,family="binomial",data=establish.dat)

## Seedling survival
sdlg.surv<-glm(surv~1, family="binomial",data=seedling.demog.dat)

##Seedling size data
recruit.size<-as.numeric(as.character(demog.dat$clumpsT10[demog.dat$seedling==1]))
hist(recruit.size)

## We use an empirical distribution
recruit.distribution<-ecdf(recruit.size)
recruit.distribution(1:45)


###3.Collect MPM parameters################################################################################

### vectorize the vital rates
HEQU.MAX <- c(NA)
##Growth
HEQU.MAX[1]<-fit.grow[,"mu.int.G"]    ## grand mean growth intercept
HEQU.MAX[2]<-fit.grow[,"beta.G"]      ## growth size slope
HEQU.MAX[5]<-fit.grow[,"theta.G"]     ##Growth overdispersion parameter
HEQU.MAX[10]<-fit.grow[,"sigma.int.G"] 
##Survival
HEQU.MAX[11]<-fit.surv[,"mu.int.S"]	   ## grand mean survival intercept
HEQU.MAX[12]<-fit.surv[,"beta.S"]		   ## survival size slope 
HEQU.MAX[20]<-fit.surv[,"sigma.int.S"]
##Flowering
HEQU.MAX[21]<-fit.flower[,"mu.int.F"]	 ## grand mean flowering intercept    
HEQU.MAX[22]<-fit.flower[,"beta.F"]	   ## flowering size slope
HEQU.MAX[30]<-fit.flower[,"sigma.int.F"]  
##Fertility
HEQU.MAX[31]<-fit.buds[,"mu.int.B"]	 ## grand mean fertility intercept  
HEQU.MAX[32]<-fit.buds[,"beta.B"]	   ## fertility slope
HEQU.MAX[35]<-heads.per.stalk           ## mean # heads per OK stalk
HEQU.MAX[36]<-seeds.per.head            ## mean # seeds per flowerhead
HEQU.MAX[40]<-fit.buds[,"sigma.int.B"]
##Stalk viability (previously abortion)
HEQU.MAX[41]<-fit.abort[,"mu.int.A"]     #grand mean, stalk viability intercept   
HEQU.MAX[50]<-fit.abort[,"sigma.int.A"]  #sd of random year effect for stalk viability
##Seeds, Recruitment, and Seedlings
HEQU.MAX[51]<-coef(establish)[1]        ## logit seed to seedling transition probability
HEQU.MAX[52]<-NA                        ## mean seedling size
HEQU.MAX[53]<-NA                        ## dispersion param for seedling size (phi of)
HEQU.MAX[54]<-coef(sdlg.surv)[1]        ## logit seedling survival probability
##Min and max sizes
HEQU.MAX[61]<-1 ## smallest size is one clump
HEQU.MAX[62]<-45 

#FUNCTION: variance-covariance matrix ########################
#This function is used in the stochastic simulation script
varCovar=function(corMat,varVec){
  out = corMat * (varVec %*% t(varVec))
  return(out)
}
varVec=HEQU.MAX[c(10,20,30,40,50)]


###4.Define demographic functions ################################################################################

#Produce the upper cumulative density for the negative binomial.
#This function is necessarly because the outer() function (line 229) would require vectorization to function
upperCumDens=function(x,params){
  out=NULL
  for(i in 1:length(x)){
    tmp=dnbinom(x=c((params[62]+1):10000),mu=exp(params[1] + params[2]*log(x[i])),
                size=params[5],log=F)
    out[i]=sum(tmp)
  }
  return(out)
}

#GROWTH FROM SIZE X TO Y. Function produces a *truncated probability density function* from each:
#1. Size at time t, (the *x* argument), and 2. Each size at time t+1 (the *y* argument).
#To compute truncation we refer to the first equation on page 280 of the Stan Modeling Language User's Guide and Reference Manual (version 2.9.0)
#(https://github.com/stan-dev/stan/releases/download/v2.9.0/stan-reference-2.9.0.pdf)
gxy<-function(x,y,upDens,params,rfx){
  
  #Numerator is the density function between the lower (1) and upper (45) bounds of the size distribution
  num <- dnbinom(x=y,mu=exp(params[1] + params[2]*log(x) + rfx[1]),size=params[5],log=F)
  #Lower cumulative density (referring to value 'y==0') 
  densLower <- dnbinom(x=c(0:(params[61]-1)),mu=exp(params[1] + params[2]*log(x)),
                  size=params[5],log=F)
  #Upper cumulative density function (referring to values '45 < y < 10001')
  denUpper  <- upDens
  #Scaling factor: the probability of values falling outside of the range '1 < y < 45'
  dens <- (1 - (densLower + denUpper))
  #Truncated density distribution
  out <- num / dens 
  return(out)
  
}

#SURVIVAL AT SIZE X. (This calculates the survival prob. of each size class)
sx<-function(x,params,rfx){
  return(invlogit(params[11] + params[12]*log(x) + rfx[2]))
}

#SURVIVAL*GROWTH
pxy<-function(x,y,upDens,params,rfx){
  sx(x,params,rfx) * gxy(x,y,upDens,params,rfx)
}

# PROBABILITY OF FLOWERING
Pfx<-function(x,params,rfx){
  return(invlogit(params[21] + params[22]*log(x) + rfx[3]))
}

#Number of STALKS
Nfx<-function(x,params,rfx){
  return(exp(params[31] + params[32]*log(x) + rfx[4]))
}

#Stalk viability (previously abortion)
Abx<-function(params,rfx){
  return(invlogit(params[41] + rfx[5]))
}

#Fertility--returns number of seedlings
Fertx<-function(x,params,rfx){
  seedlings<-Pfx(x,params,rfx)*Nfx(x,params,rfx)*Abx(params,rfx)*params[35]*params[36]*invlogit(params[51])
  return(seedlings)
}

#Seedling survival
Sdlg.Surv<-function(params){
  return(invlogit(params[54]))
}

#Seedlings size distribution (Size distribution of recruits)
Sdlg.Size<-function(y,params){
  out=c(recruit.distribution(params[61]:params[62])[1],
        recruit.distribution((params[61]+1):params[62])-recruit.distribution(params[61]:(params[62]-1)))[y]
  return(out)
}


###5. Simulation parameters#############################################################
max.yrs       <- 50000
matsize       <- 45
extra.grid    <- 1        #How many row/col add to the part of the matrix referring to sizes?
random        <- T     
allparms      <- HEQU.MAX
upDens        <- upperCumDens(c(HEQU.MAX[61]:HEQU.MAX[62]),HEQU.MAX)
varCorr       <- vCor
  
##6. Construct Matrix Projection Model (MPM)##########################################################
bigmatrix<-function(params,upDens,random=F,sigma=NULL,rand.seed=NULL){  
  
  # params: Vital rate parameters
  # upDens: upper cumulative density distribution of growth function
  # random: Use random effects? TRUE or FALSE
  # sigma: variance-covariance matrix for simulations
  # rand.seed: sequence of values to use for random number generation 
  
  matdim<-params[62]+1         ## bigmatrix dimension is max size + 1 for seedlings
  y<-c(params[61]:params[62])  #vector of sizes (from size 1 to size 45)
   
  # Random effects
  rfx <- matrix(0,1,5)  #Set random effect to 0 (in case random=F,corr=F) 
  if(random==T){
    set.seed(rand.seed)
    rfx[1,1:5] = rmvnorm(n=1, mean=rep(0,5), sigma=sigma)  
  }
 
  # Fertility matrix
  Fmat<-matrix(0,matdim,matdim)
  Fmat[1,2:matdim]<-Fertx(x=y,params=params,rfx=rfx) 
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,matdim,matdim)
  Tmat[2:matdim,2:matdim]<-t(outer(y,y,pxy,upDens,params=params,rfx=rfx)) 
  Tmat[2:matdim,1]<-Sdlg.Size(y=y,params=params) * Sdlg.Surv(params=params)
  
  # Put it all together
  # NOTE: We use "IPMmat" for consistency with other projection models.
  IPMmat<-Tmat+Fmat #sum the Tmat & Fmat to get the whole matrix
  
  return(list(IPMmat=IPMmat, Fmat=Fmat,Tmat=Tmat,rfx=rfx))
}

## Calculate deterministsic lambda 
#lambda<-Re(eigen(bigmatrix(params=HEQU.MAX,upDens=upDens,random=F)$IPMmat)$values[1]);lambda
