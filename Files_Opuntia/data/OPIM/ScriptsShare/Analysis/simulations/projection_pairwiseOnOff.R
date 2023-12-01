#File to run stochastic simulations of population projection models, 
#where each simulations is run with one of the pairwise correlations turned set to zero,
#or, in other words, "turned off"  
#'species': which of the three population projectio models should simulations refer to?
#Note that each model has a different number of pairwise correlations (10 for helianthella,
#6 for opuntia, 15 for orchis). 
projection_pairwiseOnOff=function(species){
  
  #Variables common to all simulations
  res.ipm     <- matrix(0, matsize + extra.grid, matsize + extra.grid)
  
  #Variables common to all simulations
  res.ipm       <- matrix(0, nrow=matsize + extra.grid, ncol=matsize + extra.grid)
  variances     <- rep(0,sum(upper.tri(varCorr)))  #store var(lambda_t) results here
  zeroCorInd    <- which(upper.tri(varCorr),arr.ind=T)
  
  #Progress bar 
  pb.max <- max.yrs                 #Total number of iterations in simulation
  pb <- txtProgressBar(style=3, min=0, initial=0, max=pb.max) 

  #simulate and store lambdas
  for(ind in 1:length(lambdas)){
    
    rtracker      <- rep(0,max.yrs)
    n0            <- rep(1/(matsize + extra.grid),matsize + extra.grid)
    pb.tally      <- 0        #Re-set tally count
    ptm           <- proc.time()   #track length of simulations
    
    #modify correlation
    modCorr=varCorr
    modCorr[zeroCorInd[ind,1],zeroCorInd[ind,2]]=modCorr[zeroCorInd[ind,2],zeroCorInd[ind,1]]=0
    modCorr=as.matrix(forceSymmetric(as.matrix(modCorr)))
      
    if(species=="opuntia" | species=="orchis"){
      for(g in 1:max.yrs){ #Start loop
        
        #Store kernel
        res.ipm[,]     <- bigmatrix(params=allparms, random=random, sigma=varCovar(modCorr,varVec), 
                                    lower=lower, upper=upper,rand.seed=rseed.vec[g])$IPM
        
        n0 <- res.ipm[,] %*% n0
        N  <- sum(n0)
        rtracker[g]<-log(N)
        n0 <-n0/N
        pb.tally       <- pb.tally + 1; setTxtProgressBar(pb, pb.tally) # Update the progress bar
      }
    }
    if(species=="helianthella"){
      for(g in 1:max.yrs){ #Start loop
        
        #Store kernel
        res.ipm[,]     <- bigmatrix(params=allparms, upDens=upDens, random=random, sigma=varCovar(modCorr,varVec), 
                                    rand.seed=rseed.vec[g])$IPM
        
        n0 <- res.ipm[,] %*% n0
        N  <- sum(n0)
        rtracker[g]<-log(N)
        n0 <-n0/N
        pb.tally       <- pb.tally + 1; setTxtProgressBar(pb, pb.tally) # Update the progress bar
      }
    }
    burnin=max.yrs*0.1
    variances[ind]=var(exp(rtracker[burnin:max.yrs]))
    print(proc.time() - ptm)
      
  }
  
  return(variances)
  
}