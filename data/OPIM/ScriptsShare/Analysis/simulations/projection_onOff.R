#File to run stochastic simulations of population projection models.
#'correlations': should vital rate correlations be all "on" or all "off", 
#'#              whereby off means all correlations are set to 0
#'species': run one of the three projection models; 'helianthella', 'opuntia', or 'orchis' 
projection_onOff=function(correlations,species){
  
  #Variables common to all simulations
  res.ipm     <- matrix(0, matsize + extra.grid, matsize + extra.grid)

  #Progress bar
  pb.max <- max.yrs                 #Total number of iterations in simulation
  pb <- txtProgressBar(style=3, min=0, initial=0, max=pb.max) 

  #simulate and store lambdas
  rtracker      <- rep(0,max.yrs)
  n0            <- rep(1/(matsize + extra.grid),matsize + extra.grid)
  pb.tally      <- 0        #Re-set tally count
  ptm           <- proc.time()   #track length of simulations
  
  if(species=="helianthella"){
    if(correlations == "on")  { corr=varCorr }
    if(correlations == "off") { corr=diag(5) }
    for(g in 1:max.yrs){ #Start loop
      
      #Store kernel
      res.ipm[,]     <- bigmatrix(params=allparms, upDens=upDens, random=random, 
                                  sigma=varCovar(corr,varVec), rand.seed=rseed.vec[g])$IPM
      
      n0 <- res.ipm[,] %*% n0
      N  <- sum(n0)
      rtracker[g]<-log(N)
      n0 <-n0/N
      
      pb.tally       <- pb.tally + 1; setTxtProgressBar(pb, pb.tally) # Update the progress bar
    }
  }
  if(species=="opuntia"){
    if(correlations == "on")  { corr=varCorr }
    if(correlations == "off") { corr=diag(4) }
    for(g in 1:max.yrs){ #Start loop
      
      #Store kernel
      res.ipm[,]     <- bigmatrix(params=allparms, random=random, sigma=varCovar(corr,varVec), 
                                  lower=lower, upper=upper,rand.seed=rseed.vec[g])$IPMmat
      
      n0 <- res.ipm[,] %*% n0
      N  <- sum(n0)
      rtracker[g]<-log(N)
      n0 <-n0/N
      
      pb.tally       <- pb.tally + 1; setTxtProgressBar(pb, pb.tally) # Update the progress bar
    }
  }
  if(species=="orchis"){
    if(correlations == "on")  { corr=varCorr }
    if(correlations == "off") { corr=diag(6) }
    for(g in 1:max.yrs){ #Start loop
      
      #Store kernel
      res.ipm[,]     <- bigmatrix(params=allparms, random=random, sigma=varCovar(corr,varVec), 
                                    lower=lower, upper=upper,rand.seed=rseed.vec[g])$IPM
      
      n0 <- res.ipm[,] %*% n0
      N  <- sum(n0)
      rtracker[g]<-log(N)
      n0 <-n0/N
      
      pb.tally       <- pb.tally + 1; setTxtProgressBar(pb, pb.tally) # Update the progress bar
    }
  }
  #discard initial values (to get rid of transient)
  burnin    <- round(max.yrs*0.1)
  rtracker  <- rtracker[-c(1:burnin)]

  #Finish and return
  print(proc.time() - ptm)
  return(rtracker)
 
}