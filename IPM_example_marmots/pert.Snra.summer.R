summer_pertSnra_ipm <- function(year,env,env.pert){
  
  IPMja=array(0,c(n*n.stage,n*n.stage)) # summer
  year = year
  Q.sim=env
  Q.sim.pert=env.pert
 
  
  # SURVIVAL 
  S2y <- diag(S2.fun(z,stage=1,year,Q.sim),nrow=n,ncol=n) # Survival yearlings
  S2nra <- diag(S2.fun(z,stage=2,year,Q.sim.pert),nrow=n,ncol=n) # Survival non-reproductive adults 
  S2ra <- diag(S2.fun(z,stage=3,year,Q.sim),nrow=n,ncol=n) # Survival reproductive adults 
  
  # GROWTH 
  
  # Yearlings 
  Gy <- h*t(outer(z,z,GR2.fun,stage=1,year,Q.sim)) 
  
  # Reproductive Adults 
  Gnra <- h*t(outer(z,z,GR2.fun,stage=2,year,Q.sim)) 
  
  # Non-Reproductive Adults 
  Gra <- h*t(outer(z,z,GR2.fun,stage=3,year,Q.sim)) 
  
  # Control for eviction:
  # this is equivalent to redistributing evictd sizes evenly among existing size classes 
  
  Gy=Gy/matrix(as.vector(apply(Gy,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  Gra=Gra/matrix(as.vector(apply(Gra,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  Gnra=Gnra/matrix(as.vector(apply(Gnra,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  
  # Offspring mass 
  
  OffMass <- h*t(outer(z,z,OffMass.fun,year,Q.sim)) 
  # Control for eviction:
  # this is equivalent to redistributing evictd sizes evenly among existing size classes 
  OffMass=OffMass/matrix(as.vector(apply(OffMass,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  # recruitment (as individuals that die in August still produce offspring, recruitment does not need to be multiplied by S2)
  
  Rec=(diag(R.fun(z,year,Q.sim)))/2
  
  Fkernel <- as.matrix(OffMass%*%(Rec*S2ra))
  
  # FILL IPM
  
  IPMja[(n+1):(2*n),(n+1):(2*n)]=Gy%*%S2y # Yearling stay Yearling
  IPMja[(2*n+1):(3*n),(2*n+1):(3*n)]=Gnra%*%S2nra # Non-Reproductive Adults stay Non-Reproductive Adults
  IPMja[(3*n+1):(4*n),(3*n+1):(4*n)]=Gra%*%S2ra # Reproductive Adults stay Reproductive Adults
  
  IPMja[1:n,(3*n+1):(4*n)]=Fkernel # Adults producing juveniles
  
  return(IPMja)
}