winter_pertTra_ipm <- function(year,env,env.pert){
  
  IPMaj=array(0,c(n*n.stage,n*n.stage)) # winter
  
  ### AUGUST - JUNE 
  
  year = year
  Q.sim=env
  Q.sim.pert=env.pert
  
  Sj <- diag(S1.fun(z,stage=1,year,Q.sim)) # Survival juveniles
  Sy <- diag(S1.fun(z,stage=2,year,Q.sim)) # Survival yearlings
  Snra <- diag(S1.fun(z,stage=3,year,Q.sim)) # Survival non-reproductive adults 
  Sra <- diag(S1.fun(z,stage=4,year,Q.sim)) # Survival reproductive adults 
  
  # Transition To RA or NRA 
  Ty <- diag(PR.fun(z,stage=1,year,Q.sim))
  Tnra <- diag(PR.fun(z,stage=2,year,Q.sim))
  Tra <- diag(PR.fun(z,stage=3,year,Q.sim.pert))
  
  # Growth - stage specific like for survival
  G <- h*t(outer(z,z,GR1.fun,stage=1,year,Q.sim)) 
  # Control for eviction:
  # this is equivalent to redistributing evicted sizes evenly among existing size classes 
  G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # FILL IPM
  IPMaj[(n+1):(2*n),1:n]=G%*%Sj # Juvenile to Yearling
  
  G <- h*t(outer(z,z,GR1.fun,stage=2,year,Q.sim)) 
  # Control for eviction:
  # this is equivalent to redistributing evicted sizes evenly among existing size classes 
  G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # FILL IPM
  IPMaj[(2*n+1):(3*n),(n+1):(2*n)]=G%*%(Sy*diag(1-PR.fun(z,stage=1,year,Q.sim))) # Yearling to Non-Reproductive Adult
  IPMaj[(3*n+1):(4*n),(n+1):(2*n)]=G%*%(Sy*Ty) # Yearling to Reproductive Adult
  
  
  G <- h*t(outer(z,z,GR1.fun,stage=3,year,Q.sim)) 
  G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # FILL IPM
  IPMaj[(2*n+1):(3*n),(2*n+1):(3*n)]=G%*%(Snra*diag(1-PR.fun(z,stage=2,year,Q.sim))) # Non-Reproductive Adult to Non-Reproductive Adult
  IPMaj[(3*n+1):(4*n),(2*n+1):(3*n)]=G%*%(Snra*Tnra) # Non-Reproductive Adult to Reproductive Adult
  
  G <- h*t(outer(z,z,GR1.fun,stage=4,year,Q.sim)) 
  G=G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # FILL IPM
  IPMaj[(2*n+1):(3*n),(3*n+1):(4*n)]=G%*%(Sra*diag(1-PR.fun(z,stage=3,year,Q.sim))) # Reproductive Adult to non-Reproductive Adult
  IPMaj[(3*n+1):(4*n),(3*n+1):(4*n)]=G%*%(Sra*Tra) # Reproductive Adult to Reproductive Adult
  
  return(IPMaj)
}