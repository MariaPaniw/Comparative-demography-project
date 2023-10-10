
# S0 ----
pert.S0 <- function(reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov,
                       year = 2008, HuntingLevel, i,
                       randomEff,
                       epsilon.Psi = 0,
                       epsilon.rho = 0,
                       epsilon.m0 = 0,
                       epsilon.mH = 0,
                       epsilon.mO = 0
){
  
  ## Define random effect levels based on argument "randomEff"
  if(randomEff == "match"){
    yearIdx <- year - 1997 + 1
    epsilon.Psi.use <- median(MCMC.mat[, paste0("epsilon.Psi[", yearIdx+1, "]")])
    epsilon.rho.use <- median(MCMC.mat[, paste0("epsilon.rho[", yearIdx+1, "]")])
    epsilon.m0.use <- median(MCMC.mat[, paste0("epsilon.m0[", yearIdx+1, "]")])
    epsilon.mH.use <- median(MCMC.mat[, paste0("epsilon.mH[", yearIdx, "]")])
    epsilon.mO.use <- median(MCMC.mat[, paste0("epsilon.mO[", yearIdx, "]")])
  }
  
  if(randomEff == "ignore"){
    epsilon.Psi.use <- 0
    epsilon.rho.use <- 0
    epsilon.m0.use <- 0
    epsilon.mH.use <- 0
    epsilon.mO.use <- 0
  }
  
  if(randomEff == "define"){
    epsilon.Psi.use <- epsilon.Psi
    epsilon.rho.use <- epsilon.rho
    epsilon.m0.use <- epsilon.m0
    epsilon.mH.use <- epsilon.mH
    epsilon.mO.use <- epsilon.mO
  }
  
  if(!(randomEff %in% c("match", "ignore", "define"))){
    stop("Invalid randomEff option provided. The accepted options are match, ignore, and define.")
  }
  
  
  ## Predict vital rates
  
  # Breeding probability and litter size per age class
  Psi <- rho <- rep(NA, 5)
  for(a in 2:5){
    Psi[a] <- Psi.predict(par.a = MCMC.mat[,"par.a"][i],
                          par.b = MCMC.mat[,"par.b"][i],
                          par.c = MCMC.mat[,"par.c"][i],
                          betaRC.Psi = MCMC.mat[,"betaRC.Psi"][i],
                          betaSI.Psi = MCMC.mat[,"betaSI.Psi"][i],
                          betaY.Psi = MCMC.mat[,"betaY.Psi"][i],
                          epsilon.Psi = epsilon.Psi.use,
                          age = a,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
    
    rho[a] <- rho.predict(mean.rho = MCMC.mat[,"mean.rho"][i],
                          a.eff1 = MCMC.mat[,"a.eff1"][i],
                          betaRC.rho = MCMC.mat[,"betaRC.rho"][i],
                          betaSI.rho = MCMC.mat[,"betaSI.rho"][i],
                          betaY.rho = MCMC.mat[,"betaY.rho"][i],
                          age = a,
                          epsilon.rho = epsilon.rho.use,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
  }
  
  # Denning survival
  S0 <- S0.predict(S0 = MCMC.mat[,"S0"][i],
                   betaRC.m0 = MCMC.mat[,"betaRC.m0"][i], 
                   betaSI.m0 = MCMC.mat[,"betaSI.m0"][i],
                   betaY.m0 = MCMC.mat[,"betaY.m0"][i], 
                   epsilon.m0 = epsilon.m0.use,
                   reindeerCov = pert.reindeerCov, seaIceCov = pert.seaIceCov,
                   year = year+1)
  
  # Juvenile and adult annual survival
  S <- S.predict(Mu.mH_j = MCMC.mat[,"Mu.mH[2]"][i], 
                 Mu.mH_a = MCMC.mat[,"Mu.mH[1]"][i], 
                 Mu.mO_j = MCMC.mat[,"Mu.mO[2]"][i], 
                 Mu.mO_a = MCMC.mat[,"Mu.mO[1]"][i],
                 betaRC.mO = MCMC.mat[,"betaRC.mO"][i], 
                 betaSI.mO = MCMC.mat[,"betaSI.mO"][i], 
                 betaG.mO = MCMC.mat[,"betaG.mO"][i], 
                 betaHP.mH = MCMC.mat[,"betaHP.mH"][i], 
                 betaY.mO = MCMC.mat[,"betaY.mO"][i], 
                 epsilon.mH = epsilon.mH.use, epsilon.mO = epsilon.mO.use,
                 reindeerCov = reindeerCov, seaIceCov = seaIceCov, gooseCov = gooseCov,
                 year = year, HuntingLevel = HuntingLevel)
  
  # List all age-specific vital rates
  VitalRates <- list(S0 = S0, Sj = S[1], Sa = S[2], 
                     Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                     rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Assemble projection matrix
  MPM <- assemble.arcticFoxMPM(S0 = S0, Sj = S[1], Sa = S[2], 
                               Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                               rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Return vital rates and MPM
  # return(list(VitalRates = VitalRates, MPM = MPM))
  return(MPM = MPM)
}

# Sj ----
S.predict.Sj <- function(Mu.mH_j, Mu.mH_a, Mu.mO_j, Mu.mO_a,
                            betaRC.mO, betaSI.mO, betaG.mO, 
                            betaHP.mH, betaY.mO, 
                            epsilon.mO = 0, epsilon.mH = 0,
                            reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov, pert.gooseCov,
                            year = 2008, HuntingLevel){
  
  # Determine year index
  yearIdx <- year - 1997 + 1
  
  # Make variable HPeriod based on "HuntingLevel"
  if(HuntingLevel %in% 1:2){
    HPeriod <- HuntingLevel - 1
  }else{
    stop("Invalid HuntingLevel provided. Set to either 1 (first period, higher hunting pressure) or 2 (second period, lower hunting pressure).")
  }
  
  # Calculate linear predictors on link scale (log hazard rates)
  
  # Hunting mortality
  log.mH <- c(log(Mu.mH_j) + betaHP.mH*HPeriod + epsilon.mH, # Juveniles
              log(Mu.mH_a) + betaHP.mH*HPeriod + epsilon.mH) # Adults
  
  # Natural mortality
  log.mO <- c(log(Mu.mO_j) + betaY.mO*(yearIdx-12) + betaRC.mO*pert.reindeerCov + betaG.mO*pert.gooseCov + betaSI.mO*pert.seaIceCov + epsilon.mO, # Juveniles
              log(Mu.mO_a) + betaY.mO*(yearIdx-12) + betaRC.mO*reindeerCov + betaSI.mO*seaIceCov + epsilon.mO) # Adults
  
  
  # Back-calculate & return
  mH <- exp(log.mH) # Hunting mortality hazard rate
  mO <- exp(log.mO) # Natural mortality hazard rate
  S <- exp(-(mH + mO)) # Survival probability
  names(S) <- c("Juveniles", "Adults")
  
  return(S)
}


pert.Sj <- function(reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov, pert.gooseCov,
                       year = 2008, HuntingLevel, i,
                       randomEff,
                       epsilon.Psi = 0,
                       epsilon.rho = 0,
                       epsilon.m0 = 0,
                       epsilon.mH = 0,
                       epsilon.mO = 0
){
  
  ## Define random effect levels based on argument "randomEff"
  if(randomEff == "match"){
    yearIdx <- year - 1997 + 1
    epsilon.Psi.use <- median(MCMC.mat[, paste0("epsilon.Psi[", yearIdx+1, "]")])
    epsilon.rho.use <- median(MCMC.mat[, paste0("epsilon.rho[", yearIdx+1, "]")])
    epsilon.m0.use <- median(MCMC.mat[, paste0("epsilon.m0[", yearIdx+1, "]")])
    epsilon.mH.use <- median(MCMC.mat[, paste0("epsilon.mH[", yearIdx, "]")])
    epsilon.mO.use <- median(MCMC.mat[, paste0("epsilon.mO[", yearIdx, "]")])
  }
  
  if(randomEff == "ignore"){
    epsilon.Psi.use <- 0
    epsilon.rho.use <- 0
    epsilon.m0.use <- 0
    epsilon.mH.use <- 0
    epsilon.mO.use <- 0
  }
  
  if(randomEff == "define"){
    epsilon.Psi.use <- epsilon.Psi
    epsilon.rho.use <- epsilon.rho
    epsilon.m0.use <- epsilon.m0
    epsilon.mH.use <- epsilon.mH
    epsilon.mO.use <- epsilon.mO
  }
  
  if(!(randomEff %in% c("match", "ignore", "define"))){
    stop("Invalid randomEff option provided. The accepted options are match, ignore, and define.")
  }
  
  
  ## Predict vital rates
  
  # Breeding probability and litter size per age class
  Psi <- rho <- rep(NA, 5)
  for(a in 2:5){
    Psi[a] <- Psi.predict(par.a = MCMC.mat[,"par.a"][i],
                          par.b = MCMC.mat[,"par.b"][i],
                          par.c = MCMC.mat[,"par.c"][i],
                          betaRC.Psi = MCMC.mat[,"betaRC.Psi"][i],
                          betaSI.Psi = MCMC.mat[,"betaSI.Psi"][i],
                          betaY.Psi = MCMC.mat[,"betaY.Psi"][i],
                          epsilon.Psi = epsilon.Psi.use,
                          age = a,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
    
    rho[a] <- rho.predict(mean.rho = MCMC.mat[,"mean.rho"][i],
                          a.eff1 = MCMC.mat[,"a.eff1"][i],
                          betaRC.rho = MCMC.mat[,"betaRC.rho"][i],
                          betaSI.rho = MCMC.mat[,"betaSI.rho"][i],
                          betaY.rho = MCMC.mat[,"betaY.rho"][i],
                          age = a,
                          epsilon.rho = epsilon.rho.use,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
  }
  
  # Denning survival
  S0 <- S0.predict(S0 = MCMC.mat[,"S0"][i],
                   betaRC.m0 = MCMC.mat[,"betaRC.m0"][i], 
                   betaSI.m0 = MCMC.mat[,"betaSI.m0"][i],
                   betaY.m0 = MCMC.mat[,"betaY.m0"][i], 
                   epsilon.m0 = epsilon.m0.use,
                   reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                   year = year+1)
  
  # Juvenile and adult annual survival
  S <- S.predict.Sj(Mu.mH_j = MCMC.mat[,"Mu.mH[2]"][i], 
                    Mu.mH_a = MCMC.mat[,"Mu.mH[1]"][i], 
                    Mu.mO_j = MCMC.mat[,"Mu.mO[2]"][i], 
                    Mu.mO_a = MCMC.mat[,"Mu.mO[1]"][i],
                    betaRC.mO = MCMC.mat[,"betaRC.mO"][i], 
                    betaSI.mO = MCMC.mat[,"betaSI.mO"][i], 
                    betaG.mO = MCMC.mat[,"betaG.mO"][i], 
                    betaHP.mH = MCMC.mat[,"betaHP.mH"][i], 
                    betaY.mO = MCMC.mat[,"betaY.mO"][i], 
                       epsilon.mH = epsilon.mH.use, epsilon.mO = epsilon.mO.use,
                       reindeerCov = reindeerCov, seaIceCov = seaIceCov, gooseCov = gooseCov, pert.seaIceCov = pert.seaIceCov, pert.reindeerCov = pert.reindeerCov, pert.gooseCov = pert.gooseCov,
                       year = year, HuntingLevel = HuntingLevel)
  
  # List all age-specific vital rates
  VitalRates <- list(S0 = S0, Sj = S[1], Sa = S[2], 
                     Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                     rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Assemble projection matrix
  MPM <- assemble.arcticFoxMPM(S0 = S0, Sj = S[1], Sa = S[2], 
                               Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                               rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Return vital rates and MPM
  # return(list(VitalRates = VitalRates, MPM = MPM))
  return(MPM = MPM)
}

# Sa ----
S.predict.Sa <- function(Mu.mH_j, Mu.mH_a, Mu.mO_j, Mu.mO_a,
                            betaRC.mO, betaSI.mO, betaG.mO, 
                            betaHP.mH, betaY.mO, 
                            epsilon.mO = 0, epsilon.mH = 0,
                            reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov,
                            year = 2008, HuntingLevel){
  
  # Determine year index
  yearIdx <- year - 1997 + 1
  
  # Make variable HPeriod based on "HuntingLevel"
  if(HuntingLevel %in% 1:2){
    HPeriod <- HuntingLevel - 1
  }else{
    stop("Invalid HuntingLevel provided. Set to either 1 (first period, higher hunting pressure) or 2 (second period, lower hunting pressure).")
  }
  
  # Calculate linear predictors on link scale (log hazard rates)
  
  # Hunting mortality
  log.mH <- c(log(Mu.mH_j) + betaHP.mH*HPeriod + epsilon.mH, # Juveniles
              log(Mu.mH_a) + betaHP.mH*HPeriod + epsilon.mH) # Adults
  
  # Natural mortality
  log.mO <- c(log(Mu.mO_j) + betaY.mO*(yearIdx-12) + betaRC.mO*reindeerCov + betaG.mO*gooseCov + betaSI.mO*seaIceCov + epsilon.mO, # Juveniles
              log(Mu.mO_a) + betaY.mO*(yearIdx-12) + betaRC.mO*pert.reindeerCov + betaSI.mO*pert.seaIceCov + epsilon.mO) # Adults
  
  
  # Back-calculate & return
  mH <- exp(log.mH) # Hunting mortality hazard rate
  mO <- exp(log.mO) # Natural mortality hazard rate
  S <- exp(-(mH + mO)) # Survival probability
  names(S) <- c("Juveniles", "Adults")
  
  return(S)
}


pert.Sa <- function(reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov,
                       year = 2008, HuntingLevel, i,
                       randomEff,
                       epsilon.Psi = 0,
                       epsilon.rho = 0,
                       epsilon.m0 = 0,
                       epsilon.mH = 0,
                       epsilon.mO = 0
){
  
  ## Define random effect levels based on argument "randomEff"
  if(randomEff == "match"){
    yearIdx <- year - 1997 + 1
    epsilon.Psi.use <- median(MCMC.mat[, paste0("epsilon.Psi[", yearIdx+1, "]")])
    epsilon.rho.use <- median(MCMC.mat[, paste0("epsilon.rho[", yearIdx+1, "]")])
    epsilon.m0.use <- median(MCMC.mat[, paste0("epsilon.m0[", yearIdx+1, "]")])
    epsilon.mH.use <- median(MCMC.mat[, paste0("epsilon.mH[", yearIdx, "]")])
    epsilon.mO.use <- median(MCMC.mat[, paste0("epsilon.mO[", yearIdx, "]")])
  }
  
  if(randomEff == "ignore"){
    epsilon.Psi.use <- 0
    epsilon.rho.use <- 0
    epsilon.m0.use <- 0
    epsilon.mH.use <- 0
    epsilon.mO.use <- 0
  }
  
  if(randomEff == "define"){
    epsilon.Psi.use <- epsilon.Psi
    epsilon.rho.use <- epsilon.rho
    epsilon.m0.use <- epsilon.m0
    epsilon.mH.use <- epsilon.mH
    epsilon.mO.use <- epsilon.mO
  }
  
  if(!(randomEff %in% c("match", "ignore", "define"))){
    stop("Invalid randomEff option provided. The accepted options are match, ignore, and define.")
  }
  
  
  ## Predict vital rates
  
  # Breeding probability and litter size per age class
  Psi <- rho <- rep(NA, 5)
  for(a in 2:5){
    Psi[a] <- Psi.predict(par.a = MCMC.mat[,"par.a"][i],
                          par.b = MCMC.mat[,"par.b"][i],
                          par.c = MCMC.mat[,"par.c"][i],
                          betaRC.Psi = MCMC.mat[,"betaRC.Psi"][i],
                          betaSI.Psi = MCMC.mat[,"betaSI.Psi"][i],
                          betaY.Psi = MCMC.mat[,"betaY.Psi"][i],
                          epsilon.Psi = epsilon.Psi.use,
                          age = a,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
    
    rho[a] <- rho.predict(mean.rho = MCMC.mat[,"mean.rho"][i],
                          a.eff1 = MCMC.mat[,"a.eff1"][i],
                          betaRC.rho = MCMC.mat[,"betaRC.rho"][i],
                          betaSI.rho = MCMC.mat[,"betaSI.rho"][i],
                          betaY.rho = MCMC.mat[,"betaY.rho"][i],
                          age = a,
                          epsilon.rho = epsilon.rho.use,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
  }
  
  # Denning survival
  S0 <- S0.predict(S0 = MCMC.mat[,"S0"][i],
                   betaRC.m0 = MCMC.mat[,"betaRC.m0"][i], 
                   betaSI.m0 = MCMC.mat[,"betaSI.m0"][i],
                   betaY.m0 = MCMC.mat[,"betaY.m0"][i], 
                   epsilon.m0 = epsilon.m0.use,
                   reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                   year = year+1)
  
  # Juvenile and adult annual survival
  S <- S.predict.Sa(Mu.mH_j = MCMC.mat[,"Mu.mH[2]"][i], 
                    Mu.mH_a = MCMC.mat[,"Mu.mH[1]"][i], 
                    Mu.mO_j = MCMC.mat[,"Mu.mO[2]"][i], 
                    Mu.mO_a = MCMC.mat[,"Mu.mO[1]"][i],
                    betaRC.mO = MCMC.mat[,"betaRC.mO"][i], 
                    betaSI.mO = MCMC.mat[,"betaSI.mO"][i], 
                    betaG.mO = MCMC.mat[,"betaG.mO"][i], 
                    betaHP.mH = MCMC.mat[,"betaHP.mH"][i], 
                    betaY.mO = MCMC.mat[,"betaY.mO"][i], 
                       epsilon.mH = epsilon.mH.use, epsilon.mO = epsilon.mO.use,
                       reindeerCov = reindeerCov, seaIceCov = seaIceCov, gooseCov = gooseCov, pert.seaIceCov = pert.seaIceCov, pert.reindeerCov = pert.reindeerCov,
                       year = year, HuntingLevel = HuntingLevel)
  
  # List all age-specific vital rates
  VitalRates <- list(S0 = S0, Sj = S[1], Sa = S[2], 
                     Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                     rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Assemble projection matrix
  MPM <- assemble.arcticFoxMPM(S0 = S0, Sj = S[1], Sa = S[2], 
                               Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                               rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Return vital rates and MPM
  # return(list(VitalRates = VitalRates, MPM = MPM))
  return(MPM = MPM)
}

# Psi ----
pert.Psi <- function(reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov,
                        year = 2008, HuntingLevel, i,
                        randomEff,
                        epsilon.Psi = 0,
                        epsilon.rho = 0,
                        epsilon.m0 = 0,
                        epsilon.mH = 0,
                        epsilon.mO = 0
){
  
  ## Define random effect levels based on argument "randomEff"
  if(randomEff == "match"){
    yearIdx <- year - 1997 + 1
    epsilon.Psi.use <- median(MCMC.mat[, paste0("epsilon.Psi[", yearIdx+1, "]")])
    epsilon.rho.use <- median(MCMC.mat[, paste0("epsilon.rho[", yearIdx+1, "]")])
    epsilon.m0.use <- median(MCMC.mat[, paste0("epsilon.m0[", yearIdx+1, "]")])
    epsilon.mH.use <- median(MCMC.mat[, paste0("epsilon.mH[", yearIdx, "]")])
    epsilon.mO.use <- median(MCMC.mat[, paste0("epsilon.mO[", yearIdx, "]")])
  }
  
  if(randomEff == "ignore"){
    epsilon.Psi.use <- 0
    epsilon.rho.use <- 0
    epsilon.m0.use <- 0
    epsilon.mH.use <- 0
    epsilon.mO.use <- 0
  }
  
  if(randomEff == "define"){
    epsilon.Psi.use <- epsilon.Psi
    epsilon.rho.use <- epsilon.rho
    epsilon.m0.use <- epsilon.m0
    epsilon.mH.use <- epsilon.mH
    epsilon.mO.use <- epsilon.mO
  }
  
  if(!(randomEff %in% c("match", "ignore", "define"))){
    stop("Invalid randomEff option provided. The accepted options are match, ignore, and define.")
  }
  
  
  ## Predict vital rates
  
  # Breeding probability and litter size per age class
  Psi <- rho <- rep(NA, 5)
  for(a in 2:5){
    Psi[a] <- Psi.predict(par.a = MCMC.mat[,"par.a"][i],
                          par.b = MCMC.mat[,"par.b"][i],
                          par.c = MCMC.mat[,"par.c"][i],
                          betaRC.Psi = MCMC.mat[,"betaRC.Psi"][i],
                          betaSI.Psi = MCMC.mat[,"betaSI.Psi"][i],
                          betaY.Psi = MCMC.mat[,"betaY.Psi"][i],
                          epsilon.Psi = epsilon.Psi.use,
                          age = a,
                          reindeerCov = pert.reindeerCov, seaIceCov = pert.seaIceCov,
                          year = year+1) 
    
    rho[a] <- rho.predict(mean.rho = MCMC.mat[,"mean.rho"][i],
                          a.eff1 = MCMC.mat[,"a.eff1"][i],
                          betaRC.rho = MCMC.mat[,"betaRC.rho"][i],
                          betaSI.rho = MCMC.mat[,"betaSI.rho"][i],
                          betaY.rho = MCMC.mat[,"betaY.rho"][i],
                          age = a,
                          epsilon.rho = epsilon.rho.use,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
  }
  
  # Denning survival
  S0 <- S0.predict(S0 = MCMC.mat[,"S0"][i],
                   betaRC.m0 = MCMC.mat[,"betaRC.m0"][i], 
                   betaSI.m0 = MCMC.mat[,"betaSI.m0"][i],
                   betaY.m0 = MCMC.mat[,"betaY.m0"][i], 
                   epsilon.m0 = epsilon.m0.use,
                   reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                   year = year+1)
  
  # Juvenile and adult annual survival
  S <- S.predict(Mu.mH_j = MCMC.mat[,"Mu.mH[2]"][i], 
                 Mu.mH_a = MCMC.mat[,"Mu.mH[1]"][i], 
                 Mu.mO_j = MCMC.mat[,"Mu.mO[2]"][i], 
                 Mu.mO_a = MCMC.mat[,"Mu.mO[1]"][i],
                 betaRC.mO = MCMC.mat[,"betaRC.mO"][i], 
                 betaSI.mO = MCMC.mat[,"betaSI.mO"][i], 
                 betaG.mO = MCMC.mat[,"betaG.mO"][i], 
                 betaHP.mH = MCMC.mat[,"betaHP.mH"][i], 
                 betaY.mO = MCMC.mat[,"betaY.mO"][i], 
                 epsilon.mH = epsilon.mH.use, epsilon.mO = epsilon.mO.use,
                 reindeerCov = reindeerCov, seaIceCov = seaIceCov, gooseCov = gooseCov,
                 year = year, HuntingLevel = HuntingLevel)
  
  # List all age-specific vital rates
  VitalRates <- list(S0 = S0, Sj = S[1], Sa = S[2], 
                     Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                     rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Assemble projection matrix
  MPM <- assemble.arcticFoxMPM(S0 = S0, Sj = S[1], Sa = S[2], 
                               Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                               rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Return vital rates and MPM
  # return(list(VitalRates = VitalRates, MPM = MPM))
  return(MPM = MPM)
}

# Rho ----
pert.Rho <- function(reindeerCov, seaIceCov, gooseCov, pert.seaIceCov, pert.reindeerCov,
                        year = 2008, HuntingLevel, i,
                        randomEff,
                        epsilon.Psi = 0,
                        epsilon.rho = 0,
                        epsilon.m0 = 0,
                        epsilon.mH = 0,
                        epsilon.mO = 0
){
  
  ## Define random effect levels based on argument "randomEff"
  if(randomEff == "match"){
    yearIdx <- year - 1997 + 1
    epsilon.Psi.use <- median(MCMC.mat[, paste0("epsilon.Psi[", yearIdx+1, "]")])
    epsilon.rho.use <- median(MCMC.mat[, paste0("epsilon.rho[", yearIdx+1, "]")])
    epsilon.m0.use <- median(MCMC.mat[, paste0("epsilon.m0[", yearIdx+1, "]")])
    epsilon.mH.use <- median(MCMC.mat[, paste0("epsilon.mH[", yearIdx, "]")])
    epsilon.mO.use <- median(MCMC.mat[, paste0("epsilon.mO[", yearIdx, "]")])
  }
  
  if(randomEff == "ignore"){
    epsilon.Psi.use <- 0
    epsilon.rho.use <- 0
    epsilon.m0.use <- 0
    epsilon.mH.use <- 0
    epsilon.mO.use <- 0
  }
  
  if(randomEff == "define"){
    epsilon.Psi.use <- epsilon.Psi
    epsilon.rho.use <- epsilon.rho
    epsilon.m0.use <- epsilon.m0
    epsilon.mH.use <- epsilon.mH
    epsilon.mO.use <- epsilon.mO
  }
  
  if(!(randomEff %in% c("match", "ignore", "define"))){
    stop("Invalid randomEff option provided. The accepted options are match, ignore, and define.")
  }
  
  
  ## Predict vital rates
  
  # Breeding probability and litter size per age class
  Psi <- rho <- rep(NA, 5)
  for(a in 2:5){
    Psi[a] <- Psi.predict(par.a = MCMC.mat[,"par.a"][i],
                          par.b = MCMC.mat[,"par.b"][i],
                          par.c = MCMC.mat[,"par.c"][i],
                          betaRC.Psi = MCMC.mat[,"betaRC.Psi"][i],
                          betaSI.Psi = MCMC.mat[,"betaSI.Psi"][i],
                          betaY.Psi = MCMC.mat[,"betaY.Psi"][i],
                          epsilon.Psi = epsilon.Psi.use,
                          age = a,
                          reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                          year = year+1) 
    
    rho[a] <- rho.predict(mean.rho = MCMC.mat[,"mean.rho"][i],
                          a.eff1 = MCMC.mat[,"a.eff1"][i],
                          betaRC.rho = MCMC.mat[,"betaRC.rho"][i],
                          betaSI.rho = MCMC.mat[,"betaSI.rho"][i],
                          betaY.rho = MCMC.mat[,"betaY.rho"][i],
                          age = a,
                          epsilon.rho = epsilon.rho.use,
                          reindeerCov = pert.reindeerCov, seaIceCov = pert.seaIceCov,
                          year = year+1) 
  }
  
  # Denning survival
  S0 <- S0.predict(S0 = MCMC.mat[,"S0"][i],
                   betaRC.m0 = MCMC.mat[,"betaRC.m0"][i], 
                   betaSI.m0 = MCMC.mat[,"betaSI.m0"][i],
                   betaY.m0 = MCMC.mat[,"betaY.m0"][i], 
                   epsilon.m0 = epsilon.m0.use,
                   reindeerCov = reindeerCov, seaIceCov = seaIceCov,
                   year = year+1)
  
  # Juvenile and adult annual survival
  S <- S.predict(Mu.mH_j = MCMC.mat[,"Mu.mH[2]"][i], 
                 Mu.mH_a = MCMC.mat[,"Mu.mH[1]"][i], 
                 Mu.mO_j = MCMC.mat[,"Mu.mO[2]"][i], 
                 Mu.mO_a = MCMC.mat[,"Mu.mO[1]"][i],
                 betaRC.mO = MCMC.mat[,"betaRC.mO"][i], 
                 betaSI.mO = MCMC.mat[,"betaSI.mO"][i], 
                 betaG.mO = MCMC.mat[,"betaG.mO"][i], 
                 betaHP.mH = MCMC.mat[,"betaHP.mH"][i], 
                 betaY.mO = MCMC.mat[,"betaY.mO"][i], 
                 epsilon.mH = epsilon.mH.use, epsilon.mO = epsilon.mO.use,
                 reindeerCov = reindeerCov, seaIceCov = seaIceCov, gooseCov = gooseCov,
                 year = year, HuntingLevel = HuntingLevel)
  
  # List all age-specific vital rates
  VitalRates <- list(S0 = S0, Sj = S[1], Sa = S[2], 
                     Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                     rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Assemble projection matrix
  MPM <- assemble.arcticFoxMPM(S0 = S0, Sj = S[1], Sa = S[2], 
                               Psi2 = Psi[2], Psi3 = Psi[3], Psi4 = Psi[4], Psi5 = Psi[5],
                               rho2 = rho[2], rho3 = rho[3], rho4 = rho[4], rho5 = rho[5])
  
  
  ## Return vital rates and MPM
  # return(list(VitalRates = VitalRates, MPM = MPM))
  return(MPM = MPM)
}
