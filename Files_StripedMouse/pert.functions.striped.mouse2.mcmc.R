# immature survival
pert.Si = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp, pert_food, pert_dens_F,i){
  
  Si <- S.predict('I',pert_temp,pert_food,pert_dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  #LS <- 2.65
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}


# philo survival
pert.Sp = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_food, pert_temp, pert_dens_F,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',pert_temp,pert_food,pert_dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}


# breeder survival
pert.Sb = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp, pert_food, pert_dens_F,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',pert_temp,pert_food,pert_dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}


# PsiIB
pert.PsiIB = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp, pert_food, pert_dens_F,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',pert_temp,pert_food,pert_dens_F,i) # PsiIB
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}

# PsiPB
pert.PsiPB = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp, pert_food, pert_dens_F,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',pert_temp,pert_food,pert_dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}

# Bp (breeding probability)

pert.Bp = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp_t, pert_food_t, pert_dens_F_t,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = pert_temp_t, food_t = pert_food_t, dens_F_t = pert_dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}


# Litter probability
pert.Lp = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_food_t, pert_dens_F_t,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = pert_food_t, dens_F_t = pert_dens_F_t, i) # litter prob
  LS <- ls(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # litter size
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}



# Litter size
pert.Ls = function(temp, temp_t, food, food_t, dens_F, dens_F_t, pert_temp_t, pert_food_t, pert_dens_F_t,i){
  
  Si <- S.predict('I',temp,food,dens_F,i) # immature
  Sp <- S.predict('P',temp,food,dens_F,i) # philopatric
  Sb <- S.predict('B',temp,food,dens_F,i) # breeder
  
  PsiIB <- Psi.predict('I',temp,food,dens_F,i)
  PsiPB <- Psi.predict('P',temp,food,dens_F,i)
  
  Bp <- bp(temp_t = temp_t, food_t = food_t, dens_F_t = dens_F_t, i) # breeding prob
  
  Lp <- lp(food_t = food_t, dens_F_t = dens_F_t, i) # litter prob
  LS <- ls(temp_t = pert_temp_t, food_t = pert_food_t, dens_F_t = pert_dens_F_t, i) # litter size
  
  
  A <- matrix(data = c(Si*PsiIB*Lp*LS, Si*(1-PsiIB), Si*PsiIB,
                       Sp*PsiPB*Lp*LS, Sp*(1-PsiPB), Sp*PsiPB,
                       Sb*Bp*Lp*LS,            0,       Sb),nrow=3)
  return(A)
}






