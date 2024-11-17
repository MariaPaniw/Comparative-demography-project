####### Simulating
#https://www.jamesuanhoro.com/post/2018/05/07/simulating-data-from-regression-models/
# assume the coefficients arise from a multivariate normal distribution with the estimated coefficients acting as 
# means and the variance-covariance matrix of the regression coefficients as 
# the variance-covariance matrix for the multivariate normal distribution

# VITAL RATES
#Survival# 

# Survival ~ time + age_class2 + scot_tmin + scot_pop
# 
# age_class2 = either fledgings (0 year old) or adults (pooled age class)
# 
# #Reproduction#
# H = estimated proportion of females with at least one gosling (data column "repG" - female reproduced at least one gosling (1/0))
# G = expected number of goslings per successful female (data column "Ngoslings" -  number of goslings per female with a brood)
# F = proportion of goslings fledging (data columns Nfledglings & failfled - number of fledglings per female with a brood)

library("popbio")
library("arm")
library("boot")
library("ggplot2")
library("pracma")

# Load covariates
setwd("/Users/maria/Dropbox/teaching/esin/Geese")
env=read.csv("env_covar_scaled_geese.csv",header=T,sep=';')
head(env)

year=env$year

load("RepMods_outAll_v2.rdata")

best.mod.H# Binomial
best.mod.G# Poisson
best.mod.F  # Binomial


########################
# Simulate 100 values of parameters
########################
# R 
coef.R=c(fixef(best.mod.H)[1],fixef(best.mod.H)[2],fixef(best.mod.H)[3],fixef(best.mod.H)[4])
coefs.R<- mvrnorm(n = 100, mu = coef.R, Sigma = vcov(best.mod.H))

# gos
coef.gos=c(fixef(best.mod.G)[1],fixef(best.mod.G)[2],fixef(best.mod.G)[3])
coefs.gos<- mvrnorm(n = 100, mu = coef.gos, Sigma = vcov(best.mod.G))

# F
coef.F=c(fixef(best.mod.F)[1],fixef(best.mod.F)[2],fixef(best.mod.F)[3],fixef(best.mod.F)[4])
coefs.F<- mvrnorm(n = 100, mu = coef.F, Sigma = vcov(best.mod.F))

#Survival 

load("Survival_ModOut_v2.rdata")

best.mod=results.best[[1]]$results

# I took out the time coefficients

betas=best.mod$beta[c(1,28:31),]
beta.est=betas[,1]
names(beta.est)=row.names(betas)

coef.surv=beta.est

coef.surv.f=c(coef.surv[1],coef.surv[3],coef.surv[4])
coef.surv.ad=c(coef.surv[1],coef.surv[2],coef.surv[3],coef.surv[4],coef.surv[5])

vcv.f=best.mod$beta.vcv[c(1,29,30),c(1,29,30)]
vcv.ad=best.mod$beta.vcv[c(1,28:31),c(1,28:31)]

coefs.surv.f=mvrnorm(n = 100, mu = coef.surv.f, Sigma = vcv.f)
coefs.surv.ad=mvrnorm(n = 100, mu = coef.surv.ad, Sigma = vcv.ad)

### Min, max of covariates

# Note that I am only perturbing temperatures in Scotland for now because this is the only covariate that appears in > 1 vital rate 

# But this procedure can be done for all covariates

########

max.temp.scot=max(env$scot_tmin)
min.temp.scot=min(env$scot_tmin)

# Here, I account for covariation (so what are the values of the other covariates when temp is maximum or minimum)

temp.sv.when.temp.scot.max=env$t_mjunmjul[env$scot_tmin==max.temp.scot]
temp.sv.when.temp.scot.min=env$t_mjunmjul[env$scot_tmin==min.temp.scot]

scot.pop.when.temp.scot.max=env$scot_pop[env$scot_tmin==max.temp.scot]
scot.pop.when.temp.scot.min=env$scot_pop[env$scot_tmin==min.temp.scot]

rain.when.temp.scot.max=env$helg_p_aprmay[env$scot_tmin==max.temp.scot]
rain.when.temp.scot.min=env$helg_p_aprmay[env$scot_tmin==min.temp.scot]

so.when.temp.scot.max=env$SO[env$scot_tmin==max.temp.scot]
so.when.temp.scot.min=env$SO[env$scot_tmin==min.temp.scot]

pop.when.temp.scot.max=env$pop_ad[env$scot_tmin==max.temp.scot]
pop.when.temp.scot.min=env$pop_ad[env$scot_tmin==min.temp.scot]

fox.when.temp.scot.max=env$fox[env$scot_tmin==max.temp.scot]
fox.when.temp.scot.min=env$fox[env$scot_tmin==min.temp.scot]

rain.ja.when.temp.scot.max=env$p_mjulmaug[env$scot_tmin==max.temp.scot]
rain.ja.when.temp.scot.min=env$p_mjulmaug[env$scot_tmin==min.temp.scot]


####################### For each coefficient sample, calculate lamba under different perturbations
sens.out=NULL

for(i in 1:100){
  
  ############# 1) PERTURBATION NO COVARIATION, I.E., THE OTHER COVARIATES ARE AT THEIR MEAN = 0
  
  ####  TEMPERATURE Scotland (scot_tmin) 
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,0,0,0)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(0) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,0,0)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,max.temp.scot,0)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,max.temp.scot,0,0)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,0,0,0)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(0) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,0,0)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,min.temp.scot,0)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,min.temp.scot,0,0)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaTemp.scot=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp.scot-min.temp.scot)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Mean minumum temperature in Scotland",
                                     driver="Temperature",
                                     driver.type="C",
                                     sens=deltaTemp.scot,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
####################################################################
  
############# 2) PERTURBATION NO COVARIATION, I.E., THE OTHER COVARIATES ARE AT THEIR MEAN = 0
  
  
  ####  TEMPERATURE Scotland (scot_tmin) 
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,so.when.temp.scot.max,temp.sv.when.temp.scot.max,rain.when.temp.scot.max)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,pop.when.temp.scot.max,fox.when.temp.scot.max)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(0) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,rain.ja.when.temp.scot.max,fox.when.temp.scot.max)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,max.temp.scot,scot.pop.when.temp.scot.max)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,max.temp.scot,scot.pop.when.temp.scot.max,scot.pop.when.temp.scot.max)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,so.when.temp.scot.min,temp.sv.when.temp.scot.min,rain.when.temp.scot.min)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,pop.when.temp.scot.min,fox.when.temp.scot.min)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(0) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,rain.ja.when.temp.scot.min,fox.when.temp.scot.min)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,min.temp.scot,scot.pop.when.temp.scot.min)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,min.temp.scot,scot.pop.when.temp.scot.min,scot.pop.when.temp.scot.min)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  deltaTemp.scot.cov=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp.scot-min.temp.scot)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Mean minumum temperature in Scotland",
                                     driver="Temperature",
                                     driver.type="C",
                                     sens=deltaTemp.scot.cov,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=1,
                                     sim=i))
}


sens_goose=data.frame(species="Branta leucopsis", study.doi="10.1111/gcb.14773",year.of.publication=2019,
           group="Birds",continent="Europe",driver=sens.out$driver,driver.type="C",
           stage.age="all",vital.rates="all",sens=sens.out$sens,cov=sens.out$cov,mat=2,n.vr=5,n.pam=17,dens=1,
           biotic_interactions=1,lambda.sim=0,study.length=28)

write.csv(sens_goose,"sens_goose.csv",row.names = F)

# Change density plot line colors by groups
ggplot(sens.out, aes(x=sens, color=factor(cov))) +
  geom_density()

ggplot(sens.out, aes(x=l_ratio, color=factor(cov))) +
  geom_density()

