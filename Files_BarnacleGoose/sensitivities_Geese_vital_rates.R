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

# coefficients are intercept at time 1 and age 1 (fledglings); followed by the other times and then covariates

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

### Min, max of covariates:
#########
max.temp.sv=max(env$t_mjunmjul)
min.temp.sv=min(env$t_mjunmjul)

temp.scot.when.temp.sv.max=env$scot_tmin[env$t_mjunmjul==max.temp.sv]
temp.scot.when.temp.sv.min=env$scot_tmin[env$t_mjunmjul==min.temp.sv]

scot.pop.when.temp.sv.max=env$scot_pop[env$t_mjunmjul==max.temp.sv]
scot.pop.when.temp.sv.min=env$scot_pop[env$t_mjunmjul==min.temp.sv]

rain.when.temp.sv.max=env$helg_p_aprmay[env$t_mjunmjul==max.temp.sv]
rain.when.temp.sv.min=env$helg_p_aprmay[env$t_mjunmjul==min.temp.sv]

so.when.temp.sv.max=env$SO[env$t_mjunmjul==max.temp.sv]
so.when.temp.sv.min=env$SO[env$t_mjunmjul==min.temp.sv]

pop.when.temp.sv.max=env$pop_ad[env$t_mjunmjul==max.temp.sv]
pop.when.temp.sv.min=env$pop_ad[env$t_mjunmjul==min.temp.sv]

fox.when.temp.sv.max=env$fox[env$t_mjunmjul==max.temp.sv]
fox.when.temp.sv.min=env$fox[env$t_mjunmjul==min.temp.sv]

rain.ja.when.temp.sv.max=env$p_mjulmaug[env$t_mjunmjul==max.temp.sv]
rain.ja.when.temp.sv.min=env$p_mjulmaug[env$t_mjunmjul==min.temp.sv]

########

max.temp.scot=max(env$scot_tmin)
min.temp.scot=min(env$scot_tmin)

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

############

max.rain=max(env$helg_p_aprmay)
min.rain=min(env$helg_p_aprmay)

temp.sv.when.rain.max=env$t_mjunmjul[env$helg_p_aprmay==max.rain]
temp.sv.when.rain.min=env$t_mjunmjul[env$helg_p_aprmay==min.rain]

scot.pop.when.rain.max=env$scot_pop[env$helg_p_aprmay==max.rain]
scot.pop.when.rain.min=env$scot_pop[env$helg_p_aprmay==min.rain]

temp.scot.when.rain.max=env$scot_tmin[env$helg_p_aprmay==max.rain]
temp.scot.when.rain.min=env$scot_tmin[env$helg_p_aprmay==min.rain]

so.when.rain.max=env$SO[env$helg_p_aprmay==max.rain]
so.when.rain.min=env$SO[env$helg_p_aprmay==min.rain]

pop.when.rain.max=env$pop_ad[env$helg_p_aprmay==max.rain]
pop.when.rain.min=env$pop_ad[env$helg_p_aprmay==min.rain]

fox.when.rain.max=env$fox[env$helg_p_aprmay==max.rain]
fox.when.rain.min=env$fox[env$helg_p_aprmay==min.rain]

rain.ja.when.rain.max=env$p_mjulmaug[env$helg_p_aprmay==max.rain]
rain.ja.when.rain.min=env$p_mjulmaug[env$helg_p_aprmay==min.rain]

############

max.rain.ja=max(env$p_mjulmaug)
min.rain.ja=min(env$p_mjulmaug)

temp.sv.when.rain.ja.max=env$t_mjunmjul[env$p_mjulmaug==max.rain.ja]
temp.sv.when.rain.ja.min=env$t_mjunmjul[env$p_mjulmaug==min.rain.ja]

scot.pop.when.rain.ja.max=env$scot_pop[env$p_mjulmaug==max.rain.ja]
scot.pop.when.rain.ja.min=env$scot_pop[env$p_mjulmaug==min.rain.ja]

temp.scot.when.rain.ja.max=env$scot_tmin[env$p_mjulmaug==max.rain.ja]
temp.scot.when.rain.ja.min=env$scot_tmin[env$p_mjulmaug==min.rain.ja]

so.when.rain.ja.max=env$SO[env$p_mjulmaug==max.rain.ja]
so.when.rain.ja.min=env$SO[env$p_mjulmaug==min.rain.ja]

pop.when.rain.ja.max=env$pop_ad[env$p_mjulmaug==max.rain.ja]
pop.when.rain.ja.min=env$pop_ad[env$p_mjulmaug==min.rain.ja]

fox.when.rain.ja.max=env$fox[env$p_mjulmaug==max.rain.ja]
fox.when.rain.ja.min=env$fox[env$p_mjulmaug==min.rain.ja]

rain.when.rain.ja.max=env$helg_p_aprmay[env$p_mjulmaug==max.rain.ja]
rain.when.rain.ja.min=env$helg_p_aprmay[env$p_mjulmaug==min.rain.ja]

####################### For each coefficient sample, calculate lamba under different perturbations
sens.out=NULL

for(i in 1:100){
  
  ############# 1) PERTURBATION
  
  ####  Proportion reproductive Temp
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,so.when.temp.sv.max,max.temp.sv,rain.when.temp.sv.max)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(0) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,0,0)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                             phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,so.when.temp.sv.min,min.temp.sv,rain.when.temp.sv.min)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(0) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,0,0)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaTemp.sv=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp.sv-min.temp.sv)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Mean june-july temperature on Svalbard",
                                     driver="Temperature",
                                     driver.type="C",
                                     stage.age="adult",vital.rates="prop reproductive",
                                     sens=deltaTemp.sv,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
  ####  Proportion reproductive Rain
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,so.when.rain.max,temp.sv.when.rain.max,max.rain)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(0) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,0,0)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,so.when.rain.min,temp.sv.when.rain.min,min.rain)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(0) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,0,0)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaRain=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain-min.rain)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Rainfall Helgeland",
                                     driver="Rain",
                                     driver.type="C",
                                     stage.age="adult",vital.rates="prop reproductive",
                                     sens=deltaRain,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
  #### PROPORTION OF GOSLINGS FLEDGING
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,0,0,0)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(0) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,max.rain.ja,fox.when.rain.ja.max)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,0,0,0)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(0) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,min.rain.ja,fox.when.rain.ja.min)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaRain.ja=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.rain.ja-min.rain.ja)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Rainfall July Aug",
                                     driver="Rain",
                                     driver.type="C",
                                     stage.age="adult",vital.rates="prop fledging",
                                     sens=deltaRain.ja,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
  
  ####  FELDGLING SURVIVAL
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,0,0,0)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(1) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,0,0)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,max.temp.scot,pop.when.temp.scot.max)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,0,0,0)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(1) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,0,0)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,min.temp.scot,pop.when.temp.scot.min)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,0,0,0)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaTemp.scot=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp.scot-min.temp.scot)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Mean minumum temperature in Scotland",
                                     driver="Temperature",
                                     driver.type="C",
                                     stage.age="fledgling",vital.rates="survival",
                                     sens=deltaTemp.scot,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
  ####  ADULT SURVIVAL
  
  ### Maximum
  R.sim.max=coefs.R[i,]*c(1,0,0,0)
  R.max=inv.logit(sum(R.sim.max))
  
  gos.sim.max=coefs.gos[i,]*c(1,0,0)
  gos.max=exp(sum(gos.sim.max))
  
  # I just assume exp(1) gosling for now
  F.sim.max=coefs.F[i,]*c(1,0,0,0)
  F.sim.max=inv.logit(sum(F.sim.max))
  Fp.max=F.sim.max
  
  surv.f.max=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.max=inv.logit(sum(surv.f.max))
  
  surv.ad.max=coefs.surv.ad[i,]*c(1,1,max.temp.scot,pop.when.temp.scot.max,pop.when.temp.scot.max)
  
  phi.ad.max=inv.logit(sum(surv.ad.max))
  
  mpm.max=matrix(c(0,phi.ad.max*R.max*gos.max*Fp.max*0.5, 
                   phi.f.max,phi.ad.max),nrow=2,ncol=2,byrow=T)
  
  ### Minimum
  R.sim.min=coefs.R[i,]*c(1,0,0,0)
  R.min=inv.logit(sum(R.sim.min))
  
  gos.sim.min=coefs.gos[i,]*c(1,0,0)
  gos.min=exp(sum(gos.sim.min))
  
  # I just assume exp(1) gosling for now
  F.sim.min=coefs.F[i,]*c(1,0,0,0)
  F.sim.min=inv.logit(sum(F.sim.min))
  Fp.min=F.sim.min
  
  surv.f.min=coefs.surv.f[i,]*c(1,0,0)
  
  phi.f.min=inv.logit(sum(surv.f.min))
  
  surv.ad.min=coefs.surv.ad[i,]*c(1,1,min.temp.scot,pop.when.temp.scot.min,pop.when.temp.scot.min)
  
  phi.ad.min=inv.logit(sum(surv.ad.min))
  
  mpm.min=matrix(c(0,phi.ad.min*R.min*gos.min*Fp.min*0.5, 
                   phi.f.min,phi.ad.min),nrow=2,ncol=2,byrow=T)
  
  
  deltaTemp.scot=abs(lambda(mpm.max)-lambda(mpm.min))/abs((max.temp.scot-min.temp.scot)/1)
  
  sens.out=rbind(sens.out,data.frame(spec.driver="Mean minumum temperature in Scotland",
                                     driver="Temperature",
                                     driver.type="C",
                                     stage.age="adult",vital.rates="survival",
                                     sens=deltaTemp.scot,
                                     l_ratio=abs(log(lambda(mpm.max)/lambda(mpm.min))),
                                     cov=0,
                                     sim=i))
  
}

sens_goose_vr=data.frame(study.doi="10.1111/gcb.14773",year.of.publication=2019,
                      group="Birds",species="Branta leucopsis",continent="Europe",driver=sens.out$driver,driver.type="C",
                      stage.age=sens.out$stage.age,vital.rates=sens.out$vital.rates,sens=sens.out$sens,cov=sens.out$cov,mat=2,n.vr=5,n.pam=17,dens=1,
                      biotic_interactions=1,lambda.sim=0,study.length=28)

write.csv(sens_goose_vr,"sens_goose_vital_rates.csv",row.names = F)

