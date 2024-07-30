
### PREPARATIONS ###############

library(mgcv)
library(lubridate)
library(boot)
library(Matrix)

rm(list=ls())

setwd("~/OneDrive - Universität Zürich UZH/Master Thesis/pert_analyses/meerkats")

##################### Load data ##############

load("boot.GAM")
# It’s an R object, saved as a list (150 items). Each item has different resampled datasets I used to fit the meerkats models

# where you pick one of 50-100 datasets (depending on how fast it goes), refit GAM models (as the code about does), and then, instead of loading the GAMs like you normally do, use these new Gams
# this works for all vital rates, except help.transBH and dom.transF2. I guess here I must have used averages, not GAMs. So, in these 2 cases, keep the mean GAMs that you already have (you can load just these two functions from you GAM folder instead of all of them, if you know what I mean)
load("meerkats-master/GAMs/help.transBH.rda")
load("meerkats-master/GAMs/dom.transF2.rda")

# climate data ##########
load("meerkats-master/SuppMat/clim.obs.rda")
load("meerkats-master/SuppMat/dens.obs.rda")

mu.cov=data.frame(month=clim.obs$month,
                        year=clim.obs$year,
                        tempSD=clim.obs$tempSD,
                        rainSD=clim.obs$rainSD,
                        density=dens.obs$density)

mu.cov.year=aggregate(cbind(rainSD,tempSD,density)~year,data=mu.cov,mean)
mu.cov.month=aggregate(cbind(rainSD,tempSD,density)~month,data=mu.cov,mean)

mu.cov.sd=aggregate(cbind(rainSD,tempSD,density)~month,data=mu.cov,sd)

min.rain=2015
max.rain= as.numeric(as.character(mu.cov.year$year[mu.cov.year$rainSD==max(mu.cov.year$rainSD)]))

mean.rain=0
sd.rain=1

min.temp=as.numeric(as.character(mu.cov.year$year[mu.cov.year$tempSD==min(mu.cov.year$tempSD)]))
max.temp= as.numeric(as.character(mu.cov.year$year[mu.cov.year$tempSD==max(mu.cov.year$tempSD)]))

mean.temp=0
sd.temp=1

min.density=as.numeric(as.character(mu.cov.year$year[mu.cov.year$density==min(mu.cov.year$density)]))
max.density= as.numeric(as.character(mu.cov.year$year[mu.cov.year$density==max(mu.cov.year$density)]))

mean.density=mu.cov.month$density
sd.density=sd(mu.cov$density)


##### 
# helpV2=dataOrig[dataOrig$stage=="H",]
# helpV2=helpV2[!is.na(helpV2$mass),]
# helpV2$pregCat=factor(helpV2$pregCat)
# 
# domV2=droplevels(dataOrig[dataOrig$stage=="D",])
# 
# 
# domV3=dataOrig[dataOrig$stage=="D",]
# domV3=domV3[!is.na(domV3$mass),]

###################################################

# data2=read.csv("~/data2/LH_table_FULL.csv")
# # data2=read.csv("/Users/mariapaniw/Dropbox/Meerkats/data/SuperComputerInput/LH_table_FULL.csv")
# 
# data2$date=as.Date(data2$date,format="%m/%d/%Y")
# data2$date=as.character(data2$date)
# data2=droplevels(data2[data2$date=="2016-12-15",])
# data2$date=as.Date(data2$date)
# 
# data2=data2[data2$ageM>0,]
# 
# load("~/data2/N0_2016.rda")
# # load("/Users/mariapaniw/Dropbox/Meerkats/SuppMat/N0_2016.rda")
# 
# dens0=sum(N0)/data2$home.range2M[1]

ageM.pup=c("1","2","3")
ageM.juv=c("4","5","6")
ageM.sub=c("7","8","9","10","11","12")

# define empty vectors for sensitivities
Pups.Gr.R=NULL # 1. Sensitivity of growth of pups to rain
Pups.Surv.R=NULL # 2. Sensitivity of survival of pups to rain

Juv.Gr.R=NULL # 3. Sensitivity of growth of juveniles to rain
Juv.Surv.R=NULL # 4. Sensitivity of survival of juveniles to rain

Sub.Gr.R=NULL # 5. Sensitivity of growth of subadults to rain
Sub.Surv.R=NULL # 6. Sensitivity of survival of subadults to rain

H.G.R=NULL # 7. Sensitivity of growth of helpers to rain
H.Surv.R=NULL # 8. Sensitivity of survival of helpers to rain

D.G.R=NULL # 9. Sensitivity of growth of dominants to rain
D.Surv.R=NULL # 10. Sensitivity of survival of dominants to rain

H.Rec.R=NULL # 11. Sensitivity of helper recruitment rate to rain
D.Rec.R=NULL # 12. Sensitivity of dominant recruitment rate to rain

# START LOOP ################


iterations <- 50

for(pi in 1:iterations){
  ### refit GAMs ########################
  # BOOSTRAP NUMBER 
  
  ### CREATE NEW GAMS FROM DATA
  pups=boot.GAM[[pi]]$pups
  pups2=pups[!is.na(pups$massNext),]
  
  # Survival
  pup.surv=gam(surv ~ ageM + te(density, tempSD, k = 3, bs = "cr") + te(mass,by = ageM, k = 3, bs = "cr") +
                 te(mass, density, k = 3, bs = "cr") + 
                 te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass, tempSD, k = 3, bs = "cr") +
                 te(rainSD, month, k = c(3, 4),bs = c("cr", "cc")) + te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")) ,
               data=pups,family="binomial",gamma=1.4)
  
  #Growth
  pup.gr=gam(massNext ~ ageM + s(year, bs = "re") + te(density, month, k = c(3,4), bs = c("cr", "cc")) + te(density, rainSD, k = 3, bs = "cr") + 
               te(density, tempSD, k = 3, bs = "cr") + 
               te(mass, by = ageM, k = 3, bs = "cr") + te(mass, density, k = 3, bs = "cr") + 
               te(mass, rainSD, k = 3, bs = "cr") + te(mass, tempSD, k = 3, bs = "cr") +
               te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD, tempSD, k = 3, bs = "cr") + 
               te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")) ,data=pups2,gamma=1.4)
  
  #Growth variation
  pups2$resid=log(as.numeric(residuals(pup.gr))^2)
  
  pup.gr.var=gam(resid ~ ageM + s(year, bs = "re") + te(mass, by = ageM, k = 3, bs = "cr") +
                   te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + 
                   te(mass, rainSD, k = 3, bs = "cr") + te(rainSD, month, k = c(3,4), bs = c("cr", "cc")) +
                   te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=pups2,gamma=1.4)
  
  #################### JUVENILES
  juv=boot.GAM[[pi]]$juv
  
  juv2=juv[!is.na(juv$massNext),]
  
  # Survival
  juv.surv=gam(surv ~  te(density, rainSD, k = 3, bs = "cr") +
                 te(mass, tempSD, k = 3, bs = "cr") +
                 te(rainSD, tempSD, k = 3, bs = "cr"),
               data=juv,family="binomial",gamma=1.4)
  # Growth
  juv.gr=gam(massNext ~ ageM + s(year, bs = "re") +
               te(density, month, k = c(3,4), bs = c("cr", "cc")) + te(density, rainSD, k = 3, bs = "cr") + 
               te(density, tempSD, k = 3, bs = "cr") + te(mass, by = ageM,k = 3, bs = "cr") + 
               te(mass, density, k = 3, bs = "cr") + 
               te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass,rainSD, k = 3, bs = "cr") +
               te(mass, tempSD, k = 3, bs = "cr") + 
               te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + 
               te(rainSD, tempSD, k = 3, bs = "cr") +
               te(tempSD, month, k = c(3, 4),bs = c("cr", "cc")),data=juv2,gamma=1.4)
  
  # Grwoth variation
  juv2$resid=log(as.numeric(residuals(juv.gr))^2)
  
  juv.gr.var=gam(resid ~ s(year, bs = "re") + te(density, month, k = c(3, 4),bs = c("cr", "cc")) + te(density, rainSD, k = 3, bs = "cr") + 
                   te(density, tempSD, k = 3, bs = "cr") + te(mass, density, k = 3, bs = "cr") + 
                   te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + 
                   te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=juv2,gamma=1.4)
  
  #################### SUBADULT
  subadult=boot.GAM[[pi]]$subadult
  
  subadult2=subadult[!is.na(subadult$massNext),]
  
  # Survival
  sub.surv=gam(surv ~ ageM +  te(density, month, k = c(3,4), bs = c("cr", "cc")) +
                 te(density, rainSD, k = 3, bs = "cr") + 
                 te(density, tempSD, k = 3, bs = "cr") +   te(mass, tempSD,k = 3, bs = "cr") +
                 te(rainSD, tempSD, k = 3, bs = "cr"),data=subadult,family="binomial",gamma=1.4)
  
  # Growth 
  sub.gr=gam(massNext ~ ageM + s(year, bs = "re") + te(density, month, k = c(3,4), bs = c("cr", "cc")) +
               te(density, rainSD, k = 3, bs = "cr") + 
               te(density, tempSD, k = 3, bs = "cr") + te(mass, by = ageM,k = 3, bs = "cr") + 
               te(mass, density, k = 3, bs = "cr") + 
               te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass,rainSD, k = 3, bs = "cr") +
               te(mass, tempSD, k = 3, bs = "cr") + 
               te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD,tempSD, k = 3, bs = "cr") +
               te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")) ,data=subadult2,gamma=1.4)
  
  # Growth variation
  subadult2$resid=log(as.numeric(residuals(sub.gr))^2)
  
  sub.gr.var=gam(resid ~ ageM+s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) + 
                   te(density, rainSD, k = 3, bs = "cr") + 
                   te(density, tempSD, k = 3, bs = "cr") + te(mass, by = ageM, k = 3, bs = "cr") +
                   te(mass, month, k = c(3, 4), bs = c("cr",  "cc")) + te(mass, rainSD, k = 3, bs = "cr") + te(mass, tempSD, k = 3, bs = "cr") +
                   te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(rainSD, tempSD, k = 3, bs = "cr") + te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),data=subadult2,gamma=1.4)
  
  ################### HELPER 
  
  helpS=boot.GAM[[pi]]$helpS
  help2=boot.GAM[[pi]]$help2
  
  helpE=boot.GAM[[pi]]$helpE
  
  helpHD=boot.GAM[[pi]]$helpHD
  
  helpNPH=boot.GAM[[pi]]$helpNPH
  
  helpNPD=boot.GAM[[pi]]$helpNPD
  
  helpFH=boot.GAM[[pi]]$helpFH
  
  
  helpFH2=boot.GAM[[pi]]$helpFH2
  
  helpSH=boot.GAM[[pi]]$helpSH
  
  helpSH2=boot.GAM[[pi]]$helpSH2
  
  
  helpR=boot.GAM[[pi]]$helpR
  
  # Survival
  help.surv=gam(surv~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) + 
                  te(density, rainSD, k = 3, bs = "cr") + te(density, tempSD, k = 3, bs = "cr") + te(mass, density, k = 3, bs = "cr") + 
                  te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + 
                  te(mass, rainSD, k = 3, bs = "cr") + te(mass, tempSD, k = 3, bs = "cr") + 
                  te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),
                data=helpS,family="binomial",gamma=1.4)
  
  # Growth
  help.gr=gam(massNext ~ pregCat + s(year, bs = "re") + te(density, month,k = c(3, 4), bs = c("cr", "cc")) + 
                te(density, rainSD, k = 3, bs = "cr") + te(density, tempSD, k = 3, bs = "cr") + 
                te(mass, by = pregCat, k = 3, bs = "cr") + te(mass, density, k = 3,  bs = "cr") + 
                te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + 
                te(mass, tempSD, k = 3, bs = "cr") + te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                te(rainSD, tempSD, k = 3, bs = "cr") + 
                te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")) ,data=help2,gamma=1.4)
  
  # Growth variation
  help2$resid=log(as.numeric(residuals(help.gr))^2)
  
  help.gr.var=gam(resid ~ pregCat + s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                    te(density, rainSD, k = 3, bs = "cr") + 
                    te(density, tempSD, k = 3, bs = "cr") + te(mass, by = pregCat,k = 3, bs = "cr") +
                    te(mass, density, k = 3, bs = "cr") + 
                    te(mass, rainSD, k = 3, bs = "cr") + te(rainSD, month, k = c(3,4), bs = c("cr", "cc")),data=help2,gamma=1.4)
  
  
  # Emigration
  help.emig=gam(emig ~s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                  te(density, rainSD, k = 3, bs = "cr") +
                  te(mass,month, k = c(3, 4), bs = c("cr", "cc")) + 
                  te(mass, tempSD,k = 3, bs = "cr") + te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                  te(rainSD, tempSD, k = 3, bs = "cr") + te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),
                data=helpE,family="binomial",gamma=1.4)
  
  
  # Transition to dominant
  toDom=gam(trans ~ te(density, rainSD, k = 3, bs = "cr") + 
              te(mass, density, k = 3, bs = "cr") +
              te(rainSD, month, k = c(3, 4), bs = c("cr","cc")),
            data=helpHD,family="binomial",gamma=1.4)
  
  # Transitions NP HELP (stay HELP) 
  help.transNPH=gam(trans ~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                      te(density, rainSD, k = 3, bs = "cr") + 
                      te(density, tempSD, k = 3, bs = "cr") + te(mass, month, k = c(3,4), bs = c("cr", "cc")) +
                      te(mass, tempSD, k = 3, bs = "cr") + 
                      te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD,tempSD, k = 3, bs = "cr") +
                      te(tempSD, month, k = c(3, 4),bs = c("cr", "cc")),
                    data=helpNPH,family="binomial",gamma=1.4)
  
  
  # Transitions NP HELP (go to DOM) 
  help.transNPD=gam(trans ~ te(mass, k = 3, bs = "cr") + te(month, k = 4, bs = "cc"),data=helpNPD,family="binomial",gamma=1.4)
  
  # Transitions FIRST HELP (stay HELP): go to SECOND or ABORT
  help.transFH=gam(trans ~ s(year, bs = "re") + te(mass, density, k = 3, bs = "cr") + 
                     te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD,tempSD, k = 3, bs = "cr") +
                     te(tempSD, month, k = c(3, 4),bs = c("cr", "cc")),data=helpFH,family="binomial",gamma=1.4)
  
  
  # Transitions FIRST HELP (stay HELP): if ABORT, go to FIRST or NP
  help.transFH2=gam(trans ~ te(density, month, k = c(3, 4), bs = c("cr", "cc")) + 
                      te(density, rainSD, k = 3, bs = "cr") + te(density, tempSD,k = 3, bs = "cr") +
                      te(mass, rainSD, k = 3, bs = "cr") + 
                      te(mass, tempSD, k = 3, bs = "cr") + te(rainSD, month, k = c(3,4), bs = c("cr", "cc")),
                    data=helpFH2,family="binomial",gamma=1.4)
  
  
  # Transitions SECOND HELP (stay HELP): go to BIRTH or ABORT
  help.transSH=gam(trans ~s(year, bs = "re") + te(density, month, k = c(3, 4),bs = c("cr", "cc")) + 
                     te(mass, density, k = 3, bs = "cr") + te(mass, month, k = c(3,  4), bs = c("cr", "cc")),
                   data=helpSH,family="binomial",gamma=1.4)
  
  
  # Transitions SECOND HELP (stay HELP): if ABORT, go to FIRST or NP
  help.transSH2=gam(trans ~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                      te(mass, density, k = 3, bs = "cr") + 
                      te(mass, rainSD, k = 3, bs = "cr") + te(rainSD, month, k = c(3,4), bs = c("cr", "cc")) +
                      te(tempSD, month, k = c(3, 4),  bs = c("cr", "cc")),
                    data=helpSH2,family="binomial",gamma=1.4)
  
  
  # Transitions BIRTH HELP (stay HELP): FIRST or NP
  # GET FROM NORMAL GAM OBJECT
  
  # Recruitment
  help.pups=gam(pups ~ te(density, month, k = c(3, 4), bs = c("cr","cc")) +
                  te(density, rainSD, k = 3, bs = "cr") + te(mass, month, k = c(3, 4), bs = c("cr", "cc")) +
                  te(rainSD, month,k = c(3, 4), bs = c("cr", "cc")) ,
                data=helpR,family="poisson",gamma=1.4)
  
  
  newHT=boot.GAM[[pi]]$newHT
  
  
  help.off.mass=gam(mass~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr","cc")) +
                      te(massM, month, k = c(3, 4), bs = c("cr", "cc")) + 
                      te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                      te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),data=newHT,gamma=1.4)
  # Offspring mass variation 
  newHT$resid=log(as.numeric(residuals(help.off.mass))^2)
  
  
  help.off.mass.var=gam(resid ~  te(massM, density,k = 3, bs = "cr"),data=newHT,gamma=1.4)
  
  ################### DOMINANT
  
  domS=boot.GAM[[pi]]$domS
  
  
  dom2=boot.GAM[[pi]]$dom2
  
  domNP=boot.GAM[[pi]]$domNP
  
  domF=boot.GAM[[pi]]$domF
  
  
  domSec=boot.GAM[[pi]]$domSec
  
  
  domSec2=boot.GAM[[pi]]$domSec2
  
  domB=boot.GAM[[pi]]$domB
  
  domB2=boot.GAM[[pi]]$domB2
  
  domR=boot.GAM[[pi]]$domR
  
  
  # Survival
  dom.surv=gam(surv ~ te(density, month, k = c(3, 4), bs = c("cr", "cc")) + 
                 te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass, rainSD, k = 3, bs = "cr") + 
                 te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + 
                 te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")) ,
               data=domS,family="binomial",gamma=1.4)
  
  # Growth
  dom.gr=gam(massNext ~ pregCat + s(year, bs = "re") + te(density, month,k = c(3, 4), bs = c("cr", "cc")) + 
               te(mass, by = pregCat,k = 3, bs = "cr") + te(mass, density, k = 3, bs = "cr") + 
               te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass,rainSD, k = 3, bs = "cr") +
               te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
               te(rainSD, tempSD, k = 3, bs = "cr") + 
               te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=dom2,gamma=1.4)
  
  
  # Growth variation 
  dom2$resid=log(as.numeric(residuals(dom.gr))^2)
  
  dom.gr.var=gam(resid ~ pregCat + s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(density, rainSD, k = 3, bs = "cr") + 
                   te(mass, by = pregCat, k = 3, bs = "cr") + te(mass, density,k = 3, bs = "cr") +
                   te(mass, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + 
                   te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=dom2,gamma=1.4)
  
  # Transitions NP DOM: go to FIRST or NP
  dom.transNP=gam(trans ~ s(year, bs = "re") + te(density, month, k = c(3, 4),bs = c("cr", "cc")) + 
                    te(density, rainSD, k = 3, bs = "cr") + 
                    te(mass, density, k = 3, bs = "cr") + te(rainSD, month, k = c(3,4), bs = c("cr", "cc")) +
                    te(rainSD, tempSD, k = 3, bs = "cr") + 
                    te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),
                  data=domNP,family="binomial",gamma=1.4)
  
  
  # Transitions FIRST DOM: go to SECOND or ABORT
  dom.transF=gam(trans ~ te(density, month, k = c(3, 4), bs = c("cr", "cc")) + 
                   te(density, tempSD, k = 3, bs = "cr") + te(mass, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(mass, tempSD, k = 3, bs = "cr") + 
                   te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD,tempSD, k = 3, bs = "cr") +
                   te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=domF,family="binomial",gamma=1.4)
  
  # Transitions FIRST DOM: if ABORT, go to FIRST or NP
  #GET FROM NORMAL GAMS
  
  # Transitions SECOND DOM: go to BIRTH or ABORT
  dom.transSec=gam(trans ~ s(year, bs = "re") + te(density, month, k = c(3, 4),bs = c("cr", "cc")) +
                     te(density, rainSD, k = 3, bs = "cr") + 
                     te(mass, month, k = c(3, 4), bs = c("cr", "cc")) + te(mass,tempSD, k = 3, bs = "cr") +
                     te(rainSD, month, k = c(3, 4),bs = c("cr", "cc")) + te(rainSD, tempSD, k = 3, bs = "cr") + 
                     te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=domSec,family="binomial",gamma=1.4)
  
  
  # Transitions SECOND DOM: if ABORT, go to FIRST or  NP
  dom.transSec2=gam(trans ~   te(mass, k = 3, bs = "cr") + te(tempSD, k = 3, bs = "cr"),data=domSec2,family="binomial",gamma=1.4)
  
  # Transitions BIRTH DOM: go to NP or not
  dom.transB=gam(trans ~ s(year, bs = "re") + te(density, month, k = c(3, 4),bs = c("cr", "cc")) +
                   te(density, rainSD, k = 3, bs = "cr") + 
                   te(mass, rainSD, k = 3, bs = "cr") + te(mass, tempSD, k = 3,bs = "cr") +
                   te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) +
                   te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=domB,family="binomial",gamma=1.4)
  
  # Transitions BIRTH DOM: if not to NP, go to FIRST or SECOND
  dom.transB2=gam(trans ~ s(year, bs = "re") + te(mass, month, k = c(3, 4), bs = c("cr","cc")) +
                    te(mass, rainSD, k = 3, bs = "cr") + te(mass, tempSD,k = 3, bs = "cr") +
                    te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD, tempSD, k = 3, bs = "cr") +
                    te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),data=domB2,family="binomial",gamma=1.4)
  
  
  # Recruitment
  dom.pups=gam(pups ~ te(mass, density, k = 3, bs = "cr"),data=domR,family="poisson",gamma=1.4)
  
  # Offspring mass
  
  
  newDT=boot.GAM[[pi]]$newDT
  
  dom.off.mass=gam(mass~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr","cc")) +
                     te(density, tempSD, k = 3, bs = "cr") + te(massM,density, k = 3, bs = "cr") +
                     te(massM, month, k = 4, bs = "cc") + 
                     te(massM, rainSD, k = 3, bs = "cr") + te(massM, tempSD, k = 3,bs = "cr") + 
                     te(rainSD, month, k = c(3, 4), bs = c("cr", "cc")) + te(rainSD, tempSD, k = 3, bs = "cr") +
                     te(tempSD,month, k = c(3, 4), bs = c("cr", "cc")),data=newDT,gamma=1.4)
  
  # Offspring mass variation 
  newDT$resid=log(as.numeric(residuals(dom.off.mass))^2)
  
  dom.off.mass.var=gam(resid ~ s(year, bs = "re") + te(density, month, k = c(3, 4), bs = c("cr", "cc")) +
                         te(density, tempSD, k = 3, bs = "cr") + 
                         te(massM,rainSD, k = 3, bs = "cr") + 
                         te(rainSD, month, k = c(3, 4),bs = c("cr", "cc")) + te(rainSD, tempSD, k = 3, bs = "cr") + 
                         te(tempSD, month, k = c(3, 4), bs = c("cr", "cc")),data=newDT,gamma=1.4)
  
  ###################################################################
  # Rest of code starting with ##### STEP 1: create functions
  #########################
  
  
  # pick one of 50-100 datasets, refit GAM models as above, and then, instead of loading GAMs like normally, use these new GAMs
  # this works for all vital rates, except help.transBH and dom.transF2, here keep the mean GAMs that you already have (you can load two functions from your GAM folder instead of all of them)
  
  # put together IPM ##########################
  
  # load vital rate GAMs for the two vital rates above
  # file.names <- list.files(path = "/meerkats-master/GAMs")
  
  
  
  
  
  ### EIGENANALYSIS THROUGH ITERATION
  
  get.eigen.stuff <- function(mat){ 
    sz <- dim(mat)[1]
    t.now <- runif(sz)
    t.now <- t.now/sum(t.now)
    t.next <- mat%*%t.now
    t.next <- t.next/sum(t.next)
    i <- 0
    while (sum(abs(t.next-t.now))>0.0000001){
      i <- i+1
      # print(i)
      t.now <- t.next
      t.next <- mat%*%t.now
      lambda <- sum(t.next)/sum(t.now)
      t.next <- t.next/sum(t.next)
    }
    
    return(lambda)
  }
  
  
  # STEP 1: create functions ###########
  
  # x = mass (discretized, 100 bins)
  # ageM = age in months
  # rain/tempSD = rainfall/temperature standardized deviation (variation)
  # density = population density
  # pregCat = pregnancy category 
  # massM = mass of mother
  
  f.pup.surv=function(x,ageM,month,rainSD,tempSD,density){
    
    new.data=expand.grid(mass=x,ageM=ageM.pup,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density)
    new.data$pred=predict(pup.surv,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.pup.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(mass=x,ageM=ageM.pup,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(pup.gr,newdata = new.data)
    var=2*exp(as.numeric(predict(pup.gr.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    
    return(exp(-temp2)/temp1)
  }
  
  
  f.juv.surv=function(x,density,rainSD,tempSD){
    
    new.data=expand.grid(mass=x,density=density,rainSD=rainSD,
                         tempSD=tempSD)
    new.data$pred=predict(juv.surv,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.juv.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
    new.data=expand.grid(mass=x,ageM=ageM.juv,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(juv.gr,newdata = new.data)
    var=2*exp(as.numeric(predict(juv.gr.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    
    return(exp(-temp2)/temp1)
  }
  
  f.sub.surv=function(x,ageM,month,density,rainSD,tempSD){
    
    new.data=expand.grid(mass=x,ageM=ageM.sub,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density)
    new.data$pred=predict(sub.surv,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.sub.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(mass=x,ageM=ageM.sub,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(sub.gr,newdata = new.data)
    var=2*exp(as.numeric(predict(sub.gr.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    
    return(exp(-temp2)/temp1)
  }
  
  f.help.surv=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(help.surv,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.gr=function(x,y,pregCat,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(mass=x,pregCat=pregCat,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(help.gr,newdata = new.data)
    var=2*exp(as.numeric(predict(help.gr.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    return(exp(-temp2)/temp1)
  }
  
  f.help.emig=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(help.emig,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.toDom=function(x,month,density,rainSD){
    
    new.data=expand.grid(mass=x,month=month,density=density,rainSD=rainSD)
    new.data$pred=predict(toDom,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.NPH=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(help.transNPH,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.NPD=function(x,month){
    
    new.data=expand.grid(mass=x,month=month)
    new.data$pred=predict(help.transNPD,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.FH=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(help.transFH,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.FH2=function(x,month,density,rainSD,tempSD){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density)
    new.data$pred=predict(help.transFH2,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.SH=function(x,month,density,year){
    
    new.data=expand.grid(mass=x,month=month,density=density,year=year)
    new.data$pred=predict(help.transSH,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.SH2=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(help.transSH2,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.BH=function(x,month){
    
    new.data=expand.grid(mass=x,month=month)
    new.data$pred=predict(help.transBH,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.pups=function(x,month,density,rainSD){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         density=density)
    new.data$pred=predict(help.pups,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.help.off.mass=function(x,y,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(massM=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(help.off.mass,newdata = new.data)
    
    var=2*exp(as.numeric(predict(help.off.mass.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    
    # vr=residuals(help.off.mass)^2
    # sd=sqrt(mean(vr))
    # temp1 <- sqrt(2*pi)*sd
    # temp2 <- ((y-mu)^2)/(2*sd^2)
    return(exp(-temp2)/temp1)
  }
  
  f.dom.surv=function(x,month,density,rainSD,tempSD){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density)
    new.data$pred=predict(dom.surv,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.gr=function(x,y,pregCat,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(mass=x,pregCat=pregCat,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(dom.gr,newdata = new.data)
    var=2*exp(as.numeric(predict(dom.gr.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    return(exp(-temp2)/temp1)
    
  }
  
  f.dom.NP=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(dom.transNP,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.F=function(x,month,density,rainSD,tempSD){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density)
    new.data$pred=predict(dom.transF,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.F2=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(dom.transF2,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.Sec=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(dom.transSec,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.Sec2=function(x,tempSD){
    
    new.data=expand.grid(mass=x,tempSD=tempSD)
    new.data$pred=predict(dom.transSec2,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.B=function(x,month,density,rainSD,tempSD,year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(dom.transB,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.B2=function(x,month,density,rainSD,tempSD,year=year){
    
    new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    new.data$pred=predict(dom.transB2,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.pups=function(x,density){
    
    new.data=expand.grid(mass=x,density=density)
    
    new.data$pred=predict(dom.pups,newdata = new.data,type="response")
    return(new.data)
  }
  
  f.dom.off.mass=function(x,y,month,rainSD,tempSD,density,year){
    
    new.data=expand.grid(massM=x,month=month,rainSD=rainSD,
                         tempSD=tempSD,density=density,year=year)
    mu=predict(dom.off.mass,newdata = new.data)
    
    var=2*exp(as.numeric(predict(dom.off.mass.var,newdata = new.data)))
    temp1 <- sqrt(2*pi*var)
    temp2 <- ((y-mu)^2)/(2*var)
    # vr=residuals(dom.off.mass)^2
    # sd=sqrt(mean(vr))
    # temp1 <- sqrt(2*pi)*sd
    # temp2 <- ((y-mu)^2)/(2*sd^2)
    return(exp(-temp2)/temp1)
    
  }
  
  # STEP 2: FILL IN KERNELS ##########
  #### SIMULATE VARYING DENSITY 
  minMass=3.89 # minimum observed mass
  maxMass=6.88 # maximum observed mass
  
  n.bins = 100; n.stage = 20; # number of bins and total number of life-cycle stages considered 
  b <- minMass+c(0:n.bins)*(maxMass-minMass)/n.bins 
  z <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # bin midpoint 
  h <- (maxMass - minMass)/n.bins # bin width
  
  ageM.pup.tot=c("1","2","3")
  ageM.juv.tot=c("4","5","6")
  ageM.sub.tot=c("7","8","9","10","11","12")
  
  m=n.bins
  
  # STEP 3: SENSITIVITY ANALYSES ###############
  ### 1. Sens of pup growth to rain ---------------
  lambdaR.1=NULL
  
  rain=as.character(c(min.rain,max.rain))
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.1[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  Pups.Gr.R[pi]=abs((lambdaR.1[2]-lambdaR.1[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  ### 2. Sens of pup survival to rain ###############
  lambdaR.2=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,pert.rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,pert.rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,pert.rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.2[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  Pups.Surv.R[pi]=abs((lambdaR.2[2]-lambdaR.2[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  ### 3. Sens of juvenile growth to rain ################
  lambdaR.3=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.3[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  Juv.Gr.R[pi]=abs((lambdaR.3[2]-lambdaR.3[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  
  ### 4. Sens of juvenile survival to rain ##############
  lambdaR.4=NULL
  
  
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,pert.rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,pert.rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,pert.rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.4[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  Juv.Surv.R[pi]=abs((lambdaR.4[2]-lambdaR.4[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))

  
  
  ### 5. Sens of subadults growth to rain ####################
  lambdaR.5=NULL
  
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,pert.rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.5[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  
  Sub.Gr.R[pi]=abs((lambdaR.5[2]-lambdaR.5[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
   

  
  ### 6. Sens of subadults survival to rain ###################
  lambdaR.6=NULL
  
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,pert.rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.6[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  Sub.Surv.R[pi]=abs((lambdaR.6[2]-lambdaR.6[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  ### 7. Sens of helper growth to rain #########################
  lambdaR.7=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.7[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  H.G.R[pi]=abs((lambdaR.7[2]-lambdaR.7[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  ### 8. Sens of helper survival to rain #########################
  lambdaR.8=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,pert.rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.8[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  H.Surv.R[pi]=abs((lambdaR.8[2]-lambdaR.8[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  ### 9. Sens of dominant growth to rain ##################
  lambdaR.9=NULL
  
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,pert.rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.9[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  D.G.R[pi]=abs((lambdaR.9[2]-lambdaR.9[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  ### 10. Sens of dominant survival to rain ####################
  lambdaR.10=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,pert.rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.10[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  D.Surv.R[pi]=abs((lambdaR.10[2]-lambdaR.10[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  
  ### 11. Sens of helper recruitment to rain #######################
  lambdaR.11=NULL
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,pert.rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,pert.rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.11[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  H.Rec.R[pi]=abs((lambdaR.11[2]-lambdaR.11[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  
  
  
  ### 12. Sens of dominant recruitment to rain #######################
  lambdaR.12=NULL
  
  
  for(j in 1:2){
    
    IPMmonth=array(0,c(12,n.bins*n.stage,n.bins*n.stage))
    
    
    for(i in 1:12){
      
      
      IPM=array(0,c(n.bins*n.stage,n.bins*n.stage))
      ### Get observed values for given year and month
      tempSD=mean.temp
      rainSD=mean.rain
      pert.rainSD=mu.cov$rainSD[mu.cov$year==rain[j]][i] # max or min rain
      
      #### REGULAR IPM
      
      density=mu.cov.month$density[i] # density of the month
      
      ### Age 1 Pups
      ageM.pup=ageM.pup.tot[1]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
      
      ### Age 2 Pups
      ageM.pup=ageM.pup.tot[2]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
      
      ### Age 3 Pups
      ageM.pup=ageM.pup.tot[3]
      G=h*t(outer(z,z,f.pup.gr,ageM.pup,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.pup.surv(z,ageM.pup,i,rainSD,tempSD,density)$pred)
      IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
      
      ### Age 4 Juvenile
      ageM.juv=ageM.juv.tot[1]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
      
      ### Age 5 Juvenile
      ageM.juv=ageM.juv.tot[2]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
      
      ### Age 6 Juvenile
      ageM.juv=ageM.juv.tot[3]
      
      G=h*t(outer(z,z,f.juv.gr,ageM.juv,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
      
      IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
      
      ### Age 7 Subadult
      ageM.sub=ageM.sub.tot[1]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
      
      ### Age 8 Subadult
      ageM.sub=ageM.sub.tot[2]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
      
      ### Age 9 Subadult
      ageM.sub=ageM.sub.tot[3]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
      
      ### Age 10 Subadult
      ageM.sub=ageM.sub.tot[4]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
      
      ### Age 11 Subadult
      ageM.sub=ageM.sub.tot[5]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
      
      ### Age 12 Subadult
      ageM.sub=ageM.sub.tot[6]
      
      G=h*t(outer(z,z,f.sub.gr,ageM.sub,i,rainSD,tempSD,density,rain[j]))
      # control for eviction
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      S=diag(f.sub.surv(z,ageM.sub,i,density,rainSD,tempSD)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
      
      ### HELPER 
      
      S=diag(f.help.surv(z,i,density,rainSD,tempSD,rain[j])$pred)
      Shelp1=S
      E=diag(f.help.emig(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      Ph_d=diag(f.toDom(z,i,density,rainSD)$pred)
      
      ### not pregnant -> staying not pregnant helper
      pregCat="np"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp1=G
      Thh.np.p1=diag(f.help.NPH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
      
      ### not pregnant -> pregnant helper 1 month 
      
      IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
      
      ### not pregnant -> not pregnant dominant
      Thd.p1=diag(f.help.NPD(z,i)$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
      
      ### not pregnant -> pregnant dominant 1 month
      
      IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
      
      ### pregnant helper 1 month -> pregnant helper 2 month 
      
      pregCat="first"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp2=G
      Thh.p1.p2=diag(f.help.FH(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
      
      ### pregnant helper 1 month -> pregnant helper 1 month
      Thh.p1.p1=diag(f.help.FH2(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
      
      ### pregnant helper 1 month -> non-pregnant helper 
      IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
      
      ### pregnant helper 1 month -> pregnant dominant 2 months 
      IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
      
      ### pregnant helper 2 month -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp3=G
      Thh.p2.b=diag(f.help.SH(z,i,density,rain[j])$pred)
      
      IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
      
      ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
      Thh.p2.p1=diag(f.help.SH2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
      
      ### pregnant helper 2 month -> abortion and back to non-pregnant helper
      
      IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
      
      ### pregnant helper 2 month -> birth as dominant
      
      IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
      
      ### pregnant helper 2 month -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
      
      ### helper given birth -> non-pregnant helper
      pregCat="birth"
      G=h*t(outer(z,z,f.help.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      Ghelp4=G
      Thh.b.p1=as.numeric(f.help.BH(z,i)$pred)
      
      IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
      
      ### helper given birth -> pregnant helper
      
      IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
      
      ### helper given birth -> non-pregnant dominant
      
      IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      ### helper given birth -> pregnant dominant
      
      IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
      
      
      ### DOMINANT
      S=diag(f.dom.surv(z,i,density,rainSD,tempSD)$pred)
      
      ### non-pregnant -> non-pregnant
      pregCat="np"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.np.p1=diag(f.dom.NP(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
      
      ### non-pregnant -> pregnant
      
      IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
      
      ### pregnant month 1  -> pregnant month 2
      pregCat="first"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p1.p2=diag(f.dom.F(z,i,density,rainSD,tempSD)$pred)
      
      IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
      
      ### pregnant month 1  -> pregnant month 1
      if(rain[j]=="1997"|rain[j]=="1998"){
        
        Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,"1999")$pred)
      }else{ Td.p1.p1=diag(f.dom.F2(z,i,density,rainSD,tempSD,rain[j])$pred) }
      
      
      IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
      
      ### pregnant month 1  -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
      
      ### pregnant month 2 -> birth
      pregCat="second"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.p2.b=diag(f.dom.Sec(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
      
      ### pregnant month 2 -> pregnant 1 month
      Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
      
      ### pregnant month 2 -> non-pregnant
      
      IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
      
      ### after birth -> non-pregnant
      
      pregCat="birth"
      G=h*t(outer(z,z,f.dom.gr,pregCat,i,rainSD,tempSD,density,rain[j]))
      G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      G[is.na(G)]=0
      
      Td.b.np=diag(f.dom.B(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
      
      ### after birth -> pregnant 1 month
      Td.b.p1=diag(f.dom.B2(z,i,density,rainSD,tempSD,rain[j])$pred)
      
      IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
      
      ### after birth -> pregnant 2 month
      
      IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
      
      ###### RECRUITMENT
      
      # from helper
      R=diag(as.numeric(f.help.pups(z,i,density,rainSD)$pred),n.bins)
      Rhelp=R
      D=h*t(outer(z,z,f.help.off.mass,i,rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      Dhelp=D
      IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
      
      # from dominant
      R=diag(f.dom.pups(z,density)$pred)
      D=h*t(outer(z,z,f.dom.off.mass,i,pert.rainSD,tempSD,density,rain[j]))
      D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
      D[is.na(D)]=0
      
      IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
      
      ##### FILL IN YEARLY IPM
      
      IPMmonth[i,,]=IPM
      
      
    }
    
    IPM.full=Matrix(IPMmonth[12,,])%*%Matrix(IPMmonth[11,,])%*%Matrix(IPMmonth[10,,])%*%Matrix(IPMmonth[9,,])%*%Matrix(IPMmonth[8,,])%*%Matrix(IPMmonth[7,,])%*%Matrix(IPMmonth[6,,])%*%Matrix(IPMmonth[5,,])%*%Matrix(IPMmonth[4,,])%*%Matrix(IPMmonth[3,,])%*%Matrix(IPMmonth[2,,])%*%Matrix(IPMmonth[1,,])
    
    lambdaR.12[j]=get.eigen.stuff(as.matrix(IPM.full))
    
    
  }
  
  
  D.Rec.R[pi]=abs((lambdaR.12[2]-lambdaR.12[1])/((mean(mu.cov$rainSD[mu.cov$year==rain[2]])-mean(mu.cov$rainSD[mu.cov$year==rain[1]]))/1))
  
  }







# STEP 4: SAVE OUTPUT #########################

# compile output for each stage separately and then merge at the end and save

# dominant adults
Results_domadults=data.frame(study.doi="10.1126/science.aau5905",
                   year.of.publication="2019",
                   group="Mammals",
                   species="Suricata suricatta",
                   continent="Africa",
                   driver="rain",
                   driver.type="C",
                   stage.age="dominant adult",
                   vital.rates=rep(c("growth","survival","recruitment"),each=50),
                   sens=c(D.G.R,
                          D.Surv.R,
                          D.Rec.R),
                   mat=1, # age at sexual maturity Paniw et al. 2020
                   n.vr=28, # number of vital rates with covariates
                   n.pam=1422, # number of total parameters of all these vital rates
                   dens=1, # density dependence in it?
                   biotic_interactions=0, # any biotic interactions?
                   lambda.sim=0,
                   study.length=20)

# helper adults
Results_helperadults=data.frame(study.doi="10.1126/science.aau5905",
                                year.of.publication="2019",
                                group="Mammals",
                                species="Suricata suricatta",
                                continent="Africa",
                                driver="rain",
                                driver.type="C",
                                stage.age="helper adult",
                                vital.rates=rep(c("growth","survival","recruitment"),each=50),
                                sens=c(H.G.R,
                                       H.Surv.R,
                                       H.Rec.R),
                                mat=1, # age at sexual maturity Paniw et al. 2020
                   n.vr=28, # number of vital rates with covariates
                   n.pam=1422, # number of total parameters of all these vital rates
                   dens=1, # density dependence in it?
                   biotic_interactions=0, # any biotic interactions?
                   lambda.sim=0,
                   study.length=20)

# subadults
Results_subadults=data.frame(study.doi="10.1126/science.aau5905",
                             year.of.publication="2019",
                             group="Mammals",
                             species="Suricata suricatta",
                             continent="Africa",
                             driver="rain",
                             driver.type="C",
                             stage.age="subadult",
                             vital.rates=rep(c("growth","survival"),each=50),
                             sens=c(Sub.Gr.R,
                                    Sub.Surv.R),
                             mat=1, # age at sexual maturity Paniw et al. 2020
                   n.vr=28, # number of vital rates with covariates
                   n.pam=1422, # number of total parameters of all these vital rates
                   dens=1, # density dependence in it?
                   biotic_interactions=0, # any biotic interactions?
                   lambda.sim=0,
                   study.length=20)

# pups
Results_pups=data.frame(study.doi="10.1126/science.aau5905",
                        year.of.publication="2019",
                        group="Mammals",
                        species="Suricata suricatta",
                        continent="Africa",
                        driver="rain",
                        driver.type="C",
                        stage.age="pups",
                        vital.rates=rep(c("growth","survival"),each=50),
                        sens=c(Pups.Gr.R,
                               Pups.Surv.R),
                        mat=1, # age at sexual maturity Paniw et al. 2020
                   n.vr=28, # number of vital rates with covariates
                   n.pam=1422, # number of total parameters of all these vital rates
                   dens=1, # density dependence in it?
                   biotic_interactions=0, # any biotic interactions?
                   lambda.sim=0,
                   study.length=20)

# merge all Results df above
Results_all=rbind(Results_domadults,Results_helperadults,Results_subadults,Results_pups)

write.csv(Results_all,"SensVR_Rain_Meerkats.csv",row.names = F)


