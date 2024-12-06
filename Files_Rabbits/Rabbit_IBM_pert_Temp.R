
rm(list = ls())

# PACKAGES
library(dbplyr)
library(ggplot2)
library(tidyverse)
library(extraDistr)
library(boot) # inv.logit function
library(lhs)
library(reshape)

# Record slurm objects if run on a cluster
# args=commandArgs(TRUE)
# foldNumber <- as.integer(args[1]) # foldNumber<-1
# totalFold <- as.integer(args[2]) # totalFold<-69
# setwd(args[3])

setwd("/Users/maria/Dropbox/teaching/esin/Rabbits/")

#----------------------------------------------------------
# Set vital rate parameters
#----------------------------------------------------------

params <- data.frame(AA = 9.123643,
                     AJ = 4.490206,
                     ds = 0.2180814,
                     V = 17.79108,
                     MN = 0.8221069,
                     r4 = -2.575,
                     r6 = 2.106802,
                     r9 = 4.741612,
                     dr = 4.176849,
                     Ls = 3.2,
                     Ld = 3.9,
                     Lr = 3.5,
                     # Krab=25,
                     Krab = 14)


#----------------------------------------------------------
# Set up simulation parameters
#----------------------------------------------------------

data.fine_gcm <- read.csv("data_fine_gcm.csv") 

#subset data:

data.fine_gcm=data.fine_gcm[data.fine_gcm$rcp%in%"historical"&data.fine_gcm$gcm%in%c("GISS-E2-H")&data.fine_gcm$nucleus%in%"Donana-Aljarafe",]

rcp <- "historical"

#----------------------------------------------------------
# Create vital rate functions
#----------------------------------------------------------

# survival for newborns
PSB <- function(MN){
  1-MN
}

# survival for juveniles
PSJ <- function(n,Fm){
  S1=exp(params$AJ-(params$ds * (n/params$Krab)))
  S2=1-exp(-params$V/Fm)
  S=(S1/(1+S1))*S2
  return(S^30.4)
}

# survival for adults
PSA <- function(n,Fm){
  S1=exp(params$AA-(params$ds*(n/params$Krab)))
  S2=1-exp(-params$V/Fm)
  S=(S1/(1+S1))*S2
  return(S^30.4)
}

# breeding season
PB <- function(Ta,D,DI,W){
  B1=4.542-(0.605*Ta)+(0.029*(Ta^2))-(0.006*D)-(0.017*DI)-W
  return(1/(1+exp(B1)))
}

# reproduction
PR <- function(n,rA){
  R1=rA-params$dr*(n/params$Krab)
  return(exp(R1)/(1+exp(R1)))
}


#----------------------------------------------------------
# Set up IBM objects
#----------------------------------------------------------

years = c(1:16)
months = c(1:12)

#hold rabbit

dens.rabbit=NULL


# Rabbit initial populations (can be generate as commented out or loaded)
# temp.conejos.T=data.frame(ID=1:30,
#                           edad=round(rtnorm(30,17,5,a=0,b=33),0))
# temp.conejos.T$estado=NA
# temp.conejos.T$sexo=sample(c("M","F"),size=nrow(temp.conejos.T),replace = TRUE)
# temp.conejos.T=temp.conejos.T[,c("ID","edad","estado","sexo")]
# temp.conejos.T$estado[temp.conejos.T$edad%in%(0)]="B"
# temp.conejos.T$estado[temp.conejos.T$edad%in%(1:3)]="J"
# temp.conejos.T$estado[temp.conejos.T$edad>=4]="A"

temp.conejos.T=read.csv("initial_population.csv")
temp.conejos.T$nucleus="Donana-Aljarafe"

scen=unique(data.fine_gcm$gcm)
### Get min, max, mean of covariates per GCM

#mean temperature per year:

mean.temp=aggregate(ta~yearl+gcm,mean,data=data.fine_gcm)
sd.temp=sd(mean.temp$ta)

min.temp=NULL
max.temp=NULL
mean.Fm=list(NULL)

for(x in 1: length(scen)){
  
  sub=mean.temp[mean.temp$gcm%in%scen[x],]
  
  sub2=data.fine_gcm[data.fine_gcm$gcm%in%scen[x],]
  min.temp[x]=as.numeric(as.character(sub$yearl[sub$ta==min(sub$ta)]))
  max.temp[x]=as.numeric(as.character(sub$yearl[sub$ta==max(sub$ta)]))
  mean.Fm[[x]]= round(aggregate(Fm~month,mean,data=sub2))
  mean.Fm[[x]]$WL=0
  mean.Fm[[x]]$WL[mean.Fm[[x]]$Fm==0]=(-1.592)
 
}

name.temp=c("min","max")

#----------------------------------------------------------
# IBM
#----------------------------------------------------------
repeats = c(1:100)

for (r in rcp) {
  for (s in scen) {
    for (rep in repeats) {
      
      print(paste("Running GCM ",s, ", iteration ", rep))
      
      for(p in 1:2){
        ibm.rabbit=NULL
        for(y in years){
          for(m in months){
            for(nuc in 1){
              if(m==1&y==years[1]){
                #rabbit pop
                temp.conejos <- temp.conejos.T %>% mutate(gcm=s,rcp="historical")
                maxID_rab <- max(temp.conejos$ID)
              } else { 
                if(m==1&y>years[1]) {
                  #rabbit pop
                  temp.conejos=ibm.rabbit %>% 
                    filter(nucleus%in%"Donana-Aljarafe"&gcm%in%s&rcp%in%"historical"&run%in%rep&year%in%(y-1)&month%in%12) # Antes: month%in%time.vec[time-1]
                  ibm.rabbit <- ibm.rabbit[!(ibm.rabbit$nucleus%in%"Donana-Aljarafe"&ibm.rabbit$gcm%in%s&ibm.rabbit$rcp%in%"historical"&ibm.rabbit$run%in%rep&ibm.rabbit$year%in%(y-1)&ibm.rabbit$month%in%12),]
                } else {
                  #rabbit pop
                  temp.conejos=ibm.rabbit %>% 
                    filter(nucleus%in%"Donana-Aljarafe"&gcm%in%s&rcp%in%"historical"&run%in%rep&year%in%y&month%in%(m-1))
                  ibm.rabbit <- ibm.rabbit[!(ibm.rabbit$nucleus%in%"Donana-Aljarafe"&ibm.rabbit$gcm%in%s&ibm.rabbit$rcp%in%"historical"&ibm.rabbit$run%in%rep&ibm.rabbit$year%in%y&ibm.rabbit$month%in%(m-1)),]
                }}
              
              
              if(nrow(temp.conejos) != 0){
                
                #update rabbit variables
                temp.conejos$edad <- temp.conejos$edad + 1
                temp.conejos$estado[temp.conejos$edad%in%(1:3)]="J"
                temp.conejos$estado[temp.conejos$edad>=4]="A"
                #temp.conejos$estado <- as.factor(temp.conejos$estado)
                temp.conejos$month=m
                temp.conejos$year=y
                temp.conejos$run=rep
                temp.conejos$surv=NA
                temp.conejos$offspring=NA
                temp.conejos$densi=NA
                temp.conejos$densf=NA
                temp.conejos$bs=NA
                temp.conejos$Kits=NA
                temp.conejos$Juveniles=NA
                temp.conejos$Adults=NA
                #temp.conejos$nucleus="Donana-Aljarafe"
                date=data.fine_gcm$date[data.fine_gcm$yearl %in% y & data.fine_gcm$month %in% m &
                                          data.fine_gcm$gcm %in% s & 
                                          data.fine_gcm$rcp %in% "historical" &
                                          data.fine_gcm$nucleus %in% "Donana-Aljarafe"]
                temp.conejos$date=date
                dens=30
                # dens=50
                densV2=as.numeric(nrow(temp.conejos))
                temp.conejos$densi=dens #initial density for this month
                
                # obtain variables needed for function of survival and repr success
                FmL=mean.Fm[[which(scen==s)]]$Fm[m]
                WL=mean.Fm[[which(scen==s)]]$WL[m]
                if(p==1){
                  TaL=data.fine_gcm$ta[data.fine_gcm$yearl%in%min.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"] 
                  DL=data.fine_gcm$mdl[data.fine_gcm$yearl%in%min.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"]
                  DIL=data.fine_gcm$DI[data.fine_gcm$yearl%in%min.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"]
                  
                }else{
                  TaL=data.fine_gcm$ta[data.fine_gcm$yearl%in%max.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"] 
                  DL=data.fine_gcm$mdl[data.fine_gcm$yearl%in%max.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"]
                  DIL=data.fine_gcm$DI[data.fine_gcm$yearl%in%max.temp[which(scen==s)]&data.fine_gcm$month%in%m&data.fine_gcm$gcm%in%s&data.fine_gcm$rcp%in%"historical"&data.fine_gcm$nucleus%in%"Donana-Aljarafe"]
                  
                }
                
                #pup recruitment
                temp.conejos$offspring[temp.conejos$repr%in%1]=rtpois(n=nrow(temp.conejos[temp.conejos$repr%in%1, ]),lambda=params$Ld,a=0,b=10) # no pups for each rabbit female that reproduced the previous month and survived the current month
                born.babies=as.numeric(sum(temp.conejos$offspring,na.rm = T))
                #add pups to the population
                if(born.babies>0){
                  data.babies=data.frame(ID=(maxID_rab+1):(maxID_rab+born.babies),gcm=s,rcp="historical",edad=0,
                                         estado="B",sexo=sample(c("M","F"),size=born.babies,replace = TRUE),
                                         month=m,run=rep,surv=NA,offspring=NA,repr=NA,year=y,densi=dens,densf=NA,
                                         bs=NA,date=date,nucleus="Donana-Aljarafe",
                                         Kits=Kits,Juveniles=Juveniles,Adults=Adults) #nucleus=coord_points$Nucleus[point],point=coord_points$Point[point]
                  temp.conejos <- rbind(temp.conejos,data.babies)
                  maxID_rab <- max(data.babies$ID)
                }
                
                #newborn survival
                temp.conejos$surv[temp.conejos$estado%in%"B"]=rbinom(nrow(temp.conejos[temp.conejos$estado%in%"B",]),1,PSB(params$MN))
                
                #Adult survival
                temp.conejos$surv[temp.conejos$estado%in%"A"]=rbinom(nrow(temp.conejos[temp.conejos$estado%in%"A", ]),1,PSA(dens,FmL))
                temp.conejos$surv[temp.conejos$estado%in%"J"]=rbinom(nrow(temp.conejos[temp.conejos$estado%in%"J", ]),1,PSJ(dens,FmL))
                
                #reproduction
                #is the current month inside the breeding season?
                is.BS=PB(TaL,DL,DIL,WL)
                ifelse(is.BS<0.5,(temp.conejos$bs=0),(temp.conejos$bs=1))
                temp.conejos$repr=NA
                #if it is not inside the breeding season -> no females reproduce
                #if it is inside the breeding season ->determine what females reproduce
                if(is.BS<0.5){
                  temp.conejos$repr=NA
                }else{
                  temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                      temp.conejos$edad<6&temp.conejos$surv==1]=
                    rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                             &temp.conejos$estado%in%"A"&temp.conejos$edad<6
                                             &temp.conejos$surv%in%1, ]),1,PR(dens,params$r4))
                  temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                      temp.conejos$edad%in%6:9&temp.conejos$surv==1]=
                    rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                             &temp.conejos$estado%in%"A"&temp.conejos$edad%in%6:9
                                             &temp.conejos$surv%in%1, ]),1,PR(dens,(params$r4+params$r6)))
                  temp.conejos$repr[temp.conejos$sexo%in%"F"&temp.conejos$estado%in%"A"&
                                      temp.conejos$edad>9&temp.conejos$surv==1]=
                    rbinom(nrow(temp.conejos[temp.conejos$sexo%in%"F"
                                             &temp.conejos$estado%in%"A"&temp.conejos$edad>9
                                             &temp.conejos$surv%in%1, ]),1,PR(dens,(params$r4+params$r9)))
                }
                
                # update rabbit variables
                densf <- nrow(temp.conejos[temp.conejos$surv%in%1, ])
                temp.conejos$densf= densf #final density for this month
                Kits= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"B" & temp.conejos$surv%in%1,]))
                temp.conejos$Kits=Kits
                Juveniles= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"J" & temp.conejos$surv%in%1,]))
                temp.conejos$Juveniles=Juveniles
                Adults= as.numeric(nrow(temp.conejos[temp.conejos$estado%in%"A" & temp.conejos$surv%in%1,]))
                temp.conejos$Adults=Adults
                
                # create dens.temp
                dens.temp<-data.frame(scen=s,rcp="historical",year=y,month=m,densi=dens,densf=densf,
                                      Kits=Kits,Juveniles=Juveniles,Adults=Adults,date=date,
                                      nucleus="Donana-Aljarafe",run=rep,temp.pert=name.temp[p]) 
                
                # add rabbit population to ibm.rabbit and add rabbit density to dens.rabbit
                dens.rabbit <- rbind(dens.rabbit,dens.temp) # include initial and final no of rabbits, and the no of kits, juveniles and adults in the current month (include all individuals at the beginning of the month plus newborns regardless whether they survive or not to the current month)
                ibm.rabbit <- rbind(ibm.rabbit,temp.conejos) # include the information associated to all the previous individuals (e.g. if the individual survived or not)

              }
              
              
            }
            if(sum(dens.rabbit[dens.rabbit$run %in% rep & dens.rabbit$scen %in% s & dens.rabbit$rcp %in% "historical" & dens.rabbit$year %in% y & dens.rabbit$month %in% m,"densf"]) %in% 0) break # break the loop if no lynxes nor rabbits survived in any nuclei for the current run, climate scenario, year and month
          }
          if(sum(dens.rabbit[dens.rabbit$run %in% rep & dens.rabbit$scen %in% s & dens.rabbit$rcp %in% "historical" & dens.rabbit$year %in% y & dens.rabbit$month %in% m,"densf"]) %in% 0) break # break the loop if no lynxes nor rabbits survived in any nuclei for the current run, climate scenario, year and month
        }
        
        rm(ibm.rabbit)

      }
      
    }
  }
}

write.csv(dens.rabbit,"abundance_rabbit_noCov.csv",row.names = F)

library(ggplot2)

ggplot(dens.rabbit,aes(date,densf,group = run,colour = run))+
  geom_line()+
  facet_grid(temp.pert~.)

