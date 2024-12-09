
setwd("/Users/maria/Dropbox/teaching/esin/Rabbits")

# Get covariate data:

data.fine_gcm <- read.csv("data_fine_gcm.csv") 

#subset data:

data.fine_gcm=data.fine_gcm[data.fine_gcm$rcp%in%"historical"&data.fine_gcm$gcm%in%c("GISS-E2-H")&data.fine_gcm$nucleus%in%"Donana-Aljarafe",]

scen=unique(data.fine_gcm$gcm)


########### TEMPERATURE

ab=read.csv("abundance_rabbit_cov_vr_temp.csv")

head(ab)


delta.min=NULL
delta.max=NULL


for(i in 1:100){#for each run
    
    sub=ab[ab$run==i,]
    for(j in 6:15){# pick specific seasons
      
      delta.min=rbind(delta.min,
                           data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="min"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="min"],
                                      run=i))
      delta.max=rbind(delta.max,
                           data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="max"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="max"],
                                      run=i))
    }
}



#Sensitivity
lambda.min=aggregate(delta~run,mean,data=delta.min)$delta
lambda.max=aggregate(delta~run,mean,data=delta.max)$delta

mean.temp=aggregate(ta~yearl+gcm,mean,data=data.fine_gcm)

sd.temp=sd(mean.temp$ta)
max.temp=mean.temp$ta[mean.temp$ta==max(mean.temp$ta)]
min.temp=mean.temp$ta[mean.temp$ta==min(mean.temp$ta)]

sens_vr=abs(lambda.max-lambda.min)/(abs(max.temp-min.temp)/sd.temp)


sens_rabbit_vr_temp=data.frame(species="Oryctolagus cuniculus", study.doi="10.1371/journal.pone.0048988",year.of.publication=2012,
                    group="Mammals",continent="Europe",driver="Temperature",driver.type="C",
                    stage.age="all",vital.rates="reproduction",sens=sens_vr,mat=0.33,n.vr=5,n.pam=19,dens=1,
                    biotic_interactions=0,lambda.sim=1,study.length=NA)


setwd("/Users/maria/Dropbox/teaching/esin/Rabbits")


########### FM juvenile survival

ab=read.csv("abundance_rabbit_cov_vr_Jsurv.csv")

head(ab)


delta.min=NULL
delta.max=NULL


for(i in 1:100){#for each run
  
  sub=ab[ab$run==i,]
  for(j in 6:15){# pick specific seasons
    
    delta.min=rbind(delta.min,
                    data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="min"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="min"],
                               run=i))
    delta.max=rbind(delta.max,
                    data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="max"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="max"],
                               run=i))
  }
}



#Sensitivity
lambda.min=aggregate(delta~run,mean,data=delta.min)$delta
lambda.max=aggregate(delta~run,mean,data=delta.max)$delta

max.FM.year=aggregate(Fm~yearl+gcm,max,data=data.fine_gcm)
sd.FM=sd(data.fine_gcm$Fm)

max.Fm=max(max.FM.year$Fm)
min.Fm=min(max.FM.year$Fm)

sens_vr=abs(lambda.max-lambda.min)/(abs(max.temp-min.temp)/sd.temp)


sens_rabbit_vr_fm_J=data.frame(species="Oryctolagus cuniculus", study.doi="10.1371/journal.pone.0048988",year.of.publication=2012,
                               group="Mammals",continent="Europe",driver="# dry months",driver.type="A",
                               stage.age="juvenile",vital.rates="survival",sens=sens_vr,mat=0.33,n.vr=5,n.pam=19,dens=1,
                               biotic_interactions=0,lambda.sim=1,study.length=NA)


########### FM adult survival

ab=read.csv("abundance_rabbit_cov_vr_Asurv.csv")

head(ab)


delta.min=NULL
delta.max=NULL


for(i in 1:100){#for each run
  
  sub=ab[ab$run==i,]
  for(j in 6:15){# pick specific seasons
    
    delta.min=rbind(delta.min,
                    data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="min"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="min"],
                               run=i))
    delta.max=rbind(delta.max,
                    data.frame(delta=sub$densf[sub$month==6&sub$year==(j+1)&sub$temp.pert=="max"]/sub$densf[sub$month==6&sub$year==(j)&sub$temp.pert=="max"],
                               run=i))
  }
}



#Sensitivity
lambda.min=aggregate(delta~run,mean,data=delta.min)$delta
lambda.max=aggregate(delta~run,mean,data=delta.max)$delta

max.FM.year=aggregate(Fm~yearl+gcm,max,data=data.fine_gcm)
sd.FM=sd(data.fine_gcm$Fm)

max.Fm=max(max.FM.year$Fm)
min.Fm=min(max.FM.year$Fm)

sens_vr=abs(lambda.max-lambda.min)/(abs(max.temp-min.temp)/sd.temp)


sens_rabbit_vr_fm_A=data.frame(species="Oryctolagus cuniculus", study.doi="10.1371/journal.pone.0048988",year.of.publication=2012,
                               group="Mammals",continent="Europe",driver="# dry months",driver.type="A",
                               stage.age="adult",vital.rates="survival",sens=sens_vr,mat=0.33,n.vr=5,n.pam=19,dens=1,
                               biotic_interactions=0,lambda.sim=1,study.length=NA)

sens_rabbit_vr=rbind(sens_rabbit_vr_temp,sens_rabbit_vr_fm_J,sens_rabbit_vr_fm_A)

write.csv(sens_rabbit_vr,"sens_rabbit_vr.csv",row.names = F)

