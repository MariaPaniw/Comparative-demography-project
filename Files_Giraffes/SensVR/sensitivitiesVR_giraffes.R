
rm(list=ls())

setwd("/Users/esin/Library/CloudStorage/OneDrive-UniversitätZürichUZH/Master Thesis/pert_analyses/Giraffes")

# rainfall perturbed

# survival of calves (neo) ############################

SurvNeo=read.csv("abund.rain.survneo.csv")

head(SurvNeo)

# Total abundance per season, run, and max/min rain
ab.tot = aggregate(ID~TIME.SIM+run+rain, sum, data=SurvNeo)

head(ab.tot)

# Get changes in abundace per year (for short-rain season to next); the time steps here are 4-month seasons
# discard first 12 season for transient dynamics
count=seq(1,148,length.out=50)

min.survneo=NULL
max.survneo=NULL

for(i in 1:100){#for each run
  
  sub=ab.tot[ab.tot$run==i,]
  for(j in 1:(length(count)-1)){# pick specific seasons
    
    min.survneo=rbind(min.survneo,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="min"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="min"],
                                    run=i))
    max.survneo=rbind(max.survneo,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="max"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="max"],
                                    run=i))
  }
}

#Sensitivity
lambda.min=aggregate(delta~run,mean,data=min.survneo)$delta
lambda.max=aggregate(delta~run,mean,data=max.survneo)$delta

cov.data=read.csv("abund_validation2023.csv")  # load observed rainfall anomalies

min.rain=min(cov.data$rain)
max.rain=max(cov.data$rain)

Sens_SurvNeo=abs(lambda.max-lambda.min)/(abs(max.rain-min.rain)/0.62)



# survival of juveniles ###################################

SurvJ=read.csv("abund.rain.survj.csv")

head(SurvJ)

# Total abundance per season, run, and max/min rain
ab.tot = aggregate(ID~TIME.SIM+run+rain, sum, data=SurvJ)

head(ab.tot)

# Get changes in abundace per year (for short-rain season to next); the time steps here are 4-month seasons
# discard first 12 season for transient dynamics
count=seq(1,148,length.out=50)

min.survj=NULL
max.survj=NULL

for(i in 1:100){#for each run
  
  sub=ab.tot[ab.tot$run==i,]
  for(j in 1:(length(count)-1)){# pick specific seasons
    
    min.survj=rbind(min.survj,
                      data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="min"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="min"],
                                 run=i))
    max.survj=rbind(max.survj,
                      data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="max"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="max"],
                                 run=i))
  }
}

#Sensitivity
lambda.min=aggregate(delta~run,mean,data=min.survj)$delta
lambda.max=aggregate(delta~run,mean,data=max.survj)$delta


Sens_SurvJ=abs(lambda.max-lambda.min)/(abs(max.rain-min.rain)/0.62)

  
# survival of adults #########################################

SurvA=read.csv("abund.rain.survad.csv")

head(SurvA)

# Total abundance per season, run, and max/min rain
ab.tot = aggregate(ID~TIME.SIM+run+rain, sum, data=SurvA)

head(ab.tot)

# Get changes in abundace per year (for short-rain season to next); the time steps here are 4-month seasons
# discard first 12 season for transient dynamics
count=seq(1,148,length.out=50)

min.surva=NULL
max.surva=NULL

for(i in 1:100){#for each run
  
  sub=ab.tot[ab.tot$run==i,]
  for(j in 1:(length(count)-1)){# pick specific seasons
    
    min.surva=rbind(min.surva,
                    data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="min"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="min"],
                               run=i))
    max.surva=rbind(max.surva,
                    data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="max"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="max"],
                               run=i))
  }
}

#Sensitivity
lambda.min=aggregate(delta~run,mean,data=min.surva)$delta
lambda.max=aggregate(delta~run,mean,data=max.surva)$delta


Sens_SurvA=abs(lambda.max-lambda.min)/(abs(max.rain-min.rain)/0.62)


# compare

final=rbind(data.frame(sens=Sens_SurvNeo,stage="calves"),
            data.frame(sens=Sens_SurvJ,stage="juveniles"),
            data.frame(sens=Sens_SurvA, stage="adults"))

library(ggplot2)

ggplot(final,aes(x=sens,fill=stage))+
  geom_density(alpha=0.5)


# Save output
Results=data.frame(species="Giraffa camelopardalis",
                   study.doi="10.1111/gcb.16970",
                   year.of.publication="2023",
                   group="Mammals",
                   continent="Africa",
                   driver="Rainfall",
                   driver.type="C",
                   stage.age=rep(c("calves","juveniles","adults")),
                   vital.rates="survival",
                   sens=c(Sens_SurvNeo,Sens_SurvJ,Sens_SurvA),
                   mat=6, # Bond et al. 2023
                   n.vr=3, # number of vital rates with covariates
                   n.pam=12, # number of total parameters of these vital rates
                   dens=1,
                   biotic_interactions=0,
                   lambda.sim=1, # was lambda calculated analytically (0) or using simulation (1)?
                   study.length=8
                   )

write.csv(Results,"SensVR_Giraffes.csv",row.names = F)


