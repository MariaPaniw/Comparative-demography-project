
rm(list=ls())

setwd("/Users/esinickin/Desktop/Comparative-demography-project-main/Files_Giraffes")

#rainfall perturbed, no covariation: 

rain=read.csv("abund.pert.rain.csv")

head(rain)

# Total abundance per season, run, and max/min rain
ab.tot = aggregate(ID~TIME.SIM+run+rain, sum, data=rain)

head(ab.tot)

# Get changes in abundace per year (for short-rain season to next); the time steps here are 4-month seasons
# discard first 12 season for transient dynamics
count=seq(1,148,length.out=50)

delta.min.rain=NULL
delta.max.rain=NULL

for(i in 1:100){#for each run
  
  sub=ab.tot[ab.tot$run==i,]
  for(j in 1:(length(count)-1)){# pick specific seasons
    
    delta.min.rain=rbind(delta.min.rain,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="min"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="min"],
                                    run=i))
    delta.max.rain=rbind(delta.max.rain,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="max"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="max"],
                                    run=i))
  }
}

#Sensitivity
lambda.min=aggregate(delta~run,mean,data=delta.min.rain)$delta
lambda.max=aggregate(delta~run,mean,data=delta.max.rain)$delta

cov.data=read.csv("abund_validation2023.csv")  # load observed rainfall anomalies

min.rain=min(cov.data$rain)
max.rain=max(cov.data$rain)

no.cov=abs(lambda.max-lambda.min)/(abs(max.rain-min.rain)/0.62)

# l ratios
# abs(log(lambda(mpm.max)/lambda(mpm.min)))
l_ratio.no.cov=abs(log(lambda.max/lambda.min))

####################################
#rainfall perturbed, covariation: 

rain=read.csv("abund.pert.rain.cov.csv")

head(rain)

# Total abundance per season, run, and max/min rain
ab.tot = aggregate(ID~TIME.SIM+run+rain, sum, data=rain)

head(ab.tot)

# Get changes in abundace per year (for short-rain season to next); the time steps here are 4-month seasons
# discard first 12 season for transient dynamics
count=seq(1,148,length.out=50)

delta.min.rain=NULL
delta.max.rain=NULL

for(i in 1:100){#for each run
  
  sub=ab.tot[ab.tot$run==i,]
  for(j in 1:(length(count)-1)){# pick specific seasons
    
    delta.min.rain=rbind(delta.min.rain,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="min"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="min"],
                                    run=i))
    delta.max.rain=rbind(delta.max.rain,
                         data.frame(delta=sub$ID[sub$TIME.SIM==count[j+1]&sub$rain=="max"]/sub$ID[sub$TIME.SIM==count[j]&sub$rain=="max"],
                                    run=i))
  }
}

#Sensitivity
lambda.min=aggregate(delta~run,mean,data=delta.min.rain)$delta
lambda.max=aggregate(delta~run,mean,data=delta.max.rain)$delta

cov.data=read.csv("abund_validation2023.csv")  # load observed rainfall anomalies

min.rain=min(cov.data$rain)
max.rain=max(cov.data$rain)

cov=abs(lambda.max-lambda.min)/(abs(max.rain-min.rain)/0.62)

# l ratios
# abs(log(lambda(mpm.max)/lambda(mpm.min)))
l_ratio.cov=abs(log(lambda.max/lambda.min))


# compare

final=rbind(data.frame(sens=cov,covar="yes"),
            data.frame(sens=no.cov,covar="no"))

library(ggplot2)

ggplot(final,aes(x=sens,fill=covar))+
  geom_density(alpha=0.5)


# Save output

# no.cov = sens with no cov
# cov = sens with cov

Results=data.frame(species="Giraffa camelopardalis",
                   study.doi="10.1111/gcb.16970",
                   year.of.publication="2023",
                   group="Mammals",
                   continent="Africa",
                   driver="Rainfall",
                   driver.type="C",
                   stage.age="all",
                   vital.rates="all",
                   sens=c(no.cov, cov),
                   cov=rep(c(0,1),each=100),
                   mat=6, # Bond et al. 2023
                   n.vr=3, # number of vital rates with covariates
                   n.pam=12, # number of total parameters of these vital rates
                   dens=1,
                   biotic_interactions=0,
                   lambda.sim=1, # was lambda calculated analytically (0) or using simulation (1)?
                   study.length=8,
                   l_ratio=c(l_ratio.no.cov,l_ratio.cov)
                   )

write.csv(Results, "Sens_Giraffes.csv", row.names = F)


