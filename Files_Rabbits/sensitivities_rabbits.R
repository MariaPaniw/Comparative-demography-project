
setwd("/Users/maria/Dropbox/teaching/esin/Rabbits")

# Get covariate data:

data.fine_gcm <- read.csv("data_fine_gcm.csv") 

#subset data:

data.fine_gcm=data.fine_gcm[data.fine_gcm$rcp%in%"historical"&data.fine_gcm$gcm%in%c("GISS-E2-H")&data.fine_gcm$nucleus%in%"Donana-Aljarafe",]

scen=unique(data.fine_gcm$gcm)


#temperature perturbed, no covariation: 

ab=read.csv("abundance_rabbit_noCov.csv")

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

no.cov=abs(lambda.max-lambda.min)/(abs(max.temp-min.temp)/sd.temp)

no.cov.l.rat=abs(log(lambda.max/lambda.min))

####################################
#temperature perturbed, covariation: 

ab=read.csv("abundance_rabbit_cov.csv")

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

cov=abs(lambda.max-lambda.min)/(abs(max.temp-min.temp)/sd.temp)

cov.l.rat=abs(log(lambda.max/lambda.min))

# compare

final=rbind(data.frame(sens=cov,covar=1),
            data.frame(sens=no.cov,covar=0))

final.l.rat=rbind(data.frame(sens=cov.l.rat,covar=1),
                  data.frame(sens=no.cov.l.rat,covar=0))

sens_rabbit=data.frame(species="Oryctolagus cuniculus", study.doi="10.1371/journal.pone.0048988",year.of.publication=2012,
                    group="Mammals",continent="Europe",driver="Temperature",driver.type="C",
                    stage.age="all",vital.rates="all",sens=final$sens,cov=final$covar,mat=0.33,n.vr=5,n.pam=19,dens=1,
                    biotic_interactions=0,lambda.sim=1,study.length=NA,l_ratio=final.l.rat$sens)

write.csv(sens_rabbit,"sens_rabbit.csv",row.names = F)



library(ggplot2)

ggplot(final,aes(x=sens,fill=as.factor(covar)))+
  geom_density(alpha=0.5)

ggplot(final.l.rat,aes(x=sens,fill=as.factor(covar)))+
  geom_density(alpha=0.5)

#### Plot convergence:

library(ggplot2)

delta.max$year=rep(c(6:15),length(unique(delta.max$run)))
delta.max$run=factor(delta.max$run)


ggplot(delta.max,aes(x=year,y=delta,col=run)) +
  geom_line(alpha=0.5)+guides(col="none")+
  labs(y=expression(paste("Annual growth rate (",N[t+1]/N[t],")" )),
       x="Year")+
  xlim(6,15.7)+
  theme_bw(base_size=18)+
  theme(legend.title = element_blank(),
        # legend.text = element_text(size=20),
        legend.position = c(0.8,0.87),
        strip.background = element_blank())

ggsave(filename="rabbits_covergence_Annual_growth.pdf",width=7, height=5)



