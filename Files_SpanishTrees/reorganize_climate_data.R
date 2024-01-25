# Climate data of several tree species of spain
# 08.08.2023
# Esin Ickin

# clean environment
rm(list=ls())

# load R packages
library(tidyverse)

# set wd
setwd("~/Files_SpanishTrees")

# read climate data (scenario A2)
mean_2010 <- read.table(file="data/clima/MEAN_2010.csv",header=T,sep=";",dec=".")

mean_A2_2020 <- read.table(file="data/clima/MEAN_A2_2020.csv",header=T,sep=";",dec=".")

mean_A2_2030 <- read.table(file="data/clima/MEAN_A2_2030.csv",header=T,sep=";",dec=".")

mean_A2_2040 <- read.table(file="data/clima/MEAN_A2_2040.csv",header=T,sep=";",dec=".")

mean_A2_2050 <- read.table(file="data/clima/MEAN_A2_2050.csv",header=T,sep=";",dec=".")

mean_A2_2060 <- read.table(file="data/clima/MEAN_A2_2060.csv",header=T,sep=";",dec=".")

mean_A2_2070 <- read.table(file="data/clima/MEAN_A2_2070.csv",header=T,sep=";",dec=".")

mean_A2_2080 <- read.table(file="data/clima/MEAN_A2_2080.csv",header=T,sep=";",dec=".")

mean_A2_2090 <- read.table(file="data/clima/MEAN_A2_2090.csv",header=T,sep=";",dec=".")

# merge the climate data of each year
climate_all_years <- rbind(mean_2010,mean_A2_2020,mean_A2_2030,mean_A2_2040,mean_A2_2050,mean_A2_2060,mean_A2_2070,mean_A2_2080,mean_A2_2090)

# take mean, standard deviation, max, min of temperature, precipitation, temperature anomaly, precipitation anomaly per ID (patch) across time
summary_climate <- climate_all_years %>%
  group_by(ID) %>%
  summarize(mean_temp = mean(temp),
            sd_temp = sd(temp),
            max_temp = max(temp),
            min_temp = min(temp),
            mean_precip=mean(precip),
            sd_precip=sd(precip),
            max_precip=max(precip),
            min_precip=min(precip),
            mean_anom_temp=mean(anom_temp),
            sd_anom_temp=sd(anom_temp),
            max_anom_temp=max(anom_temp),
            min_anom_temp=min(anom_temp),
            mean_anom_pl=mean(anom_pl),
            sd_anom_pl=sd(anom_pl),
            max_anom_pl=max(anom_pl),
            min_anom_pl=min(anom_pl))

# calculate covariation
# which is the value of a covariate when another covariate is at its maximum: 
# Precip_when_Temp_max=cov$Precip[which(cov$Temp)==max(cov$Temp)][1]
# for each cell
# across years

climate_covariation <- climate_all_years %>%
  group_by(ID) %>%
  summarise(Precip_when_Temp_max = precip[which.max(temp)],
            Precip_when_Temp_min = precip[which.min(temp)],
            AnomTemp_when_Temp_max = anom_temp[which.max(temp)],
            AnomTemp_when_Temp_min = anom_temp[which.min(temp)],
            AnomPrecip_when_Temp_max = anom_pl[which.max(temp)],
            AnomPrecip_when_Temp_min = anom_pl[which.min(temp)],
            
            Temp_when_Precip_max = temp[which.max(precip)],
            Temp_when_Precip_min = temp[which.min(precip)],
            AnomTemp_when_Precip_max = anom_temp[which.max(precip)],
            AnomTemp_when_Precip_min = anom_temp[which.min(precip)],
            AnomPrecip_when_Precip_max = anom_pl[which.max(precip)],
            AnomPrecip_when_Precip_min = anom_pl[which.min(precip)],
            
            Precip_when_AnomTemp_max = precip[which.max(anom_temp)],
            Precip_when_AnomTemp_min = precip[which.min(anom_temp)],
            Temp_when_AnomTemp_max = temp[which.max(anom_temp)],
            Temp_when_AnomTemp_min = temp[which.min(anom_temp)],
            AnomPrecip_when_AnomTemp_max = anom_pl[which.max(anom_temp)],
            AnomPrecip_when_AnomTemp_min = anom_pl[which.min(anom_temp)],

            Precip_when_AnomPrecip_max = precip[which.max(anom_pl)],
            Precip_when_AnomPrecip_min = precip[which.min(anom_pl)],
            Temp_when_AnomPrecip_max = temp[which.max(anom_pl)],
            Temp_when_AnomPrecip_min = temp[which.min(anom_pl)],
            AnomTemp_when_AnomPrecip_max = anom_temp[which.max(anom_pl)],
            AnomTemp_when_AnomPrecip_min = anom_temp[which.min(anom_pl)]            
            )

# join summary climate and climate covariation
summary_climate=right_join(summary_climate,climate_covariation)
write.csv(summary_climate,"ClimateCovariates.csv",row.names = F)


