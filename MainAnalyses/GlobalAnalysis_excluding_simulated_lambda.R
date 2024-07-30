###
# This is a script containing the global analysis
# but only including lambda that were analytically calculated
# thus exluding simulated lambda

# Author of this script: Esin Ickin
# Date: 30.07.2024
###


# 0. Prepare session ############################################

rm(list=ls())

# set wd
setwd("/Users/esinickin/Desktop/MainAnalyses")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)

# 1. Load data #########################################################

# make sure the directory is correct

# here we removed species where lambda was calculated using simulations and λ = (Nt+1/Nt)
# for that see all_studies.xlsx

df=read.csv("AllSens.csv")

df=df[df$lambda.sim==0,]

length(levels(factor(df$species))) # 21
levels(factor(df$species))

#write.csv(df,"NonSimSpecies.csv",row.names = F)
#df=read.csv("NonSimSpecies.csv")

### 2.1 EDA & edit data ##################
str(df)

summary(df)

#hist(df$n.vr)
#hist(df$mat)
#hist(df$n.pam)


sum(is.na(df$sens))
sum(is.infinite(df$sens))

# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform mat
df$mat=log(df$mat)
# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# 3. GLMM: Sens to all Climate Variables ##################################################################

climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))

climate_df=filter(climate_df,!driver %in% c("Winterlength"))

# remove sens == 0 if there is 
min(climate_df$sens)
#climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any

climate_df$driver=factor(climate_df$driver)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

levels(factor(climate_df$driver))
levels(factor(climate_df$driver.type))
levels(factor(climate_df$species))


m1 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m1)



### 3.1 Plot #################

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(climate_df$mat), max(climate_df$mat), length.out = 1000), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))

levels(pred.df$dens) <- c("No density effects","Density effects",NA)

# with data points added
# but less
# one data point per species
# maybe the mean sens per species?
mean.climate_df <- climate_df %>%
  group_by(species,mat,cov,dens,group) %>%
  summarise(sens = mean(sens))

levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density effects","Density effects",NA)
levels(mean.climate_df$dens) <- c("No density effects","Density effects",NA)

library(ggrepel)


plot.clim <- ggplot(pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=3) +
  
  geom_jitter(data = mean.climate_df, aes(x  = mat, y = sens, color = cov), alpha = 0.8, size = 6, width = 0, height = 0) +
  
  geom_text_repel(data = mean.climate_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 2,
             box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
             max.overlaps = 100) +
  
  
  facet_grid(dens ~ .) +
  scale_fill_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) + 
  scale_color_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) +
  labs(
    x = "Log age at sexual maturity ( in years)",
    y = "Scaled population growth sensitivities (|S|)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 22), axis.text = element_text(size = 20),
        strip.text = element_text(size = 20),legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.ticks = element_line(color = "black"),
        legend.position = "bottom") +
  geom_rect(
    data = unique(pred.df[c("dens")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim
ggsave("SensAllClim_NoSims.png", plot.clim, width = 10, height = 8, dpi = 300)
