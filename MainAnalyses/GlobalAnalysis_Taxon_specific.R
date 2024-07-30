###
# This is a script containing the main global analysis 
# But we fitted the GLMM for each taxa separately

# Author of this script: Esin Ickin
# Date: 30.07.2024
###


# 0. Prepare session ############################################

rm(list=ls())

# set wd
setwd("~/Documents/Master Thesis/pert_analyses")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)

# 1. Load data ###########################################

# load df
df=read.csv("AllSens.csv")
length(levels(factor(df$species))) # 36 species

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

# log-transform age at sex. maturity
df$mat=log(df$mat)
# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)


climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","Winterlength")) # remove SAM and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
#climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any


# 3. Mammals ####################################################

## 3.1 filter df ######################

# filter for mammals
mammals_df=climate_df[climate_df$group=="Mammals",]

mammals_df$driver=factor(mammals_df$driver)
mammals_df$dens=factor(mammals_df$dens)
mammals_df$biotic_interactions=factor(mammals_df$biotic_interactions)
mammals_df$cov=factor(mammals_df$cov) #is there covariation 0/1
mammals_df$study.doi=factor(mammals_df$study.doi)
mammals_df$species=factor(mammals_df$species)

# check species, drivers, etc.
levels(mammals_df$species)
levels(mammals_df$driver)

## 3.2 fit GLMMM ########################

m1 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|species) , family = Gamma(link="log"), data = mammals_df)

summary(m1)

## 3.3 Plot #################

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(mammals_df$mat), max(mammals_df$mat), length.out = 1000), cov = levels(mammals_df$cov),dens=levels(mammals_df$dens))))

levels(pred.df$dens) <- c("No density effects","Density effects",NA)

plot.clim.mammals <- ggplot(pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=0.8) +
  facet_grid(dens ~ .) +
  scale_fill_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) + 
  scale_color_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) +
  labs(
    x = "Log generation time (year)",
    y = "Scaled population growth sensitivities (|S|)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18),
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

plot.clim.mammals
ggsave("SensClimMammals.png", plot.clim.mammals, width = 10, height = 8, dpi = 300)


# 4. Birds ####################################################

## 4.1 filter df ######################

# filter for birds
birds_df=climate_df[climate_df$group=="Birds",]

birds_df$driver=factor(birds_df$driver)
birds_df$dens=factor(birds_df$dens)
birds_df$biotic_interactions=factor(birds_df$biotic_interactions)
birds_df$cov=factor(birds_df$cov) #is there covariation 0/1
birds_df$study.doi=factor(birds_df$study.doi)
birds_df$species=factor(birds_df$species)

# check species, drivers, etc.
levels(birds_df$species)
levels(birds_df$driver)

## 4.2 fit GLMMM ########################

m2 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|species) , family = Gamma(link="log"), data = birds_df)

summary(m2)

## 4.3 Plot #################

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m2,xlevels=list(mat=seq(min(birds_df$mat), max(birds_df$mat), length.out = 1000), cov = levels(birds_df$cov),dens=levels(birds_df$dens))))

levels(pred.df$dens) <- c("No density effects","Density effects",NA)

plot.clim.birds <- ggplot(pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=0.8) +
  facet_grid(dens ~ .) +
  scale_fill_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) + 
  scale_color_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) +
  labs(
    x = "Log generation time (year)",
    y = "Scaled population growth sensitivities (|S|)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18),
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

plot.clim.birds

ggsave("SensClimBirds.png", plot.clim.birds, width = 10, height = 8, dpi = 300)


# 5. Plants ####################################################

## 5.1 filter df ######################

# filter for plants
plants_df=climate_df[climate_df$group=="Plants",]

plants_df$driver=factor(plants_df$driver)
plants_df$dens=factor(plants_df$dens)
plants_df$biotic_interactions=factor(plants_df$biotic_interactions)
plants_df$cov=factor(plants_df$cov) #is there covariation 0/1
plants_df$study.doi=factor(plants_df$study.doi)
plants_df$species=factor(plants_df$species)

# check species, drivers, etc.
levels(plants_df$species)
levels(plants_df$driver)

## 4.2 fit GLMMM ########################

m3 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|species) , family = Gamma(link="log"), data = plants_df)

summary(m3)

length(levels(factor(plants_df[plants_df$dens==0,]$species)))

length(levels(factor(plants_df[plants_df$dens==1,]$species)))

## 4.3 Plot #################

# NEW PLOT
# with data points added as mean sensitivities of each species
mean.plants_df <- plants_df %>%
  group_by(species,mat,cov,dens,group) %>%
  summarise(sens = mean(sens))

climate_df$dens=factor(climate_df$dens)
levels(climate_df$dens)
levels(mean.plants_df$dens)

levels(climate_df$dens) <- c("No density effects","Density effects",NA)
levels(mean.plants_df$dens) <- c("No density effects","Density effects",NA)

library(ggrepel)



# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m3,xlevels=list(mat=seq(min(plants_df$mat), max(plants_df$mat), length.out = 1000), cov = levels(plants_df$cov),dens=levels(plants_df$dens))))

levels(pred.df$dens) <- c("No density effects","Density effects",NA)

plot.clim.plants <- ggplot(pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=3) +

  geom_jitter(data = mean.plants_df, aes(x = mat, y = sens, color = cov), alpha = 0.8, size = 6, width = 0, height = 0) +
  
  # if you want labels included: 
  
# geom_text_repel(data = mean.plants_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
 #            box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #        max.overlaps = 100) +
  
  
  facet_grid(dens ~ .) +
  scale_fill_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) + 
  scale_color_manual(name = "Covariation:",values = c("#b7b3e6","#029356"), labels = c("No","Yes")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18),
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

plot.clim.plants


ggsave("SensClimPlants.png", plot.clim.plants, width = 10, height = 8, dpi = 300)

