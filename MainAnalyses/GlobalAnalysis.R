###
# This is a script containing two main analyses
# GLMM 1: Sensitivities to all climatic variables (the global analysis)
# GLMM 2: Sensitivities to temperature and rain

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

# 1. Load data #########################################################

# We can import each species' sensitivities separately here or alternatively import the already compiled dataframe with all the species

# Alternatively import all sensitivitesdata
df=read.csv("AllSens.csv")


# Species's sensitivites:

# make sure the directory is correct for each species

# Mammals
AFox=read.csv("ArcticFox/Sens_ArcticFox_MCMC.csv") # [Nater et al. 2021]
StripedMouse=read.csv("StripedMouse/Sens_StripedMouse_MCMC.csv") # [Nater et al. 2018]
Giraffe=read.csv("Giraffes/Sens_Giraffes.csv") # [Bond et al. 2023]
MouseLemur=read.csv("mouse lemurs/Sens_MouseLemurs.csv") # [Ozgul et al. 2023]
Reindeer=read.csv("reindeer/SensReindeer.csv") # [Hansen & Gamelon et al. 2019]
Marmot=read.csv("marmots/final_results/SENS_MARMOTS_MCMC.csv") # [Paniw et al. 2020]
Meerkat=read.csv("meerkats/Sens_Meerkats.csv") # [Paniw et al. 2019]


# Birds
Petrel=read.csv("SkuaPetrel/Output_Sens_MCMC/SensPetrel_MCMC.csv") # [Quéroué et al. 2019]
MPenguin=read.csv("MagellanicPenguin/Sens_MPenguins.csv") # [Clark-Wolf et al. 2022]
Certhia=read.csv("birds/output/Certhia_Sens.csv") # [Malchow et al. 2023]
Linaria=read.csv("birds/output/Linaria_Sens.csv") # [Malchow et al. 2023]
Lophophanes=read.csv("birds/output/Lophophanes_Sens.csv") # [Malchow et al. 2023]
PrunellaCollaris=read.csv("birds/output/PrunellaCollaris_Sens.csv") # [Malchow et al. 2023]
PrunellaModularis=read.csv("birds/output/PrunellaModularis_Sens.csv") # [Malchow et al. 2023]
Pyrrhula=read.csv("birds/output/Pyrrhula_Sens.csv") # [Malchow et al. 2023]
Sitta=read.csv("birds/output/Sitta_Sens.csv") # [Malchow et al. 2023]
Turdus=read.csv("birds/output/Turdus_Sens.csv") # [Malchow et al. 2023]
Dipper=read.csv("dipper/Sens_Dipper_Resampling.csv") # [Gamelon et al. 2017]
Albatross=read.csv("Albatross/Sens_Albatross.csv") # [Jenouvrier et al. 2018]
Goose=read.csv("BarnacleGoose/sens_goose.csv") # [Layton-Matthews et al. 2020]
Jay=read.csv("SiberianJay/sens_jays.csv") # [Layton-Matthews et al. 2018]


# Plants
Cistus=read.csv("shrubs/Output/Sens_Cistus.csv") # [Paniw et al. 2023]
Halimium=read.csv("shrubs/Output/Sens_Halimium.csv") # [Paniw et al. 2023]
Protea=read.csv("Protea/Sens_Protea.csv") # [Merow et al. 2014]
Opuntia=read.csv("OpuntiaImbricata/Sens_OpuntiaImbricata.csv") # [Evers et al. 2020]
Draco=read.csv("Dracocephalum/Sensitivities_Dracocephalum.csv") # [Evers et al. in review]

Trees=read.csv("SpainTrees/SensSpainTrees.csv") # 11 tree species [David Garcia-Callejas et al. 2016]

DewyP=read.csv("/Users/esinickin/Desktop/DewyPines/Sens_DewyPines.csv") # [Conquet et al. in prep]



# 2. Compile one big data file ###########################################

df=rbind(AFox,StripedMouse,Giraffe,MouseLemur,Reindeer,Marmot,Meerkat,
           Petrel,MPenguin,Certhia,Linaria,Lophophanes,PrunellaCollaris,PrunellaModularis,Pyrrhula,Sitta,Turdus,Dipper,
           Cistus,Halimium,Protea,Opuntia,DewyP,Albatross,Goose,Jay,Trees,Draco)
# 
length(levels(factor(df$species))) # 36

#write.csv(df,"AllSens.csv",row.names = F)



### 2.1 EDA & edit data ##################

str(df)
summary(df)

#hist(df$sens)
#hist(df$n.vr)
#hist(df$mat)
#hist(df$n.pam)


sum(is.na(df$sens))
sum(is.infinite(df$sens))

# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform age at sexual maturity
df$mat=log(df$mat)

# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# 3. GLMM: Sensitivities to all climatic variables ##################################################################

# filter df so that we only include sensitivities to climatic variables
climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","Winterlength")) # remove SAM and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see species Marmota flaviventris by Paniw et al. 2020)

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

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))

# GLMM 1
m1 <- glmer(sens ~ cov *dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m1)

# DHARMa ####################
library(DHARMa)

# Simulate residuals
simulationOutput <- simulateResiduals(fittedModel = m1, n = 250)

# Residual vs. Fitted Values Plot
png("residuals_vs_fitted.png", width = 8, height = 6, units = "in", res = 300)
plotResiduals(simulationOutput)
dev.off()

# Q-Q Plot
png("qq_plot.png", width = 8, height = 6, units = "in", res = 300)
plotQQunif(simulationOutput)
dev.off()

# Residuals vs. Key Predictors
png("residuals_vs_cov.png", width = 8, height = 6, units = "in", res = 300)
plotResiduals(simulationOutput, form = climate_df$cov)
dev.off()

png("residuals_vs_dens.png", width = 8, height = 6, units = "in", res = 300)
plotResiduals(simulationOutput, form = climate_df$dens)
dev.off()

# Overall Residual Plot
png("overall_residual_plot.png", width = 8, height = 6, units = "in", res = 300)
plot(simulationOutput)
dev.off()

### 3.1 GLMM: Test effect of study length ###############

# here we want to quickly check whether study length has an effect

m1.2 <- glmer(sens ~ cov *dens + mat + log(study.length) + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m1.2)

mean(climate_df$study.length)


### 3.2 Plot #################

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(climate_df$mat), max(climate_df$mat), length.out = 1000), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))

levels(pred.df$dens) <- c("No density effects","Density effects",NA)

# add data points that show mean sensitivities of each species
mean.climate_df <- climate_df %>%
  group_by(species,mat,cov,dens,group) %>%
  summarise(sens = mean(sens))

levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density effects","Density effects",NA)
levels(mean.climate_df$dens) <- c("No density effects","Density effects",NA)

levels(pred.df$cov) <- c("Observed env change","Simplified env change",NA)
levels(climate_df$cov) <- c("Observed env change","Simplified env change",NA)
levels(mean.climate_df$cov) <- c("Observed env change","Simplified env change",NA)

library(ggrepel) # to add labels

# plot
plot.clim <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = dens), alpha = 0.2) +
  geom_line(data = pred.df, aes(col = dens), linewidth = 2) +
  facet_grid(cov ~ .) +
  
  geom_jitter(data = mean.climate_df, aes(x = mat, y = sens, color = dens), alpha = 0.8, size = 6, width = 0, height = 0) +
  
  # can add labels of species if needed: 
  
  # geom_text_repel(data = mean.climate_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
  #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #             max.overlaps = 100) +
  
  scale_fill_manual(name = "Density feedbacks:", values = c("#b7b3e6", "#029356"), labels = c("No", "Yes")) + 
  scale_color_manual(name = "Density feedbacks:", values = c("#b7b3e6", "#029356"), labels = c("No", "Yes")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)"
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  geom_rect(
    data = unique(pred.df[c("cov")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim

ggsave("SensAllClimSpp.png", plot.clim, width = 10, height = 8, dpi = 300)


summary(m1)

### NOTE MARIA : FROME HERE ON, CODE NEEDS TO BE CHANGED TO ADJUST FOR NEW PLOTTING

# 4. GLMM: Sens to Temperature vs Rain ##################################################################

# split climate further into temperature and rain

levels(factor(climate_df$driver))

# group the drivers into rain and temperature
climate_df$driver[climate_df$driver %in% c("Pbr", "Pwn", "rain","fallR","Precipitation","precipitation","Rain","Rainfall","ROS")]="rain"

climate_df$driver[climate_df$driver %in% c("lagged temperature","prevwinterT","sea ice","SST","SSTA_b","SSTA_b_l","SSTA_m","SSTA_m_l","summerT","Tat","Tbr","temperature","Temperature","Twn","SST1","SST2","SST3","SSTG")]="temperature"

climate_df=filter(climate_df, !driver %in% c("SAM","Q","PET","Winterlength"))
# remove SAM (a driver that influences climate like rain and temperature, so not really a climatic driver)
# remove Q because now I can't sort it into rain or temperature
# remove winterlength again for the same reasons as above

# remove sens == 0 if there is 
min(climate_df$sens)
#climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution

# from chr to factors
climate_df$driver=factor(climate_df$driver)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species) 

length(levels(climate_df$species)) # the only species removed is marmot because we could not sort Q into rain or temperature category

length(levels(climate_df$species))
length(levels(factor((climate_df[climate_df$dens=="No density effects",]$species))))

length(levels(factor((climate_df[climate_df$dens=="Density effects",]$species))))


levels(climate_df$dens)

# GLMM:
m2 <- glmer(sens ~ cov*dens + cov*driver + dens*driver + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m2)


#### diagnostic plots ################

# Simulate residuals
simulationOutput2 <- simulateResiduals(fittedModel = m2, n = 250)

# overall residual plot
png("overall_residual_plot2.png", width = 8, height = 6, units = "in", res = 300)
plot(simulationOutput2)
dev.off()

### 4.1 Plot #############
# plot sens~mat

pred.df2 <- data.frame(Effect(c("mat","cov","driver","dens"), m2,xlevels=list(mat=seq(min(climate_df$mat), max(climate_df$mat), length.out = 1000), cov = levels(climate_df$cov),driver = levels(climate_df$driver),dens=levels(climate_df$dens))))


levels(pred.df2$driver) <- c("Rain","Temperature",NA)
levels(pred.df2$dens) <- c("No density effects","Density effects",NA)

levels(climate_df$driver) <- c("Rain","Temperature",NA)
levels(climate_df$dens) <- c("No density effects","Density effects",NA)

# to add data points that show mean sensitivities of each species
mean.climate_df <- climate_df %>%
  group_by(species,mat,cov,dens) %>%
  summarise(sens = mean(sens))

library(ggrepel)

# plot
plot.raintemp <- ggplot(pred.df2, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=2.5) +
  facet_grid(dens ~ driver) +
  
  geom_jitter(data = mean.climate_df, aes(x = mat, y = sens, color = cov), alpha = 0.8, size = 5) + 
  
  # to add species labels
 #geom_text_repel(data = mean.climate_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3,
  #              box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
   #            max.overlaps = 1000) +

  scale_fill_manual(name = "Covariation:",values = c("#ecbeb3","#029356"), labels = c("No","Yes")) + 
  scale_color_manual(name = "Covariation:",values = c("#ecbeb3","#029356"), labels = c("No","Yes")) +
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
    data = unique(pred.df2[c("dens", "driver")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.raintemp

 ggsave("SensRainTempSpp.png", plot.raintemp, width = 10, height = 8, dpi = 300)


# 6. SI #########################
# some plots for the supporting information
# showing the difference in sensitivities without and with covariation
 
# for each taxon

df.no.cov=filter(climate_df,cov=="0")
df.cov=filter(climate_df,cov=="1")

new.df=df.no.cov
colnames(new.df)[colnames(new.df)=="sens"] <- "sens.no.cov"

new.df$sens.cov <- df.cov$sens
new.df$difference <- new.df$sens.no.cov - new.df$sens.cov
hist(new.df$difference)

levels(new.df$driver)=c("Rain","Temperature")
levels(new.df$dens)=c("No","Yes")


m.violinmix <-  ggplot(new.df[new.df$group=="Mammals",], aes(x = species, y = difference, color = dens, fill = dens)) +
  scale_fill_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  scale_color_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  
  geom_violin(aes(fill = dens,fill = after_scale(colorspace::lighten(fill, .5))), alpha = .5, size = 1.5, width = .5,bw=.1, scale = "width") +
  
  geom_jitter(aes(color=dens), alpha = .5, size = 2, position = position_dodge2(width = .45),show.legend = T) +
  
  stat_summary(geom = "point", fun = median, shape = 23, size = 4, color = "black", stroke = 1, position = position_dodge2(width = .5), show.legend = F) +
  facet_grid(driver ~ .) +
  labs(
    x = "Species",
    y = expression(italic("Sno cov - Scov"))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 15),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "white", size = 0),
    strip.placement = "outside",
    legend.position = "bottom",
    axis.ticks = element_line(color = "black")) +
  guides(fill = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_rect(
    data = unique(new.df[c("driver")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )


m.violinmix

ggsave("SensMammals.png", m.violinmix, width = 10, height = 8, dpi = 300)


p.violinmix <-  ggplot(new.df[new.df$group=="Plants",], aes(x = species, y = difference, color = dens, fill = dens)) +
  scale_fill_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  scale_color_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  
  geom_violin(aes(fill = dens,fill = after_scale(colorspace::lighten(fill, .5))), alpha = .5, size = 1.5, width = .5,bw=0.1, scale = "width") +
  
  geom_jitter(aes(color=dens), alpha = .5, size = 2, position = position_dodge2(width = .45),show.legend = T) +
  
  stat_summary(geom = "point", fun = median, shape = 23, size = 4, color = "black", stroke = 1, position = position_dodge2(width = .5), show.legend = F) +
  facet_grid(driver ~ .) +
  labs(
    x = "Species",
    y = expression(italic("Sno cov - Scov"))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "white", size = 0),
    strip.placement = "outside",
    legend.position = "bottom",
    axis.ticks = element_line(color = "black")) +
  guides(fill = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_rect(
    data = unique(new.df[c("driver")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

p.violinmix

ggsave("SensPlants.png", p.violinmix, width = 10, height = 8, dpi = 300)


b.violinmix <-  ggplot(new.df[new.df$group=="Birds",], aes(x = species, y = difference, color = dens, fill = dens)) +
  scale_fill_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  scale_color_viridis_d(name = "Density:", option = "D", begin = .2, end = .6) +
  
  geom_violin(aes(fill = dens,fill = after_scale(colorspace::lighten(fill, .5))), alpha = .5, size = 1.5, width = .5,bw=0.1, scale = "width") +
  
  geom_jitter(aes(color=dens), alpha = .5, size = 2, position = position_dodge2(width = .45),show.legend = T) +
  
  stat_summary(geom = "point", fun = median, shape = 23, size = 4, color = "black", stroke = 1, position = position_dodge2(width = .5), show.legend = F) +
  facet_grid(driver ~ .) +
  labs(
    x = "Species",
    y = expression(italic("Sno cov - Scov"))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 15),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "white", size = 0),
    strip.placement = "outside",
    legend.position = "bottom",
    axis.ticks = element_line(color = "black")) +
  guides(fill = "none") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_rect(
    data = unique(new.df[c("driver")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

b.violinmix

ggsave("SensBirds.png", b.violinmix, width = 10, height = 8, dpi = 300)


