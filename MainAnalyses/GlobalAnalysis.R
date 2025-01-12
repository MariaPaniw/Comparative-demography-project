###
# This is a script containing the main analyses:
# 1. Load all data
# 2. Prepare data for analyse
# 3. GLMM 1: Sensitivities to all climatic variables (the global analysis)
# 4. GLMM 2: Mean sensitivities per species
# 5. GLMM 3: Sensitivities for plants only
# 6. GLMM 4: Sensitivities to temperature and rain
# 7. GLMM 5: Sensitivities removing simulated lambdas (form IBMs and integrated population models)
# 8. GLMM 6: Log response ratios
# 9. SI plots

# Author of this script: Esin Ickin
# Date: 30.07.2024
###


# 0. Prepare session ############################################

rm(list=ls())

# set wd
setwd("~/Downloads")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)
library(MuMIn)

# 1. Load data #########################################################

# We can import each species' sensitivities separately here or alternatively import the already compiled dataframe with all the species

# Alternatively import all sensitivitesdata
df=read.csv("AllSens.csv")

# Species's sensitivites:

# make sure the directory is correct for each species

# Mammals
AFox=read.csv("Files_ArcticFox/Sens_ArcticFox_MCMC.csv") # [Nater et al. 2021]
StripedMouse=read.csv("Files_StripedMouse/Sens_StripedMouse_MCMC.csv") # [Nater et al. 2018]
Giraffe=read.csv("Files_Giraffes/Sens_Giraffes.csv") # [Bond et al. 2023]
MouseLemur=read.csv("Files_MouseLemur/Sens_MouseLemurs.csv") # [Ozgul et al. 2023]
Reindeer=read.csv("Files_Reindeer/SensReindeer.csv") # [Hansen & Gamelon et al. 2019]
Marmot=read.csv("Files_Marmot/final_results/SENS_MARMOTS_MCMC.csv") # [Paniw et al. 2020]
Meerkat=read.csv("Files_Meerkat/Sens_Meerkats.csv") # [Paniw et al. 2019]
Rabbits=read.csv("Files_Rabbits/sens_raabit.csv")# [Tablado et al. 2012]


# Birds
Petrel=read.csv("Files_SkuaPetrel/Output_Sens_MCMC/SensPetrel_MCMC.csv") # [Quéroué et al. 2019]
MPenguin=read.csv("Files_MagellanicPenguin/Sens_MPenguins.csv") # [Clark-Wolf et al. 2022]
Certhia=read.csv("Files_SwissBirds/output/Certhia_Sens.csv") # [Malchow et al. 2023]
Linaria=read.csv("Files_SwissBirds/output/Linaria_Sens.csv") # [Malchow et al. 2023]
Lophophanes=read.csv("Files_SwissBirds/output/Lophophanes_Sens.csv") # [Malchow et al. 2023]
PrunellaCollaris=read.csv("Files_SwissBirds/output/PrunellaCollaris_Sens.csv") # [Malchow et al. 2023]
PrunellaModularis=read.csv("Files_SwissBirds/output/PrunellaModularis_Sens.csv") # [Malchow et al. 2023]
Pyrrhula=read.csv("Files_SwissBirds/output/Pyrrhula_Sens.csv") # [Malchow et al. 2023]
Sitta=read.csv("Files_SwissBirds/output/Sitta_Sens.csv") # [Malchow et al. 2023]
Turdus=read.csv("Files_SwissBirds/output/Turdus_Sens.csv") # [Malchow et al. 2023]
Dipper=read.csv("Files_Dipper/Sens_Dipper_Resampling.csv") # [Gamelon et al. 2017]
Albatross=read.csv("Files_Albatross/Sens_Albatross.csv") # [Jenouvrier et al. 2018]
Goose=read.csv("Files_BarnacleGoose/sens_goose.csv") # [Layton-Matthews et al. 2020]
Jay=read.csv("Files_SiberianJay/sens_jays.csv") # [Layton-Matthews et al. 2018]
EmperorPanguin=read.csv("Files_Emperor_penguin/sens_emperor_peng.csv") # [Jenouvrier et al. 2012)


# Plants
Cistus=read.csv("Files_Shrubs/Sens_Cistus.csv") # [Paniw et al. 2023]
Halimium=read.csv("Files_Shrubs/Sens_Halimium.csv") # [Paniw et al. 2023]
Protea=read.csv("Files_Protea/Sens_Protea.csv") # [Merow et al. 2014]
Opuntia=read.csv("Files_Opuntia/Sens_OpuntiaImbricata.csv") # [Evers et al. 2020]
Draco=read.csv("Dracocephalum/Sensitivities_Dracocephalum.csv") # [Evers et al. in review]
Trees=read.csv("Files_SpainTrees/SensSpainTrees.csv") # 11 tree species [David Garcia-Callejas et al. 2016]
DewyP=read.csv("DewyPines/Sens_DewyPines.csv") # [Conquet et al. 2024]
Lavandula=read.csv("Files_Lavandula/sens_lavandula.csv") # [Sanchez-Mejia et al. in prep]



# 2. Compile one big data file ###########################################

df=rbind(AFox,StripedMouse,Giraffe,MouseLemur,Reindeer,Marmot,Meerkat,Rabbits,
           Petrel,MPenguin,Certhia,Linaria,Lophophanes,PrunellaCollaris,PrunellaModularis,Pyrrhula,Sitta,Turdus,Dipper,
           Cistus,Halimium,Protea,Opuntia,DewyP,Albatross,Goose,Jay,Trees,Draco,EmperorPanguin,Lavandula)
# 
length(levels(factor(df$species))) # 41

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
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any

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

r.squaredGLMM(m1)

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

levels(climate_df$dens) <- c("No density dependence","Density dependence")
levels(mean.climate_df$dens) <- c("No density dependence","Density dependence")

levels(pred.df$cov) <- c("Yes","No",NA)
levels(climate_df$cov) <- c("Yes","No")
levels(mean.climate_df$cov) <- c("Yes","No")

library(ggrepel) # to add labels

# plot
plot.clim <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.3) +
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  facet_grid(dens ~ .) +
  
  geom_jitter(data = mean.climate_df, aes(x = mat, y = sens, color = cov), alpha = 0.7, size = 5, width = 0, height = 0) +
  
  # can add labels of species if needed: 
  
  # geom_text_repel(data = mean.climate_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
  #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #             max.overlaps = 100) +
  scale_fill_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) + 
  scale_color_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)"
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  geom_rect(
    data = unique(pred.df[c("dens")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim

ggsave("SensAllClimSpp.png", plot.clim, width = 10, height = 8, dpi = 300)

# 4. Mean sensitivities per species ################################# 

#### 4.1 Removing outliers
sub=droplevels(climate_df[climate_df$mat>=(-1)&climate_df$mat<2,])

mean.df = aggregate(sens~species+mat+cov+dens+group, mean, data=sub)

m2=glmer(sens ~ dens + mat + (1|group) , family = Gamma(link="log"), data = mean.df[mean.df$cov=="Yes",])

summary(m2)

r.squaredGLMM(m2)

a= ggplot(mean.df[mean.df$cov=="Yes",], aes(x=dens, y=sens)) +
  
  geom_boxplot(width=0.3, fill="grey90") + 
  geom_jitter(data=mean.df[mean.df$cov=="Yes",],aes(col=group),size=4,width=0.2,alpha=0.6) +
  # stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", show.legend = F, 
  #              position = position_dodge(.1),col="red") +
  scale_colour_manual(values=c("blue4","darkorchid","darkorange"))+
  ylab("Scaled population growth sensitivities (|S|)")+
  xlab("")+
  theme_bw(base_size=18)+theme(panel.grid = element_blank())+
  theme(panel.grid.major = element_blank(),
        legend.title = element_blank(),
        # legend.text = element_text(size=20),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour="black"))

### MEAN SENS PER SPECIES

#### 4.2 All species 
mean.df = aggregate(sens~species+mat+cov+dens+group, mean, data=climate_df)

m2.2=glmer(sens ~ dens + mat + (1|group) , family = Gamma(link="log"), data = mean.df[mean.df$cov=="Yes",])

summary(m2.2)

r.squaredGLMM(m2.2)

b=ggplot(mean.df[mean.df$cov=="Yes",], aes(x=dens, y=sens)) +
  
  geom_boxplot(width=0.3, fill="grey90") + 
  geom_jitter(data=mean.df[mean.df$cov=="Yes",],aes(col=group),size=4,width=0.2,alpha=0.6) +
  # stat_summary(fun.data = "mean_cl_boot", geom = "pointrange", show.legend = F, 
  #              position = position_dodge(.1),col="red") +
  scale_colour_manual(values=c("blue4","darkorchid","darkorange"))+
  ylab("Scaled population growth sensitivities (|S|)")+
  xlab("")+
  theme_bw(base_size=18)+theme(panel.grid = element_blank())+
  theme(panel.grid.major = element_blank(),
        legend.title = element_blank(),
        # legend.text = element_text(size=20),
        legend.position = c(0.8,0.87),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour="black"))

library(patchwork)

x= b + a +plot_annotation(tag_levels = 'A')

ggsave(x,filename="mean_per_species.pdf",width=11, height=6)

# 5. GLMM: Plants only ##################################################################

sub=climate_df[climate_df$group%in%"Plants",]

length(unique(sub$species))

m1 <- glmer(sens ~ cov*dens + mat + n.vr + par.per.vr + (1|species) , family = Gamma(link="log"), data =climate_df[climate_df$group%in%"Plants",])

summary(m1)

r.squaredGLMM(m1)

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(sub$mat), max(sub$mat), length.out = 1000), cov = levels(sub$cov),dens=levels(sub$dens))))

levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)

mean.sub = aggregate(sens~species+mat+cov+dens+group, mean, data=sub)

levels(sub$dens)
levels(mean.sub$dens)

levels(sub$dens) <- c("No density dependence","Density dependence")
levels(mean.sub$dens) <- c("No density dependence","Density dependence")

library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(sub$cov) <- c("Yes","No")
levels(mean.sub$cov) <- c("Yes","No")

# plot 
plot.clim.plants <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.3) +
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  facet_grid(dens ~ .) +
  
  geom_jitter(data = mean.sub, aes(x = mat, y = sens, color = cov), alpha = 0.7, size = 5, width = 0, height = 0) +
  
  # can add labels of species if needed: 
  
  # geom_text_repel(data = mean.climate_df, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
  #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #             max.overlaps = 100) +
  scale_fill_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) + 
  scale_color_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)"
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  geom_rect(
    data = unique(pred.df[c("dens")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim.plants

ggsave(plot.clim.plants,filename="main_plants.pdf",width=8, height=7)

# 6. GLMM: Sens to Temperature vs Rain ##################################################################

# split climate further into temperature and rain

levels(factor(climate_df$driver))

# group the drivers into rain and temperature
climate_df$driver[climate_df$driver %in% c("Pbr", "Pwn", "rain","fallR","Precipitation","precipitation","Rain","Rainfall","ROS","Snow")]="rain"

climate_df$driver[climate_df$driver %in% c("lagged temperature","prevwinterT","sea ice","SST","SSTA_b","SSTA_b_l","SSTA_m","SSTA_m_l","summerT","Tat","Tbr","temperature","Temperature","Twn","SST1","SST2","SST3","SSTG")]="temperature"

climate_df_sub=filter(climate_df, !driver %in% c("SAM","Q","PET","Winterlength"))
# remove SAM (a driver that influences climate like rain and temperature, so not really a climatic driver)
# remove Q because now I can't sort it into rain or temperature
# remove winterlength again for the same reasons as above

# remove sens == 0 if there is 
min(climate_df_sub$sens)
climate_df_sub=filter(climate_df_sub,sens!=0) # need to remove zeros because of gamma distribution

# from chr to factors
climate_df_sub$driver=factor(climate_df_sub$driver)
climate_df_sub$dens=factor(climate_df_sub$dens)
climate_df_sub$biotic_interactions=factor(climate_df_sub$biotic_interactions)
climate_df_sub$cov=factor(climate_df_sub$cov) #is there covariation 0/1
climate_df_sub$study.doi=factor(climate_df_sub$study.doi)
climate_df_sub$group=factor(climate_df_sub$group) # plants, birds, mammals
climate_df_sub$species=factor(climate_df_sub$species) 

length(levels(climate_df_sub$species)) # the only species removed is marmot because we could not sort Q into rain or temperature category

length(levels(climate_df_sub$species))
length(levels(factor((climate_df_sub[climate_df_sub$dens=="No density dependence",]$species))))
length(levels(factor((climate_df_sub[climate_df_sub$dens=="Density dependence",]$species))))

levels(climate_df_sub$dens)

# GLMM:
m2 <- glmer(sens ~ cov*dens + cov*driver + dens*driver + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df_sub)

summary(m2)

r.squaredGLMM(m2)

### 4.1 Plot #############
# plot sens~mat

pred.df2 <- data.frame(Effect(c("mat","cov","driver","dens"), m2,xlevels=list(mat=seq(min(climate_df_sub$mat), max(climate_df_sub$mat), length.out = 1000), cov = levels(climate_df_sub$cov),driver = levels(climate_df_sub$driver),dens=levels(climate_df_sub$dens))))

levels(pred.df2$driver) <- c("Rain","Temperature",NA)
levels(pred.df2$dens) <- c("No density dependence","Density dependence",NA)

levels(climate_df_sub$driver) <- c("Rain","Temperature",NA)
levels(climate_df_sub$dens) <- c("No density dependence","Density dependence",NA)

# to add data points that show mean sensitivities of each species
mean.climate_df_sub = aggregate(sens~species+mat+cov+dens+group+driver, mean, data=climate_df_sub)

library(ggrepel)

# plot
plot.raintemp <- ggplot(pred.df2, aes(x = mat, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.2) +
  geom_line(aes(col = cov),linewidth=2.5) +
  facet_grid(dens ~ driver) +
  
  geom_jitter(data = mean.climate_df_sub, aes(x = mat, y = sens, color = cov), alpha = 0.8, size = 5) + 
  
  # # to add species labels
  # geom_text_repel(data = mean.climate_df_sub, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3,
  #              box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #            max.overlaps = 1000) +
  # 
  scale_fill_manual(name = "Environmental covariation:",values = c("#ecbeb3","#029356"), labels = c("Yes","No")) + 
  scale_color_manual(name = "Environmental covariation:",values = c("#ecbeb3","#029356"), labels = c("Yes","No")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)") + 
  theme_minimal() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18),
        strip.text = element_text(size = 19),legend.text = element_text(size = 18),
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

ggsave(plot.raintemp,filename="main_plot_rain_temp.pdf",width=10, height=9)

# 7. GLMM Without studies that simulated lambda ################################################

sub=climate_df[climate_df$lambda.sim=="0",] # results also hold if you use simulated only 

m3 <- glmer(sens ~ cov*dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df[climate_df$lambda.sim=="0",])

summary(m3)

r.squaredGLMM(m3)

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m3,xlevels=list(mat=seq(min(sub$mat), max(sub$mat), length.out = 1000), cov = levels(sub$cov),dens=levels(sub$dens))))

levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)

# add data points that show mean sensitivities of each species
# mean.sub <- sub %>%
#   group_by(species,mat,cov,dens,group) %>% summarise(sens = mean(sens))

mean.sub = aggregate(sens~species+mat+cov+dens+group, mean, data=sub)

levels(sub$dens)
levels(mean.sub$dens)

levels(sub$dens) <- c("No density dependence","Density dependence")
levels(mean.sub$dens) <- c("No density dependence","Density dependence")

library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(sub$cov) <- c("Yes","No")
levels(mean.sub$cov) <- c("Yes","No")

# plot
plot.clim.sub <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.3) +
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  facet_grid(dens ~ .) +
  
  geom_jitter(data = mean.sub, aes(x = mat, y = sens, color = cov), alpha = 0.7, size = 5, width = 0, height = 0) +
  
  # can add labels of species if needed: 
  
  # geom_text_repel(data = mean.sub, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
  #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #             max.overlaps = 100) +
  scale_fill_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) + 
  scale_color_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Scaled population growth sensitivities (|S|)"
  ) + 
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  geom_rect(
    data = unique(pred.df[c("dens")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim.sub

ggsave(plot.clim.sub,filename="main_plot_no_simulated_Lambda.pdf",width=8, height=7)

# 8. GLMM Log Response ratios ################################################

climate_df=filter(climate_df,l_ratio!=0) # need to remove zeros because of gamma distribution if there are any

# Change order of facto levels

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))

m4 <- glmer(l_ratio ~ cov*dens + mat+ n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m4)

r.squaredGLMM(m4)

##### Note that I remove the smallest age at sexual maturity for plotting

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m4,xlevels=list(mat=seq(-0.4, max(climate_df$mat), length.out = 100), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))

levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)

# add data points that show mean sensitivities of each species

mean.climate_df = aggregate(l_ratio~species+mat+cov+dens+group, mean, data=climate_df[climate_df$mat>(-1),])

levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density dependence","Density dependence")
levels(mean.climate_df$dens) <- c("No density dependence","Density dependence")

library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(climate_df$cov) <- c("Yes","No")
levels(mean.climate_df$cov) <- c("Yes","No")

# plot
plot.clim.lr <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
  geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.3) +
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  facet_grid(dens ~ .) +
  
  geom_jitter(data = mean.climate_df, aes(x = mat, y = l_ratio, color = cov), alpha = 0.7, size = 5, width = 0, height = 0) +
  
  # can add labels of species if needed: 
  
  # geom_text_repel(data = mean.climate_df, aes(x = mat, y = l_ratio, label = species,fontface = "italic"), size = 3.5,
  #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
  #             max.overlaps = 100) +
  scale_fill_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) + 
  scale_color_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
  labs(
    x = "Log age at sexual maturity (in years)",
    y = "Log response ratio (|L|)"
  ) + 
  # ylim(0,2)+
  theme_minimal() +
  theme(
    axis.title = element_text(size = 20), 
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 19),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom"
  ) +
  geom_rect(
    data = unique(pred.df[c("dens")]),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    color = "black",
    fill = NA,
    inherit.aes = FALSE
  )

plot.clim.lr

ggsave(plot.clim.lr,filename="main_plot_L_ratio.pdf",width=8, height=7)


# 9. SI #########################
# some plots for the supporting information
# showing the difference in sensitivities without and with covariation
 
# filter df so that we only include sensitivities to climatic variables
climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","Winterlength")) # remove SAM and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see species Marmota flaviventris by Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any

climate_df$driver=factor(climate_df$driver)
climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))


# for each taxon

levels(climate_df$driver)
# group the drivers into rain and temperature
climate_df$driver[climate_df$driver %in% c("Pbr", "Pwn", "rain","fallR","Precipitation","precipitation","Rain","Rainfall","ROS","Snow")]="rain"

climate_df$driver[climate_df$driver %in% c("lagged temperature","prevwinterT","sea ice","SST","SSTA_b","SSTA_b_l","SSTA_m","SSTA_m_l","summerT","Tat","Tbr","temperature","Temperature","Twn","SST1","SST2","SST3","SSTG")]="temperature"

levels(factor(climate_df$driver))
climate_df=filter(climate_df,driver!="Q")
climate_df$driver=factor(climate_df$driver)

df.no.cov=filter(climate_df,cov=="0")
df.cov=filter(climate_df,cov=="1")

new.df=df.no.cov

colnames(new.df)[colnames(new.df)=="sens"] <- "sens.no.cov"

new.df$sens.cov <- df.cov$sens
new.df$difference <- new.df$sens.no.cov - new.df$sens.cov
hist(new.df$difference)

View(new.df)

levels(new.df$driver)

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


