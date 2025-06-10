###########################################################

# This script contains the main analysis but 
# with populations as random effect

#############################################################
rm(list=ls())

# set wd
 setwd("/Users/maria/Dropbox/teaching/esin/")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)
library(MuMIn)


######################## POPULATION 
 
df=read.csv("AllSens.csv")
length(unique(df$species))

df$population=1

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

climate_df$species_pop=factor(paste(climate_df$species,climate_df$population,sep="_"))

climate_df$species=factor(climate_df$species)
climate_df$population=factor(climate_df$population)

unique(climate_df$species_pop)


# GLMM 1

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))

m1 <- glmer(sens ~ cov *dens +mat + n.vr + par.per.vr + (1+cov|group/species/population) , family = Gamma(link="log"), data = climate_df)

summary(m1)

0.1226584/(0.1226584+0.2228329+1.3131606 +0.4168371+0.8591424)
0.2228329/(0.1226584+0.2228329+1.3131606 +0.4168371+0.8591424)
1.3131606 /(0.1226584+0.2228329+1.3131606 +0.4168371+0.8591424)
0.4168371/(0.1226584+0.2228329+1.3131606 +0.4168371+0.8591424)
0.8591424/(0.1226584+0.2228329+1.3131606 +0.4168371+0.8591424)

r.squaredGLMM(m1)


# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(climate_df$mat), max(climate_df$mat), length.out = 1000), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))

levels(pred.df$dens) <- c("No density feedbacks","Density feedbacks",NA)

# add data points that show mean sensitivities of each species
# mean.climate_df <- climate_df %>%
#   group_by(species,mat,cov,dens,group) %>% summarise(sens = mean(sens))

mean.climate_df = aggregate(sens~species+mat+cov+dens+group, mean, data=climate_df)

levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density feedbacks","Density feedbacks",NA)
levels(mean.climate_df$dens) <- c("No density feedbacks","Density feedbacks",NA)

library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(climate_df$cov) <- c("Yes","No",NA)
levels(mean.climate_df$cov) <- c("Yes","No",NA)

# plot V"
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
    strip.text = element_text(size = 20),
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

ggsave(plot.clim,filename="main_plot_population.pdf",width=8, height=7)

