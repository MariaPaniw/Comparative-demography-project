###########################################################

# This script contains the main analysis

#############################################################
rm(list=ls())

# set wd
 setwd("/Users/esinickin/Desktop/PNAS NEXUS")

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(viridis)
library(effects)
library(MuMIn)

df=read.csv("/Users/esinickin/Documents/Master Thesis/MainAnalyses/SciAdv/AllSens.csv") # 41

head(df)
unique(df$species)

### EDA & edit data ##################


# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform generation time
df$mat=log(df$mat)
# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# GLMM: Sens to all Climate Variables ##################################################################
# climate and not splitting by rain/temp
# cov*dens + mat + par + random effects

# filter df so that we only include climatic covariates (what do you think?)
climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","PET","Winterlength")) # remove SAM, PET, and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$sens)
climate_df=filter(climate_df,sens!=0) # need to remove zeros because of gamma distribution if there are any


climate_df$driver=factor(climate_df$driver)

climate_df$densV2=climate_df$dens
climate_df$densV2[climate_df$biotic_interactions==1]=1

climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

climate_df$lambda.sim=factor(climate_df$lambda.sim)


# aggregate sensitivities modeled per lambda type and whether density was included in models to check distribution

sum=aggregate(sens~lambda.sim+dens,function(x) length(unique(x)), data=climate_df)

ggplot(sum, aes(fill=lambda.sim, y=sens, x=dens)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Density effects")+
  ylab("# of |S| values")


# Change order of facto levels

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))

#cov.num=read.csv("cov.considered.csv")

m1 <- glmer(sens ~ cov*dens + mat+ n.vr +par.per.vr+ (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)


summary(m1) 

# new:
1.738/(1.738+0.241+0.767)
0.241/(1.738+0.241+0.767)
0.767/(1.738+0.241+0.767)

r.squaredGLMM(m1)

### CHECK EFFECT OF STUDY LENGTH

m1.1 <- glmer(sens ~ cov*dens + mat+ n.vr +log(study.length)+ par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)


summary(m1.1)


### CHECK IF ADDING BIOTIC INTERACTIONS DOES SOMETHING

climate_df$densV2=factor(climate_df$densV2)

m1.2 <- glmer(sens ~ cov*densV2 + mat+ n.vr +log(study.length)+ par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)

summary(m1.2)

##### NOTE THAT I REMOVE THE SMALLEST AGE AT SEXUAL MATURITY (THE RABBIT) FROM THE PLOT BECAUSE OTHERWISE THE PLOT LOOKS AWFUL

# Main Plot #######################
# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(-0.4, max(climate_df$mat), length.out = 100), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))


levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)


# add data points that show mean sensitivities of each species
# and error bars for variation of the mean


mean.climate_df <- climate_df %>%
  filter(mat > -1) %>%
  group_by(species, mat, cov, dens, group) %>%
  summarise(
    sens_mean = mean(sens),
    sens_sd = sd(sens),
    sens_se = sd(sens) / sqrt(n()), # standard error
    sens_lower = sens_mean - sens_se,
    sens_upper = sens_mean + sens_se,
    sens_margin = (qt(0.975,df=(n()-1)) * sd(sens)) / (sqrt(n())), # calculate 95 CI for visualization purposes only to capture the variation of the mean sensitivities
    lower_interval = mean(sens) - sens_margin,
    upper_interval = mean(sens) + sens_margin,
    .groups = "drop"
  )


levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density dependence","Density dependence")
levels(mean.climate_df$dens) <- c("No density dependence","Density dependence")


library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(climate_df$cov) <- c("Yes","No")
levels(mean.climate_df$cov) <- c("Yes","No")



# species to label
species_labels=c('Dracocephalum austriacum','Lavandula stoechas','Vulpes lagopus','Marmota flaviventer','Protea repens','Pinus pinaster','Fagus sylvatica','Suricata suricatta','Cinclus cinclus','Branta leucopsis','Cistus libanotis','Halimium halimifolium','Linaria cannabina')

# sample size label
levels(pred.df$dens)
annotation_df <- data.frame(
  dens = c("No density dependence", "Density dependence"),  
  mat = c(Inf, Inf),        # right side of x-axis
  sens = c(Inf, Inf),       # top of y-axis
  label = c("n[species] == 28", "n[species] == 13")
)
annotation_df$dens=factor(annotation_df$dens, levels = c("No density dependence", "Density dependence"))

# figure label
panel_labels <- data.frame(
  dens = c("No density dependence", "Density dependence"),  
  mat = -Inf,     # left side
  sens = Inf,     # top side
  label = c("A", "B")
)
panel_labels$dens=factor(panel_labels$dens, levels = c("No density dependence", "Density dependence"))


# MAIN FIGURE PLOT
set.seed(123)  # For reproducibility
mean.climate_df$mat_jit <- jitter(mean.climate_df$mat, amount = 0.1)

# also only want to have points that are labeled
label_points <- filter(mean.climate_df, species %in% species_labels)

plot.clim <- 
  
  ggplot(data = pred.df, aes(x = mat, y = fit)) +
  
  geom_ribbon(
    data = pred.df, 
    aes(ymin = lower, ymax = upper, fill = cov), 
    alpha = 0.3
  ) +
  
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  
  facet_grid(dens ~ .) +
  
  geom_point(
    data = mean.climate_df, 
    aes(x = mat_jit, y = sens_mean, color = cov), 
    alpha = 0.9, 
    size = 4) +
  
  # just for visualization purposes
  geom_errorbar(
    data = mean.climate_df,
    aes(x = mat_jit, ymin = lower_interval, ymax = upper_interval, color = cov),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 1,
    alpha = 0.8
  ) +
  
  scale_fill_manual(name = "Environmental covariation:", values = c( "#197f64","#df9e41"), labels = c("Yes","No")) + 
  
  scale_color_manual(name = "Environmental covariation:", values = c( "#197f64","#df9e41"), labels = c("Yes","No")) +
  
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
  ) +
  
  # comment out if you don't want labels:
  geom_text_repel(
    data = subset(mean.climate_df, species %in% species_labels),
    aes(x = mat_jit, y = sens_mean, label = species, fontface = "italic"),
    size = 4,
    box.padding = unit(0.6, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.color = "black",
    segment.size = 0.6,
    force = 1,
    max.overlaps = 100,
    nudge_x = 0.2,
    nudge_y = 0.3
  ) +
  geom_text(
    data = annotation_df,
    aes(x = mat, y = sens, label = label),
    hjust = 1.1, vjust = 1.2,
    parse = TRUE,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = panel_labels,
    aes(x = mat, y = sens, label = label),
    hjust = -0.5, vjust = 1.2,
    size = 6,
    fontface = "bold",
    inherit.aes = FALSE
  )

plot.clim


# save as tiff file for publication
ggsave(
  filename = "Fig1_06_26/main_plot_empty_Fig1_06_26.tiff",
  plot = plot.clim,
  device = "tiff",
  dpi = 600,
  width = 10,
  height = 8,
  units = "in",
  compression = "lzw"
)


# MAIN PLOT ALL SPECIES
set.seed(123)  # For reproducibility
mean.climate_df$mat_jit <- jitter(mean.climate_df$mat, amount = 0.1)

# also only want to have points that are labeled
label_points <- filter(mean.climate_df, species %in% species_labels)

all.plot.clim <- 
  
  ggplot(data = pred.df, aes(x = mat, y = fit)) +
  
  geom_ribbon(
    data = pred.df, 
    aes(ymin = lower, ymax = upper, fill = cov), 
    alpha = 0.3
  ) +
  
  geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
  
  facet_grid(dens ~ .) +
  
  # for species with label
  geom_point(
    data = mean.climate_df, 
    aes(x = mat_jit, y = sens_mean, color = cov), 
    alpha = 1, 
    size = 4.5) +
  # for species with label
  geom_errorbar(
    data = mean.climate_df,
    aes(x = mat_jit, ymin = sens_lower, ymax = sens_upper, color = cov),
    inherit.aes = FALSE,
    width = 0.15,
    linewidth = 1,
    alpha = 0.9
  ) +
  

  scale_fill_manual(name = "Environmental covariation:", values = c( "#197f64","#df9e41"), labels = c("Yes","No")) + 
  
  scale_color_manual(name = "Environmental covariation:", values = c( "#197f64","#df9e41"), labels = c("Yes","No")) +
  
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
  ) +
  geom_text_repel(
    data = mean.climate_df,
    aes(x = mat_jit, y = sens_mean, label = species, fontface = "italic"),
    size = 4,
    box.padding = unit(0.6, "lines"),
    point.padding = unit(0.4, "lines"),
    segment.color = "black",
    segment.size = 0.6,
    force = 1,
    max.overlaps = 100,
    nudge_x = 0.2,
    nudge_y = 0.3
  ) +
  geom_text(
    data = annotation_df,
    aes(x = mat, y = sens, label = label),
    hjust = 1.1, vjust = 1.2,
    parse = TRUE,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = panel_labels,
    aes(x = mat, y = sens, label = label),
    hjust = -0.5, vjust = 1.2,
    size = 6,
    fontface = "bold",
    inherit.aes = FALSE
  )

all.plot.clim

#ggsave(plot.clim,filename="PNASNEXUS_PLOTS/main_plot_Fig1.png",width=8, height=7)


######################## REMOVE OUTLIERS - I LEFT THIS OUT IN THE END

# sub=droplevels(climate_df[climate_df$mat>=(-1)&climate_df$mat<2,])
# 
# length(unique(sub$species))
# length(unique(sub$species[sub$dens%in%"Density dependence"]))
# 
# m1 <- glmer(sens ~ cov*dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df[climate_df$mat>(-1)&climate_df$mat<2,])
# 
# 
# summary(m1)
# 
# 1.359/c(1.359+5.104e-01+8.750e-01)
# 5.104e-01/c(1.359+5.104e-01+8.750e-01)
# 8.750e-01/c(1.359+5.104e-01+8.750e-01)
# 
# r.squaredGLMM(m1)
# 
# 
# # plot sens~mat
# pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(sub$mat), max(sub$mat), length.out = 1000), cov = levels(sub$cov),dens=levels(sub$dens))))
# 
# levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)
# 
# # add data points that show mean sensitivities of each species
# # mean.sub <- sub %>%
# #   group_by(species,mat,cov,dens,group) %>% summarise(sens = mean(sens))
# 
# mean.sub = aggregate(sens~species+mat+cov+dens+group, mean, data=sub)
# 
# levels(sub$dens)
# levels(mean.sub$dens)
# 
# levels(sub$dens) <- c("No density dependence","Density dependence")
# levels(mean.sub$dens) <- c("No density dependence","Density dependence")
# 
# library(ggrepel) # to add labels
# 
# levels(pred.df$cov) <- c("Yes","No",NA)
# levels(sub$cov) <- c("Yes","No")
# levels(mean.sub$cov) <- c("Yes","No")
# 
# # plot V"
# plot.clim.sub <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
#   geom_ribbon(data = pred.df, aes(ymin = lower, ymax = upper, fill = cov), alpha = 0.3) +
#   geom_line(data = pred.df, aes(col = cov), linewidth = 1.2) +
#   facet_grid(dens ~ .) +
# 
#   geom_jitter(data = mean.sub, aes(x = mat, y = sens, color = cov), alpha = 0.7, size = 5, width = 0, height = 0) +
# 
#   # can add labels of species if needed:
# 
#   # geom_text_repel(data = mean.sub, aes(x = mat, y = sens, label = species,fontface = "italic"), size = 3.5,
#   #               box.padding = unit(0.4, "lines"), point.padding = unit(0.2, "lines"),
#   #             max.overlaps = 100) +
#   scale_fill_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
#   scale_color_manual(name = "Environmental covariation:", values = c( "#b7b3e6","#029356"), labels = c("Yes","No")) +
#   labs(
#     x = "Log age at sexual maturity (in years)",
#     y = "Scaled population growth sensitivities (|S|)"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 20),
#     axis.text = element_text(size = 18),
#     strip.text = element_text(size = 19),
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 20),
#     axis.ticks = element_line(color = "black"),
#     legend.position = "bottom"
#   ) +
#   geom_rect(
#     data = unique(pred.df[c("dens")]),
#     aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
#     color = "black",
#     fill = NA,
#     inherit.aes = FALSE
#   )
# 
# plot.clim.sub
# 
# ggsave(plot.clim.sub,filename="main_no_outlier_Mat.pdf",width=8, height=7)

### MEAN SENS PER SPECIES SUBSET 

sub=droplevels(climate_df[climate_df$mat>=(-1)&climate_df$mat<2,])

mean.df = aggregate(sens~species+mat+cov+dens+group, mean, data=sub)

m2=glmer(sens ~ dens + mat + (1|group), 
         family = Gamma(link="log"),
         data = mean.df[mean.df$cov=="Yes",])

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

#library(patchwork)

#x= b + a + plot_annotation(tag_levels = 'A')
#x

library(cowplot)
x=plot_grid(b,a, labels = c("A","B",label_size = 18))
x
ggsave(x,filename="mean_per_species.pdf",width=11, height=6)


######################## PLANTS ONLY

sub=climate_df[climate_df$group%in%"Plants",]

length(unique(sub$species))

m1 <- glmer(sens ~ cov*dens + mat + n.vr + par.per.vr + (1|species) , family = Gamma(link="log"), data =climate_df[climate_df$group%in%"Plants",])


summary(m1)

0.6289 /c(0.6289 +0.6738)
0.6738/c(0.6289 +0.6738)


r.squaredGLMM(m1)


# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(sub$mat), max(sub$mat), length.out = 1000), cov = levels(sub$cov),dens=levels(sub$dens))))

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

# plot V"
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



######################## LOG RATIOS

df=read.csv("AllSens_new.csv")

head(df)

df=df[!df$species%in%"Rangifer tarandus",]

reindeer=read.csv("/Users/maria/Dropbox/teaching/esin/reindeer/SensReindeer.csv")

df=rbind(df,reindeer)

df=df[!df$species%in%"Branta leucopsis",]

goose=read.csv("/Users/maria/Dropbox/teaching/esin/Geese/sens_goose_full.csv")

df=rbind(df,goose)

rabbits= read.csv("/Users/maria/Dropbox/teaching/esin/Rabbits/sens_rabbit.csv")

df=rbind(df,rabbits)

lav= read.csv("/Users/maria/Dropbox/teaching/esin/Lavandula/sens_lavandula.csv")

df=rbind(df,lav)

emp= read.csv("/Users/maria/Dropbox/teaching/esin/emperor penguin/sens_emperor_peng.csv")

df=rbind(df,emp)

# add new column with mean number of parameters per vital rate
df$par.per.vr=df$n.pam/df$n.vr

# log-transform generation time
df$mat=log(df$mat)
# log-transform also number of vital rates and parameters per vital rate
df$n.vr=log(df$n.vr)
df$par.per.vr=log(df$par.per.vr)

# GLMM: Log ratios ##################################################################
# climate and not splitting by rain/temp
# cov*dens + mat + par + random effects

# filter df so that we only include climatic covariates (what do you think?)
climate_df=filter(df, driver.type == "C")

levels(factor(climate_df$driver))
levels(factor(climate_df$species))

climate_df=filter(climate_df,!driver %in% c("SAM","PET","Winterlength")) # remove SAM, PET, and winterlength since they are not really climatic drivers like temperatures or percipitation, however we keep Q because it's a composite measure (for more see Paniw et al. 2020)

# remove sens == 0 if there is 
min(climate_df$l_ratio)
climate_df=filter(climate_df,l_ratio!=0) # need to remove zeros because of gamma distribution if there are any


climate_df$driver=factor(climate_df$driver)

# climate_df$dens[climate_df$biotic_interactions==1]=1

climate_df$dens=factor(climate_df$dens)
climate_df$biotic_interactions=factor(climate_df$biotic_interactions)
climate_df$cov=factor(climate_df$cov) #is there covariation 0/1
climate_df$study.doi=factor(climate_df$study.doi)
climate_df$group=factor(climate_df$group) # plants, birds, mammals
climate_df$species=factor(climate_df$species)

climate_df$lambda.sim=factor(climate_df$lambda.sim)



# Change order of facto levels

climate_df$cov=factor(climate_df$cov,levels=c("1","0"))


m1 <- glmer(l_ratio ~ cov*dens + mat+ n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df)


# AIC(m0,m1)

summary(m1)

1.398/(1.398+6.116e-01+9.661e-01)
6.116e-01/(1.398+6.116e-01+9.661e-01)
9.661e-01/(1.398+6.116e-01+9.661e-01)

r.squaredGLMM(m1)

##### NOTE THAT I REMOVE THE SMALLEST AGE AT SEXUAL MATURITY (THE RABBIT) FROM THE PLOT BECAUSE OTHERWISE THE PLOT LOOKS AWFUL

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(-0.4, max(climate_df$mat), length.out = 100), cov = levels(climate_df$cov),dens=levels(climate_df$dens))))

levels(pred.df$dens) <- c("No density dependence","Density dependence",NA)


# add data points that show mean sensitivities of each species
# mean.climate_df <- climate_df %>%
#   group_by(species,mat,cov,dens,group) %>% summarise(sens = mean(sens))

mean.climate_df = aggregate(l_ratio~species+mat+cov+dens+group, mean, data=climate_df[climate_df$mat>(-1),])

levels(climate_df$dens)
levels(mean.climate_df$dens)

levels(climate_df$dens) <- c("No density dependence","Density dependence")
levels(mean.climate_df$dens) <- c("No density dependence","Density dependence")

# mean.climate_df[mean.climate_df$l_ratio>3,]

library(ggrepel) # to add labels

levels(pred.df$cov) <- c("Yes","No",NA)
levels(climate_df$cov) <- c("Yes","No")
levels(mean.climate_df$cov) <- c("Yes","No")

# plot
plot.clim <- ggplot(data = pred.df, aes(x = mat, y = fit)) +
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

plot.clim

ggsave(plot.clim,filename="main_plot_L_ratio.pdf",width=8, height=7)

####################### DIFFERENCES RAIN TEMP

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
#climate_df_sub=filter(climate_df_sub,sens!=0) # need to remove zeros because of gamma distribution

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

1.279/(1.279+4.872e-01+8.458e-01)
4.872e-01/(1.279+4.872e-01+8.458e-01)
8.458e-01/(1.279+4.872e-01+8.458e-01)

r.squaredGLMM(m2)

### Plots #############
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


######################## REMOVE SIMULATED

sub=climate_df[climate_df$lambda.sim=="0",] # results also hold if you use simulated only 

m1 <- glmer(sens ~ cov*dens + mat + n.vr + par.per.vr + (1+cov|group/species) , family = Gamma(link="log"), data = climate_df[climate_df$lambda.sim=="0",])

summary(m1)

r.squaredGLMM(m1)

9.018e-01/(9.018e-01+5.170e-01+9.076e-01)
5.170e-01/(9.018e-01+5.170e-01+9.076e-01)
9.076e-01/(9.018e-01+5.170e-01+9.076e-01)

# plot sens~mat
pred.df <- data.frame(Effect(c("mat","cov","dens"), m1,xlevels=list(mat=seq(min(sub$mat), max(sub$mat), length.out = 1000), cov = levels(sub$cov),dens=levels(sub$dens))))

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

# plot V"
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


# SI #########################

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
#climate_df_sub=filter(climate_df_sub,sens!=0) # need to remove zeros because of gamma distribution

# from chr to factors
climate_df_sub$driver=factor(climate_df_sub$driver)
climate_df_sub$dens=factor(climate_df_sub$dens)
climate_df_sub$biotic_interactions=factor(climate_df_sub$biotic_interactions)
climate_df_sub$cov=factor(climate_df_sub$cov) #is there covariation 0/1
climate_df_sub$study.doi=factor(climate_df_sub$study.doi)
climate_df_sub$group=factor(climate_df_sub$group) # plants, birds, mammals
climate_df_sub$species=factor(climate_df_sub$species) 

length(levels(climate_df_sub$species)) # the only species removed is marmot because we could not sort Q into rain or temperature category

levels(climate_df_sub$driver)


df.no.cov=filter(climate_df_sub,cov=="0")
df.cov=filter(climate_df_sub,cov=="1")

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
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(angle = 52, hjust = 1, face = "italic"),
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

ggsave("SensBirds.png", b.violinmix, width = 12, height = 8, dpi = 300)


