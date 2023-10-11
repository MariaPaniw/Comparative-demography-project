################################################################################
#### Covariate Extraction
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)
library(raster)
library(coda)
library(bayesplot)
library(BayesianTools)
library(Hmisc)
library(xtable)
library(pROC)


# Set working directory
setwd("~/SwissBirds")

# Read in species counts and covariates
load("MHB_aggregated_counts.Rdata")
load("aggrVariables_zTrafod.Rdata")

# Rename
dat <- MHB_counts_sp
env <- yearly_envStacks_standard
rm(MHB_counts_sp, yearly_envStacks_standard)

# The sightings are by blocks. The cell and block IDs match between the raster objects and the species counts.

################################################################################
#### Data Cleaning
################################################################################
# Let's reorganize and put everything into a single
# tibble/dataframe
dat_clean <- lapply(1:length(dat), function(x) {

  # Species we are currently looking at
  species <- names(dat)[x]

  # Identify the sampling locations (cells within blocks)
  locations <- dat[[x]]$sample_cells_abd %>%
    setNames(c("Cell", "Block")) %>%
    nest(Cells = -Block) %>%
    mutate(Cells = map(Cells, function(x) {
        pull(x, Cell)
    }))

  # Identify counts within each block
  counts <- dat[[x]]$MHB_counts

  # Now combine the two. Each row consists of the counts across multiple
  # years within each block. We also keep track of the nested "Cells", so
  # we later know from which cells to extract covariate values.
  all <- tibble(Species = species, locations, counts)
  return(all)

}) %>% do.call(rbind, .) %>% arrange(Species, Block)

# We can transform this to a tidy dataframe and to binary
# presence/absence data
dat_clean <- dat_clean %>%
  pivot_longer(X1999:X2019, names_to = "Year", values_to = "Count") %>%
  mutate(Year = as.numeric(substr(Year, start = 2, stop = 5))) %>%
  mutate(Presence = if_else(Count == 0 | is.na(Count), F, T)) %>%
  arrange(Species, Year, Block)

# Let's also clean the environmental covariates
env_clean <- tibble(
    Year = as.numeric(substr(names(env), start = 6, stop = 9))
  , Maps = lapply(env, function(x) {x})
)

# Look at them
print(dat_clean)
print(env_clean)

# Now we can join the two
dat <- left_join(dat_clean, env_clean, by = "Year")

# Take a look again
print(dat)

################################################################################
#### Covariate Extraction
################################################################################
# Extract covariates
pb <- txtProgressBar(min = 1, max = nrow(dat), style = 3)
covariates <- lapply(1:nrow(dat), function(i) {
  covars <- raster::extract(dat$Maps[[i]], dat$Cells[[i]])
  covars <- apply(covars, 2, mean, na.rm = T)
  setTxtProgressBar(pb, i)
  return(covars)
}) %>% bind_rows()

# Combine with our data
dat <- tibble(dat, covariates)
print(dat)

# Now you can delete all that you don't need anymore
dat_final <- dat %>%
  dplyr::select(-c(Cells, Maps))
print(dat_final)

# remove NAs
dat_final=na.omit(dat_final)
# remove Presence=FALSE
dat_final=filter(dat_final,Presence == TRUE)

# species 2
pyrrhula.env <- filter(dat_final, Species == "Pyrrhula.pyrrhula")  %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))
write.csv(pyrrhula.env,"PyrrhulaCovariates.csv")

# species 3
lophophanes.env <- filter(dat_final, Species == "Lophophanes.cristatus") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(lophophanes.env,"LophophanesCovariates.csv")

# species 4
certhia.env <- filter(dat_final, Species == "Certhia.familiaris") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(certhia.env,"CerthiaCovariates.csv")

# species 5
sitta.env <- filter(dat_final, Species == "Sitta.europaea") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(sitta.env,"SittaCovariates.csv")

# species 7
prunellamod.env <- filter(dat_final, Species == "Prunella.modularis") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(prunellamod.env,"PrunellaModularisCovariates.csv")


# species 8
linaria.env <- filter(dat_final, Species == "Linaria.cannabina") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(linaria.env,"LinariaCovariates.csv")

# species 9
turdus.env <- filter(dat_final, Species == "Turdus.torquatus") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(turdus.env,"TurdusCovariates.csv")

# species 10
prunellacoll.env <- filter(dat_final, Species == "Prunella.collaris") %>%
  group_by(Year) %>%
  summarise(breeding_tmp = mean(breeding_tmp),
            breeding_pr = mean(breeding_pr),
            postbreed_tmp = mean(postbreed_tmp),
            winter_tmin = mean(winter_tmin),
            winter_pr = mean(winter_pr))

write.csv(prunellacoll.env,"PrunellaCollarisCovariates.csv")


plot(pyrrhula.env$breeding_tmp, type = "l", col = "blue", xlab = "X Axis Label", ylab = "Breeding temperature", main = "Breeding temperature of 8 species")
lines(lophophanes.env$breeding_tmp, col = "green")
lines(linaria.env$breeding_tmp, col = "red")
lines(prunellacoll.env$breeding_tmp, col = "purple")
lines(prunellamod.env$breeding_tmp, col = "orange")
lines(sitta.env$breeding_tmp, col = "yellow")
lines(turdus.env$breeding_tmp, col = "pink")
lines(certhia.env$breeding_tmp, col = "brown")
legend("topleft", legend = c("Pyrrhula", "Lophophanes", "Linaria", "Prunella c.", "Prunella m.", "Sitta", "Turdus", "Certhia"),
       col = c("blue", "green", "red", "purple", "orange", "yellow", "pink", "brown"), lty = 1,cex=0.5)

# MCMC POSTERIOR DISTRIBUTION EXTRACTION ------------
# MCMC files
# species names
species_names <- c("Chaffinch", 
                   "Bullfinch",
                   "Crested tit",
                   "Eurasian woodcreeper",
                   "Eurasian nuthatch",
                   "Goldcrest",
                   "Dunnock",
                   "Common linnet",
                   "Ring ouzel",
                   "Alpine accentor")

# parameter names
par_names_tex <- c("Dens-dep. $b^{-1}$","Fec. $\\rho_0$","Juv. surv. $s_j$","Ad. surv. $s_a$","Emig. prob. $e_1$","Disp. dist. $\\overline d$","Dispersion $\\nu$","\\beta_{fc,p1}","\\beta_{fc,p2}","\\beta_{fc,t1}","\\beta_{fc,t2}","\\beta_{jS,p1}","\\beta_{jS,p2}","\\beta_{jS,t1}","\\beta_{jSt2}","\\beta_{aS,p1}","\\beta_{aS,p2}","\\beta_{aS,t1}")
par_names     <- c("DensDep",           "Fecund",        "juvSurv",         "adSurv",         "EmigProb",         "DispDist",                 "GPsize",           "Fcb.pr1",      "Fcb.pr2",      "Fcb.tm1",      "Fcb.tm2",      "jSb.pr1",      "jSb.pr2",      "jSb.tm1",      "jSb.tm2",     "aSb.pr1",      "aSb.pr2",      "aSb.tm"       )


## get info about the MCMC to identify the correct file names to read:

samplr <- "CA_mcmc" # likelihood type
iter <- 6e+04       # number of iterations in chain
specs <- spec0 <- 2:10       # species IDs
spec0[9] <- 0               # species IDs as used in batch numbers
batches <- matrix(rep((500+spec0),each=3) + (4:6)*10, nrow=3)					#--- Changed batch nr. 
thinningrate <- 100  # rate of thinning

sel_specs <- c(2,3,4,5,7,8,9,10) # -> 6 not converged

calc_quantiles <- c(0.10,0.20,0.25,0.50,0.75,0.80,0.90,0.025,0.975) 					
#--- added additional quantiles here
use_climquants <- paste0(calc_quantiles[c(1,3,4,5,7)]*100,"%")					#--- we use 0.1 and 0.9 
post_samplesize <- 396															#--- sample size to 396
pred_resolution <- 0.01


load("results/CA_out/CA_comb_Spec2_it60000_542-562.Rdata")
sp2 <- sampler_list # for species 2

load("results/CA_out/CA_comb_Spec3_it60000_543-563.Rdata")
sp3 <- sampler_list # for species 3

load("results/CA_out/CA_comb_Spec4_it60000_544-564.Rdata")
sp4 <- sampler_list # etc.

load("results/CA_out/CA_comb_Spec5_it60000_545-565.Rdata")
sp5 <- sampler_list

load("results/CA_out/CA_comb_Spec7_it60000_547-567.Rdata")
sp7 <- sampler_list

load("results/CA_out/CA_comb_Spec8_it60000_548-568.Rdata")
sp8 <- sampler_list

load("results/CA_out/CA_comb_Spec9_it60000_549-569.Rdata")
sp9 <- sampler_list

load("results/CA_out/CA_comb_Spec10_it60000_540-560.Rdata")
sp10 <- sampler_list

sp  <- c(2,3,4,5,7,8,9,10) # species ID
length(sp)
post_sample_list <- list()

for( i in sp){
  file_name <- paste0("sp",i)
  # convert to BayesianTools objects
  BT_smplst <- convertCoda(get(file_name), names = par_names)
  
  # remove loaded object again
  rm(file_name)
  
  # select needed parameters, discard burn-in (= first half of chain) & draw posterior sample
  parSel <- c(1:4,8:18)
  post_sample <- getSample(BT_smplst, parametersOnly = T, coda = F, start = 1+iter/2, whichParameters = parSel, numSamples = post_samplesize) #-- changed from 3 to 2 as we now have 3*60.000 iterations
  post_sample_name <- paste0("post_sample",i) 
  post_sample_list[[post_sample_name]] <- as.data.frame(post_sample) # list of dataframes of the 8 species
  
}

saveRDS(post_sample_list,"Post_sample_list.rds")


