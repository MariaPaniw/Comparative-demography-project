############################################################################
#
# This script processes the results of the dewy-pine projections under
# four scenarios for each covariate: 
# 
# (1) Focal = min value, other = mean
# (2) Focal = max value, other = mean
# (3) Focal = min value, other = observed value when focal is min
# (4) Focal = max value, other = observed value when focal is min
#
# Author: Eva Conquet
#
###########################################################################

###########################################################################
#
# 1. House keeping and loading libraries and data ----
#
###########################################################################

## 1.1. House keeping ----
# -------------------

rm(list = ls())

# set wd
# ...

## 1.2. Loading libraries ----
# -----------------------

library(ggplot2)
library(viridis)


## 1.3. Loading data ----
# ------------------

for(i in 1:length(list.files(path = "Output/", pattern = "IBM_Natural"))){
  
  load(paste0("Output/", list.files(path = "Output/", pattern = "IBM_Natural"))[i])
}




###########################################################################
#
# 2. Putting lambdas into a table ----
#
###########################################################################

nb_proj = nrow(ibm_vertedero_summerT_min_other_cov_mean[[1]]$log_lambda) * length(ibm_vertedero_summerT_min_other_cov_mean)
nb_years = ncol(ibm_vertedero_summerT_min_other_cov_mean[[1]]$log_lambda)

# Creating the empty table
results_df = expand.grid(population = c("SierraCarboneraY5", "SierraRetinY5", "Vertedero"),
                         projection = seq(1, nb_proj),
                         year = seq(1, nb_years), 
                         focal_cov = c("summerT", "prevwinterT", "fallR", "prevfallR", 
                                       "dens", "size", "TSF"),
                         focal_min_max = c("min", "max"),
                         other_mean_obs = c("mean", "obs"))


# Fill in data
unique(paste(results_df$focal_cov[which(results_df$population == "SierraCarboneraY5")], 
             results_df$focal_min_max[which(results_df$population == "SierraCarboneraY5")], 
             results_df$other_mean_obs[which(results_df$population == "SierraCarboneraY5")]))

results_df$lambda = NA

results_df$lambda[which(results_df$population == "SierraCarboneraY5")] = 
  c(c(ibm_sierracarboneray5_summerT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_dens_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_size_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_summerT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevwinterT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_fallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_prevfallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_dens_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_size_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierracarboneray5_TSF_max_other_cov_obs[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "SierraRetinY5")] = 
  c(c(ibm_sierraretiny5_summerT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_summerT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_fallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_dens_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_size_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_TSF_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_summerT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevwinterT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_fallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_prevfallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_dens_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_size_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_sierraretiny5_TSF_max_other_cov_obs[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "Vertedero")] = 
  c(c(ibm_vertedero_summerT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_mean[[1]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_mean[[2]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_mean[[3]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_mean[[4]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_summerT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_prevwinterT_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_fallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_prevfallR_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_dens_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_size_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_TSF_min_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_summerT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_prevwinterT_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_fallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_prevfallR_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_dens_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_size_max_other_cov_obs[[4]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_obs[[1]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_obs[[2]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_obs[[3]]$log_lambda),
    c(ibm_vertedero_TSF_max_other_cov_obs[[4]]$log_lambda))


# Transform log lambdas into lambdas and set infinite values to 0. 
# Infinite values are caused by the lack of aboveground individuals (the projections
# continue because there are still seeds in the seedbank).
results_df$lambda = exp(results_df$lambda)
results_df$lambda[which(is.nan(results_df$lambda))] = 0


# Get mean lambda per projection for each scenario
results_df = aggregate(lambda ~ population + focal_cov + focal_min_max + other_mean_obs + projection, 
                       data = results_df,
                       mean, na.rm = T)

# Format data
results_df$focal_full = NA
results_df$focal_full[which(results_df$focal_cov == "summerT")] = "Next summer mean max.\ndaily temperature"
results_df$focal_full[which(results_df$focal_cov == "prevwinterT")] = "Previous winter mean max.\ndaily temperature"
results_df$focal_full[which(results_df$focal_cov == "fallR")] = "Next fall cumulative\nrainfall"
results_df$focal_full[which(results_df$focal_cov == "prevfallR")] = "Previous fall cumulative\nrainfall"
results_df$focal_full[which(results_df$focal_cov == "dens")] = "Aboveground density of\nlarge individuals"
results_df$focal_full[which(results_df$focal_cov == "size")] = "Size"
results_df$focal_full[which(results_df$focal_cov == "TSF")] = "Time Since Fire"

results_df$scenario = paste0("Focal covariate = ", results_df$focal_full, " (", results_df$focal_min_max, " value) ",
                            "\nOther covariates = ", results_df$other_mean_obs, " value")
results_df$focal_scenario = paste0(results_df$focal_full, " = ", results_df$focal_min_max, " value")
results_df$other_scenario = paste0(results_df$other_mean_obs, " value")




###########################################################################
#
# 3. Plotting results ----
#
###########################################################################

png(filename = "DewyPine_Lambdas.png", 
     width = 25,
     height = 20,
     units = "cm",
     bg = "white",
     res = 600)

ggplot(results_df, aes(x = focal_scenario, y = lambda, fill = other_scenario)) +
  geom_boxplot(alpha = 0.6) +
  stat_summary(fun = mean, position = position_dodge(0.75)) +
  scale_fill_manual(name = "Other covariates",
                      values = viridis(2, begin = 0.4, end = 0.8, option = "B"),
                      labels = c("Mean value", "Observed value at\nfocal covariate value")) +
  xlab("Focal covariate value") +
  ylab("Population growth rate (\u03BB)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))


###########################################################################
#
# 3. Save results ----
#
###########################################################################
write.csv(results_df, "results_df.csv", row.names = F)



dev.off()