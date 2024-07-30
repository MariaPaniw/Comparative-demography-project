############################################################################
#
# This script processes the results of the dewy-pine projections under
# two scenarios for each focal covariate related to each focal vital rate 
# (reproductive and non-reproductive survival, and reproduction probability): 
# 
# (1) Focal covariate = min value, 
#     other covariates of focal vital rate = observed when focal covariate = min,
#     covariates of other vital rates = mean
# (2) Focal covariate = max value, 
#     other covariates of focal vital rate = observed when focal covariate = max,
#     covariates of other vital rates = mean
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


## 1.2. Loading libraries ----
# -----------------------

library(ggplot2)
library(viridis)


## 1.3. Loading data ----
# ------------------

for(i in 1:length(list.files(path = "Output/", pattern = "IBM_[^N]")[-grep("other", list.files(path = "Output/", pattern = "IBM_[^N]"))][-grep("Progress", list.files(path = "Output/", pattern = "IBM_[^N]")[-grep("other", list.files(path = "Output/", pattern = "IBM_[^N]"))])])){
  
  load(paste0("Output/", list.files(path = "Output/", pattern = "IBM_[^N]")[-grep("other", list.files(path = "Output/", pattern = "IBM_[^N]"))][-grep("Progress", list.files(path = "Output/", pattern = "IBM_[^N]")[-grep("other", list.files(path = "Output/", pattern = "IBM_[^N]"))])])[i])
}




###########################################################################
#
# 2. Putting lambdas into a table ----
#
###########################################################################

nb_proj = nrow(ibm_vertedero_flowering_dens_max[[1]]$log_lambda) * length(ibm_vertedero_flowering_dens_max)
nb_years = ncol(ibm_vertedero_flowering_dens_max[[1]]$log_lambda)

# Creating the empty table
results_df = expand.grid(population = c("Bujeo", "Prisoneros", "Retin", "MonteraTorero", "SCarbDist", 
                                        "SierraCarboneraY5", "SierraRetinY5", "Vertedero"),
                         projection = seq(1, nb_proj),
                         year = seq(1, nb_years), 
                         focal_vr = c("non_repro_survival", "repro_survival", "flowering"),
                         focal_cov = c("summerT", "prevwinterT", "fallR", "prevfallR", 
                                       "prevwinterR", "dens"),
                         focal_min_max = c("min", "max"))

results_df = results_df[-which(results_df$focal_vr %in% c("non_repro_survival",
                                                          "repro_survival") &
                               results_df$focal_cov %in% c("prevwinterT",
                                                           "prevfallR",
                                                           "prevwinterR")), ]
results_df = results_df[-which(results_df$focal_vr == "flowering" &
                                 results_df$focal_cov %in% c("summerT",
                                                             "fallR")), ]


# Fill in data
unique(paste(results_df$focal_vr[which(results_df$population == "Bujeo")],
             results_df$focal_cov[which(results_df$population == "Bujeo")], 
             results_df$focal_min_max[which(results_df$population == "Bujeo")]))

results_df$lambda = NA

results_df$lambda[which(results_df$population == "Bujeo")] = 
  c(c(ibm_bujeo_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_bujeo_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_bujeo_flowering_prevwinterR_min[[1]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_min[[2]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_min[[3]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_min[[4]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_bujeo_flowering_dens_min[[1]]$log_lambda),
    c(ibm_bujeo_flowering_dens_min[[2]]$log_lambda),
    c(ibm_bujeo_flowering_dens_min[[3]]$log_lambda),
    c(ibm_bujeo_flowering_dens_min[[4]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_summerT_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_bujeo_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_fallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_bujeo_flowering_prevwinterR_max[[1]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_max[[2]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_max[[3]]$log_lambda),
    c(ibm_bujeo_flowering_prevwinterR_max[[4]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_bujeo_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_bujeo_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_bujeo_flowering_dens_max[[1]]$log_lambda),
    c(ibm_bujeo_flowering_dens_max[[2]]$log_lambda),
    c(ibm_bujeo_flowering_dens_max[[3]]$log_lambda),
    c(ibm_bujeo_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "Prisoneros")] = 
  c(c(ibm_prisoneros_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_prisoneros_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_prisoneros_flowering_prevwinterR_min[[1]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_min[[2]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_min[[3]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_min[[4]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_min[[1]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_min[[2]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_min[[3]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_min[[4]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_summerT_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_prisoneros_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_fallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_prisoneros_flowering_prevwinterR_max[[1]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_max[[2]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_max[[3]]$log_lambda),
    c(ibm_prisoneros_flowering_prevwinterR_max[[4]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_prisoneros_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_prisoneros_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_max[[1]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_max[[2]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_max[[3]]$log_lambda),
    c(ibm_prisoneros_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "Retin")] = 
  c(c(ibm_retin_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_retin_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_retin_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_retin_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_retin_rSurv_summerT_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_retin_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_retin_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_retin_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_retin_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_retin_rSurv_fallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_retin_flowering_prevwinterR_min[[1]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_min[[2]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_min[[3]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_min[[4]]$log_lambda),
    c(ibm_retin_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_retin_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_retin_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_retin_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_retin_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_retin_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_retin_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_retin_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_retin_flowering_dens_min[[1]]$log_lambda),
    c(ibm_retin_flowering_dens_min[[2]]$log_lambda),
    c(ibm_retin_flowering_dens_min[[3]]$log_lambda),
    c(ibm_retin_flowering_dens_min[[4]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_retin_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_retin_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_retin_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_retin_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_retin_rSurv_summerT_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_retin_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_retin_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_retin_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_retin_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_retin_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_retin_rSurv_fallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_retin_flowering_prevwinterR_max[[1]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_max[[2]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_max[[3]]$log_lambda),
    c(ibm_retin_flowering_prevwinterR_max[[4]]$log_lambda),
    c(ibm_retin_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_retin_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_retin_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_retin_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_retin_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_retin_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_retin_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_retin_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_retin_flowering_dens_max[[1]]$log_lambda),
    c(ibm_retin_flowering_dens_max[[2]]$log_lambda),
    c(ibm_retin_flowering_dens_max[[3]]$log_lambda),
    c(ibm_retin_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "MonteraTorero")] = 
  c(c(ibm_monteratorero_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_monteratorero_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_monteratorero_flowering_prevwinterR_min[[1]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_min[[2]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_min[[3]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_min[[4]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_min[[1]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_min[[2]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_min[[3]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_min[[4]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_summerT_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_monteratorero_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_fallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_monteratorero_flowering_prevwinterR_max[[1]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_max[[2]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_max[[3]]$log_lambda),
    c(ibm_monteratorero_flowering_prevwinterR_max[[4]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_monteratorero_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_monteratorero_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_max[[1]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_max[[2]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_max[[3]]$log_lambda),
    c(ibm_monteratorero_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "SCarbDist")] = 
  c(c(ibm_scarbdist_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_scarbdist_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_scarbdist_flowering_prevwinterR_min[[1]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_min[[2]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_min[[3]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_min[[4]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_min[[1]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_min[[2]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_min[[3]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_min[[4]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_summerT_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_scarbdist_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_fallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_scarbdist_flowering_prevwinterR_max[[1]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_max[[2]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_max[[3]]$log_lambda),
    c(ibm_scarbdist_flowering_prevwinterR_max[[4]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_scarbdist_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_scarbdist_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_max[[1]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_max[[2]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_max[[3]]$log_lambda),
    c(ibm_scarbdist_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "SierraCarboneraY5")] = 
  c(c(ibm_sierracarboneray5_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_sierracarboneray5_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_min[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_min[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_min[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_min[[4]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_summerT_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevwinterT_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_fallR_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_prevfallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_sierracarboneray5_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_max[[1]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_max[[2]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_max[[3]]$log_lambda),
    c(ibm_sierracarboneray5_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "SierraRetinY5")] = 
  c(c(ibm_sierraretiny5_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_sierraretiny5_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_min[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_min[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_min[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_min[[4]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_summerT_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevwinterT_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_fallR_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_prevfallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_sierraretiny5_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_max[[1]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_max[[2]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_max[[3]]$log_lambda),
    c(ibm_sierraretiny5_flowering_dens_max[[4]]$log_lambda))

results_df$lambda[which(results_df$population == "Vertedero")] = 
  c(c(ibm_vertedero_nrSurv_summerT_min[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_min[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_min[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_min[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_min[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_min[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_min[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_min[[4]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_min[[1]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_min[[2]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_min[[3]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_min[[4]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_min[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_min[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_min[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_min[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_min[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_min[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_min[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_min[[4]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_min[[1]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_min[[2]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_min[[3]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_min[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_vertedero_nrSurv_dens_min[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_min[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_min[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_min[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_min[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_min[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_min[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_min[[4]]$log_lambda),
    c(ibm_vertedero_flowering_dens_min[[1]]$log_lambda),
    c(ibm_vertedero_flowering_dens_min[[2]]$log_lambda),
    c(ibm_vertedero_flowering_dens_min[[3]]$log_lambda),
    c(ibm_vertedero_flowering_dens_min[[4]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_max[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_max[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_max[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_summerT_max[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_max[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_max[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_max[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_summerT_max[[4]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_max[[1]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_max[[2]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_max[[3]]$log_lambda),
    c(ibm_vertedero_flowering_prevwinterT_max[[4]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_max[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_max[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_max[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_fallR_max[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_max[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_max[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_max[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_fallR_max[[4]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_max[[1]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_max[[2]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_max[[3]]$log_lambda),
    c(ibm_vertedero_flowering_prevfallR_max[[4]]$log_lambda),
    rep(NA, nb_proj * nb_years),
    c(ibm_vertedero_nrSurv_dens_max[[1]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_max[[2]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_max[[3]]$log_lambda),
    c(ibm_vertedero_nrSurv_dens_max[[4]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_max[[1]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_max[[2]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_max[[3]]$log_lambda),
    c(ibm_vertedero_rSurv_dens_max[[4]]$log_lambda),
    c(ibm_vertedero_flowering_dens_max[[1]]$log_lambda),
    c(ibm_vertedero_flowering_dens_max[[2]]$log_lambda),
    c(ibm_vertedero_flowering_dens_max[[3]]$log_lambda),
    c(ibm_vertedero_flowering_dens_max[[4]]$log_lambda))


# Remove first 20 years
results_df = results_df[-which(results_df$year %in% seq(1, 20)), ]


# Transform log lambdas into lambdas and set infinite values to 0. 
# Infinite values are caused by the lack of aboveground individuals (the projections
# continue because there are still seeds in the seedbank).
results_df$lambda = exp(results_df$lambda)
# results_df$lambda[which(is.nan(results_df$lambda) | is.infinite(results_df$lambda))] = 0


# Get mean lambda per projection for each scenario
results_df = aggregate(lambda ~ population + focal_vr + focal_cov + focal_min_max + projection, 
                       data = results_df,
                       mean, na.rm = T)

# Format data
results_df$focal_vr_full = NA
results_df$focal_vr_full[which(results_df$focal_vr == "non_repro_survival")] = "Non-reproductive\nsurvival"
results_df$focal_vr_full[which(results_df$focal_vr == "repro_survival")] = "Reproductive\nsurvival"
results_df$focal_vr_full[which(results_df$focal_vr == "flowering")] = "Reproduction"

results_df$focal_cov_full = NA
results_df$focal_cov_full[which(results_df$focal_cov == "summerT")] = "Next summer mean max.\ndaily temperature"
results_df$focal_cov_full[which(results_df$focal_cov == "prevwinterT")] = "Previous winter mean max.\ndaily temperature"
results_df$focal_cov_full[which(results_df$focal_cov == "fallR")] = "Next fall cumulative\nrainfall"
results_df$focal_cov_full[which(results_df$focal_cov == "prevfallR")] = "Previous fall cumulative\nrainfall"
results_df$focal_cov_full[which(results_df$focal_cov == "prevwinterR")] = "Previous winter cumulative\nrainfall"
results_df$focal_cov_full[which(results_df$focal_cov == "dens")] = "Aboveground density of\nlarge individuals"

results_df$focal_cov_cat = NA
results_df$focal_cov_cat[which(results_df$focal_cov %in% c("summerT", "prevwinterT"))] = "Temperature"
results_df$focal_cov_cat[which(results_df$focal_cov %in% c("fallR", "prevfallR", "prevwinterR"))] = "Rainfall"
results_df$focal_cov_cat[which(results_df$focal_cov == "dens")] = "Density"




###########################################################################
#
# 3. Plotting results ----
#
###########################################################################

png(filename = "DewyPine_Lambdas_VitalRates.png", 
     width = 25,
     height = 20,
     units = "cm",
     bg = "white",
     res = 600)

ggplot(results_df, aes(x = focal_vr_full, y = lambda, fill = focal_cov_cat)) +
  geom_boxplot(alpha = 0.6) +
  # stat_summary(fun = mean, position = position_dodge(0.75)) +
  scale_fill_manual(name = "Focal covariate",
                    values = viridis(3, begin = 0.4, end = 0.8, option = "B")) +
  xlab("Focal vital-rate") +
  ylab("Population growth rate (\u03BB)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(),
        axis.title = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.title = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"))

dev.off()




###########################################################################
#
# 2. Save results ----
#
###########################################################################

write.csv(results_df, "VR_results_df.csv", row.names = F)
