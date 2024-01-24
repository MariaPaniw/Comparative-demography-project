
# The African Striped Mouse

Corresponding author: Chloé R. Nater  
Study DOI:

## Species

The short-lived African striped mouse (_Rhabdomys pumilio_) lives in the dry regions of South Africa and has a generation time of less than one year (Mallarino et al. 2018; Nater et al. 2018). 

## Model

Nater et al. (2018) built a female-only stage-structured matrix population model with vital rates linked to temperature, food availability, and population density.
The author provided the climate and population data and the R code. For the vital rates survival and maturation custom prediction functions were defined. For the remaining vital rates (breeding probability, litter probability, and litter size) R's custom "predict()" function was used.
The functions use z-standardized covariates mentioned above and the lagged form of each. The population matrix was then assembled using the vital rate functions.
We computed the scaled sensitivities by perturbing each covariate 100 times, resampling the parameters of the vital rate functions from the provided distributions to get the uncertainties around the estimated sensitivities. We calculated the sensitivities of all vital rates and also per vital rate.

## Files Overview

Main code including model and analysis:
- SensitivityAnalysis_StripedMouse_MCMC.Rmd

Inputs:
- 170503_allmodels_final.RData
- pert.functions.striped.mouse2.mcmc.R

Outputs:
- Sens_StripedMouse_MCMC.csv
- Sens_per_VR_StripedMouse_MCMC.csv
