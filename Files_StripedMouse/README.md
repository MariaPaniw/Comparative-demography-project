
# The African Striped Mouse

Corresponding author: Chloé R. Nater
Study DOI:

## Species

The short-lived African striped mouse (Rhabdomys pumilio) lives in the dry regions of South Africa and has a generation time of less than one year (Mallarino et al. 2018; Nater et al. 2018). 

## Model

Nater et al. (2018) built a female-only stage-structured matrix population model with vital rates linked to temperature, food availability, and population density.
The climate and population data, and the R code was provided by the author. For the vital rates survival and maturation custom predictions functions were defined. For the remaining vital rates (breeding probability, litter probability, and litter size) R's custom "predict()" function was used.
The functions use z-standardized covariates mentioned above as well as the lagged form of each. The population matrix was then assembled using the vital rate functions.
We computed the scaled sensitivities by perturbing each covariate 100 times, resampling the parameters of the vital rate functions from the provided distributions to get the uncertainties around the estimated sensitivities. 

## Files Overview

Main code including model and analysis:
