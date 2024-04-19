
# The Common Sugarbush

Corresponding author: Cory Merow  
Study DOI: [10.1111/ecog.00839](https://doi.org/10.1111/ecog.00839)  

## Species

The common sugarbush (_Protea repens_) is a shrub found throughout the Mediterranean climate of the Cape Floristic Region in South Africa.

## Model

Merow et al. (2014) linked individual survival to factors like size, winter temperature, and precipitation, constructing then an integral projection model. All data and code were accessible online. However, to obtain the posterior samples we ran the regression models again.  

We calculated the sensitivity of lambda (with and without covariation) to temperature and precipitation 100 times, each time resampling parameters from the simulated distributions.
We computed two types of sensitivities, perturbing all vital rates and perturbing each vital rate.

## Files Overview

Main code (including model and sensitivity analysis):
- ProteaSensAnalyses_MCMC.rmd

Inputs:
- PR_adultsurv.post
- PR_flowering.post
- PR_germ.post
- PR_growth.post
- PR_offspring.post
- PR_seed.post
- PR_seedhead.post
- PR_seedlingsurv.post
- Protea_IPM_Shared_Data.Rdata

Outputs:
- Sens_Protea.csv (sensitivities)
- SensVR_Protea.csv (sensitivities per vital rate)

