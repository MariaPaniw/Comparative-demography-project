
# The Dewy Pine

Corresponding author: Eva Conquet  
Study DOI: In Preparation

## Species

The dewy pine (_Drosophyllum lusitanicum_) is a carnivorous subshrub that is endemic to the Southwest Iberian Peninsula and Northern Morocco. It is located in the heathlands, characterized by frequent fire events and high biodiversity (Correia & Freitas 2002; Garrido et al. 2003; Paniw et al. 2017). The dewy pine has a generation time of one year (Paniw et al. 2017).

## Model

Conquet et al. used GAMs to estimate survival and flowering probability, growth and seedling size, and the number of flowers per individual using individual demographic data from 2016 until 2021 from eight populations in Southern Spain. The vital rates are linked to five covariates: temperature, rainfall, population density, size, and time since fire. They then built an individual-based model (IBM). The asymptotic population growth rate λ was calculated by projecting the population 50 years into the future, calculating λ (Nt+1/Nt), and taking the average of the last 25 years to discard transient dynamics. Sensitivities to each covariate were computed twice, first assuming covariation among covariates, and second, assuming average values for all covariates when perturbing one. This was repeated 100 times to get the uncertainties around the sensitivities. We also calculated vital rate specific sensitivities.

## Files Overview

Models:
1) For covariation vs no covariation:
 - Projections_Covariation_IBM_PopulationName.R: IBMs, "PopulationName" is a placeholder for any of the eight populations
 - Projections_Covariation_ResultsProcessing.R: Processing results of the IBMs
 - Covariation_Sensitivities.R: Calculating the results from the processing
 
 2) For vital-rate specific sensitivities:
 - Projections_VitalRates_IBM_PopulationName.R: IBMs, "PopulationName" is a placeholder for any of the eight populations
 - Projections_VitalRates_ResultsProcessing.R: Processing the results of the IBMs
 - VR_Sensitivities.R: Calculating the results from the processing

Final outputs: 
1) Sens_DewyPines.csv
2) SensVR_DewyPines.csv
