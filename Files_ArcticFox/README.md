# The Arctic Fox

Corresponding author: Chloé R. Nater  
Study DOI: [10.1002/ecs2.3546](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.3546)

## Species
The Arctic Fox (_Vulpus lagopus_) in Svalbard Norway is a generalist and an important apex predator and scavenger without any natural predators or competitors.
They have a generation time of four years and thus are considered to be long-lived (Hemphill et al. 2020).

## Model
Nater et al. (2021) used a Bayesian integrated population model (IPM). They linked the vital rates to three covariates:
Sea ice extent, availability of reindeer carcasses, and goose population size. 
The first three steps of the model (see SensAnalysis_ArcticFox_MCMC.Rmd) were conducted by the corresponding author.
Functions to compute vital rates (pregnancy rate, litter size, denning survival, juvenile, and adult survival) based on specific values of standardized and de-trended environmental covariates were defined.
These vital rates were then assembled into a projection matrix.
We computed the matrices under two hunting scenarios (low vs high pressure) for each covariate, calculated the average, and then the scaled sensitivities.
We repeated our sensitivity analysis 100 times, each time resampling from the MCMC posterior distributions to obtain the associated uncertainties.
We calculated two types of sensitivities. First, we perturbed all vital rates and the second time we did perturbations per vital rate.

## Files Overview
Main code including vital rate models, population model, and sensitivity analysis: 
- SensAnalysis_ArcticFox_MCMC.Rmd

Inputs needed to run the code above: 
- ArcticFoxIPM_EnvCov.rds
- ArcticFox_IPM_PosteriorSamples.RData
- ArcticFoxIPM_PosteriorSummary.rds

Outputs of the main code:
- Sens_ArcticFox_MCMC.csv (Sensitivities of all vital rates)
- Sens_per_VR_ArcticFox_MCMC.csv (Sensitivities per vital rate)
