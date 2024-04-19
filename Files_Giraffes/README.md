
# The Masai Giraffe

Corresponding authors: Monica L. Bond, Derek E. Lee  
Study DOI: [10.1111/gcb.16970](https://doi.org/10.1111/gcb.16970)  

## Species

The study population of Masai giraffes (_Giraffa camelopardalis tippelskirchi_ [subspp.] or _G. tippelkirchi_ [spp.]) is located in northern Tanzania.


## Model

Bond et al. (2023) developed a stochastic, socially structured individual-based model (IBM) using long-term demographic data. The two variables used in the model are rain and population density. We ran the simulation 100 times to get the sensitivities to rain with and without covariation. Then, we computed the vital-rate-specific sensitivities, again 100 times, by perturbing the covariates in each vital rate, namely survival of calves, survival of juveniles, and survival of adults.

## Files Overview

Sensitivity of all vital rates:
Inputs for the main code:
- abund_validation2023.csv
- InitPopGiraffe.csv
  
Main code (including model):
- giraffe_IBM_pertRain.R
- giraffe_IBM_pertRainCov.R

Outputs of the main codes:
- abund.pert.rain.csv
- abund.pert.rain.cov.csv
- use these outputs and "abund_validation2023.csv" in this script to calculate the sensitivities:
  - sensitivities_giraffes.R
  - this will give the final output: Sens_Giraffes.csv

Sensitivities per vital rate (all in folder SensVR, except for the first two inputs):
Inputs for the main code:
- abund_validation2023.csv
- InitPopGiraffe.csv

Main code (including model):
- giraffe_IBM_Rain_SurvNeo.R (survival of newborns)
- giraffe_IBM_Rain_SurvJuv.R (survival of juveniles)
- giraffe_IBM_Rain_SurvAd.R (survival of adults)

Outputs of the main codes:
- abund.rain.survneo.csv
- abund.rain.survj.csv
- abund.rain.survad.csv
- use these outputs and "abund_validation2023.csv" in this script to calculate the sensitivities:
  - sensitivitiesVR_giraffes.R
  - this will give the final output: SensVR_Giraffes.csv



