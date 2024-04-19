
# Meerkats

Corresponding author: Maria Paniw  
Study DOI: [10.1126/science.aau5905](https://www.science.org/doi/10.1126/science.aau5905)

## Species

Meerkats (_Suricata suricatta_) are small social mammals living up to 12 years (Bateman et al. 2013; Clutton-Brock & Manser 2016).  
The study population is located in the Kuruman River Reserve in South Africa.

## Model

The covariates used in the model were rainfall, temperature, and population density. Paniw et al. (2019) used the most parsimonious generalized additive models (GAMs) of vital rates (survival, growth, reproduction, stage transitions, and emigration) to construct a density-dependent, environment-specific, mass-stage-classified integral projection model.  
The data and code were provided by the author.  
We calculated the sensitivities of all vital rates 100 times, each time resampling from the parameter distributions. Sensitivities per vital rates were also calculated in a similar manner.

## File Overview

Main code including model and analysis:
- Sens_Meerkats_all_VR.R (perturbing rain and temperature in all vital rates at the same time)
- SensVR_Rain_Meerkats.R (perturbing rain in single vital rates)
- SensVR_Temp_Meerkats.R (perturbing temperature in single vital rates)

Inputs:
- boot.GAM
- meerkats-master (folder)

Outputs:
- Sens_Meerkats.csv (output of Sens_Meerkats_all_VR.R)
- SensVR_Rain_Meerkats.csv (output of SensVR_Rain_Meerkats.R)
- SensVR_Temp_Meerkats.csv (output of SensVR_Temp_Meerkats.R)
