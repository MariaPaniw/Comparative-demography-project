
# The Pontic Dragonhead

Corresponding authors: Sanne M. Evers, Tomáš Dostálek, Zuzana Münzbergová  
Study DOI: in Preparation

## Species

_Dracocephalum austriacum_ is from the family Lamiaceae. The four study populations are located in the Bohemian Karst in Central Europe. Evers et al. analyzed individual-level data from 2003 - 2019.

## Model

The authors modeled five vital rates, namely survival probability, change in size, flowering probability, seed production probability, and the number of seeds produced, as being size-dependent on the log-transformed number of stems. They used functional linear models to link climatic drivers to demography. The climatic and abiotic drivers were potential evapotranspiration (PET), precipitation, and shading. They then built an integral projection model. We calculated the sensitivities to PET and precipitation. We ran the IPM through all the locations and years and then averaged it. This was repeated 100 times to get the uncertainties around our estimates.

## Files Overview

Main code (including model):
- SensNoCov.r (sensitivities without covariation) (1)
- SensCov.r (sensitivities with covariation) (2)

Inputs:
- R/functions_ipmr.R
- R/functions_GAM.R
- results/rds/VR_FLM.rds
- results/rds/state_independent_VR.rds
- data/CHELSA_data.csv
- data/Dracocephalum_with_vital_rates.csv

Outputs:
- SensNoCov.csv (1)
- SensCov.csv (2)
- load these two csv files above into the script "Sensitivities.R" to get the final output:
  - Sensitivities_Dracocephalum.csv
