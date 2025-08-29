
# The Pontic Dragonhead

Corresponding authors: Sanne M. Evers, Tomáš Dostálek, Zuzana Münzbergová  
Study DOI: 10.5281/zenodo.16942068

## Species

_Dracocephalum austriacum_ is from the family Lamiaceae. The four study populations are located in the Bohemian Karst in Central Europe. Evers et al. analyzed individual-level data from 2003 - 2019.

## Model

The authors modeled five vital rates, namely survival probability, change in size, flowering probability, seed production probability, and the number of seeds produced, as being size-dependent on the log-transformed number of stems. They used functional linear models to link climatic drivers to demography. The climatic and abiotic drivers were potential evapotranspiration (PET), precipitation, and temperature. They then built an integral projection model. We calculated the sensitivities to temperature and precipitation.

## Files Overview

Main code (including model):
- new_ipm_stoch_analysis.R

Inputs:
- R/functions_ipmr.R
- R/functions_GAM.R
- results/rds/VR_FLM.rds
- results/rds/state_independent_VR.rds
- data/CHELSA_data.csv
- data/Dracocephalum_with_vital_rates.csv

Outputs:
- Sensitivities_Dracocephalum.csv
