
# The European Rabbit

Corresponding authors: Ana Morales González and Sanne Evers  
Study DOI: Re-analysis of [10.1371/journal.pone.0048988](https://doi.org/10.1371/journal.pone.0048988) 

## Species

The European rabbit (_Oryctolagus cuniculus_) is native to the Iberian Peninsula. While rabbit populations have become inveasive across the world, their numbers decrease across their home ranges. Rabbits are short-lived and are able to reproduce at 4 months of age, and therefore effects of environmental drivers on populaitons are strongly mediated by reproductions The study by Tablado et al. [(2012)](https://doi.org/10.1371/journal.pone.0048988) compiled information for stage-specific reproductive and survival rates across the range of rabbit occurrences and fitted an individual-based model. Here, Morales-González and Evers use the model and supplement it with data on climate and other biotic parameters from the Doñana region to project population dynamics 1996-2006.


## Model

A stochastic, stage-structured individual-based model (IBM) was developed here, as per Tablado et al. (2012). The main climatic varaible used in the model is mean monthly temperature from which (in combination with informatino on rainfall) other indices (such as food availablity and number of dry months) are derived. Density-dependence is included in all vital-rate models. 


## Files Overview

Sensitivity of all vital rates:
Inputs for the main code:
- initial_population.csv
- data_fine_gcm.csv (data on all environmental drivers based on global climate-circulation models)
  
Main code:
- rabbit_IBM_pertTemp.R (global perturbations of temperatures fixing other covariates at their averages)
- rabbit_IBM_pertTemp_cov.R (global perturbations of temperatures allowing other covariates to vary with changes in temperatures)
- rabbit_IBM_pertTemp_vr.R (vital-rate specific perturbations)

Outputs of the main codes:
- abundance_rabbit_cov.csv
- abundance_rabbit_cov.csv
- abundance_rabbit_vr_Asurv.csv
- abundance_rabbit_vr_Jsurv.csv
- use these outputs to calculate the sensitivities:
  - sensitivities_rabbits.R
  - sensitivities_rabbits_vr.R
  - this will give the final output: sens_rabbit.csv, sens_rabbit_vr.csv

