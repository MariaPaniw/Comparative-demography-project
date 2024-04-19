
# The Svalbard Reindeer

Corresponding authors: Marlène Gamelon, Brage B. Hansen 
Study DOI: [10.1038/s41467-019-09332-5](https://doi.org/10.1038/s41467-019-09332-5)  

## Species

The study population of the wild Svalbard reindeer (_Rangifer tarandus_) is located in central Spitsbergen, Svalbard, Norway (Hansen et al. 2019).

## Model

The authors built GLMMs, linking climate and density to survival and fecundity. The two climatic variables were rain-on-snow (ROS) and winter length. For more information see the corresponding study (link above). We rebuilt the vital rate models and a female-only MPM containing six age classes. The sensitivity of lambda (with and without covariation) was calculated 100 times, each time resampling parameters from the simulated distributions.
We computed two types of sensitivities, perturbing all vital rates and perturbing each vital rate.

## Files Overview

Main code (including model and sensitivity analysis):
- SensAnalysis_Reindeer.rmd

Inputs:
- coeffic_lmer_withintercept.Rdata
- coefficF_lmer_withintercept.Rdata
- data_reindeer.txt

Outputs:
- SensReindeer.csv (sensitivities)
- SensReindeer_VR.csv (sensitivities per vital rate)

