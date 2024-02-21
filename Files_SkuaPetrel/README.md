
# Halobaena caerulea and Catharacta lönnbergi

Corresponding author: Maud Quéroué  
Study DOI: [10.1002/ecm.1459](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1459)

## Species

Skuas and petrels were studied on Mayes Island in the Southern Ocean where the two species breed during the austral summer (Quéroué et al. 2021).
The blue petrel (Halobaena caerulea) is a small long-lived seabird with a generation time of four years (Cherel et al. 2002). 
The brown skua (Catharacta lönnbergi), preying on the petrel, is a medium-sized long-lived seabird also with a generation time of four years (Mougeot et al. 1998).

## Model

Quéroué et al. (2021) built a multispecies integrated population model. They incorporated the effects of regional and global climatic conditions and intra- and interspecific density on vital rates.
The climate and population data were obtained from the author's GitHub repository: https://github.com/maudqueroue/MultispeciesIPM_SkuaPetrel. The code for the vital rate models & population model and regression coefficients were provided by the author.
In the original model, the standardized environmental covariates were specific for each species & vital rate. Since they all fell within similar ranges and the focus was not on the within-year fluctuations of these covariates but instead on specific values (maximum, minimum, and covariation),
we decided to aggregate the average time series for these covariates: Southern annular mode, intraspecific density, interspecific density, sea surface temperature anomalies, chlorophyll a concentration.
Whereas the latter two covariates were only considered in the models for the petrel. The two stages in the model were non-breeders and breeders with the vital rates survival, reproduction, and breeding probability. 
λ was calculated by running the model for 20 years, calculating λ (Nt+1/Nt), and averaging it. The model could not be run for more than 20 years. The first 10 λ were treated as “burn-in” and were omitted. For uncertainties, we repeated everything 10 times. 
We computed sensitivities with and without covariation.

## Files Overview

Main code including model and analysis:
- SkuaPetrelSensAnalyses_MCMC.Rmd

Inputs:
- Data/Counts (folder)
- mcmc_env_cov_reg_coef.Rdata

Output:
- SensPetrel_MCMC.csv
- SensSkua_MCMC.csv
