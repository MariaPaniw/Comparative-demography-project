
# Swiss Birds

Corresponding author: Anne Malchow  
Study DOI: [10.1098/rstb.2022.0194](https://royalsocietypublishing.org/doi/10.1098/rstb.2022.0194)

## Species

1. Eurasian bullfinch (_Pyrrhula pyrrhula_)
2. European crested tit (_Lophophanes cristatus_)
3. Eurasian treecreeper (_Certhia familiaris_)
4. Eurasian nuthatch (_Sitta europaea_)
5. Dunnock (_Prunella modularis_)
6. Common linnet (_Linaria cannabina_)
7. Ring ouzel (_Turdus torquatus_)
8. Alpine accentor (_Prunella collaris_)

## Model

Malchow et al. (2023) developed spatially explicit, process-based models for eight Swiss breeding bird populations using distribution and abundance data.
They chose bird species with a one-year generation time and shared common traits. They built a female-only, two-stage population model consisting of the stages "juvenile" and "breeding adult" and three vital rates: juvenile and adult survival and fecundity.
Each vital rate is linked to relevant climatic predictors. In the fecundity model, mean temperature and total precipitation during the breeding season were the climatic predictors. In the juvenile survival model, it was the mean temperature during autumn and total precipitation during winter. 
For adult survival, the predictors were also the total precipitation during winter and the minimum temperature during winter. The regression coefficients were provided by the author. The climate and population data and the code were obtained from the author's GitHub repository: https://github.com/UP-macroecology/Malchow_DemogEnv_2022.

In the main code, the vital rate models were first defined. Then, species-specific maximum and minimum covariates and the covariation were specified. The covariates are standardized with a mean of zero and a standard deviation of one. Then, we assembled the vital rate functions into a matrix in a function. To compute the sensitivities, we set the habitat suitability index to 100 for all species. We calculated the sensitivities of all vital rates, once without and once with covariation among covariates. This was repeated 100 times, each time resampling parameters from the posterior distributions. We also calculated the sensitivities per vital rate.

## Files Overview

Main code(s) including model and analysis (species-specific) are in the code folder:
- Certhia_MCMC.Rmd
- Linaria_MCMC.Rmd
- Lophophanes_MCMC.Rmd
- PrunellaCollaris_MCMC.Rmd
- PrunellaModularis_MCMC.Rmd
- Pyrrhula_MCMC.Rmd
- Sitta_MCMC.Rmd
- Turdus_MCMC.Rmd

Inputs:
- Post_sample_list.rds
Species-specific climate data:
- CerthiaCovariates.csv
- LinariaCovariates.csv
- LophophanesCovariates.csv
- PrunellaCollarisCovariates.csv
- PrunellaModularisCovariates.csv
- PyrrhulaCovariates.csv
- SittaCovariates.csv
- TurdusCovariates.csv

Outputs:
- are in the output folder
- two csv files for each species, one containing the "usual" sensitivities and one vital-rate specific sensitivities

Species-specific climate data were extracted in this script: CovariateExtraction_Birds.R and for this, these were the inputs: 
- MHB_aggregated_counts.Rdata
- aggrVariables_zTrafod.Rdata
- CA_comb_Spec10_it60000_540-560.Rdata
- CA_comb_Spec2_it60000_542-562.Rdata
- CA_comb_Spec3_it60000_543-563.Rdata
- CA_comb_Spec4_it60000_544-564.Rdata
- CA_comb_Spec5_it60000_545-565.Rdata
- CA_comb_Spec7_it60000_547-567.Rdata
- CA_comb_Spec8_it60000_548-568.Rdata
- CA_comb_Spec9_it60000_549-569.Rdata

Outputs:
Outputs of each species are in the folder "Output". Each csv file starts with the species name same as above. Sensitivities of all vital rates and sensitivities per vital rate are in two different files for each species. Example: Certhia_Sens.csv and Certhia_Sens_per_VR.csv
