
# The Dewy Pine

Corresponding author: Eva Conquet  
Study DOI: In Preparation

## Species

The dewy pine (_Drosophyllum lusitanicum_) is a carnivorous subshrub that is endemic to the Southwest Iberian Peninsula and Northern Morocco. It is located in the heathlands, characterized by frequent fire events and high biodiversity (Correia & Freitas 2002; Garrido et al. 2003; Paniw et al. 2017). The dewy pine has a generation time of one year (Paniw et al. 2017).

## Model

Conquet et al. used GAMs to estimate survival and flowering probability, growth and seedling size, and the number of flowers per individual using individual demographic data from 2016 until 2021 from eight populations in Southern Spain. The vital rates are linked to five covariates: temperature, rainfall, population density, size, and time since fire. They then built an individual-based model (IBM). The asymptotic population growth rate λ was calculated by projecting the population 10 years into the future, calculating λ (Nt+1/Nt), and averaging it. The model could not be run for more than 10 years. Sensitivities to each covariate were computed twice, first assuming covariation among covariates, and second, assuming average values for all covariates when perturbing one. This was repeated 100 times to get the uncertainties around the sensitivities.

## Files Overview

Model:
1) Projections_IBM_Natural_Vertedero.R: Includes the IBM and projections of the Vertedero population
2) Projections_IBM_Natural_SierraRetinY5.R: Includes the IBM and projections of the Sierra Retin Y5 population
3) Projections_IBM_Natural_SierraCarboneraY5.R: Includes the IBM and projections of the Sierra Carbonera Y5 population
4) Projections_ResultsProcessing.R: Processes the results (stored in the output folder) of the projections
Input: "Data" folder
Output: 1) Output of the first script is stored in the created "Output" folder which is then needed to run the second script and this 2) yields the output: results_df.csv

Sensitivities:  
- Calculate_Sensitivities.R: Calculates the sensitivities  
For this to work, results_df.csv and the folder "Data" are needed as inputs
Final output: Sens_DewyPines.csv
