
# Spanish Trees

Corresponding author: David García-Callejas  
Study DOI: [10.1093/jpe/rtw081](https://academic.oup.com/jpe/article/10/5/731/3062498)

## Species

García-Callejas et al. (2016) built spatially explicit integral projection models of several tree species across Spain. These species include:
- _Fagus sylvatica_
- _Pinus halepensis_
- _Pinus nigra_
- _Pinus pinaster_
- _Pinus pinea_
- _Pinus sylvestris_
- _Pinus uncinata_
- _Quercus faginea_
- _Quercus ilex_
- _Quercus robur/petraea_
- _Quercus suber_

## Model

The vital rates in the model are sapling recruitment, sapling transition to the adult stage, survival, and growth of adult trees. These are the environmental covariates included in the vital rate models: temperature, precipitation, and their anomalies.  
We calculated the sensitivities to temperature and precipitation by projecting the populations 90 years under different scenarios into the future and then approximating the long-term population growth rate λ with (Nt+1/Nt) and averaging it. To get the scaled sensitivities to temperature and precipitation (Morris et al. 2020) we need to calculate λ under different scenarios. These scenarios were done in separate R scripts and the scaled sensitivities were then calculated in another script. 

## Files Overview

Each scenario (which covariate perturbed and whether covariation was included or not) has its own R script with the obvious naming. For example, if temperature is perturbed (i.e. set to its maximum value) and other covariates are kept at their mean values, then the R script name would be "MaxTempNoCov.R". And when temperature is set to its minimum value it would be "MinTempNoCov.R". And with the results of these two, we can calculate the scaled sensitivities to temperature of each species. (Scaled Sensitivity to T = (lambda(MaxTempNoCov)-lambda(MinTempNoCov))/((MaxTemp-MinTemp)/SD), (Morris et al. 2020))

Main R scripts including IPMs and projections:
- MaxTempCov.R
- MinTempCov.R
- MaxPrecipCov.R
- MinPrecipCov.R
  
- MaxTempNoCov.R
- MinTempNoCov.R
- MaxPrecipNoCov.R
- MinPrecipNoCov.R

These are the inputs the main scripts above need:
- files in the "R" folder
- files in the "data" folder and sub-folders
- ClimateCovariates.csv which was generated with this script: reorganize_climate_data.R
    - and for this, we need the inputs in the folder "data/clima"

The main R scripts will generate folders and subfolders to store the results of each species and simulation run. For example, the results of the R script MaxTempCov.R will go into the folder MaxTempCov/run1 when it's the first run. When running the model again, change the subfolder to run2, etc. The output is stored in the files named "results.spi.csv with i being the ID number of a species (i.e., MaxTempCov/run1/results.sp1.csv)

The sensitivities to temperature and precipitation of each species are then calculated in the script "sensitivities_trees.R" with the results above. This generates this file: "SensSpainTrees.csv"



