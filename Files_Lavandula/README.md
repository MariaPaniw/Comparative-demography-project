# Lavandula stoechas

Corresponding author: T.Sanchez-Mejía  
Study DOI: Sanchez-Mejía et al., in prep; [10.1098/rspb.2022.1494](https://royalsocietypublishing.org/doi/10.1098/rspb.2022.1494)

## Species

French lavender (_Lavandula stoechas_) is a common species in Mediterranean shurblands, including Doñana Naitonal Park, where the species was studied (abundance data since 2007; individual demogrpahic data since 2019).  

## Model
Paniw et al. (2023) built a metapopulation model including this species. Sanchez-Mejía then expanded on this model with adiditional invidual-based data on the species. Here, we constructed an integral projection model (IPM) based on individual plant size (instead of using three discrete size categories); and vital rate were fitted as functions of different seasonal climatic varaibles. Climatic variables (1980-2024) were obtained from http://icts.ebd.csic.es/es/datos-meteorologicos and https://meteorologia-palacio.icts-donana.es/hydroyear/. Plant density was also used as a covararite, but affected only reproduction rates as an offset. The vital-rate models included a ramdom site effect, and we used a representative site to project population dynamics (averaging over site did not affect results).     
The R code of the model, climate and population data, and posterior samples were provided by the corresponding authors. We computed the scaled sensitivities for all vital rates and for each vital rate (VR). This was done 100 times, each time resampling parameters using parametric bootstrapping.

## Files Overview

Main codes: sens_lavandula.R (for global sensitivities) and sens_lavandula_vr.R (for vital-rate specific sensitivities)

Inputs:
- Lavandula_GAM.Rdata (vital-rate models parameterized using generalized additive models)
- Lavandula_data.Rdata (all the climatic and individual plant monitoring data necessary to run the main codes)

Outputs:
- sens_lavandula.csv
- sens_lavandula_vr.csv



