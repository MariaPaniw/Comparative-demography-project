
# The Black-Browed Albatross

Corresponding author: Stéphanie Jenouvrier  
Study DOI: [10.1111/1365-2656.12827](https://doi.org/10.1111/1365-2656.12827)  

## Species

The black-browed albatross (_Thalassarche melanophris_) is a long-lived seabird. The study population is located at Kerguelen Island in the colony of Cañon des Sourcils Noirs (Jenouvrier et al. 2018). The capture-mark-recapture program has been going on since 1979 and at least 200 breeding pairs were monitored annually.

## Model

The authors linked sea surface temperatures (SSTs) to vital rates and built a matrix population model. To obtain the uncertainties around our estimates we extracted the standard errors of some parameters from Tbl. S2.4b of a previous study (Desprez et al. 2018). The sensitivity of lambda (with and without covariation) to SSTs was calculated 100 times, each time resampling parameters from the simulated distributions.

## Files Overview

Main code:
- MAIN_BBA.R

Inputs:
- invlogit.R
- invlogitG.R
- parameter_ENVSTO.R
- popmat.R
- SST.mat

Output:
- Sens_Albatross.csv

