
# The Gray Mouse Lemur

Corresponding author: Arpat Ozgul  
Study DOI: [10.1073/pnas.2214244120](https://www.pnas.org/doi/10.1073/pnas.2214244120)

## Species

The gray mouse lemur (_Microcebus murinus_) is a small short-lived lemur with a generation time of one year (Eberle & Kappeler 2004; Schliehe-Diecks et al. 2012). The study population is located in the Kirindy forest.

## Model

Ozgul et al. (2023) modeled stage-specific vital rates using GLMs and incorporated thee covariates: Temperature, precipitation, and density. They built a two-stage and two-sex matrix population model.  
The climate and population data were provided by the author. We obtained the structure of the vital rates model, regression coefficients, and their standard errors from Table 1 of the paper. Then we rebuilt the MPM based on the annual life cycle illustrated in Figure 6 of that paper.  
We simulated the distributions of the regression coefficients. Then, we computed the scaled sensitivities by repeatedly perturbing the covariates 100 times and resampling parameters from the distributions. 

## File Overview

Main code including models and analysis:
- SensitivityAnalysis_MouseLemur.rmd

Input:
- covariates.csv

Output:
- Sens_MouseLemurs.csv
- SensVR_MouseLemurs.csv
