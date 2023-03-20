This readme.md file was generated on 2023-03-20 by Esin Ickin

# Comparative demography project

## Introduction
This project focuses on the population-level consequences of climate effects in vital-rate models. The main hypothesis is that when we do not consider the nuanced ways in which climatic variation affects individuals across the life cycle, we substantially misrepresent  climate-change effects on population fitness.
To test this hypothesis across different species of plants and animals, we are collecting studies in which vital rates were fit as a function of climatic and/or biotic factors, and a population model was developed from which fitness (λ) can be computed. From this, changes in fitness can be calculated when covariates in vital-rate models are perturbed, either separately or considering covariation with other covariates. Such sensitivities can then be compared across species. 

Ultimately, the goal is to develop this repository to make all data and analyses for all species freely accessible. 

## Methodological information
We are building integral and matrix population models from vital rates that are parameterized as a function of climatic and/or biotic variables. We then quantify how sensitive λ is to changes in these variables. 

## File overview
You can find two R code templates, one is for matrix population models, and the other is for integral projection models. Both are an example of how to parameterize vital-rate functions using coefficients and covariate inputs to then build a population model and compute sensitivities. Both templates are provided as an R Markdown file, an R script, and pdf.

## R code
The code is divided into 5 parts.

#### 1) Vital rates
The vital rates are defined as functions and we need the coefficients from vital-rate models. Alternatively, it is possible to have (G)LM(M) models as an object and then just use the “predict” function.

#### 2) Covariates
These are the input data for the vital rate functions.

#### 3) Population model
Here, we use the vital rate functions to construct a population model that can give us the population growth rate (λ).

#### 4) Scaled sensitivity analyses
Here, we calculate scaled sensitivities according to Morris et al. 2020 (DOI: https://doi.org/10.1073/pnas.1918363117)

#### 5) Sensitivity analyses at equilibrium dynamics
Here, we perform “classic” sensitivity analyses (see Paniw et al. 2019; DOI: 10.1126/science.aau5905)

## Contact
For further information, please don’t hesitate to contact us.

Esin Ickin
MSc Student
Population Ecology Research Group
https://www.popecol.org/

University of Zurich | Department of Evolutionary Biology and Environmental Studies
Winterthurerstr. 190
8057 Zürich - Switzerland
Office: Y19-M-88
esin.ickin@uzh.ch
https://www.ieu.uzh.ch


