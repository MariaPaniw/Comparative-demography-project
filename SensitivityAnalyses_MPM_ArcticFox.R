### 13-03-2023 - Esin Ickin ##

### This code is an example on how to parameterize vital-rate functions using coefficients and covariate inputs to then build a matrix population model and compute sensitivities.


### 1) Vital rates ##############################

# The example is from a paper on gray mouse lemurs: Ozgul et al. (in press)
# In the paper, we parameterized five vital rates as a function of population density, rainfall, and temperature.

# The vital rates are defined as functions, and we need the coefficients from vital rate models (here GLMs). We would appreciate a similar input on your part. Alternatively, it possible to have GLM models as an object, and then just use the "predict" function, instead of building out own functions. 



### 2) Covariates ###########################################

# There are the input data for the vital rate functions. Here, we have a time series of data. This is the best format because it allows us to calculate not only means and variances, but also covariances, the latter being important to calculate scaled sensitivities. 



### 3) Population model ######################################

# Here, we use the vital rate function to construct an annual population model that can give us the population growth rate (lambda).

# In the following example, the MPM is constructed with mean covariate values.
# In the perturbations, the covariate values are changed.



### 4) Scaled sensitivity analyses ##################################

# Here, we calculate scaled sensitivities, according to Morris et al. 2020 (DOI: https://doi.org/10.1073/pnas.1918363117)

# Note that this is a step that Esin will implement in her MS thesis. With the information given in 1-3, we should be able to run these analyses.



### 5) Sensitivity analyses at equilibrium dynamics #############################

# Here, we perform "classic" sensitivity analyses (see Paniw et al. 2019; DOI: 10.1126/science.aau5905)






