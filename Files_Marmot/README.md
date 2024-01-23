
# The Yellow-Bellied Marmot

Corresponding author: Maria Paniw  
Study DOI: [10.1111/ele.13459](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13459)  

## Species

Yellow-bellied marmots (_Marmota flaviventer_) are large rodents that can experience strong seasonal fluctuations in environmental conditions. The generation of this species is two years (Paniw et al. 2020).  

## Model

Paniw et al. (2020) built a seasonal stage-, mass- and environmental-specific integral projection model that accounts for seasonal demographic covariation using a latent climatic variable Q that depicts a measure of environmental quality.
Random year variaiton was sconsidered as a separate covariate. The vital rates survival, transition to the reproductive stage, growth, offspring mass, and recruitment were defined as functions and the coefficients were obtained from mean posterior values from Bayesian GLMs.  
The stages in the model were: juvenile, yearling, non-reproductive adult, and reproductive adult.  
We repeated the computations of the sensitivities 100 times and sampled parameters from the MCMC posterior distributions to obtain the uncertainties around our estimates.  
We calculated sensitivities without covariation to Q assuming no year effect. To get the sensitivities with covariation, we added the covariation among the covariates, i.e. the year value was set to the year when Q was at its maximum or minimum.
To get the sensitivities per vital rate, we created a separate R file containing the perturbation functions and repeated each analysis again.

## File Overview

Main code including model and sensitivity analysis:
- SensitivityAnalysis_Marmots_MCMC.Rmd

Inputs:
- perturbation_func_marmots_MCMC (folder)
- outputBay1_latent3.rda

Output:
- AllSens_Marmots_MCMC.csv
