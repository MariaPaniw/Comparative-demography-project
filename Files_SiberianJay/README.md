
# The Siberian Jay

Corresponding author: Kate Layton-Matthews
Study DOI: [10.1007/s00442-018-4100-z](https://doi.org/10.1007/s00442-018-4100-z)

## Species

We used 15 years of data from a long-term study on Siberian jays, located near Arvidsjaur, northern Sweden (65°40′ N, 19°10′ E). In total, 4341 sightings from 1166 individuals were recorded during the study period. Scots pine (Pinus sylvestris) and Norway spruce (Picea abies) are the common tree species at the study site, where the latter is selectively removed in commercially exploited forests. The southern area of the study site is heavily managed, involving thinning, harvesting, and re-planting in 80–120 year cycles. Forests in the northern area of the study site have not been managed for at least 200 years, and thus, are structurally more diverse and contain trees of all age classes.

## Model

We built a periodic matrix population model. The covariates used in the vital-rate models were mean winter snow depth (December–March), average temperature during the breeding season (April–May), and population density. The sensitivity of lambda (with and without covariation) to was calculated 100 times, each time resampling parameters from the simulated distributions. We calcualted sensitivities to climatic varaibles across all vital rates or perturbing each vital rate separately. We did this for a managed and a natural population separately.

## Files Overview

Main code:
- sensitivities_Jays.R
- sensitivities_Jays_vital_rates.R

Inputs:
- covariates_jays.csv: all environmental covariates
- Survival_ModOuts_Jays.rdatav: survival models
- Rep_ModOuts_Jays.rdata: reproductive models

Output:
- sens_jays.csv
- sens_jays_pop.csv
- sens_jays_vital_rates.csv



