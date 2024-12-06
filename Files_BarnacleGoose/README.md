
# The Barnacle Goose

Corresponding author: Kate Layton-Matthews
Study DOI: [10.1111/gcb.14773](https://doi.org/10.1111/gcb.14773)

## Species

The Svalbard barnacle goose population is migratory and overwinters at Solway Firth, Scotland, before flying to Svalbard for breeding in summer. The barnacle goose reaches sexual maturity at the age of 2 years. Individual mark–recapture data from both sexes were collected from the nesting islands and coastal area around Ny-Ålesund from 1990-2017 (for this study).

## Model

We built a matrix population model comprised of 2 states, fledglings and adults. Vital rates were modeled as a function of various seasonal temperature and rainfall measures at breeding and overwindering grounds, as well as populaiton density and predation pressure were mean daily minimum temperatures October-March in Scotland and in April-May in Helgeland, mean precipitation in April-May in Helgeland, the flyway population size at the wintering grounds in Scotland, spring onset, adult numbers in Svalbard, and fox predation. The sensitivity of lambda (with and without covariation) to was calculated 100 times, each time resampling parameters from the simulated distributions. We calcualted sensitivities to climatic varaibles across all vital rates or perturbing each vital rate separately.

## Files Overview

Main code:
- sensitivities_Geese_full.R
- - sensitivities_Geese_vital_rate.R

Inputs:
- env_covar_scaled_geese.csv: all environemtnal covariates
- RepMods_outAll_v2.rdata: reproduction models
- Survival_ModOut_v2.rdata: survival models

Output:
- sens_goose_full.csv
- sens_goose_vital_rate.csv


