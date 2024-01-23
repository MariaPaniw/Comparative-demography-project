
# The Tree Cholla Cactus

Corresponding authors: Sanne M. Evers & Aldo Compagnoni  
Study DOI: [10.1111/gcb.15519](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.15519) & [10.1002/ecm.1228](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1228)

## Species

The tree cholla cactus (_Opuntia imbricata_) is from the family of Cactaceae nad has a generation time of at least nine years (Ohm & Miller 2014).
The study population is located at the Sevilleta National Wildlife Refuge in New Mexico, USA (Miller et al. 2009; Ohm & Miller 2014). 
## Model

The integral projection model is based on the R script from Aldo Compagnoni (Compagnoni et al. 2016). Further, two vital rate models are from the analysis of Sanne M. Evers (Evers et al. 2021).  
The R code and inputs were provided by the author (Evers). We calculated the scaled sensitivities by peturbing different temperature anomaly variables and repeated the analysis 10 times.

## Files Overview

Main code including model and analysis:
- OPIM_IPM_Resampling.R

Inputs:
- data (folder)
- results (folder)

Outputs:
- Sens_w_resampling.csv
- Use this output in the script "ReOrganize.R" under the section "OPUNTIA IMBRICATA" to reorganize it and get SENS_OPIM.csv
