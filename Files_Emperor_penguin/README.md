# The Emperor penguin

Corresponding author: Stéphanie Jenouvrier  
Study DOI: [10.1111/j.1365-2486.2012.02744.x](https://doi.org/10.1111/j.1365-2486.2012.02744.x)
For details on constructing the periodic matrix population model, see also DOI: [10.1086/652436](https://doi.org/10.1086/652436)

## Species
The emperor penguin (_Aptenodytes forsteri_) The population model here is based on monitoring of breeding emperor penguins at Dumont D’Urville, Terre Adélie, in Antarctica. The colony has been monitored every year, during the breeding season (March–December), from 1962 onwards. Various components of the sea ice environment affect different parts of the emperor penguin life cycle during different seasons. Emperor penguins breed on motionless sea ice (i.e. fast ice) during the Antarctic winter. They arrive at the colony sometime in late March to early April while sea ice is thickening, and leave the colony in late December before the ice breaks up. The colony site is usually far from the ocean, and during the breeding season emperor penguins travel long distances to feed in ice-free areas, such as polynyas, within the sea ice cover.

Emperor Penguins become sexually mature at three years of age.

## Model
Jenouvrier et al. (2012) used a sex- (males and females) and stage-structured (pre-breeders, breeding pairs, non-breeder) periodic (seasonal) matrix population model following Jenouvrier et al. [2010](https://doi.org/10.1086/652436). The climatic covariates in vital-rate models were proportional anomalies in sea-ice concentration (SIC), relative to the mean from 1979 to 2007 in the pre-breeding, laying, incubating, and rearing seasons. The code presented here uses parameters from the vital-rate models (different functions) loaded to the main script to construct the periodic MPM. Some of the parameters are obtained from model averageing, as explained in Jenouvrier et al. (2012). Breeding in this model is dependent on male/female ratios in the populations, which are generated internally in the model. Due to this, we used randomly sampled initial population sizes, and then project the annual population dynamics for 1000 years until the sex ratios and the population growth rates, λ, converged to a stable distribution. We calculated λ as the dominant eigenvalue of the annual product of the seasonal matrix population model in year 1000. Deu to the strcuture of the model, we could not perform sensitivity analyses by perturbing only specific vital rates.    

## Files Overview
Main code including vital rate models, population model, and sensitivity analysis: 
- SensAnalysis_EmperorPenguins.R

Inputs needed to run the code above: 
- SIC.xlsx

Outputs of the main code:
- Sens_emperor_peng.csv
