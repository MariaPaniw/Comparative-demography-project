
# The Magellanic Penguin

Corresponding authors: T.J. Clark-Wolf, P. Dee Boersma, Briana Abrahams  
Study DOI: [10.1073/pnas.2209821120](https://doi.org/10.1073/pnas.2209821120)  

## Species

The study population of _Spheniscus magellanicus_ is at the Punta Tombo colony. These penguins breed after the age of four, sometimes even skipping breeding due to poor body conditions.

## Model

Clark-Wolf and colleagues built a pre-breeding, three-stage, female-only integrated population model- The authors assumed a simple stage structure (juveniles, adults, and immigrants). They linked climatic and oceanographic variables to vital rates. For more details see the corresponding study (link above). The code of the IPM and data were obtained from the author's [GitHub repository](https://github.com/teejclark/Press_Pulse). We then computed the sensitivities to average monthly max temperature, average monthly precipitation, and sea surface temperature anomalies. The analyses were repeated 100 times, to get the uncertainties around our estimates. We also calculated the sensitivities per vital rate.

## Files Overview

Main code (including model and sensitivity analysis):
- first, you need to run the script CRData.R
- SensAnalysis_MagellanicPenguins.R (all sensitivities)
- SensVR_MagellanicPenguins.R (sensitivities per vital rate)

Inputs:
- all files in the folder "Data"

Outputs:
- Sens_MPenguins.csv (all sensitivities)
- SensVR_MPenguins.csv (sensitivities per vital rate)

