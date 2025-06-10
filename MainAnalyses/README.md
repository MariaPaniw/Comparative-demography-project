These R scripts here contain the main analyses done for the publication. The input of these scripts are "AllSens.csv" and for VR_Analysis it's "AllSensVR.csv". These input files contain the calculated sensitivities of all species.
main_plot_Fig1.tiff is the main output of GlobalAnalysis.R
main_plot_Fig2.tiff is the main output of VR_Analysis.R

- GlobalAnalysis.R: Main global analysis (i.e., climatic drivers perturbed across all vital rates) which includes all the GLMMs to assess sensitivities of species to climatic drivers.
- VR_Analysis.R: Includes GLMM to assess vital-rate specific sensitivities of each species, where applicable.
- GlobalAnalysis_without_RhabdomysPumilio.R: This script contains all the main analyses but without the African striped mouse. We did this to thest whether the exclusion of this species which was a monthly MPM and not an annual makes a difference.
- Checkin_GammaDist.R: Checking that gamma distribution is the best fit.
- checking_sensitivity_single_populations.R: This script contains the main analysis but only on single-population species
- checkin_sensitivity_population_randomeffect.R: This script contains the main analysis but has population as random effect
