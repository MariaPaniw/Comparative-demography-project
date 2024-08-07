These R scripts here contain the main analyses done for the publication. The input of these scripts are "AllSens.csv" and for VR_Analysis it's "AllSensVR.csv". These input files contain the calculated sensitivities of all species.

- GlobalAnalysis.R: Main global analysis which includes the GLMM to assess sensitivities of species to climatic drivers and the GLMM to assess sensitivities of species to rain vs temperature.
- VR_Analysis.R: Includes GLMM to assess vital-rate specific sensitivities of each species.
- GlobalAnalysis_Taxon_specific.R: Includes the main GLMM but for mammals, birds, and plants separately fitted.
- GlobalAnalysis_excluding_simulated_lambda.R: Main analysis but without the species where lambda was calculated using simulations. So here's there are only species where lambda was analytically calculated.
- GlobalAnalysis_without_RhabdomysPumilio.R: This script contains all the main analyses but without the African striped mouse. We did this to thest whether the exclusion of this species which was a monthly MPM and not an annual makes a difference.
- Checkin_GammaDist.R: Checking that gamma distribution is the best fit.
