This readme.md file was generated on 2023-03-20 by Esin Ickin

🌎🐧🦒🌿🌍

**GENERAL INFORMATION**
1. Title of Database/Repository: Comparative demography project
2. Corresponding Author Information:
- First author:
   - Name: Esin Ickin; Institution: Department of Evolutionary Biology and Environmental Studies, University of Zurich, Switzerland; Email: ickin.esin@gmail.com

- Senior author (Species: Marmot, Shrubs, Meerkat):
   - Name: Maria Paniw; Institution: Doñana Biological Station (EBD-CSIC), Spain; Email: maria.paniw@ebd.csic.es
 
3. List of all co-authors and their contributions can be found in "ListofAuthors.docx"
     
**SHARING/ACCESS INFORMATION**
1. Licenses/restrictions placed on the data and script: Please contact the authors if you want to use the data and scripts for further analyses and to prevent data misuse.
2. Sources for data:
- Arctic Fox: [10.1002/ecs2.3546](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.3546)
- Striped Mouse: [10.1111/1365-2656.12888](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12888)
- Marmot: [10.1111/ele.13459](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13459)
- Shrubs: [10.1098/rspb.2022.1494](https://royalsocietypublishing.org/doi/10.1098/rspb.2022.1494)
- Petrel: [10.1002/ecm.1459](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1459)
- Swiss Birds: [10.1098/rstb.2022.0194 ](https://royalsocietypublishing.org/doi/10.1098/rstb.2022.0194)
- Dipper: [10.1126/sciadv.1602298](https://www.science.org/doi/10.1126/sciadv.1602298)
- Meerkat: [10.1126/science.aau5905](https://www.science.org/doi/10.1126/science.aau5905)
- Mouse Lemur: [10.1073/pnas.2214244120](https://www.pnas.org/doi/10.1073/pnas.2214244120)
- Opuntia: [10.1111/gcb.15519](https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.15519) & [10.1002/ecm.1228](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1228)
- Spanish Trees: [10.1093/jpe/rtw081](https://academic.oup.com/jpe/article/10/5/731/3062498)
- Dewy Pine: In Preparation
- Dracocephalum: In Preparation
- Albatross: [10.1111/1365-2656.12827](https://doi.org/10.1111/1365-2656.12827)
- Giraffes: [10.1111/gcb.16970](https://doi.org/10.1111/gcb.16970)
- Magellanic Penguin: [10.1073/pnas.2209821120](https://doi.org/10.1073/pnas.2209821120)
- Protea repens: [10.1111/ecog.00839](https://doi.org/10.1111/ecog.00839)
- Reindeer: [10.1038/s41467-019-09332-5](https://doi.org/10.1038/s41467-019-09332-5)
- Barnacle goose: [10.1111/gcb.14773](https://doi.org/10.1111/gcb.14773)
- Siberian jay: [10.1007/s00442-018-4100-z](https://doi.org/10.1007/s00442-018-4100-z)

**DATA & FILE OVERVIEW**
1. Introduction: This project focuses on the population-level consequences of climate effects in vital-rate models. The main hypothesis is that when we do not consider the nuanced ways in which climatic variation affects individuals across the life cycle, we substantially misrepresent climate-change effects on population fitness. To test this hypothesis across different species of plants and animals, we are collecting studies in which vital rates were fit as a function of climatic and/or biotic factors, and a population model was developed from which fitness (λ) can be computed. From this, changes in fitness can be calculated when covariates in vital-rate models are perturbed, either separately or considering covariation with other covariates. Such sensitivities can then be compared across species. Ultimately, the goal is to develop this repository to make all data and analyses for all species freely accessible.  
2. Methodological Information: Models include: Matrix population models (example: Dipper), integral projection models (example: Marmot), integrated population models (similar to MPM), and individual-based models (example: see word file or giraffes). The idea is to standardize them so that we can run sensitivity analyses.
3. File overview: There is one R script or R Markdown file in each species- or study-specific folder that contains the population models and sensitivity analyses. All the input files are in the same folder, as well as the output files. The R code is divided into 4 parts: 1) Vital rates: This part contains the vital rate models, which vary among studies. 2) Covariates: These are the input data for the vital rate functions. 3) Population model: Here, the vital rates functions are used to construct a population model that can give us the population growth aka fitness (λ). 4) Scaled Sensitivity analyses: Sensitivities (1. assuming no covariation 2. assuming covariation) are calculated according to Morris et al. 2020 ([https://doi.org/10.1073/pnas.1918363117](https://doi.org/10.1073/pnas.1918363117))
4. Main Analysis: This is the folder containing the scripts for the main analyses and input files for these scripts. We fitted GLMM to assess sensitivities of populations to climatic drivers. For more information see README.md file in this folder.
   
