# APcomp
Repository for the R-based scripts used in Akabane et al. (2025)
“Vegetation and fire regimes in the Neotropics over the last 21,000 years”, Earth System Dynamics

This repository contains scripts for:


1 - Data retrieval from the Neotoma database using the neotoma2 R package. (1_Get_NeotomaData.R)
Select a region, download all pollen data from the region of interest, create folders to store data, save metadata (location, site id...) and classification of taxa (ecological groups), so these can be checked.

2 - Data processing, including the calculation of arboreal pollen (AP) percentages from the downloaded datasets. (2_Neotoma_data_treatment.R)
Calculations involve TRSH, UPHE, and PALM.
Calculate arboreal pollen percentages based on TRSH+PALM.
Calculate herb pollen percentages based on UPHE.
Calculate mangrove pollen percentages based on MANG.

3 - Composite curve construction for defined subregions, using tools provided by the paleofire R package. (3_Paleofire_script.R)


Single_single script includes all the described steps.


The scripts primarily address AP calculations and composite construction, as workflows for charcoal influx data have already been well established in previous studies (e.g., Blarquez et al., 2014). In contrast, equivalent workflows for AP are rarely implemented and, to my knowledge, not yet openly available. Therefore, Step 2 bridges data downloaded in Neotoma format with the procedures commonly applied to charcoal data.

Datasets obtained from PANGAEA or manually extracted from publications are not included in the automated scripts but are also made available in this repository.

This repository also includes the main transformed and processed datasets used to generate the graphs in Akabane et al. (2025):

It contains data on charcoal influx (CHAR) and arboreal pollen (TRSH) percentages, along with their corresponding references and site locations.
Please refer to RefSupTable.csv for additional information (reference, location, site elevation, source database) and for checking the codes used for the charcoal influx files (CHAR).




Suggestions for improvement are welcome and may be incorporated into future releases.
Please feel free to contact me if you have any questions.

e-mail: thomask.akabane@gmail.com

References: 
Blarquez, O., Vannière, B., Marlon, J. R., Daniau, A. L., Power, M. J., Brewer, S., and Bartlein, P. J.: Paleofire: An R package to analyse sedimentary charcoal records from the Global Charcoal Database to reconstruct past biomass burning, Comput. Geosci., 72, 255–261, https://doi.org/10.1016/j.cageo.2014.07.020, 2014.
