# APcomp
Repository for the R-based scripts used in Akabane et al. (2025)
“Vegetation and fire regimes in the Neotropics over the last 21,000 years”, Earth System Dynamics, DOI: 10.5194/esd-16-1887-2025

This repository includes scripts for:

1. Data retrieval from the Neotoma database

Using the neotoma2 R package (1_Get_NeotomaData.R).
This step selects the region of interest, downloads all available pollen datasets, creates the directory structure, and saves metadata (e.g. site location, site ID) and taxonomic classifications (ecological groups) for subsequent verification.

2. Data processing and calculation of arboreal pollen percentages

(2_Neotoma_data_treatment.R)
This step includes the calculation of arboreal pollen (AP) percentages from the downloaded datasets. Calculations involve the pollen groups TRSH, UPHE, and PALM:
Arboreal pollen percentages based on TRSH + PALM
Herb pollen percentages based on UPHE
Mangrove pollen percentages based on MANG

3. Composite curve construction for defined subregions
Using tools provided by the paleofire R package (3_Paleofire_script.R).

<strong> ✔️ Note: A single, consolidated script including all steps described above is also provided and reflects the most up-to-date version of the workflow. ✔️ </strong>


The scripts primarily address AP calculations and composite construction, as workflows for charcoal influx data have already been well established in previous studies (e.g., Blarquez et al., 2014). In contrast, equivalent workflows for AP are rarely implemented and, to my knowledge, not yet openly available. Therefore, Step 2 bridges data downloaded in Neotoma format with the procedures commonly applied to charcoal data.

<p style="color:#b22222;">
<strong>⚠️ Note:</strong> The <em>paleofire</em> R package may encounter installation issues due to the retirement of the <strong>rgdal</strong> package.
To ensure reproducibility, the scripts corresponding to the <em>paleofire</em> functions used in this study are provided as supplementary material.
Users should refer to <strong>Blarquez et al. (2014)</strong> for the original methodological framework underlying these functions.
</p>

Datasets obtained from PANGAEA or manually extracted from publications are not included in the automated scripts but are also made available in this repository.

This repository also includes the main transformed and processed datasets used to generate the graphs in Akabane et al. (2025):

It contains data on charcoal influx (CHAR) and arboreal pollen (TRSH) percentages, along with their corresponding references and site locations.
Please refer to RefSupTable.csv for additional information (reference, location, site elevation, source database) and for checking the codes used for the charcoal influx files (CHAR).


Suggestions for improvement are welcome and may be incorporated into future releases.
Please feel free to contact me if you have any questions.

e-mail: thomask.akabane@gmail.com

References: 
Blarquez, O., Vannière, B., Marlon, J. R., Daniau, A. L., Power, M. J., Brewer, S., and Bartlein, P. J.: Paleofire: An R package to analyse sedimentary charcoal records from the Global Charcoal Database to reconstruct past biomass burning, Comput. Geosci., 72, 255–261, https://doi.org/10.1016/j.cageo.2014.07.020, 2014.
