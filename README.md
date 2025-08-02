# APcomp
Repository for the R-based script used in Akabane et al. (in review) – Vegetation and fire regimes in the Neotropics over the last 21,000 years. DOI: https://doi.org/10.5194/egusphere-2025-1424

This repository will be updated alongside the publication of the manuscript.

It contains raw data on charcoal influx and arboreal pollen (AP) percentages, along with their corresponding references and site locations.
Included here are the scripts used in the study, organized into three main steps:

1 - Data retrieval from the Neotoma database using the neotoma2 R package.

2 - Data processing, including the calculation of arboreal pollen percentages from the downloaded datasets.

3 - Composite curve construction for defined subregions, using tools provided by the paleofire R package.

The scripts here focus primarily on AP calculations and composites, as workflows for charcoal influx have already been well established in previous studies (e.g., Blarquez et al., 2014). In contrast, similar workflows for AP are rarely implemented and, to my knowledge, not yet openly available. Thus, Step 2 provides the link between data downloaded in Neotoma format and its preparation for use with the same methods typically applied to charcoal data.
Additional datasets obtained from PANGAEA or manually extracted from publications are not included in the automated script, but they are also made available in this repository.

Suggestions for improvement are welcome, as the current version of the scripts is not yet fully organized.
A clearer description and improved instructions for the scripts will be added in the near future.

Current link for the manuscript: https://egusphere.copernicus.org/preprints/2025/egusphere-2025-1424/


References:
Blarquez, O., Vannière, B., Marlon, J. R., Daniau, A. L., Power, M. J., Brewer, S., and Bartlein, P. J.: Paleofire: An R package to analyse sedimentary charcoal records from the Global Charcoal Database to reconstruct past biomass burning, Comput. Geosci., 72, 255–261, https://doi.org/10.1016/j.cageo.2014.07.020, 2014.
