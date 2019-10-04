# megafaunalFrugivore

## Author & Reference
**Author:** Jun Ying Lim

**Contact:** junyinglim@gmail.com

**Reference:** Lim, J.Y., Svenning, J.-C., Göldel, B., Faurby, S. & Kissling, W.D. Past and future extinctions shape the body size - fruit size relationship between palms and mammalian frugivores.

## System requirements
To run this code you will need a current installation of R and some R packages. This includes both the analysis of data and the plotting of the figures in the manuscript. 

Code was tested on R version 3.5.1 and the following packages (and versions accordingly): `stringr` (v1.4.0), `scales` (v1.0.0), `cowplot` (v1.0.0), `reshape2` (v1.4.3), `MuMIn` (v1.43.6), `viridis` (v0.5.1), `wesanderson` (v0.3.6), `RColorBrewer` (v1.1-2), `ggplot2` (v3.2.1), `rgdal` (v1.4-4), `sp` (v1.3-1), `relaimpo` (v2.2-3), `car` (v3.0-2), `care` (v1.1.10), `spatialreg` (v1.1-3), `spdep` (v1.1-2)

## Installation guide:
Installation of all neccessary packages may be performed by running the following line of code in R.
`install.packages(c("rgdal, "relaimpo", "care", "ggplot2", "MuMIn", "plyr", "car", "spdep", "spatialreg", "wesanderson", "reshape2", "scales", "cowplot", "RColorBrewer", "stringr"), dependencies = TRUE)`


## Instructions for use:
To perform a full analysis, the following R scripts have to be run in sequence. A description of 

**Step 1: Compute body size and fruit size at the scale of botanical countries (`generateTDWGdata.R`)**
* Generate summary statistics (i.e., maximum and median body size and fruit size for each botanical country)
* Imputation of trait values for palm species without data
* Implement probabilistic simulations of future extinction

**Step 2: Implement model averaging (`analyzeData.R`)**
* Performs model averaging of ordinary least squares (OLS) and spatial autoregressive (SAR) linear models of fruit size at both global and regional scales
* Projects fruit size under future scenarios of defaunation
* Uses convenience functions in `ggmodavg.R`

**Step 3: Plot figures (`figures.R`)**
* Produces all figures in the manuscript

## Data:
Raw data is included (in the `data` folder) to allow for reproducibility of our results

* Palm checklist (presence-absence) data at scale of botanical countries (`palms_in_tdwg3.csv`)
* Palm fruit size trait data (`PalmTraits_10.csv`)
* Botanical country climatic variables (`TDWG_Environment_AllData_2019Feb.csv`)
* Mammal body size and IUCN Red List status (`Phylacine_Trait_data.csv`; see Faurby et al. "PHYLACINE 1.2: The phylogenetic atlas of mammal macroecology." Ecology 99.11 (2018): 2626)
* Mammal checklist (presence-absence) data at the scale of botanical countries for both current and present-natural scenarios (`mammal_curr_occ.csv`, `mammal_presnat_occ.csv`)







