
[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# An Integrated Open Population Distance Sampling Model

## What is in this repository?
This repository contains code for a workflow analysing line transect data
using an integrated distance sampling model (IDSM). 
The model utilize age-structured survey data and auxiliary
data from marked individuals to jointly estimate population dynamics and
the temporal variation in the demographic rates (recruitment rate and
survival probability) that determine the temporal transition. The model
was spesifically written for data collected through the [Norwegian
monitoring program for tetraonid
birds](https://honsefugl.nina.no/Innsyn/en) (mainly Willow Ptarmigan
*Lagopus lagopus*), but can be used for other systems that collect
age-structured line transect distance sampling data.

The model itself is written and implemented in NIMBLE 
(see [here](https://r-nimble.org/) for more information about NIMBLE for R). 

All the NIMBLE code is found in the "NIMBLE code" folder. 

Additional R functions used to download and wrangle the data, simulate data, 
prepare data in correct format, setting up and running the model, as well 
as plotting and quality control of results is contained in the R folder. 
Refer to each function's roxygen documentation for details. 

The complete workflows for analysing simulated data, real data from the 
Lierne area only, and data from all areas for which public data on willow
ptarmigan in Norway are available are provided in the following master
scripts: 

- "Analysis_SimData.R" for a single simulated dataset
- "Analysis_SimData_Replicates.R" for multiple simulated datasets including replicate runs
- "Analysis_RealData_LierneVest.R" for real data from Lierna area only
- "Analysis_RealData.R" for real data from all areas
- "_targets.R" for real data from all area, organized in a targets pipeline (more detailed information will be provided in the v2.0 update of the repository)
  
## Additional dependencies
There are a three dependencies that need to be manually installed to
run the workflow. 
First, you need to install NIMBLE (follow instructions given here: [https://r-nimble.org/download](https://r-nimble.org/download)). 
Second, the analysis uses code from the nimbleDistance package (https://github.com/scrogster/nimbleDistance). 
to estimate the half normal detection distribution. 
Finally, we use the LivingNorwayR package (https://livingnorway.github.io/LivingNorwayR/) to download and
wrangle the line transect distance sampling survey data.

## Contact
Erlend Nilsen: erlend.nilsen@nina.no

Chloé Nater: chloe.nater@nina.no
