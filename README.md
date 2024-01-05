---
editor_options: 
  markdown: 
    wrap: 72
---

[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

# An Integrated Open Population Distance Sampling Model

This repository contain code for an integrated distance sampling
approach. The model utilize age-structured survey data and auxiliary
data from marked individuals to jointly estimate population dynamics and
the temporal variation in the demographic rates (recruitment rate and
survival probability) that determine the temporal transition. The model
was spesifically written for data collected through the [Norwegian
monitoring program for tetraonid
birds](https://honsefugl.nina.no/Innsyn/en) (mainly Willow Ptarmigan
*Lagopus lagopus*), but can be used for other systems that collect
age-structured line transect distance sampling data.

The code is written in NIMBLE (see [here](https://r-nimble.org/) for
more information about NINBLE for R). NIMBLE is built in R, but compiles
the models and algorithms using C++. Using NIMBLE is in general
substantially faster than obtaining the posterior samples from BUGS code
based on traditional MCMC sampling.

All the NIMBLE code is found in the "NIMBLE code" folder. Additional R
code used to download and wrangle the data, prepare data in correct
format for the modelling, as well as code for plotting and quality
control is contained in the R folder. The folder also contain functions
that were used to simulate data that can be used to assess the model
robustness. The functions are documented using roxygen-skeletons.

There are a few dependencies that might need to be manually installed to
run the model code. First, you need to install NIMBLE (follow
instructions given [here](https://r-nimble.org/download)). Second, we
use code from the nimbleDistance package
([here](https://r-nimble.org/download)) to estimate the half normal
detection function. Finally, we use the LivingNorwayR-package (located
[here](https://livingnorway.github.io/LivingNorwayR/)) to download and
wrangle the line transect distance sampling survey data.
