
# SETUP #
#-------#

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set switches

# Re-simulate data
resimulate <- FALSE

# Recruitment per adult or per adult female
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE


## Set scaling parameter
scale1 <- 1000

# SIMULATE DATA #
#---------------#

if(resimulate){
  source("DataSimulation_Full.R")
}else{
  AllSimData <- readRDS("SimData_Full.rds")
}


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- list(
  nim.data = list(
    sumR_obs = AllSimData$Rep.data$sumR_obs,
    sumAd_obs = AllSimData$Rep.data$sumAd_obs,
    y = AllSimData$DS.data$d,
    L = AllSimData$DS.data$L,
    N_a_line_year = AllSimData$DS.data$DS.count,
    Survs1 = AllSimData$RT.data$Survs1,
    Survs2 = AllSimData$RT.data$Survs2
  ),
  
  nim.constants = list(
    N_years = AllSimData$SimParams$Tmax,
    Tmin.RT = AllSimData$SimParams$Tmin.RT,
    Tmax.RT = AllSimData$SimParams$Tmax.RT,
    W = AllSimData$SimParams$W,
    scale1 = scale1,
    A = colSums(AllSimData$DS.data$L)*(AllSimData$SimParams$W/scale1)*2,
    N_obs = length(AllSimData$DS.data$d),
    Year_obs = AllSimData$DS.data$d_year,
    N_sites = AllSimData$SimParams$Jmax,
    sumR_obs_year = AllSimData$Rep.data$sumR_obs_year,
    N_sumR_obs = AllSimData$Rep.data$N_sumR_obs,
    N_ageC = AllSimData$SimParams$Amax
  )
)


# MODEL SETUP #
#-------------#

# Original version (zeroes-trick)
# model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM.R",
#                           customDist = FALSE,
#                           nim.data = input_data$nim.data,
#                           nim.constants = input_data$nim.constants,
#                           testRun = FALSE, initVals.seed = 0)
  
# Updated version (nimbleDistance::dHN)
model_setup <- setupModel(modelCode.path = "NIMBLE Code/RypeIDSM_dHN.R",
                          customDist = TRUE,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          testRun = FALSE, initVals.seed = 0)

# MODEL (TEST) RUN #
#------------------#
t.start <- Sys.time()
IDSM.out <- nimbleMCMC(code = model_setup$modelCode, 
                       data = input_data$nim.data, 
                       constants = input_data$nim.constants,
                       inits = model_setup$initVals, 
                       monitors = model_setup$modelParams,
                       nchains = model_setup$mcmcParams$nchains, 
                       niter = model_setup$mcmcParams$niter, 
                       nburnin = model_setup$mcmcParams$nburn, 
                       thin = model_setup$mcmcParams$nthin, 
                       samplesAsCodaMCMC = TRUE, 
                       setSeed = 0)
Sys.time() - t.start

#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t15.rds')
#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t30.rds')
saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t15_MuR.rds')
#saveRDS(IDSM.out, file = 'rypeIDSM_dHN_simData_t30_MuR.rds')

