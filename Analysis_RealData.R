
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

# (Re-)downloading data
# downloadData <- FALSE
downloadData <- TRUE

# Recruitment per adult or per adult female
R_perF <- TRUE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- FALSE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData){
  Rype_arkiv <- downloadLN(version = 1.6, save = TRUE)
}else{
  stop("downloadData = FALSE not supported yet. There is an issue with encoding when using LivingNorwayR::initializeDwCArchive() that needs to be resolved first.")
  #Rype_arkiv <- initializeDwCArchive("data/Rype_arkiv.zip")
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities and time period of interest
localities <- "Lierne Fjellst. Vest"
minYear <- 2015
maxYear <- 2020

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_LineTrans(DwC_archive = Rype_arkiv, 
                                 localities = localities,
                                 minYear = minYear, maxYear = maxYear)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
d_cmr <- wrangleData_CMR()


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = sumR.Level,
                               dataVSconstants = TRUE,
                               save = TRUE)


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
                          testRun = TRUE, initVals.seed = 0)

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

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_realData_Lierne.rds')

