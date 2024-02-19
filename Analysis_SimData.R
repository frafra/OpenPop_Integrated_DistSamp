
# GENERAL SETUP #
#---------------#

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
resimulate <- TRUE

# Recruitment per adult or per adult female
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Time variation in survival
survVarT <- FALSE

# Rodent covariate on reproduction
fitRodentCov <- FALSE

# Addition of dummy dimension for running multi-area setup
addDummyDim <- FALSE

# Random effects shared across areas
if(survVarT & addDummyDim){
  shareRE <- FALSE
}else{
  shareRE <- TRUE
}



# SET SIMUALATION PARAMETERS #
#----------------------------#

# General simulation parameters
#---

mySeed <- 0

Amax <- 2 # Number of age classes
Tmax <- 15 # Number of years
Jmax <- 50 # Number of sites/transect lines


# Vital rate parameters
#---

## Annual survival
Mu.S <- 0.35 # Average annual survival probability
if(survVarT){
  sigmaT.S <- 0.8 # SD of random year variation in survival
}else{
  sigmaT.S <- 0 # SD of random year variation in survival
}

sigmaJ.S <- 0 # SD of random site variation in survival

## Reproduction
Mu.R <- 2 # Average number of chicks in August
sigmaT.R <- 0.4 # SD of random year variation in number of chicks
sigmaJ.R <- 0 # SD of random site variation in number of chicks

## Juvenile summer survival
#Mu.sJ <- 0.2 # Average summer survival of chicks
#sigmaT.sJ <- 0 # SD of random year variation in chick survival
#sigmaJ.sJ <- 0 # SD of random site variation in survival


# Population parameters
#---

# Initial population numbers per site
N1_juv_limits <- c(3, 8)

# Average group size
avg_Gsize <- 5.6


# Data & observation parameters 
#---

## Line-transect distance sampling
min.Tlength <- 1000 # Minimum transect length
max.Tlength <- 1000  # Maximum transect length

W <- 200 # Truncation distance (max. distance at which observation is possible)

Mu.dd <- 75 # Average width parameter for half-normal detection function
sigmaT.dd <- 0.3 # SD of random year variation in detection probability
sigmaJ.dd <- 0 # SD of random line variation in detection probability

## Known-fate radio-telemetry
Tmin.RT <- 5 # First year for which radio-telemetry data has been collected
Tmax.RT <- 10 # Last year for which radio-telemetry data has been collected

# Average number of individuals fitted with transmitters each year
nind.avg.RT <- 30

# ## Nest survey
# Tmin.NS <- 1 # First year for which nest survey data has been collected
# Tmax.NS <- 15 # Last year for which nest survey data has been collected
# 
# # Average number of nests monitored each year
# nind.avg.NS <- 40



# SIMULATE DATA #
#---------------#

if(resimulate){
  AllSimData <- assembleSimData(Amax = Amax, Tmax = Tmax, Jmax = Jmax,
                                avg_Gsize = avg_Gsize, 
                                Mu.S = Mu.S, sigmaT.S = sigmaT.S, sigmaJ.S = sigmaJ.S,
                                Mu.R = Mu.R, sigmaT.R = sigmaT.R, sigmaJ.R = sigmaJ.R,
                                Mu.dd = Mu.dd, sigmaT.dd = sigmaT.dd, sigmaJ.dd = sigmaJ.dd,
                                W = W, min.Tlength = min.Tlength, max.Tlength = max.Tlength,
                                nind.avg.RT = nind.avg.RT, 
                                Tmin.RT = Tmin.RT, Tmax.RT = Tmax.RT,
                                seed = mySeed, 
                                R_perF = R_perF,
                                R_parent_drop0 = R_parent_drop0,
                                stochasticSim = TRUE,
                                plotPopSim = TRUE,
                                save = TRUE)
}else{
  AllSimData <- readRDS("SimData_Full.rds")
}


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

input_data <- prepareInputData_Sim(SimData = AllSimData,
                                   addDummyDim = addDummyDim)


# MODEL SETUP #
#-------------#

## Determine correct code path
modelCode.path <- selectCodePath(shareRE = shareRE,
                                 survVarT = survVarT,
                                 addDummyDim = addDummyDim)

## Set up model
model_setup <- setupModel(modelCode.path = modelCode.path,
                          customDist = TRUE,
                          R_perF = R_perF,
                          shareRE = shareRE, 
                          survVarT = survVarT, 
                          addDummyDim = addDummyDim,
                          fitRodentCov = fitRodentCov,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          niter = 500000, nthin = 5, nburn = 300000, nchains = 4,
                          testRun = TRUE,
                          initVals.seed = mySeed)

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
                       setSeed = mySeed)
Sys.time() - t.start

if(survVarT){
  if(addDummyDim){
    saveRDS(IDSM.out, file = paste0("rypeIDSM_dHN_multiArea_survT_simData_s", Jmax, "_t", Tmax, ".rds"))
  }else{
    saveRDS(IDSM.out, file = paste0("rypeIDSM_dHN_survT_simData_s", Jmax, "_t", Tmax, ".rds"))
  }
}else{
  if(addDummyDim){
    saveRDS(IDSM.out, file = paste0("rypeIDSM_dHN_multiArea_simData_s", Jmax, "_t", Tmax, ".rds"))
  }else{
    saveRDS(IDSM.out, file = paste0("rypeIDSM_dHN_simData_s", Jmax, "_t", Tmax, ".rds"))
  }
}




# MCMC TRACE PLOTS #
#------------------#

plotMCMCTraces(mcmc.out = IDSM.out,
               fitRodentCov = fitRodentCov)



# MODEL CHECKS AGAINST SIMULATION PARAMETERS #
#--------------------------------------------#

if(!addDummyDim){
  plotSimCheck(SimData = AllSimData,
               mcmc.out = IDSM.out)
}else{
  stop("Function plotSimCheck() is not currently adapted for running for the multi-area setup (dummy dimension not supported)")
}


# MODEL COMPARISON #
#------------------#

## Multi-area vs. single-area setup
mcmc.out <- readRDS("rypeIDSM_dHN_multiArea_simData_s50_t15.rds")
NodeNames <- dimnames(mcmc.out[[1]])[[2]]
NodeNames_drop <- gsub("[1, ", "[", NodeNames, fixed = TRUE)
for(i in 1:model_setup$mcmcParams$nchains){
  dimnames(mcmc.out[[i]])[[2]] <- NodeNames_drop
}
saveRDS(mcmc.out, "rypeIDSM_dHN_multiArea_simData_s50_t15_dimDrop.rds")

modelComp <- plotModelComparison(modelPaths = c("rypeIDSM_dHN_multiArea_simData_s50_t15_dimDrop.rds", 
                                                "rypeIDSM_dHN_simData_s50_t15.rds"), 
                                 modelChars = c("Multi-area setup", "Single-area setup"), 
                                 N_sites = input_data$nim.constants$N_sites, N_years = input_data$nim.constants$N_years,
                                 plotPath = "Plots/ModelComp_MultiVSSingleSetup",
                                 returnData = FALSE)

