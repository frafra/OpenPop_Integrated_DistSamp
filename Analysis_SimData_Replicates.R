
# GENERAL SETUP #
#---------------#

## Set number of datasets to simulate
N_datasets <- 10

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


# SIMULATE REPLICATE DATASETS #
#-----------------------------#

## Make directory (if not present)
if(!dir.exists("simData")){
  dir.create("simData")
}

## Select seeds randomly
seed.list <- sample(1:1000, size = N_datasets, replace = FALSE)
#seed.list <- c(120, 199, 216, 365, 541, 596, 706, 794, 919, 956)

for(i in 1:N_datasets){
  
  ## Set seed randomly
  mySeed <- seed.list[i]
  
  ## Simulate dataset
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
                                save = FALSE)
  
  ## Assemble input data object
  input_data <- prepareInputData_Sim(SimData = AllSimData,
                                     addDummyDim = addDummyDim)
  
  ## Save dataset and input data object with custom name
  saveRDS(AllSimData, file = paste0("simData/AllSimData_seed", mySeed, ".rds"))
  saveRDS(input_data, file = paste0("simData/inputData_seed", mySeed, ".rds"))
}

## Save list of seeds
saveRDS(seed.list, file = "simData/seedList.rds")


#*******************************************************************************

# RUN MODEL SEVERAL TIMES PER DATASET #
#-------------------------------------#

library(tidyverse)
library(tidybayes)

## Retrieve simulation seeds
simSeed.list <- readRDS("simData/seedList.rds")

## Make directories (if not present)
if(!dir.exists("simModelFits")){
  dir.create("simModelFits")
}

if(!dir.exists("simModelFits_sum")){
  dir.create("simModelFits_sum")
}

## Set number of run replicates
N_runs <- 3
runSeed.list <- list()

for(i in 1:length(simSeed.list)){
  
  ## Select seeds randomly (and store)
  runSeeds <- sample(1:100, size = N_runs, replace = FALSE)
  runSeed.list[[i]] <- runSeeds
  names(runSeed.list)[i] <- paste0("simSeed_", simSeed.list[i])
  
  ## Load corresponding dataset
  input_data <- readRDS(paste0("simData/inputData_seed", simSeed.list[i], ".rds"))
  
  for(k in 1:length(runSeeds)){
    
    ## Set run seed
    mySeed <- runSeeds[k]
    
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
                              #niter = 500000, nthin = 5, nburn = 300000, nchains = 3,
                              testRun = FALSE,
                              initVals.seed = mySeed)
    
    ## Run model
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
    
    ## Save full posteriors (incl. seed information)
    saveRDS(list(samples = IDSM.out,
                 simSeed = simSeed.list[i],
                 runSeed = runSeeds[k]),
            file = paste0("simModelFits/IDSMsamples_simSeed", simSeed.list[i], "_runSeed", runSeeds[k], ".rds"))
    
    ## Summarise posteriors to minimum necessary for plotting sim checks
    
    # Recruitment parameters
    R_year <- IDSM.out %>% tidybayes::spread_draws(R_year[year])
    Mu_R <- IDSM.out %>% tidybayes::spread_draws(Mu.R) %>% mutate(lab_code = "Mu.R")
    sigmaT_R <- IDSM.out %>% tidybayes::spread_draws(sigmaT.R) %>% mutate(lab_code = "sigmaT.R")
    
    # Survival parameters
    Mu_S1 <- IDSM.out %>% tidybayes::spread_draws(Mu.S1) %>% mutate(Surv = "S1") %>% rename(S = Mu.S1) %>% select(S, Surv)
    Mu_S <- IDSM.out %>% tidybayes::spread_draws(Mu.S) %>% mutate(Surv = "S") %>% rename(S = Mu.S) %>% select(S, Surv)
    Mu_S_data <-  tibble(S = Mu_S$S/Mu_S1$S, Surv = "S2") %>% bind_rows(., Mu_S1, Mu_S)
    
    # Detection parameters
    Mu_dd <- IDSM.out %>% tidybayes::spread_draws(mu.dd) %>% mutate(lab_code = "mu.dd")
    sigmaT_dd <- IDSM.out %>% tidybayes::spread_draws(sigmaT.dd) %>% mutate(lab_code = "sigmaT.dd")
    esw_year <- IDSM.out %>% tidybayes::spread_draws(esw[year])
    p_year <- IDSM.out %>% tidybayes::spread_draws(p[year])
    
    # Population sizes
    N_tot <- IDSM.out %>% tidybayes::spread_draws(N_tot_exp[year]) 
    
    # Population densities
    A_temp <- apply(input_data$nim.data$L, 2, sum) * input_data$nim.constants$W*2 / (1000 *1000)
    Density_year <- IDSM.out %>% tidybayes::spread_draws(N_tot_exp[year]) %>% 
      dplyr::mutate(density = (N_tot_exp/A_temp))
    
    ## Collate and save summarized posteriors
    sumPost <- list(
      sum.post = list(
        R_year = R_year, Mu_R = Mu_R, sigmaT_R = sigmaT_R,
        Mu_S_data = Mu_S_data,
        Mu_dd = Mu_dd, sigmaT_dd = sigmaT_dd, esw_year = esw_year, p_year = p_year,
        N_tot = N_tot, Density_year = Density_year
      ),
      simSeed = simSeed.list[i],
      runSeed = runSeeds[k])
    
    saveRDS(sumPost, file = paste0("simModelFits_sum/IDSMsampleSum_simSeed", simSeed.list[i], "_runSeed", runSeeds[k], ".rds"))
    
    ## Remove samples and free up disk space
    rm(IDSM.out)
    gc()
  }
  
}

## Save complete seed information
saveRDS(runSeed.list, file = "simModelFits_sum/seedInfo.rds")


#*******************************************************************************

# PLOT COMPARISON OF MODEL ESTIMATES AND SIMULATED DATA #
#-------------------------------------------------------#

plotSimCheck_replicates(thin = 100)
plotSimCheck_replicates("Temps", thin = 100)
plotSimCheck_replicates("Zissou1", thin = 100)

plotSimCheckMetrics_replicates("Temps", thin = 1)

