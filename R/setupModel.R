#' Set up model, data, and initial values for running MCMC
#'
#' @param modelCode.path string. Relative path to the model file to be used
#' @param customDist logical. If TRUE, uses custom half-normal distribution from
#' nimbleDistance package.  
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param niter integer. Number of MCMC iterations (default = 25000)
#' @param nthin integer. Thinning factor (default = 5)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 5000)
#' @param nchains integer. Number of chains to run.
#' @param testRun logical. If TRUE, sets up for a test run with 10 iterations,
#' no thinning, and no burn-in (default = FALSE)
#' @param initVals.seed integer. Seed to use for inital value simulation.
#'
#' @return list of list containing all components necessary for running model 
#' with `nimble::nimbleMCMC()`
#' @export
#'
#' @examples

setupModel <- function(modelCode.path, customDist,
                       nim.data, nim.constants,
                       niter = 50000, nthin = 10, nburn = 10000, nchains = 3,
                       testRun = FALSE, initVals.seed){

  
  ## Catch mismatches between model code name and distribution settings
  if((customDist & (!(grepl('dHN', modelCode.path, fixed = TRUE)) & !(grepl('dHR', modelCode.path, fixed = TRUE)))) |
     (!customDist & (grepl('dHN', modelCode.path, fixed = TRUE) & grepl('dHR', modelCode.path, fixed = TRUE)))){
    stop('Mismatch between model code name and distribution settings. Check inputs for modelCode.path and customDist.')
  }
  
  ## Load model code
  require('nimble')
  if(customDist){require('nimbleDistance')}
  source(modelCode.path)
  
  ## Set parameters to monitor
  params <- c("esw", "p", #"D",
              "R_year", "Mu.R", "h.Mu.R", "h.sigma.R", "sigmaT.R",
              "sigma", "mu.dd", "sigmaT.dd",
              "Density", "N_exp", "N_tot_exp",
              "Mu.D1", "sigma.D",
              "S", "Mu.S1", "Mu.S2", "h.Mu.S1", "h.Mu.S2", "h.sigma.S1", "h.sigma.S2",
              "ratio.JA1")
  
  if(grepl('dHR', modelCode.path, fixed = TRUE)){
    params <- c(params, "b")
  }
  
  ## Simulate initial values
  set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    initVals[[c]] <- simulateInits(nim.data = nim.data, nim.constants = nim.constants)
  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 50
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelCode = rypeIDSM,
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(niter = niter, nthin = nthin, 
                      nburn = nburn, nchains = nchains)
  )
   
  return(setup) 
}
