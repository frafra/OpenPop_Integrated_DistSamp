runModel <- function(parallelMCMC,
                     model_setup, MCMC.seeds, input_data){
  
  
  if(!parallelMCMC){
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
                           setSeed = MCMC.seeds)
    
  }else{
    
    ## Add toggles to constants
    input_data$nim.constants$fitRodentCov <- fitRodentCov
    input_data$nim.constants$survVarT <- survVarT
    input_data$nim.constants$R_perF <- R_perF
    
    ## Set up cluster
    this_cluster <- makeCluster(model_setup$mcmcParams$nchains)
    #clusterEvalQ(this_cluster, library(nimble))
    #clusterEvalQ(this_cluster, library(nimbleDistance))
    
    ## Collect chain-specific information
    per_chain_info <- vector("list", model_setup$mcmcParams$nchains)
    for(i in 1:model_setup$mcmcParams$nchains){
      per_chain_info[[i]] <- list(mySeed = MCMC.seeds[i],
                                  inits = model_setup$initVals[[i]])
    }
    
    ## Run chains in parallel
    IDSM.out <- parLapply(cl = this_cluster, 
                          X = per_chain_info, 
                          fun = runMCMC_allcode, 
                          model_setup = model_setup,
                          input_data = input_data)
    
    stopCluster(this_cluster)
    
  }
  
  ## Return posterior samples
  return(IDSM.out)
}