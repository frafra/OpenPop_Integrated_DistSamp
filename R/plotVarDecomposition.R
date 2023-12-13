#' Decompose variance by component and plot proportions
#'
#' @param mcmc.out n mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas in the analysis.
#' @param N_years integer. Number of years in the analysis.
#' @param fitRodentCov logical. Indicates whether (TRUE) or not (FALSE) variance
#' due to effect of a rodent covariate should be considered.
#' @param RodentOcc_data numeric vector containing the rodent occupancy data
#' used for fitting the model. Optional argument that needs to be provided if
#' fitRodenCov = TRUE.
#'
#' @return a list containing a dataframe of posterior distributions for proportions
#' of variance explained by different components, as well as the path to the 
#' plot. 
#' @export
#'
#' @examples

plotVarDecomposition <- function(mcmc.out, N_areas, N_years, fitRodentCov, RodentOcc_data = NULL){
  
  ## Check rodent data is provided if needed
  if(fitRodentCov & missing("RodentOcc_data")){
    stop("RodentOcc_data has to be provided if fitRodentCov = TRUE.")
  }
  
  ## Convert posterior samples into matrix
  sam.mat <- as.matrix(mcmc.out)
  nsamples <- nrow(sam.mat)
  
  ## Extract posterior distributions for random variance components
  # Random year variation
  varR_year <- sam.mat[, "sigmaT.R"]^2
  varS_year <- sam.mat[, "sigmaT.S"]^2
  varDet_year <- sam.mat[, "sigmaT.dd"]^2
  
  # Random area variation
  varR_area <- sam.mat[, "h.sigma.R"]^2
  varS_area <- sam.mat[, "h.sigma.S"]^2
  varDet_area <- sam.mat[, "h.sigma.dd"]^2
  
  # Residual variation
  varR_res <- sam.mat[, "sigmaR.R"]^2
  varS_res <- sam.mat[, "sigmaR.S"]^2
  varDet_res <- sam.mat[, "sigmaR.dd"]^2
  
  
  ## Calculate posterior distribution for variance due to area-specific rodent effect
  if(fitRodentCov){
    
    varR_rodent <- rep(NA, nsamples)
    
    for(i in 1:nsamples){
      
      rodent_Eff <- c()
      
      for(x in 1:N_areas){
        
        betaR.R <- sam.mat[i, paste0("betaR.R[", x, "]")]
        
        RodentOcc <- RodentOcc_data[x, ]
        for(t in 1:N_years){
          if(is.na(RodentOcc[t])){
            RodentOcc[t] <- sam.mat[i, paste0("RodentOcc[", x, ", ", t, "]")]
          }
        }
        
        rodent_Eff <- cbind(rodent_Eff, betaR.R*RodentOcc)
      }
      
      varR_rodent[i] <- var(as.vector(rodent_Eff))
    }
    
  }else{
    
    varR_rodent <- rep(0, nsamples)
  }
  
  
  ## Calculate proportions of variance explained by different components
  
  # Recruitment
  propVar.R <- data.frame(
    Parameter = "Recruitment",
    VarComponent = rep(c("Year", "Area", "RodentOcc", "Residual"), each = nsamples),
    PropVar = c(varR_year / (varR_year + varR_area + varR_rodent + varR_res),
                varR_area / (varR_year + varR_area + varR_rodent + varR_res),
                varR_rodent / (varR_year + varR_area + varR_rodent + varR_res),
                varR_res / (varR_year + varR_area + varR_rodent + varR_res))
  )
  
  # Survival
  propVar.S <- data.frame(
    Parameter = "Survival",
    VarComponent = rep(c("Year", "Area", "Residual"), each = nsamples),
    PropVar = c(varS_year / (varS_year + varS_area + varS_res),
                varS_area / (varS_year + varS_area + varS_res),
                varS_res / (varS_year + varS_area + varS_res))
  )
  
  # Detection
  propVar.Det <- data.frame(
    Parameter = "Detection",
    VarComponent = rep(c("Year", "Area", "Residual"), each = nsamples),
    PropVar = c(varDet_year / (varDet_year + varDet_area + varDet_res),
                varDet_area / (varDet_year + varDet_area + varDet_res),
                varDet_res / (varDet_year + varDet_area + varDet_res))
  )
  
  # Combined
  propVar.all <- rbind(propVar.R, propVar.S, propVar.Det)
  propVar.all$Parameter <- factor(propVar.all$Parameter, levels = c("Recruitment", "Survival", "Detection"))
  
  if(!fitRodentCov){
    propVar.all <- subset(propVar.all, VarComponent != "RodentOcc")
  }
  
  ## Make plotting directory if it does not exist yet
  ifelse(!dir.exists("Plots/VarDecomposition"), dir.create("Plots/VarDecomposition"), FALSE)
  
  ## Plot variance decomposition for all parameters
  pdf("Plots/VarDecomposition/VarDecomposition_overall.pdf", width = 6, height = 6)
  print(
    ggplot(propVar.all) + 
      geom_density(aes(x = PropVar, color = VarComponent, fill = VarComponent), alpha = 0.5) + 
      facet_wrap(~Parameter, ncol = 1, scales = "free_y") + 
      ylab("Density") + xlab("Proportion of total variance") +
      #paletteer::scale_fill_paletteer_d("nbapalettes::grizzlies_00s") + 
      #paletteer::scale_color_paletteer_d("nbapalettes::grizzlies_00s") + 
      scale_fill_manual(values = c("#3B99B1FF", "#A2C194FF", "#EABB22FF", "#E87200FF")) + 
      scale_color_manual(values = c("#3B99B1FF", "#A2C194FF", "#EABB22FF", "#E87200FF")) + 
      theme_classic()
  )
  dev.off()
  
  ## Return data and plot path
  return(list(propVar.post = propVar.all, 
              plotPath = "Plots/VarDecomposition/VarDecomposition_overall.pdf"))
}

