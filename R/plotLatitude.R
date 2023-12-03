plotLatitude <- function(mcmc.out, 
                         N_areas, area_names, area_coord, N_sites, 
                         min_years, max_years, minYear, maxYear,
                         fitRodentCov){
  
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Set year index ranges for plotting shared data collection period only
  minYearIdx_shared <- max(min_years)
  maxYearIdx_shared <- min(max_years)
  
  ## Determine area-specific year range
  area_yearIdxs <- (1:(maxYear-minYear+1))
  area_years <- area_yearIdxs + (minYear - 1)
  
  ## Calculate area- and year-specific population growth rates
  for(i in 1:N_areas){
    for(t in 1:length(area_years)){
      
      # Summarize annual average population densities (per area, averaged over lines)
      popDens_juv_t0 <- out.mat[, paste0("Density[",  i, ", 1, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      popDens_ad_t0 <- out.mat[, paste0("Density[",  i, ", 2, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      popDens_mean_t0 <- rowMeans(popDens_juv_t0 + popDens_ad_t0)
      
      if(t < length(area_years)){
        popDens_juv_t1 <- out.mat[, paste0("Density[",  i, ", 1, ", 1:N_sites[i], ", ", area_yearIdxs[t+1], "]")]
        popDens_ad_t1 <- out.mat[, paste0("Density[",  i, ", 2, ", 1:N_sites[i], ", ", area_yearIdxs[t+1], "]")]
        popDens_mean_t1 <- rowMeans(popDens_juv_t1 + popDens_ad_t1)
        
        # Calculate annual average population growth rate
        lambda_t0t1 <- popDens_mean_t1 / popDens_mean_t0
      }
      
      # Merge average density back into posterior samples
      out.mat <- cbind(out.mat, popDens_mean_t0)
      colnames(out.mat)[dim(out.mat)[2]] <- paste0("meanDens[",i, ", ", area_yearIdxs[t], "]" )
      
      
      # Merge population growth rate back into posterior samples
      if(t < length(area_years)){
        out.mat <- cbind(out.mat, lambda_t0t1)
        colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambda[",i, ", ", area_yearIdxs[t], "]" )
      }
    }
    
    ## Calculate average area-specific average densities and population growth rate over different time periods
    
    # Entire time period - average density
    lambda_mean_tot <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[1:(length(area_years)-1)], "]")])
    out.mat <- cbind(out.mat, lambda_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_tot[", i, "]")
    
    # Overlapping time period (all areas) - average density
    lambda_mean_shared <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", minYearIdx_shared:(maxYearIdx_shared-1), "]")])
    out.mat <- cbind(out.mat, lambda_mean_shared)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_shared[", i, "]")
    
    # Entire time period - population growth rate
    lambda_mean_tot <- rowMeans(out.mat[, paste0("lambda[", i, ", ", area_yearIdxs[1:(length(area_years)-1)], "]")])
    out.mat <- cbind(out.mat, lambda_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambdaAvg_tot[", i, "]")
    
    # Overlapping time period (all areas) - population growth rate
    lambda_mean_shared <- rowMeans(out.mat[, paste0("lambda[", i, ", ", minYearIdx_shared:(maxYearIdx_shared-1), "]")])
    out.mat <- cbind(out.mat, lambda_mean_shared)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambdaAvg_shared[", i, "]")
  }
  
  ## Summarize posteriors for relevant parameters
  
  # Prepare matrices for storage of results
  rRep.sum <- pSurv.sum <- betaR.sum <- data.frame()
  popDens.sum <- lambda.sum <- data.frame()
  
  for(i in 1:N_areas){
    
    # Summarize average reproductive rates
    rRep_name <- paste0("Mu.R[",  i, "]")
    rRep_add <- data.frame(Area = area_names[i],
                           Median = median(out.mat[, rRep_name]),
                           lCI = unname(quantile(out.mat[, rRep_name], probs = 0.025)),
                           uCI = unname(quantile(out.mat[, rRep_name], probs = 0.975)),
                           Mean = mean(out.mat[, rRep_name]),
                           SD = sd(out.mat[, rRep_name])
    )
    rRep.sum <- rbind(rRep.sum, rRep_add)
    
    # Summarize average annual survival rates
    pSurv_name <- paste0("Mu.S[",  i, "]")
    pSurv_add <- data.frame(Area = area_names[i],
                            Median = median(out.mat[, pSurv_name]),
                            lCI = unname(quantile(out.mat[, pSurv_name], probs = 0.025)),
                            uCI = unname(quantile(out.mat[, pSurv_name], probs = 0.975)),
                            Mean = mean(out.mat[, pSurv_name]),
                            SD = sd(out.mat[, pSurv_name])
    )
    pSurv.sum <- rbind(pSurv.sum, pSurv_add)
    
    # Summarize rodent effects
    if(fitRodentCov){
      betaR_name <- paste0("betaR.R[",  i, "]")
      betaR_add <- data.frame(Area = area_names[i],
                              Median = median(out.mat[, betaR_name]),
                              lCI = unname(quantile(out.mat[, betaR_name], probs = 0.025)),
                              uCI = unname(quantile(out.mat[, betaR_name], probs = 0.975)),
                              Mean = mean(out.mat[, betaR_name]),
                              SD = sd(out.mat[, betaR_name])
      )
      betaR.sum <- rbind(betaR.sum, betaR_add)
    }
    
    
    # Summarize average population densities
    Dens_names <- c(paste0("densAvg_tot[",  i, "]"), paste0("densAvg_shared[",  i, "]"))
    Dens_tot_add <- data.frame(Area = area_names[i],
                               SummaryPeriod = paste0(minYear, "-", maxYear),
                               Median = median(out.mat[, Dens_names[1]]),
                               lCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.025)),
                               uCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.975)),
                               Mean = mean(out.mat[, Dens_names[1]]),
                               SD = sd(out.mat[, Dens_names[1]])
    )
    Dens_shared_add <- data.frame(Area = area_names[i],
                                  SummaryPeriod = paste0(minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1),
                                  Median = median(out.mat[, Dens_names[2]]),
                                  lCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.025)),
                                  uCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.975)),
                                  Mean = mean(out.mat[, Dens_names[2]]),
                                  SD = sd(out.mat[, Dens_names[2]])
    )
    popDens.sum <- rbind(popDens.sum, Dens_tot_add, Dens_shared_add)
    
    # Summarize average population growth rates
    lambda_names <- c(paste0("lambdaAvg_tot[",  i, "]"), paste0("lambdaAvg_shared[",  i, "]"))
    lambda_tot_add <- data.frame(Area = area_names[i],
                                 SummaryPeriod = paste0(minYear, "-", maxYear),
                                 Median = median(out.mat[, lambda_names[1]]),
                                 lCI = unname(quantile(out.mat[, lambda_names[1]], probs = 0.025)),
                                 uCI = unname(quantile(out.mat[, lambda_names[1]], probs = 0.975)),
                                 Mean = mean(out.mat[, lambda_names[1]]),
                                 SD = sd(out.mat[, lambda_names[1]])
    )
    lambda_shared_add <- data.frame(Area = area_names[i],
                                    SummaryPeriod = paste0(minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1),
                                    Median = median(out.mat[, lambda_names[2]]),
                                    lCI = unname(quantile(out.mat[, lambda_names[2]], probs = 0.025)),
                                    uCI = unname(quantile(out.mat[, lambda_names[2]], probs = 0.975)),
                                    Mean = mean(out.mat[, lambda_names[2]]),
                                    SD = sd(out.mat[, lambda_names[2]])
    )
    lambda.sum <- rbind(lambda.sum, lambda_tot_add, lambda_shared_add)
  } 
  
  ## Add coordinate information to summarized posterior data
  area_coord <- area_coord %>%
    dplyr::rename(Area = spatialUnit)
  
  # Average reproductive rates
  rRep.sum <- rRep.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  # Average annual survival 
  pSurv.sum <- pSurv.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  # Rodent effects 
  if(fitRodentCov){
    betaR.sum <- betaR.sum %>%
      dplyr::left_join(area_coord, by = "Area")
  }
  
  # Average population densities
  popDens.sum <- popDens.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  # Average population growth rates
  lambda.sum <- lambda.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  
  ## Plot estimates as a function of longitude (and latitude)
  ggplot(rRep.sum) + 
    geom_pointrange(aes(x = latitudeAvg, y = Median, ymin = lCI, ymax = uCI, colour = longitudeAvg), size = 0.3, fatten = 4) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  
}