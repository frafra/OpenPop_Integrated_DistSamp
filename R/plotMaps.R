
plotMaps <- function(mcmc.out, mapNM,
                     N_areas, area_names, N_sites, 
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
  

  ## Add estimates to map objects
  
  # Average reproductive rates
  mapNM.rRep <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., rRep.sum, by = "Area") 
  
  mapNM.rRep@data <- as.data.frame(matchData)
  
  # Average annual survival 
  mapNM.pSurv <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., pSurv.sum, by = "Area") 
  
  mapNM.pSurv@data <- as.data.frame(matchData)
  
  # Rodent effects 
  if(fitRodentCov){
    mapNM.betaR <- mapNM
    matchData <- as_tibble(mapNM@data) %>%
      dplyr::left_join(., betaR.sum, by = "Area") 
    
    mapNM.betaR@data <- as.data.frame(matchData)
  }
  
  # Average population densities
  mapNM.popDens1 <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., subset(popDens.sum, SummaryPeriod == paste0(minYear, "-", maxYear)), by = "Area") 
  
  mapNM.popDens1@data <- as.data.frame(matchData)
  
  mapNM.popDens2 <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., subset(popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear)), by = "Area") 
  
  mapNM.popDens2@data <- as.data.frame(matchData)
  
  # Average population growth rates
  mapNM.lambda1 <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., subset(lambda.sum, SummaryPeriod == paste0(minYear, "-", maxYear)), by = "Area") 
  
  mapNM.lambda1@data <- as.data.frame(matchData)
  
  mapNM.lambda2 <- mapNM
  matchData <- as_tibble(mapNM@data) %>%
    dplyr::left_join(., subset(lambda.sum, SummaryPeriod != paste0(minYear, "-", maxYear)), by = "Area") 
  
  mapNM.lambda2@data <- as.data.frame(matchData)
  
  ## Plot maps and print to pdf
  
  # Average reproductive rates
  pdf("Plots/AreaMaps/Avg_rRep_Map.pdf", width = 5, height = 6)
  
  print(
    tmap::tm_shape(mapNM.rRep) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80")
  )
  
  print(
    tmap::tm_shape(mapNM.rRep) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  
  dev.off()
  
  # Average annual survival 
  pdf("Plots/AreaMaps/Avg_pSurv_Map.pdf", width = 5, height = 6)
  
  print(
    tmap::tm_shape(mapNM.pSurv) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80")
  )
  
  print(
    tmap::tm_shape(mapNM.pSurv) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  
  dev.off()
  
  # Rodent effects
  if(fitRodentCov){
    pdf("Plots/AreaMaps/betaR_Map.pdf", width = 5, height = 6)
    
    print(
      tmap::tm_shape(mapNM.betaR) + tmap::tm_polygons("Median", palette = colorspace::divergingx_hcl(10, palette = "Zissou1"), style = "cont", colorNA = "grey80")
    )
    
    print(
      tmap::tm_shape(mapNM.betaR) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
    )
    
    dev.off()
  }
 
  
  # Average population densities
  pdf("Plots/AreaMaps/Avg_popDens_Map.pdf", width = 5, height = 6)
  
  print(
    tmap::tm_shape(mapNM.popDens1) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80") #+
    #tmap::tm_layout(title = paste0("Average population density (", minYear, "-", maxYear, "), Median"))
  )
  
  print(
    tmap::tm_shape(mapNM.popDens1) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80") #+ 
    #tmap::tm_layout(title = paste0("Average population density (", minYear, "-", maxYear, "), SD"))
  )
  
  print(
    tmap::tm_shape(mapNM.popDens2) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80") #+
    #tmap::tm_layout(title = paste0("Average population density (", minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1, "), Median"))
  )
  
  print(
    tmap::tm_shape(mapNM.popDens2) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80") #+ 
    #tmap::tm_layout(title = paste0("Average population density (", minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1, "), SD"))
  )
  
  dev.off()
  
  
  # Average population growth rates
  pdf("Plots/AreaMaps/Avg_lambda_Map.pdf", width = 5, height = 6)
  
  print(
    tmap::tm_shape(mapNM.lambda1) + tmap::tm_polygons("Median", palette = colorspace::divergingx_hcl(10, palette = "PiYG"), midpoint = 1, style = "cont", colorNA = "grey80") #+
    #tmap::tm_layout(title = paste0("Average population growth rate (", minYear, "-", maxYear, "), Median"))
  )
  
  print(
    tmap::tm_shape(mapNM.lambda1) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80") #+ 
    #tmap::tm_layout(title = paste0("Average population growth rate (", minYear, "-", maxYear, "), SD"))
  )
  
  print(
    tmap::tm_shape(mapNM.lambda2) + tmap::tm_polygons("Median", palette = colorspace::divergingx_hcl(10, palette = "PiYG"), midpoint = 1, style = "cont", colorNA = "grey80") #+
    #tmap::tm_layout(title = paste0("Average population growth rate (", minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1, "), Median"))
  )
  
  print(
    tmap::tm_shape(mapNM.lambda2) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80") #+ 
    #tmap::tm_layout(title = paste0("Average population growth rate (", minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1, "), SD"))
  )
  
  dev.off()
  
  ## Write and return plot paths
  plot.paths <- c("Plots/AreaMaps/Avg_rRep_Map.pdf", 
                  "Plots/AreaMaps/Avg_pSurv_Map.pdf",
                  "Plots/AreaMaps/Avg_popDens_Map.pdf",
                  "Plots/AreaMaps/Avg_lambda_Map.pdf")
  
  if(fitRodentCov){
    plot.paths <- c(plot.paths, "Plots/AreaMaps/betaR_Map.pdf")
  }
  
  return(plot.paths)
}