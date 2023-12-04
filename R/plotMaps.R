#' Plot posterior summaries for vital rate and population parameters on map
#'
#' @param PostSum.list a list containing 5 dataframes with summarised posteriors for several
#' parameters by area: average recruitment (rRep.sum), average survival 
#' (pSurv.sum), rodent effect slope (betaR.sum, = NA when fitRodentCov = FALSE), 
#' average population density (popDens.sum), and population growth rate 
#' (lambda.sum). 
#' @param mapNM sf / dataframe object containing map of Norwegian municipalities. 
#' @param minYear integer. First year in the analysis.
#' @param maxYear integer. Last year in the analysis.
#' @param fitRodentCov logical. If TRUE, makes plot for rodent effect slope in 
#' addition to other plots. 
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/AreaMaps.
#' @export
#'
#' @examples

plotMaps <- function(PostSum.list, mapNM,
                     minYear, maxYear,
                     fitRodentCov){
  
  ## Add estimates to map objects
  
  # Average reproductive rates
  mapNM.rRep <- mapNM
  mapNM.rRep <- mapNM.rRep %>%
    dplyr::left_join(., PostSum.list$rRep.sum, by = "Area") 
  
  # Average annual survival 
  mapNM.pSurv <- mapNM
  mapNM.pSurv <- mapNM.pSurv %>%
    dplyr::left_join(., PostSum.list$pSurv.sum, by = "Area") 
  
  # Rodent effects 
  if(fitRodentCov){
    mapNM.betaR <- mapNM
    mapNM.betaR <- mapNM.betaR %>%
      dplyr::left_join(., PostSum.list$betaR.sum, by = "Area") 
  }
  
  # Average population densities
  mapNM.popDens1 <- mapNM
  mapNM.popDens1 <- mapNM.popDens1 %>%
    dplyr::left_join(., subset(PostSum.list$popDens.sum, SummaryPeriod == paste0(minYear, "-", maxYear)), by = "Area") 
  
  mapNM.popDens2 <- mapNM
  mapNM.popDens2 <- mapNM.popDens2 %>%
    dplyr::left_join(., subset(PostSum.list$popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear)), by = "Area") 

  # Average population growth rates
  mapNM.lambda1 <- mapNM
  mapNM.lambda1 <- mapNM.lambda1 %>%
    dplyr::left_join(., subset(PostSum.list$lambda.sum, SummaryPeriod == paste0(minYear, "-", maxYear)), by = "Area") 

  mapNM.lambda2 <- mapNM
  mapNM.lambda2 <- mapNM.lambda2 %>%
    dplyr::left_join(., subset(PostSum.list$lambda.sum, SummaryPeriod != paste0(minYear, "-", maxYear)), by = "Area") 
  
  
  ## Make plotting directory if it does not exist yet
  ifelse(!dir.exists("Plots/AreaMaps"), dir.create("Plots/AreaMaps"), FALSE)
  
  
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
