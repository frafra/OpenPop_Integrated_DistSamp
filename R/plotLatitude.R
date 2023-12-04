#' Plot posterior summaries for vital rate and population parameters as a function of latitude
#'
#' @param PostSum.list a list containing 5 dataframes with summarised posteriors for several
#' parameters by area: average recruitment (rRep.sum), average survival 
#' (pSurv.sum), rodent effect slope (betaR.sum, = NA when fitRodentCov = FALSE), 
#' average population density (popDens.sum), and population growth rate 
#' (lambda.sum). 
#' @param area_coord a tibble containing coordinate data for each area in analysis.
#' @param minYear integer. First year in the analysis.
#' @param maxYear integer. Last year in the analysis.
#' @param fitRodentCov logical. If TRUE, makes plot for rodent effect slope in 
#' addition to other plots. 
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/Latitude.
#' @export
#'
#' @examples

plotLatitude <- function(PostSum.list, 
                         area_coord,
                         minYear, maxYear,
                         fitRodentCov){
  
  ## Add coordinate information to summarized posterior data and plot
  area_coord <- area_coord %>%
    dplyr::rename(Area = spatialUnit,
                  Longitude = longitudeAvg,
                  Latitude = latitudeAvg)
  
  # Average reproductive rates
  rRep.sum <- PostSum.list$rRep.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_rRep <- ggplot(rRep.sum) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    ylab("Estimate") + 
    theme_classic()
  
  # Average annual survival 
  pSurv.sum <- PostSum.list$pSurv.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_pSurv <- ggplot(pSurv.sum) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    ylab("Estimate") + 
    theme_classic()
  
  # Rodent effects 
  if(fitRodentCov){
    betaR.sum <- PostSum.list$betaR.sum %>%
      dplyr::left_join(area_coord, by = "Area")
    
    p_betaR <- ggplot(betaR.sum) + 
      geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
      ylab("Estimate") + 
      paletteer::scale_color_paletteer_c("grDevices::Temps") + 
      theme_classic()
  }
  
  # Average population densities
  popDens.sum <- PostSum.list$popDens.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_popDens1 <- ggplot(subset(popDens.sum, SummaryPeriod == paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  p_popDens2 <- ggplot(subset(popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  # Average population growth rates
  lambda.sum <- PostSum.list$lambda.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_lambda1 <- ggplot(subset(lambda.sum, SummaryPeriod == paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  p_lambda2 <- ggplot(subset(lambda.sum, SummaryPeriod != paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0.5, fatten = 4, alpha = 0.5) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  
  ## Make plotting directory if it does not exist yet
  ifelse(!dir.exists("Plots/Latitude"), dir.create("Plots/Latitude"), FALSE)
  
  ## Assemble plots (vital rate parameters)
  if(fitRodentCov){
    pdf("Plots/Latitude/VitalRateParams_Latitude.pdf", width = 8, height = 6)
    print(
      ggpubr::ggarrange(p_rRep + ggtitle("Average recruitment rate") + theme(axis.title.x = element_blank()), 
                        p_pSurv + ggtitle("Average survival rate") + theme(axis.title.x = element_blank()),  
                        p_betaR + ggtitle("Rodent effect on recruitment"),
                        heights = c(1, 1, 1.2), ncol = 1, 
                        legend = "right", common.legend = TRUE)
    )
    dev.off()
  }else{
    pdf("Plots/Latitude/VitalRateParams_Latitude.pdf", width = 8, height = 4)
    print(
      ggpubr::ggarrange(p_rRep + ggtitle("Average recruitment rate") + theme(axis.title.x = element_blank()), 
                        p_pSurv + ggtitle("Average survival rate"),  
                        heights = c(1, 1.2), ncol = 1, 
                        legend = "right", common.legend = TRUE)
    )
    dev.off()
  }

  
  
  ## Assemble plots (population-level parameters)
  pdf("Plots/Latitude/PopParams1_Latitude.pdf", width = 8, height = 4)
  print(
    ggpubr::ggarrange(p_popDens1 + ggtitle(paste0("Average population density (", minYear, "-", maxYear, ")")) + theme(axis.title.x = element_blank()),
                      p_lambda1 + ggtitle(paste0("Average population growth rate (", minYear, "-", maxYear, ")")),
                      heights = c(1, 1.1), ncol = 1, 
                      legend = "right", common.legend = TRUE)
  )
  dev.off()
  
  
  shared_sumPeriod <- subset(popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear))$SummaryPeriod[1]
  
  pdf("Plots/Latitude/PopParams2_Latitude.pdf", width = 8, height = 4)
  print(ggpubr::ggarrange(p_popDens2 + ggtitle(paste0("Average population density (", shared_sumPeriod, ")")) + theme(axis.title.x = element_blank()),
                          p_lambda2 + ggtitle(paste0("Average population growth rate (", shared_sumPeriod, ")")),
                          heights = c(1, 1.1), ncol = 1, 
                          legend = "right", common.legend = TRUE)
  )
  dev.off()
  
  
  ## Write and return plot paths
  plot.paths <- c("Plots/Latitude/VitalRateParams_Latitude.pdf", 
                  "Plots/Latitude/PopParams1_Latitude.pdf",
                  "Plots/Latitude/PopParams2_Latitude.pdf")
  
  return(plot.paths)
  
}