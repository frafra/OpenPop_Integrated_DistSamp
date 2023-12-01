#' Compare posterior densities of major parameters for two or more models
#'
#' @param modelPaths vector of strings. Paths to .rds files for each models' posterior samples
#' @param modelChars vector of strings. Characteristics of each model (identifier)
#' @param N_sites integer. Number of sites to plot comparison for
#' @param N_years integer. Number of years to plot comparison for
#' @param plotPath string. Directory into which to save plots
#' @param returnData logical. If TRUE, return collated data from all models (default = FALSE)
#' 
#' @return
#' @export
#'
#' @examples

plotModelComparison <- function(modelPaths, modelChars, N_sites, N_years, plotPath, returnData = FALSE){
  
  require(coda)
  require(ggplot2)
  require(viridis)
  
  if(length(modelPaths) != length(modelChars)){
    stop("Unequal number of model paths and model characteristics provided.")
  }
  
  ## Set number of models to compare
  nMod <- length(modelPaths)
  
  ## Collate posterior data in a list
  data.list <- list()
  
  for(n in 1:nMod){
    out.sam <- as.matrix(readRDS(modelPaths[n]))
    
    # Calculate average density
    for(t in 1:N_years){
      popDens_juv <- out.sam[, paste0("Density[1, 1, ", 1:N_sites, ", ", t, "]")]
      popDens_ad <- out.sam[, paste0("Density[1, 2, ", 1:N_sites, ", ", t, "]")]
      
      popDens_mean <- rowMeans(popDens_juv + popDens_ad)
      rm(popDens_juv, popDens_ad)
      out.sam <- cbind(out.sam, popDens_mean)
      colnames(out.sam)[ncol(out.sam)] <- paste0("AvgDens[", t, "]")
    }
    
    out.data <- reshape2::melt(out.sam)
    rm(out.sam)
    colnames(out.data) <- c('Sample', 'Parameter', 'Value')
    out.data$Model <- modelChars[n]
    data.list[[n]] <- out.data
    rm(out.data)
  }
  
  ## Combine model data into single data frame
  data.all <- dplyr::bind_rows(data.list)
  rm(data.list)
  
  ## Make plotting directory if necessary
  dir.create(plotPath, showWarnings = FALSE)
  
  ## Plot overlapping posterior densities for different parameter groups
  
  # "Main" parameters
  mains <- c('Mu.D1[1]', 'sigma.D[1]', 'ratio.JA1[1]',
             'Mu.R[1]', 'sigmaT.R',
             'Mu.S[1]', 'sigmaT.S[1]', 'Mu.S1', 'epsT.S1.prop[1]')
  
  pdf(paste0(plotPath, '/ModelComp_Mains.pdf'), width = 8, height = 5)
  print(
    ggplot(subset(data.all, Parameter %in% mains)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Detection parameters
  detects <- c('mu.dd[1]', 'sigmaT.dd', paste0('esw[1, ', 1:N_years, ']'), paste0('p[1, ', 1:N_years, ']'))
  
  pdf(paste0(plotPath, '/ModelComp_Detects.pdf'), width = 8, height = 5)
  print(
    ggplot(subset(data.all, Parameter %in% detects)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  
  # Annual Recruitment
  R_year <- paste0('R_year[1, ', 1:N_years, ']')
  
  pdf(paste0(plotPath, '/ModelComp_R_year.pdf'), width = 8, height = 5)
  print(
    ggplot(subset(data.all, Parameter %in% R_year)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  
  # Site- and age-specific population sizes
  pdf(paste0(plotPath, '/ModelComp_N_a_exp.pdf'), width = 8, height = 5)
  for(j in 1:N_sites){
    print(
      ggplot(subset(data.all, Parameter %in% c(paste0('N_exp[1, 1, ', j, ', ', 1:N_years, ']'), 
                                               paste0('N_exp[1, 2, ', j, ', ', 1:N_years, ']')) &
                      Value != 0)) + 
        geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
        facet_wrap(~Parameter, scales = 'free') +
        scale_fill_viridis(discrete = T) + 
        scale_color_viridis(discrete = T) + 
        ggtitle(paste0('Site ', j)) + 
        theme_bw() + theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  
  # Site- and age-specific densities
  pdf(paste0(plotPath, '/ModelComp_Density_a.pdf'), width = 8, height = 5)
  for(j in 1:N_sites){
    print(
      ggplot(subset(data.all, Parameter %in% c(paste0('Density[1, 1, ', j, ', ', 1:N_years, ']'), 
                                                paste0('Density[1, 2, ', j, ', ', 1:N_years, ']')) &
                      Value != 0)) + 
        geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
        facet_wrap(~Parameter, scales = 'free') +
        scale_fill_viridis(discrete = T) + 
        scale_color_viridis(discrete = T) + 
        ggtitle(paste0('Site ', j)) + 
        theme_bw() + theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  
  # Average densities
  data.dens.sum <- data.all %>%
    dplyr::filter(Parameter %in% paste0("AvgDens[", 1:N_years, "]")) %>%
    dplyr::mutate(YearIdx = readr::parse_number(as.character(Parameter))) %>%
    dplyr::group_by(YearIdx, Model) %>%
    dplyr::summarise(median = median(Value),
                     lCI = quantile(Value, 0.025),
                     uCI = quantile(Value, 0.975),
                     .groups = "keep")
  
  pdf(paste0(plotPath, '/ModelComp_AvgDensity.pdf'), width = 8, height = 5)
  print(
    ggplot(data.dens.sum) + 
      geom_line(aes(x = YearIdx, y = median, color = Model)) + 
      geom_ribbon(aes(x = YearIdx, ymin = lCI, ymax = uCI, fill = Model), alpha = 0.5) + 
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) +
      scale_x_continuous(breaks = 1:max(data.dens.sum$YearIdx)) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  
  ## Return all data
  if(returnData){
    return(data.all)
  }else{
    return()
  }

}
