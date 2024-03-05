#' Visualize comparison of model estimates with true parameter values from data simulation
#'
#' @param PlotColors string specifying color palette to use for plots. 
#' Currently supports "customRainbow" (default), "Temps", and "Zissou1".
#' @param thin integer. Thinning rate to use on samples for plotting. 
#' @return a vector of png plot names. The plots can be found in Plots/TimeSeries.
#' 
#' @export 
#'
#' @examples


plotSimCheck_replicates <- function(plotColors = "customRainbow", thin = 1) {
  
  require(coda)
  require(tidyverse)
  require(tidybayes)
  require(ggforce)
  require(see)
  require(cowplot)
  
  if(!(plotColors %in% c("customRainbow", "Temps", "Zissou1"))){
    stop("Invalid plotColors. Currently supported are customRainbow, Temps, and Zissou1.")
  }
  
  # Data aggregation #
  #------------------#
  
  ## Load information on simulation and run seeds
  simSeeds <- readRDS("simData/seedList.rds")
  runSeeds <- readRDS("simModelFits/seedInfo.rds")
  
  ## Set up tibbles to collate simulated data and model posteriors
  N_simData <- D_simData <-  R_simData <- det_simData <- simParams <- tibble()
  R_year <- Mu_R <- sigmaT_R <- Mu_S <- tibble()
  Mu_dd <- sigmaT_dd <- esw_year <- p_year <- tibble()
  N_tot <- Density_year <- tibble()
  
  for(i in 1:length(simSeeds)){
    
    ## Read in simulated dataset
    SimData <- readRDS(paste0("simData/AllSimData_seed", simSeeds[i], ".rds"))
    
    ## Extract data on baseline parameters
    simParams_temp <- tibble(Mu_R = SimData$SimParams$Mu.R,
                             sigmaT_R = SimData$SimParams$sigmaT.R,
                             Mu_S = SimData$SimParams$Mu.S,
                             Mu_S1 = sqrt(SimData$SimParams$Mu.S),
                             Mu_S2 = sqrt(SimData$SimParams$Mu.S),
                             Mu_dd = SimData$SimParams$Mu.dd,
                             sigmaT_dd = SimData$SimParams$sigmaT.dd,
                             dataSetID = as.factor(i),
                             simSeed = as.factor(simSeeds[i]))
    simParams <- rbind(simParams, simParams_temp)
    
    ## Extract data on annual recruitment
    Na_temp <- apply(SimData$N.data, c(2,3), sum) 
    R_simData_temp <- as_tibble(t(Na_temp)) %>%
      dplyr::mutate(realizedR = V1 / V2, year = seq(1:SimData$SimParams$Tmax)) %>%
      dplyr::select(year, realizedR) %>%
      dplyr::mutate(predictedR = colMeans(SimData$VR.list$R),
                    dataSetID = as.factor(i),
                    simSeed = as.factor(simSeeds[i]))
    R_simData <- rbind(R_simData, R_simData_temp)
    
    ## Extract data on population size
    N_data_temp <- tibble(year = seq(1:SimData$SimParams$Tmax), 
                          N_tot = apply(SimData$N.data, 3, sum),
                          N_juv = colSums(SimData$N.data[,1,]),
                          N_ad = colSums(SimData$N.data[,2,]),
                          dataSetID = as.factor(i),
                          simSeed = as.factor(simSeeds[i])) 
    N_simData <- rbind(N_simData, N_data_temp)
    
    ## Extract data on density
    A_temp <- apply(SimData$DS.data$L, 2, sum) * SimData$SimParams$W*2 / (1000 *1000)
    D_temp <- tibble(year = seq(1:SimData$SimParams$Tmax), 
                     Mean.D = N_data_temp$N_tot / A_temp,
                     dataSetID = as.factor(i),
                     simSeed = as.factor(simSeeds[i])) 
    D_simData <- rbind(D_simData, D_temp)
    
    ## Extract data on detection
    sigma_temp <- apply(SimData$DS.data$sigma, 2, mean) # Works only if sigma is the same for all lines
    det_simData_temp <- tibble(year = seq(1:SimData$SimParams$Tmax), 
                               sigma = sigma_temp,
                               esw = sqrt(3.141593*(sigma_temp^2)/2),
                               p = ifelse(sqrt(3.141593*(sigma_temp^2)/2) > SimData$SimParams$W, 1, sqrt(3.141593*(sigma_temp^2)/2)/SimData$SimParams$W),
                               dataSetID = as.factor(i),
                               simSeed = as.factor(simSeeds[i])) 
    det_simData <- rbind(det_simData, det_simData_temp)
    
    
    for(k in 1:length(runSeeds[[i]])){
      
      ## Read in posterior samples
      postSam <- readRDS(paste0("simModelFits_sum/IDSMsampleSum_simSeed", simSeeds[i], "_runSeed", runSeeds[[i]][k], ".rds"))$sum.post
      
      ## List seed information
      seedInfo <- tibble(dataSetID = as.factor(i),
                         simSeed = as.factor(simSeeds[i]),
                         runID = as.factor(k),
                         runSeed = as.factor(runSeeds[[i]][k]))
      
      ## Set indices for thinning samples
      N_samples <- nrow(postSam$Mu_R)
      thin.idx <- seq(1, N_samples, by = thin)
      
      N_samples_t <- nrow(postSam$R_year)
      thin.idx_t <- seq(1, N_samples_t, by = thin)
      
      ## Append data
      R_year <- rbind(R_year, cbind(postSam$R_year[thin.idx_t,], seedInfo))
      Mu_R <- rbind(Mu_R, cbind(postSam$Mu_R[thin.idx,], seedInfo))
      sigmaT_R <- rbind(sigmaT_R, cbind(postSam$sigmaT_R[thin.idx,], seedInfo))
      Mu_S <- rbind(Mu_S, cbind(postSam$Mu_S_data[seq(1, nrow(postSam$Mu_S_data), by = thin),], seedInfo))
      
      Mu_dd_temp <- postSam$Mu_dd %>%
        dplyr::mutate(mu.dd = exp(mu.dd),
                      lab_code = "Mu.dd") %>%
        dplyr::rename(Mu.dd = mu.dd)
      
      Mu_dd <- rbind(Mu_dd, cbind(Mu_dd_temp[thin.idx,], seedInfo))
      sigmaT_dd <- rbind(sigmaT_dd, cbind(postSam$sigmaT_dd[thin.idx,], seedInfo))
      esw_year <- rbind(esw_year, cbind(postSam$esw_year[thin.idx_t,], seedInfo))
      p_year <- rbind(p_year, cbind(postSam$p_year[thin.idx_t,], seedInfo))
      
      N_tot <- rbind(N_tot, cbind(postSam$N_tot[thin.idx_t,], seedInfo))
      Density_year <- rbind(Density_year, cbind(postSam$Density_year[thin.idx_t,], seedInfo))
    }
  }
  
  
  # Set up for plotting #
  #---------------------#
  
  ## Plotting directory
  if(!dir.exists("Plots")){
    dir.create("Plots")
    dir.create("Plots/SimCheck_replicates")
  }
  
  if(!dir.exists("Plots/SimCheck_replicates")){
    dir.create("Plots/SimCheck_replicates")
  }
  
  ## List of plot names (required for targets integration)
  plot.paths <- c()
  
  ## Custom color palette
  source("ColorPalettes_Custom.R")
  
  ## Set plot colors
  if(plotColors == "customRainbow"){
    plot.cols <- custom_palettes("darkRainbow", n = length(simSeeds), type = "continuous")
  }
  if(plotColors == "Temps"){
    plot.cols <- paletteer::paletteer_c("grDevices::Temps", length(simSeeds))
  }
  if(plotColors == "Zissou1"){
    plot.cols <- hcl.colors(length(simSeeds), palette = "Zissou1")
  }
  
  # Plotting vital rates (all replicates) #
  #---------------------------------------#
  
  
  ## Plot average survival probabilities (Mu.S)
  p1 <- ggplot(data = Mu_S, aes(x = Surv, y = S)) +
    geom_violinhalf(aes(color = simSeed, linetype = runID), fill = NA, position = "identity", scale = "width")  +
    geom_point(aes(x = "S", y = simParams$Mu_S[1]), 
               col = "black", size = 3) +
    geom_point(aes(x = "S1", y = simParams$Mu_S1[1]), 
               col = "black", size = 3) + 
    geom_point(aes(x = "S2", y = simParams$Mu_S2[1]), 
               col = "black", size = 3) + 
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(Mu_S$runID))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3",  linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("Survival probability") +
    xlab("Survival component")
  
  
  ## Plot average recruitment rate (Mu.R)
  p2 <- ggplot(data = Mu_R, aes(x = lab_code, y = Mu.R)) +
    geom_violinhalf(aes(linetype = runSeed, color = simSeed), fill = NA, position = "identity")  +
    geom_point(aes(x = "Mu.R", y = simParams$Mu_R[1]), 
               col = "black", size = 3) +
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(Mu_R$runSeed))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("Recruitment rate") +
    xlab("")
  
  ## Plot temporal variation in recruitment rate (sigma.R)
  p3 <- ggplot(data = sigmaT_R, aes(x = lab_code, y = sigmaT.R)) +
    geom_violinhalf(aes(linetype = runSeed, color = simSeed), fill = NA, position = "identity")  +
    geom_point(aes(x = "sigmaT.R", y = simParams$sigmaT_R[1]), 
               col = "black", size = 3) +
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(sigmaT_R$runSeed))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("Year RE SD") +
    xlab("")
  
  ## Plot annual recruitment rates (R_year)
  p4 <- ggplot(data = R_year, aes(x = as.factor(year), y = R_year)) +
    geom_violinhalf(aes(fill = simSeed, linetype = runID), alpha = 1/3, color = NA, position = "identity", scale = "width")  +
    geom_point(data = R_simData, aes(x = as.factor(year), y = predictedR, fill = simSeed), col = "black", size = 2, shape = 21) +
    #geom_point(data = R_simData, aes(x = as.factor(year), y = realizedR, fill = simSeed), col = "grey40", size = 2, shape = 23) +
    scale_fill_manual(values = plot.cols) + 
    scale_color_manual(values = plot.cols) + 
    scale_linetype_manual(values = rep("solid", nlevels(R_year$runID))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("R") +
    xlab("Year")
  
  ## Putting the plots together in a multipanel plot
  p_a <- plot_grid(p1, p2, p3, nrow = 1)
  p_out <- plot_grid(p_a, p4, nrow = 2, label_size = 10)
  
  ## Plot to pdf
  pdf(paste0("Plots/SimCheck_replicates/SimCheck_VRs_", plotColors, ".pdf"), width = 12, height = 7.5)
  suppressWarnings(
    print(p_out)
  )
  dev.off()
  
  ## Plot to png
  png(paste0("Plots/SimCheck_replicates/SimCheck_VRs_", plotColors, ".png"), width = 12, height = 7.5, units = "in", res = 300)
  suppressWarnings(
    print(p_out)
  )
  dev.off()
  
  
  # Plotting detection parameters (all replicates) #
  #------------------------------------------------#
  
  ## Mean half normal scale parameter (mu.dd)
  p1 <- ggplot(data = Mu_dd, aes(x = lab_code, y = Mu.dd)) +
    geom_violinhalf(aes(linetype = runID, color = simSeed), fill = NA, position = "identity")  +
    geom_point(aes(x = "Mu.dd", y = simParams$Mu_dd[1]), 
               col = "black", size = 3) +
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(Mu_dd$runID))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("HN scale parameter") +
    xlab("")
  
  ## Temporal variation in scale parameter (sigmaT.dd) - i.e. sd of random var.
  p2 <- ggplot(data = sigmaT_dd, aes(x = lab_code, y = sigmaT.dd)) +
    geom_violinhalf(aes(linetype = runID, color = simSeed), fill = NA, position = "identity")  +
    geom_point(aes(x = "sigmaT.dd", y = simParams$sigmaT_dd[1]), 
               col = "black", size = 3) +
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(sigmaT_dd$runID))) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("sigmaT.dd") +
    xlab("")
  
  ## Effective strip width across years
  p4 <- ggplot(data = esw_year, aes(x = as.factor(year), y = esw)) +
    geom_violinhalf(aes(fill = simSeed, linetype = runID), alpha = 1/3, color = NA, position = "identity", scale = "width")  +
    scale_color_manual(values = plot.cols) + 
    scale_fill_manual(values = plot.cols) + 
    scale_linetype_manual(values = rep("solid", nlevels(esw_year$runID))) + 
    geom_point(data = det_simData, aes(x = as.factor(year), y = esw, fill = simSeed), col = "black", size = 2, shape = 21) +
    ylim(min(esw_year$esw), 300) + 
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("Eff. strip width") +
    xlab("Year")
  
  ## Detection probability across years
  p5 <- ggplot(data = p_year, aes(x = as.factor(year), y = p)) +
    geom_violinhalf(aes(fill = simSeed, linetype = runID), alpha = 1/3, color = NA, position = "identity", scale = "width")  +
    scale_color_manual(values = plot.cols) + 
    scale_fill_manual(values = plot.cols) + 
    scale_linetype_manual(values = rep("solid", nlevels(p_year$runID))) +
    geom_point(data = det_simData, aes(x = as.factor(year), y = p, fill = simSeed), col = "black", size = 2, shape = 21) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          text = element_text(size = 10), legend.position = "none") +
    ylab("Detection probability") +
    xlab("Year")
  
  ## Assemble plots
  p3 <- plot_grid(p1, p2, nrow = 1)
  p_out <- plot_grid(p3, p4, p5, nrow = 3)
  
  ## Plot to pdf
  pdf(paste0("Plots/SimCheck_replicates/SimCheck_Detects_", plotColors, ".pdf"), width = 10, height = 7.5)
  suppressWarnings(
    print(p_out)
  )
  dev.off()
  
  ## Plot to png
  png(paste0("Plots/SimCheck_replicates/SimCheck_Detects_", plotColors, ".png"), width = 10, height = 7.5, units = "in", res = 300)
  suppressWarnings(
    print(p_out)
  )
  dev.off()
  
  
  for(i in 1:length(simSeeds)){
    
    ## Plot annual population size 
    N_tot_sub <- subset(N_tot, dataSetID == i)
    N_simData_sub <- subset(N_simData, dataSetID == i)
    
    pN <- ggplot(data = N_tot_sub, aes(x = as.factor(year), y = N_tot_exp)) +
      geom_violinhalf(aes(linetype = runID), color = NA, fill = plot.cols[i], alpha = 0.25, position = "identity", scale = "width")  +
      geom_point(data = N_simData_sub, aes(x = as.factor(year), y = N_tot), color = "black", size = 2) +
      scale_linetype_manual(values = rep("solid", nlevels(N_tot_sub$runID))) + 
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
            text = element_text(size = 10), legend.position = "none") +
      ylab("Population size") +
      xlab("Year")
    
    ## Plot annual population densities 
    Density_year_sub <- subset(Density_year, dataSetID == i)
    D_simData_sub <- subset(D_simData, dataSetID == i)
    
    pD <- ggplot(data = Density_year_sub , aes(x = as.factor(year), y = density)) +
      geom_violinhalf(aes(linetype = runID), color = NA, fill = plot.cols[i], alpha = 0.25, position = "identity", scale = "width")  +
      geom_point(data = D_simData_sub, aes(x = as.factor(year), y = Mean.D), color = "black", size = 2) +
      scale_linetype_manual(values = rep("solid", nlevels(Density_year_sub$runID))) + 
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
            text = element_text(size = 10), legend.position = "none") +
      ylab("Average population density") +
      xlab("Year")
    
    ## Plot annual recruitment rates 
    R_year_sub <- subset(R_year, dataSetID == i)
    R_simData_sub <- subset(R_simData, dataSetID == i)
    
    
    pR <- ggplot(data = R_year_sub, aes(x = as.factor(year), y = R_year)) +
      geom_violinhalf(aes(linetype = runID), color = NA, fill = plot.cols[i], alpha = 0.25, position = "identity", scale = "width")  +
      geom_point(data = R_simData_sub, aes(x = as.factor(year), y = predictedR), color = "black", size = 2) +
      geom_point(data = R_simData_sub, aes(x = as.factor(year), y = realizedR), color = "grey20", size = 3, alpha = 0.25) +
      scale_linetype_manual(values = rep("solid", nlevels(Mu_S$runID))) + 
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
            text = element_text(size = 10), legend.position = "none") +
      ylab("R") +
      xlab("Year")
    
    ## Plot effective strip width across years
    esw_year_sub <- subset(esw_year, dataSetID == i)
    det_simData_sub <- subset(det_simData, dataSetID == i)
    
    pE <- ggplot(data = esw_year_sub, aes(x = as.factor(year), y = esw)) +
      geom_violinhalf(aes(linetype = runID), color = NA, fill = plot.cols[i], alpha = 0.25, position = "identity", scale = "width")  +
      scale_color_manual(values = plot.cols) + 
      scale_linetype_manual(values = rep("solid", nlevels(esw_year$runID))) + 
      geom_point(data = det_simData_sub, aes(x = as.factor(year), y = esw), color = "black", size = 2) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
            text = element_text(size = 10), legend.position = "none") +
      ylab("Eff. strip width") +
      xlab("Year")
    
    ## Plot detection probability across years
    p_year_sub <- subset(p_year, dataSetID == i)
    
    pP <- ggplot(data = p_year_sub, aes(x = as.factor(year), y = p)) +
      geom_violinhalf(aes(linetype = runID), color = NA, fill = plot.cols[i], alpha = 0.25, position = "identity", scale = "width")  +
      scale_color_manual(values = plot.cols) + 
      scale_linetype_manual(values = rep("solid", nlevels(p_year$runID))) +
      geom_point(data = det_simData_sub, aes(x = as.factor(year), y = p), color = "black", size = 2) +
      theme_cowplot() +
      theme(panel.grid.major.y = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
            text = element_text(size = 10), legend.position = "none") +
      ylab("Detection probability") +
      xlab("Year")
    
    
    ## Assemble plots
    pAll <- plot_grid(pN, pD, pR, pE, pP, ncol = 1)
    
    ## Plot to pdf
    pdf(paste0("Plots/SimCheck_replicates/SimCheck_tAll_simSeed", simSeeds[i], "_", plotColors, ".pdf"), width = 10, height = 14)
    suppressWarnings(
      print(pAll)
    )
    dev.off()
    
    ## Plot to png
    png(paste0("Plots/SimCheck_replicates/SimCheck_tAll_simSeed", simSeeds[i], "_", plotColors, ".png"), width = 10, height = 14, units = "in", res = 300)
    suppressWarnings(
      print(pAll)
    )
    dev.off()
    
  }
  
}