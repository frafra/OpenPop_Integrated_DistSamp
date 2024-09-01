#' Visualize comparison of model estimates with true parameter values from data simulation
#'
#' @param PlotColors string specifying color palette to use for plots. 
#' Currently supports "customRainbow" (default), "Temps", and "Zissou1".
#' @param thin integer. Thinning rate to use on samples for plotting. Defaults 
#' to 1 (no thinning) unless specified. 
#' @param saveMetrics logical. If TRUE (default), performance metrics are saved
#' as RDS files in the working directory. 
#' @return a vector of png plot names. The plots can be found in Plots/TimeSeries.
#' 
#' @export 
#'
#' @examples


plotSimCheckMetrics_replicates <- function(plotColors = "customRainbow", thin = 1, saveMetrics = TRUE) {
  
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
  runSeeds <- readRDS("simModelFits_sum/seedInfo.rds")
  
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
  
  ## Add true values to output data
  Mu_R$trueValue <- simParams$Mu_R[1]
  sigmaT_R$trueValue <- simParams$sigmaT_R[1]
  Mu_dd$trueValue <- simParams$Mu_dd[1]
  sigmaT_dd$trueValue <- simParams$sigmaT_dd[1]
  Mu_S <- Mu_S %>%
    dplyr::mutate(trueValue = dplyr::case_when(Surv == "S" ~ simParams$Mu_S[1],
                                               Surv == "S1" ~ simParams$Mu_S1[1],
                                               Surv == "S2" ~ simParams$Mu_S2[1]))
  
  R_year <- R_year %>%
    dplyr::left_join(R_simData, by = c("year", "dataSetID", "simSeed")) %>%
    dplyr::rename(trueValue = realizedR, trueValue_mean = predictedR)
  esw_year <- esw_year %>%
    dplyr::left_join(det_simData[, c("year", "esw", "dataSetID", "simSeed")], by = c("year", "dataSetID", "simSeed")) %>%
    dplyr::rename(esw = esw.x, trueValue = esw.y)
  p_year <- p_year %>%
    dplyr::left_join(det_simData[, c("year", "p", "dataSetID", "simSeed")], by = c("year", "dataSetID", "simSeed")) %>%
    dplyr::rename(p = p.x, trueValue = p.y)
  N_tot <- N_tot %>%
    dplyr::left_join(N_simData[, c("year", "N_tot", "dataSetID", "simSeed")], by = c("year", "dataSetID", "simSeed")) %>%
    dplyr::rename(trueValue = N_tot)
  Density_year <- Density_year[, c("year", "density", ".chain", ".iteration", ".draw", "dataSetID", "simSeed", "runID", "runSeed")] %>%
    dplyr::left_join(D_simData, by = c("year", "dataSetID", "simSeed")) %>%
    dplyr::rename(trueValue = Mean.D)
  
  ## Collate data on single and time-dependent parameters
  
  # Single parameters
  c.Mu_R <- Mu_R %>% dplyr::rename(estimate = Mu.R)
  c.sigmaT_R <- sigmaT_R %>% dplyr::rename(estimate = sigmaT.R)
  c.Mu_dd <- Mu_dd %>% dplyr::rename(estimate = Mu.dd)
  c.sigmaT_dd <- sigmaT_dd %>% dplyr::rename(estimate = sigmaT.dd)
  c.S <- Mu_S %>% dplyr::rename(estimate = S, lab_code = Surv)
  
  out.single <- dplyr::bind_rows(c.Mu_R, c.sigmaT_R,
                                 c.Mu_dd, c.sigmaT_dd,
                                 c.S)
  
  # Time-dependent parameters
  c.R_year <- R_year %>% dplyr::rename(estimate = R_year)
  c.esw_year <- esw_year %>% dplyr::rename(estimate = esw)
  c.p_year <- p_year %>% dplyr::rename(estimate = p)
  c.N_tot <- N_tot %>% dplyr::rename(estimate = N_tot_exp)
  c.Density_year <- Density_year %>% dplyr::rename(estimate = density)
  
  out.annual <- dplyr::bind_rows(c.R_year, c.esw_year, c.p_year, 
                                 c.N_tot, c.Density_year, 
                                 .id = "id") %>%
    dplyr::mutate(Parameter = dplyr::case_when(id == "1" ~ "R",
                                               id == "2" ~ "esw",
                                               id == "3" ~ "p", 
                                               id == "4" ~ "N", 
                                               id == "5" ~ "D"))
  
  # Calculate summary metrics #
  #---------------------------#

  # Single parameters
  metrics.single <- out.single %>%
    dplyr::group_by(lab_code, dataSetID, runID) %>%
    dplyr::mutate(diff = estimate - trueValue,
                  est_above_true = ifelse(estimate > trueValue, 1, 0)) %>%
    dplyr::summarise(Bayes_p = mean(est_above_true),
                     MSD = sum(diff)/n(),
                     RMSD = sqrt(sum(diff^2)/n()),
                     .groups = "keep")
  
  # Time-dependent parameters
  metrics.annual <- out.annual %>%
    dplyr::group_by(Parameter, year, dataSetID, runID) %>%
    dplyr::mutate(diff = estimate - trueValue,
                  est_above_true = ifelse(estimate > trueValue, 1, 0)) %>%
    dplyr::summarise(Bayes_p = mean(est_above_true),
                     MSD = sum(diff)/n(),
                     RMSD = sqrt(sum(diff^2)/n()),
                     .groups = "keep")
                       
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
  
  ##########################
  # SINGLE PARAMETER PLOTS #
  ##########################
  
  ## Add parameter labels to single parameter data
  out.single <- out.single %>%
    dplyr::mutate(ParamName = dplyr::case_when(lab_code == "Mu.R" ~ "mu[R]",
                                               lab_code == "sigmaT.R" ~ "sigma[R]",
                                               lab_code == "Mu.dd" ~ "mu[sigma]",
                                               lab_code == "sigmaT.dd" ~ "sigma[sigma]",
                                               lab_code %in% c("S", "S1", "S2") ~ lab_code)) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("mu[R]", "sigma[R]", "S", "S1", "S2", "mu[sigma]", "sigma[sigma]")))
  
  metrics.single <- metrics.single %>%
    dplyr::mutate(ParamName = dplyr::case_when(lab_code == "Mu.R" ~ "mu[R]",
                                               lab_code == "sigmaT.R" ~ "sigma[R]",
                                               lab_code == "Mu.dd" ~ "mu[sigma]",
                                               lab_code == "sigmaT.dd" ~ "sigma[sigma]",
                                               lab_code %in% c("S", "S1", "S2") ~ lab_code)) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("mu[R]", "sigma[R]", "S", "S1", "S2", "mu[sigma]", "sigma[sigma]")))
  
  ## Calculate metric means
  metricMeans.single <- metrics.single %>%
    dplyr::group_by(ParamName) %>%
    dplyr::summarise(Bayes_p_mean = mean(Bayes_p), 
                     MDS_mean = mean(MSD), 
                     RMSD_mean = mean(RMSD), 
                     .groups = "keep")
  
  ## Plot posterior densities for single parameters
  p_a <- ggplot(data = out.single, aes(x = estimate)) +
    geom_density(aes(linetype = runSeed, color = simSeed), fill = NA)  +
    geom_vline(aes(xintercept = trueValue), col = "black", linetype = "dashed") +
    scale_color_manual(values = alpha(plot.cols, 0.75)) + 
    scale_linetype_manual(values = rep("solid", nlevels(Mu_R$runSeed))) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("Estimate") +
    ylab("")
  
  ## Plot Bayesian p-values for single parameters
  p_b <- ggplot(data = metrics.single, aes(x = Bayes_p)) + 
    geom_histogram(aes(fill = dataSetID), position = "stack", bins = 20) + 
    geom_vline(aes(xintercept = 0.5), col = "black", linetype = "dashed") + 
    geom_vline(data = metricMeans.single, aes(xintercept = Bayes_p_mean), color = "#B065E6") + 
    geom_text(data = metricMeans.single, aes(label = round(Bayes_p_mean, digits = 3), x = Bayes_p_mean + 0.10, y = Inf, vjust = 1.5), color = "#B065E6") + 
    scale_fill_manual(values = alpha(plot.cols, 0.75)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free_y", axes = "all_x", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("Bayesian p-value") +
    ylab("")
 
  ## Plot root mean square deviation for single parameters
  p_c <- ggplot(data = metrics.single, aes(x = RMSD)) + 
    geom_histogram(aes(fill = dataSetID), position = "stack",bins = 20) + 
    geom_vline(data = metricMeans.single, aes(xintercept = RMSD_mean), color = "#B065E6") + 
    geom_text(data = metricMeans.single, aes(label = round(RMSD_mean, digits = 3), x = RMSD_mean * 1.2, y = Inf, vjust = 1.5), color = "#B065E6") +
    scale_fill_manual(values = alpha(plot.cols, 0.75)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("RMSD") +
    ylab("")
  
  ## Combine plots
  p_comb <- plot_grid(p_a, p_b, p_c, nrow = 1)
  
  ## Plot to pdf
  pdf(paste0("Plots/SimCheck_replicates/SimCheckMetrics_singleParams_", plotColors, ".pdf"), width = 9, height = 14)
  suppressWarnings(
    print(p_comb)
  )
  dev.off()
  
  ## Plot to png
  png(paste0("Plots/SimCheck_replicates/SimCheckMetrics_singleParams_", plotColors, ".png"), width = 9, height = 14, units = "in", res = 300)
  suppressWarnings(
    print(p_comb)
  )
  dev.off()
  
  
  ##################################
  # TIME-DEPENDENT PARAMETER PLOTS #
  ##################################
  
  ## Add parameter labels to annual parameter data
  out.annual <- out.annual %>%
    dplyr::mutate(ParamName = dplyr::case_when(Parameter == "R" ~ "R[t]",
                                               Parameter == "esw" ~ "esw[t]",
                                               Parameter == "p" ~ "p[t]",
                                               Parameter == "N" ~ "N[t]",
                                               Parameter == "D" ~ "D[t]")) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("R[t]", "esw[t]", "p[t]", "N[t]", "D[t]")))
  
  metrics.annual <- metrics.annual %>%
    dplyr::mutate(ParamName = dplyr::case_when(Parameter == "R" ~ "R[t]",
                                               Parameter == "esw" ~ "esw[t]",
                                               Parameter == "p" ~ "p[t]",
                                               Parameter == "N" ~ "N[t]",
                                               Parameter == "D" ~ "D[t]")) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("R[t]", "esw[t]", "p[t]", "N[t]", "D[t]")))
  
  
  ## Drop derived parameters for simples representation
  out.annual <- out.annual %>%
    dplyr::filter(Parameter %in% c("R", "D", "p"))
  
  metrics.annual <- metrics.annual %>%
    dplyr::filter(Parameter %in% c("R", "D", "p"))
  
  ## Calculate metric means
  metricMeans.annual <- metrics.annual %>%
    dplyr::group_by(ParamName) %>%
    dplyr::summarise(Bayes_p_mean = mean(Bayes_p), 
                     MDS_mean = mean(MSD), 
                     RMSD_mean = mean(RMSD), 
                     .groups = "keep")
  
  ## Summarise posteriors into medians
  out.annual.median <- out.annual %>%
    dplyr::group_by(dataSetID, runID, year, trueValue, trueValue_mean, Parameter, ParamName) %>%
    dplyr::summarise(medianEst = median(estimate), .groups = "keep")
  
  ## Calculate average slopes per parameter
  slopes.annual <- data.frame()
  for(i in unique(out.annual$ParamName)){
    
    sub.data <- subset(out.annual, ParamName == i)
    sub.data2 <- subset(out.annual.median, ParamName == i)
    
    linModel <- lm(sub.data$estimate ~ sub.data$trueValue)
    linModel2 <- lm(sub.data2$medianEst ~ sub.data2$trueValue)
    
    param.data <- data.frame(
      ParamName = i, 
      Slope = unname(round(linModel$coefficients[2], digits = 3)),
      Intercept = unname(round(linModel$coefficients[1], digits = 3)),
      Slope_median = unname(round(linModel2$coefficients[2], digits = 3)),
      Intercept_median = unname(round(linModel2$coefficients[1], digits = 3))
    )
    
    slopes.annual <- rbind(slopes.annual, param.data)
  }

  metricMeans.annual <- metricMeans.annual %>%
    dplyr::left_join(slopes.annual, by = "ParamName")
  
  ## Plot estimated versus true values for annual parameters (posterior samples)
  p_a <- ggplot(data = out.annual, aes(x = trueValue, y = estimate)) + 
    geom_point(aes(color = dataSetID, shape = runID)) + 
    geom_abline(aes(slope = 1, intercept = 0), col = "black", linetype = "dashed") + 
    geom_abline(data = metricMeans.annual, aes(slope = Slope, intercept = Intercept), col = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = paste0("a = ", Intercept), x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5), color = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = paste0("b = ", Slope), x = -Inf, y = Inf, vjust = 2.5, hjust = -0.5), color = "#B065E6") + 
    scale_color_manual(values = alpha(plot.cols, 0.25)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("True value") +
    ylab("Estimated value")
  
  ## Plot estimated versus true values for annual parameters (median)
  p_a2 <- ggplot(data = out.annual.median, aes(x = trueValue, y = medianEst)) + 
    geom_point(aes(color = dataSetID, shape = runID)) + 
    geom_abline(aes(slope = 1, intercept = 0), col = "black", linetype = "dashed") + 
    geom_abline(data = metricMeans.annual, aes(slope = Slope_median, intercept = Intercept_median), col = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = paste0("a = ", Intercept_median), x = -Inf, y = Inf, vjust = 1.5, hjust = -0.5), color = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = paste0("b = ", Slope_median), x = -Inf, y = Inf, vjust = 2.5, hjust = -0.5), color = "#B065E6") + 
    scale_color_manual(values = alpha(plot.cols, 0.5)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("True value") +
    ylab("Estimated value (median)")
  
  ## Plot Bayesian p-values for annual parameters
  p_b <- ggplot(data = metrics.annual, aes(x = Bayes_p)) + 
    geom_histogram(aes(fill = dataSetID), position = "stack", bins = 20) + 
    geom_vline(aes(xintercept = 0.5), col = "black", linetype = "dashed") + 
    geom_vline(data = metricMeans.annual, aes(xintercept = Bayes_p_mean), color = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = round(Bayes_p_mean, digits = 3), x = Bayes_p_mean + 0.10, y = Inf, vjust = 1.5), color = "#B065E6") + 
    scale_fill_manual(values = alpha(plot.cols, 0.75)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free_y", axes = "all_x", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("Bayesian p-value") +
    ylab("")
  
  ## Plot root mean square deviation for single parameters
  p_c <- ggplot(data = metrics.annual, aes(x = RMSD)) + 
    geom_histogram(aes(fill = dataSetID), position = "stack",bins = 20) + 
    geom_vline(data = metricMeans.annual, aes(xintercept = RMSD_mean), color = "#B065E6") + 
    geom_text(data = metricMeans.annual, aes(label = round(RMSD_mean, digits = 3), x = RMSD_mean * 1.2, y = Inf, vjust = 1.5), color = "#B065E6") +
    scale_fill_manual(values = alpha(plot.cols, 0.75)) + 
    facet_wrap(~ParamName, labeller = label_parsed, scales = "free", ncol = 1) + 
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "#8ccde3", linewidth = 0.25, linetype = 2), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          text = element_text(size = 10), legend.position = "none") +
    xlab("RMSD") +
    ylab("")
  
  ## Combine plots
  p_comb <- plot_grid(p_a, p_b, p_c, nrow = 1, rel_widths = c(1.33, 1, 1))
  p_comb2 <- plot_grid(p_a2, p_b, p_c, nrow = 1, rel_widths = c(1.33, 1, 1))
  
  ## Plot to pdf
  pdf(paste0("Plots/SimCheck_replicates/SimCheckMetrics_annualParams_", plotColors, ".pdf"), width = 9, height = 7)
  suppressWarnings(
    print(p_comb)
  )
  suppressWarnings(
    print(p_comb2)
  )
  dev.off()
  
  ## Plot to png
  png(paste0("Plots/SimCheck_replicates/SimCheckMetrics_annualParams_", plotColors, ".png"), width = 9, height = 7, units = "in", res = 300)
  suppressWarnings(
    print(p_comb)
  )
  dev.off()

  png(paste0("Plots/SimCheck_replicates/SimCheckMetrics_annualParams_median_", plotColors, ".png"), width = 9, height = 7, units = "in", res = 300)
  suppressWarnings(
    print(p_comb2)
  )
  dev.off()
  
  ## Save metrics for reporting
  if(saveMetrics){
    
    metrics.list <- list(metrics.single = metrics.single,
                         metricMeans.single = metricMeans.single,
                         metrics.annual = metrics.annual,
                         metricMeans.annual = metricMeans.annual)
    
    saveRDS(metrics.list, file = "SimChecks_PerformanceMetrics.rds")
  }
}