#' Prepare line transect and known fate CMR data for integrated analysis
#'
#' @param d_trans tibble containing information on transects (events). Output of
#' wrangleData_LineTrans(). 
#' @param d_obs tibble containing information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.). Output of wrangleData_LineTrans(). 
#' @param d_cmr list with 2 elements. Surv1 and Surv2 are matrices of individuals 
#' released (column 1) and known to have survived (column 2) in each year (row)
#' for season 1 and season 2, respectively. Output of wrangleData_CMR().
#' @param augment logical. If TRUE, augments line transect data with N_aug 
#' (see below) entries. Default = FALSE.
#' @param N_aug integer. The number of dummy observations to add to data if
#' augment = TRUE (default = 1000).
#' @param dataVSconstants logical. If TRUE (default) returns a list of 2 lists
#' containing data and constants for analysis with Nimble. If FALSE, returns a
#' list containing all data and constants. 
#' @param save logical. If TRUE (default) saves prepared data in working 
#' directory as .rds file.
#'
#' @return A list or list of lists, depending on argument `dataVSconstants` 
#' (see above).
#' @export
#'
#' @examples

prepareInputData <- function(d_trans, d_obs, d_cmr, augment = FALSE, N_aug = 1000, dataVSconstants = TRUE, save = TRUE){
  
  # Constants #
  #-----------#
  
  ## Numbers of years and sites
  N_sites <- n_distinct(d_trans$locationID)
  N_years <- n_distinct(d_trans$Year)
  
  ## Number of age classes
  N_ageC <- 2
  
  ## Truncation distance
  W <- 200
  
  ## Scaling parameter
  scale1 <- 1000
  
  
  # Transect characteristics #
  #--------------------------#
  
  ## Tansect lengths
  TransLen <- d_trans %>% 
    dplyr::select(locationID, Year, sampleSizeValue) %>%
    dplyr::mutate(sampleSizeValue = sampleSizeValue/scale1) %>%
    reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
    arrange(locationID)
  
  L <- TransLen %>% select(-locationID) %>% as.matrix() 
  colnames(L) <- NULL
  
  ## Total covered area 
  A <- TransLen %>% 
    dplyr::select(-locationID) 
  colnames(A) <- NULL
  A <- colSums(A)*(W/scale1)*2
  
  
  # Observation distance from transect #
  #------------------------------------#
  
  temp_dist <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::select(Year, DistanceToTransectLine) %>%
    dplyr::mutate(Year2 = Year - (min(Year)) + 1)
  
  ## Distance to transect line
  y <- temp_dist$DistanceToTransectLine
  
  ## Observation year
  Year_obs <- temp_dist$Year2
  
  ## Number of observations
  N_obs <- length(y)
  
  ## Vector of 0's of same length as y
  zeros_dist <- rep(0, length(y))
  
  ## Augmented version of the data (if applicable)
  if(augment){
    aug <- c(rep(1, N_obs), rep(0, N_aug)) # Augmented inds. have y = 0 by definition
    y_aug <- c(y, rep(NA, N_aug)) # Value of distance are missing for the augmented observation
    y <- y_aug
  }
  
  # Number of birds/line (pooled age classes) #
  #-------------------------------------------#
  
  temp <- TransLen %>% select(locationID)
  
  TaksObs <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
    dplyr::mutate(cs = unknownJuvenile+unknownunknown+FemaleAdult+MaleAdult) %>%
    reshape2::dcast(locationID~Year, value.var = "cs", sum) %>%
    dplyr::right_join(., temp, by = c("locationID" = "locationID")) %>%
    replace(., is.na(.), 0) %>%
    dplyr::arrange(locationID)
  
  N_line_year <- TaksObs %>% 
    dplyr::select(-locationID) %>% 
    as.matrix() 
  colnames(N_line_year) <- NULL
  
  
  # Number of birds/line (by age class) #
  #-------------------------------------#
  
  ## Juveniles (& unknowns)
  TaksObs_J <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
    dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
    reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
    dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
    replace(., is.na(.), 0) %>%
    dplyr::arrange(locationID)
  
  N_J_line_year <- TaksObs_J %>% 
    dplyr::select(-locationID) %>% 
    as.matrix() 
  colnames(N_J_line_year) <- NULL
  
  
  ## Adults 
  TaksObs_A <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
    dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
    reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
    dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
    replace(., is.na(.), 0) %>%
    dplyr::arrange(locationID)
  
  N_A_line_year <- TaksObs_A %>% 
    dplyr::select(-locationID) %>% 
    as.matrix() 
  colnames(N_A_line_year) <- NULL
  
  
  ## Check and combine in array
  if(!(all(N_J_line_year + N_A_line_year == N_line_year))){
    warning("Number of observed adults and juveniles does not add up correctly. Double-check data.")
  }
  
  N_a_line_year <- array(NA, dim = c(2, nrow(N_line_year), ncol(N_line_year)))
  N_a_line_year[1,,] <- N_J_line_year
  N_a_line_year[2,,] <- N_A_line_year
  

  # Recruitment #
  #-------------#
  
  temp_Rec <- d_obs %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
    dplyr::mutate(R = unknownJuvenile + unknownunknown) %>%
    dplyr::mutate(Maletemp = unknownJuvenile + unknownunknown + FemaleAdult, MaleIndeks = if_else(Maletemp == 0, 1, 0)) %>%
    dplyr::select(Year, R, MaleIndeks) %>%
    dplyr::mutate(Year2 = Year - (min(Year)) + 1) %>%
    dplyr::filter(MaleIndeks == 0) # --> drop all observations of only males
  
  R_obs <- temp_Rec$R
  R_obs_year <- temp_Rec$Year2
  N_R_obs <- length(R_obs)
  
  
  # Data assembly #
  #---------------#
  
  ## Assembling all data in a list
  input.data <- list(
    R_obs = R_obs, # Observed numbers of recruits
    R_obs_year = R_obs_year, # Year of observed numbers of recruits
    N_R_obs = N_R_obs, # Total number of observations of numbers of recruits
    
    y = y, # Distance to transect line for each individual observation
    zeros_dist = zeros_dist, # Vector of 0's of same length as y
    Year_obs = Year_obs, # Year of each observation
    N_obs = N_obs, # Total number of observations
    N_aug = ifelse(augment, N_aug, NA), # Number of augmented observations
    aug = ifelse(augment, aug, NA), # Dummy index to separate real and augmented data
    
    N_line_year = N_line_year, # Number of birds observed per site per year
    N_a_line_year = N_a_line_year, # Number of birds observed per ageclass per site per year
    L = L, # Transect length per site and year
    
    N_years = N_years, # Number of years with data
    N_sites = N_sites, # Total number of monitored sites
    
    A = A, # Total covered area per year
    W = W, # Truncation distance
    scale1 = scale1, # Scaling parameter
    N_ageC = N_ageC, # Number of age classes
    
    Survs1 = d_cmr$Survs1, # Season 1 releases & survivors
    Survs2 = d_cmr$Survs2 # Season 2 releases & survivors
  )
  
  ## Assembling Nimble data
  nim.data <- list(R_obs = input.data$R_obs, y = input.data$y, 
                   zeros.dist = input.data$zeros_dist, L = input.data$L,
                   aug = input.data$aug,
                   N_line_year = input.data$N_line_year, 
                   N_a_line_year = input.data$N_a_line_year, 
                   A = input.data$A,
                   Survs1 = d_cmr$Survs1, Survs2 = d_cmr$Survs2)
  
  ## Assembling Nimble constants
  nim.constants <- list(N_years = input.data$N_years, W = input.data$W, scale1 = scale1,
                        N_aug = input.data$N_aug,
                        N_obs = input.data$N_obs, Year_obs = input.data$Year_obs,
                        N_sites = input.data$N_sites, 
                        R_obs_year = input.data$R_obs_year, N_R_obs = input.data$N_R_obs,
                        N_ageC = N_ageC)
  
  ## Make final data list to return
  if(dataVSconstants){
    rype.data <- list(nim.data = nim.data,
                      nim.constants = nim.constants)
  }else{
    rype.data <- input.data
  }
  
  ## Optional: save data as .rds
  if(save){
    saveRDS(rype.data, file = "RypeData_forIM.rds")
  }
  
  ## Return data
  return(rype.data)
}



