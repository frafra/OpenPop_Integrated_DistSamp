
modM3.code <- nimbleCode({
  
  ####################################################
  #### Distance sampling half normal detection function; 
  ## Data; 
  
  # y <- vector with distances to transect line; 
  # N_years <- number of years
  # N_obs <- number of observations
  # Year_obs <- vector with year for each observation
  
  
  # priors for distance model; 
  
  pi <- 3.141593
  
  #random year effect for distance sampling model; 
  ## Priors for hyper-parameters
  mu.dd ~ dunif(-10, 100)
  sigma.dd ~ dunif(0, 20)
  
  for(t in 1:N_years){
    eps.dd[t] ~ dnorm(0, sd = sigma.dd) 
  }
  
  
  
  ########################################################
  for(t in 1:N_years){
    log(sigma[t]) <- mu.dd + eps.dd[t]
    sigma2[t] <- sigma[t]*sigma[t]
    
    # effective strip width
    esw[t] <- sqrt(pi * sigma2[t] / 2) 
    p[t] <- esw[t] / W
    #f0[t] <- 1/esw[t] #assuming all detected on the line
  }
  
  ########################################################   
  for (i in 1:N_obs){ 
    # LIKELIHOOD
    # using zeros trick
    y[i] ~ dunif(0, W) 
    L.f0[i] <- exp(-y[i]*y[i] / (2*sigma2[Year_obs[i]])) * 1/esw[Year_obs[i]] #y are the distances
    nlogL.f0[i] <-  -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
  }
  
  ## Hierarchical node: esw[t]; 
  ##
  ###################################################
  ## Random effects model for R (i.e. )
  ## Data: 
  
  # R_obs <- vector with number of chicks / obs [0 - 12]
  # N_Years <- number of years in time series 
  # N_obs <- number of observations 
  # Year_obs <- vector with years for each observation. 
  
  ## Priors; 
  for (t in 1:N_years){
    eps.R[t] ~ dnorm(0, sd = sigma.R)
  }
  
  mu.R  ~ dunif(-5, 5)
  sigma.R ~ dunif(0, 15)
  
  ## Likelihood;
  for (i in 1:N_R_obs){
    
    R_obs[i] ~ dpois(R_exp[i])
    # Linear predictor
    R_exp[i] <- exp(mu.R + eps.R[R_obs_year[i]])   
  }
  
  ### Derive parameter: annual predictions for R
  for (t in 1:N_years){
    R_year[t] <- exp(mu.R + eps.R[t])
  }
  
  ## Hierarchical node: R_year
  
  ########################################################################
  ### MODEL FOR Density in year 1:
  ### Simple random effects model 
  # Data; 
  
  # N_sites <- number of sites
  # N_line_year <- number of birds pr. line 
  # L <- length of transect lines
  
  
  ## Priors; 
  
  for(j in 1:N_sites){
    eps.D1[j] ~ dnorm(0, sd = sigma.D)
  }
  
  mu.D1 ~ dunif(-3, 30)
  sigma.D ~ dunif(0, 20)
  
  # State model
  for (j in 1:N_sites){
    N_exp[j, 1] ~ dpois(Density[j, 1]*L[j, 1]*(W/scale1)*2)      ## Expected number of birds
    Density[j, 1] <- exp(mu.D1 + eps.D1[j])             ## random effects model for spatial variation in density for year 1
    ## Detection model year 1
    N_line_year[j, 1] ~ dpois(p[1]* N_exp[j,1])
  }
  #####################################################
  ## Model for survival; 
  
  S ~ dunif(0.25, 0.65)
  
  
  #####################################################    
  ### Model for year 2 - n.years; 
  ### post-breeding census
  
  for(j in 1:N_sites){
    for(t in 2:N_years){
      
      ## Process model
      Density[j, t] <- (Density[j, t-1] * S) + (Density[j, t-1]*S*R_year[t]/2) 
      N_exp[j, t] <- Density[j, t]*L[j, t]*(W/scale1)*2
      
      ## Detection model year 2 - T
      N_line_year[j, t] ~ dpois(p[t]*N_exp[j, t])
      
      ## Pop structure - derived
      D1[j, t] <- Density[j, t-1] * S
      D2[j, t] <- Density[j, t-1] * (S*R_year[t]/2)
      
    }}
  
  ####################################################
  ## Observation model
  ## P is estimated in distance sampling component - based on 
  ## distance to transect line data; 
  
  ####################################################
  ### Derived parameters; Nt and Dt
  
  for (t in 1:N_years){
    N_tot_exp[t] <- sum(N_exp[1:N_sites, t])    ## Summing up expected number of birds in covered area; 
    D[t] <- N_tot_exp[t] / A[t]       ## Deriving density as N/A     
  }
  
})