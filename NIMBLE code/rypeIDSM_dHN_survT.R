rypeIDSM <- nimbleCode({
  
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
  sigmaT.dd ~ dunif(0, 20)
  
  for(t in 1:N_years){
    epsT.dd[t] ~ dnorm(0, sd = sigmaT.dd) 
  }
  
  
  
  ########################################################
  for(t in 1:N_years){
    log(sigma[t]) <- mu.dd + epsT.dd[t]
    sigma2[t] <- sigma[t] * sigma[t]
    
    # effective strip width
    esw[t] <- sqrt(pi * sigma2[t] / 2) 
    p[t] <- min(esw[t], W) / W
  }
  
  ########################################################   
  for (i in 1:N_obs){ 
    # LIKELIHOOD
    # using nimbleDistance::dHN
    y[i] ~ dHN(sigma = sigma[Year_obs[i]], Xmax = W, point = 0)
  }
  

  ###################################################
  ## Random effects model for R (i.e. )
  ## Data: 
  
  # R_obs <- vector with number of chicks / obs [0 - 12]
  # N_Years <- number of years in time series 
  # N_obs <- number of observations 
  # Year_obs <- vector with years for each observation. 
  
  ## Priors; 
  for (t in 1:N_years){
    epsT.R[t] ~ dnorm(0, sd = sigmaT.R)
  }
  
  Mu.R  ~ dunif(0, 10)
  sigmaT.R ~ dunif(0, 5)
  
  ## Constraints;
  R_year[1:N_years] <- exp(log(Mu.R) + epsT.R[1:N_years])
  
  ## Likelihood;
  for (i in 1:N_sumR_obs){
    
    sumR_obs[i] ~ dpois(R_year[sumR_obs_year[i]]*sumAd_obs[i])
  }
  
  
  
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
  
  Mu.D1 ~ dunif(0, 10)
  sigma.D ~ dunif(0, 20)
  
  ratio.JA1 ~ dunif(0, 1)
  
  ## State model
  for (j in 1:N_sites){
    
    N_exp[1, j, 1] ~ dpois(Density[1, j, 1]*L[j, 1]*W*2) 
    N_exp[2, j, 1] ~ dpois(Density[2, j, 1]*L[j, 1]*W*2) 
    
    Density[1, j, 1] <- exp(log(Mu.D1) + eps.D1[j])*ratio.JA1             ## random effects model for spatial variation in density for year 1
    Density[2, j, 1] <- exp(log(Mu.D1) + eps.D1[j])*(1-ratio.JA1)
    
    ## Detection model year 1
    for(x in 1:N_ageC){
      N_a_line_year[x, j, 1] ~ dpois(p[1]*N_exp[x, j, 1])
    }
    
    #N_line_year[j, 1] ~ dpois(p[1]* sum(N_exp[1:N_ageC, j, 1]))
  }
  
  #####################################################
  ## Model for survival; 
  
  ## Priors
  Mu.S ~ dunif(0, 1)
  Mu.S1 ~ dunif(0.25, 0.9)
  
  sigmaT.S ~ dunif(0, 5)
  epsT.S1.prop ~ dunif(0, 1) # Proportion of random year effect that will be allocated to season 1
  
  ## Constraints
  logit(S[1:N_years]) <- logit(Mu.S) + epsT.S[1:N_years]
  
  S1[1:N_years] <- logit(Mu.S1) + epsT.S1.prop*epsT.S[1:N_years]
  S2[1:N_years] <- S[1:N_years]/S1[1:N_years]
  
  for(t in 1:N_years){
    epsT.S[t] ~ dnorm(0, sd = sigmaT.S) # Temporal RE
  }
  
  ## Data likelihoods
  for (t in 1:N_years_RT){
    
    Survs1[t, 2] ~ dbinom(S1[year_Survs[t]], Survs1[t, 1])
    Survs2[t, 2] ~ dbinom(S2[year_Survs[t]], Survs2[t, 1])
    
  }
  
  
  #####################################################    
  ### Model for year 2 - n.years; 
  ### post-breeding census
  
  for(j in 1:N_sites){
    for(t in 2:N_years){
      
      ## Process model
      Density[2, j, t] <- sum(Density[1:N_ageC, j, t-1])*S[t-1] 
      
      if(R_perF){
        Density[1, j, t] <- (Density[2, j, t]/2)*R_year[t] 
      }else{
        Density[1, j, t] <- Density[2, j, t]*R_year[t]
      }
      
      N_exp[1:N_ageC, j, t] <- Density[1:N_ageC, j, t]*L[j, t]*W*2
      
      ## Detection model year 2 - T
      for(x in 1:N_ageC){
        N_a_line_year[x, j, t] ~ dpois(p[t]*N_exp[x, j, t])
      }
      
      #N_line_year[j, t] ~ dpois(p[t]*sum(N_exp[1:N_ageC, j, t]))
      
    }
  }
  
  ####################################################
  ## Observation model
  ## P is estimated in distance sampling component - based on 
  ## distance to transect line data; 
  
  ####################################################
  ### Derived parameters; Nt and Dt
  
  for (t in 1:N_years){
    N_tot_exp[t] <- sum(N_exp[1, 1:N_sites, t] + N_exp[2, 1:N_sites, t])    ## Summing up expected number of birds in covered area; 
    #D[t] <- N_tot_exp[t] / A[t]       ## Deriving density as N/A     
  }
  
})
