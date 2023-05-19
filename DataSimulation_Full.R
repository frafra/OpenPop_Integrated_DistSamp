
mySeed <- 0
set.seed(mySeed)


#########################
# SIMULATION PARAMETERS #
#########################


# General simulation parameters #
#-------------------------------#

Amax <- 2 # Number of age classes
Tmax <- 15 # Number of years
#Tmax <- 30 # Number of years
Jmax <- 50 # Number of sites/transect lines


# Vital rate parameters #
#-----------------------#

## Annual survival
Mu.S <- 0.35 # Average annual survival probability
sigmaT.S <- 0 # SD of random year variation in survival
sigmaJ.S <- 0 # SD of random site variation in survival

## Reproduction
Mu.R <- 2 # Average number of chicks in August
sigmaT.R <- 0.5 # SD of random year variation in number of chicks
sigmaJ.R <- 0 # SD of random site variation in number of chicks

## Juvenile summer survival
#Mu.sJ <- 0.2 # Average summer survival of chicks
#sigmaT.sJ <- 0 # SD of random year variation in chick survival
#sigmaJ.sJ <- 0 # SD of random site variation in survival


# Population parameters #
#-----------------------#

# Initial population numbers per site
N1 <- matrix(NA, ncol = Amax, nrow = Jmax)
for(j in 1:Jmax){
  N1[j,2] <- round(runif(1, 3, 8)) # Number of adults
  N1[j,1] <- rpois(1, N1[j,2]*Mu.R) # Number of juveniles
}

# Average group size
avg_Gsize <- 5.6


# Data & observation parameters #
#-------------------------------#

## Line-transect distance sampling
min.Tlength <- 1000 # Minimum transect length
max.Tlength <- 1000  # Maximum transect length

W <- 200 # Truncation distance (max. distance at which observation is possible)
  
#pi <- 3.141593 # Approximation of pi

Mu.dd <- 75 # Average width parameter for half-normal detection function
sigmaT.dd <- 0 # SD of random year variation in detection probability
sigmaJ.dd <- 0 # SD of random line variation in detection probability

## Known-fate radio-telemetry
Tmin.RT <- 5 # First year for which radio-telemetry data has been collected
Tmax.RT <- 10 # Last year for which radio-telemetry data has been collected

# Average number of individuals fitted with transmitters each year
nind.avg.RT <- 30


# ## Nest survey
# Tmin.NS <- 1 # First year for which nest survey data has been collected
# Tmax.NS <- 15 # Last year for which nest survey data has been collected
# 
# # Average number of nests monitored each year
# nind.avg.NS <- 40


#########################
# VITAL RATE SIMULATION #
#########################

## Function to simulate year- and site-specific vital rates
simulate.VRs <- function(Tmax, Jmax,
                         Mu.S, Mu.R, 
                         sigmaT.S, sigmaT.R, 
                         sigmaJ.S, sigmaJ.R){
  
  # Sample random year effects for all vital rates
  epsilonT.S <- rnorm(Tmax, mean = 0, sd = sigmaT.S)
  epsilonT.R <- rnorm(Tmax, mean = 0, sd = sigmaT.R)
  #epsilonT.sJ <- rnorm(Tmax, mean = 0, sd = sigmaT.sJ)
  
  # Sample random site effects for all vital rates
  epsilonJ.S <- rnorm(Jmax, mean = 0, sd = sigmaJ.S)
  epsilonJ.R <- rnorm(Jmax, mean = 0, sd = sigmaJ.R)
 #epsilonJ.sJ <- rnorm(Jmax, mean = 0, sd = sigmaJ.sJ)
  
  # Calculate year- and site-specific vital rates
  S <- R <- sJ <- matrix(NA, nrow = Jmax, ncol = Tmax)
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      
      S[j,t] <- plogis(qlogis(Mu.S) + epsilonT.S[t] + epsilonJ.S[j])
      R[j,t] <- exp(log(Mu.R) + epsilonT.R[t] + epsilonJ.R[j])
      #sJ[j,t] <- plogis(qlogis(Mu.sJ) + epsilonT.sJ[t] + epsilonJ.sJ[j])
      
    }
  }
  
  # Arrange vital rates in a list
  VR.list <- list(S = S, R = R)
  
}

## Simulate vital rates
VR.list <- simulate.VRs(Tmax, Jmax,
                        Mu.S, Mu.R, 
                        sigmaT.S, sigmaT.R,
                        sigmaJ.S, sigmaJ.R)


#########################
# POPULATION SIMULATION #
#########################

## Function for simulating temporal population dynamics in each site
simulate.pop <- function(Amax, Tmax, Jmax, VR.list, N1, stochastic = TRUE, plot = FALSE){
  
  # Prepare arrays for storing population projections
  N <- array(NA, dim = c(Jmax, Amax, Tmax))
  #Chicks <- matrix(NA, nrow = Jmax, ncol = Tmax)
  
  # Initialize population projection
  N[,,1] <- N1
  
  # Print simulation info
  if(stochastic){
    message('Simulating population dynamics with demographic stochasticity...')
  }else{
    message('Simulating deterministic population dynamics...')
  }
  
  # Simulate survival & reproduction over time in each site
  for(j in 1:Jmax){
    for(t in 1:(Tmax-1)){
      
      if(stochastic){
    
        # - Survivors
        N[j,2,t+1] <- rbinom(1, size = sum(N[j,1:Amax,t]), prob = VR.list$S[j,t])
      
        # - Reproduction of survivors
        N[j,1,t+1] <- rpois(1, lambda = N[j,2,t+1]*(VR.list$R[j,t+1]))
      
        # - Young-of-the-year recruiting in the end of summer
       # N[j,1,t+1] <- rbinom(1, size = Chicks[j,t+1], prob = VR.list$sJ[j,t+1])
    
    
      }else{
    
        # - Survivors
        N[j,2,t+1] <- sum(N[j,1:Amax,t])*VR.list$S[j,t]
      
        # - Reproduction of survivors
        N[j,1,t+1] <- N[j,2,t+1]*(VR.list$R[j,t+1])
      
        # - Young-of-the-year recruiting in the end of summer
        #N[j,1,t+1] <- Chicks[j,t+1]*VR.list$sJ[j,t+1]
      }
    }
  }
  
  # Plot trajectories
  #par(mfrow = c(1,2))
  plot(colSums(N[,1,]), type = 'l', col = 'red', lty = 'dashed', 
       ylab = 'Number', xlab = 'Time', ylim = c(0, max(colSums(N[,1,]), na.rm = T)))
  lines(colSums(N[,2,]), col = 'blue')
  legend('topleft', legend = c('Recruits', 'Adults'), 
               col = c('red', 'blue'),
               lty = c('dashed', 'solid'), 
               cex = 0.8, y.intersp = 0.1, bty = 'n')
  
  matplot(apply(N, 1, colSums), type = 'l', ylab = 'Number per site', xlab = 'Time')
  
  # Arrange simulated numbers in a list & return
  SimData <- list(N = N)
  return(SimData)
}

## Simulate population trajectories for all sites
SimData <- simulate.pop(Amax, Tmax, Jmax, VR.list, N1, stochastic = TRUE, plot = TRUE)


################################
# GROUP AGGREGATION SIMULATION #
################################

## Function for simulating group aggregation of population
simulate.groups <- function(Jmax, Tmax, N.age, avg_Gsize, discard0 = TRUE){
  
  # Set up group data frame and matrix
  G.age <- data.frame()
  G.no <- matrix(0, nrow = Jmax, ncol = Tmax)
  
  for(j in 1:Jmax){
    for(t in 1:Tmax){
      
      if(sum(N.age[j, 1:2, t]) > 0){
        # Set number of groups to simulate
        G.noTot <- sum(N.age[j, 1:2, t])/avg_Gsize
        G.noTot <- ifelse(G.noTot < 1, 1, round(G.noTot))
        
        # Distribute juveniles among groups
        G.juv <- rmultinom(1, size = N.age[j, 1, t], prob = rep(1/G.noTot, G.noTot))
        
        # Distribute adults among groups
        G.ad <- rmultinom(1, size = N.age[j, 2, t], prob = rep(1/G.noTot, G.noTot))
        
        # Remove any potential 0 observations
        if(discard0){
          obs0 <- which(G.juv + G.ad == 0)
          if(length(obs0) > 0){
            G.juv <- G.juv[-obs0]
            G.ad <- G.ad[-obs0]
          }
        }
        
        
        # Collate in data frame
        G.data <- data.frame(site = j,
                             year = t,
                             no_juv = G.juv,
                             no_ad = G.ad)
        G.age <- rbind(G.age, G.data)
        
        # Set number of groups for site-year
        G.no[j, t] <- nrow(G.data)
      }
    }
  }    
  
  return(G.age)
}

## Simulate group aggregation of population
G.age <- simulate.groups(Jmax, Tmax, N.age = SimData$N, avg_Gsize, discard0 = TRUE)

#################################
# LINE TRANSECT DATA SIMULATION #
#################################

## Function for simulating hierarchical distance sampling data
simulate.HDSdata <- function(Jmax, Tmax, G.age,
                             Mu.dd, sigmaT.dd, sigmaJ.dd,
                             W, min.Tlength, max.Tlength,
                             discard0){
  
  # Set site-specific transect lengths
  # NOTE: These are not currently influencing the data simulation in any way
  L <- matrix(runif(Jmax*Tmax, min.Tlength, max.Tlength), nrow = Jmax, ncol = Tmax)
  
  # Sample random year and site effects on sigma parameter
  epsilonT <- rnorm(Tmax, mean = 0, sd = sigmaT.dd)
  epsilonJ <- rnorm(Jmax, mean = 0, sd = sigmaJ.dd)
  
  # Calculate year- and site-specific sigma parameter
  Sigma <- exp(log(Mu.dd) + outer(epsilonJ, epsilonT, FUN = '+'))
  
  # Copy group data
  data <- G.age
  
  # Sample distance from transect line for each group (assuming uniform distribution)
  data$distance <- round(runif(nrow(data), 0, W))
  
  # Calculate group detection probabilities p based on distances d & simulate observation process
  data$p <- NA
  data$obs <- NA
  for(i in 1:nrow(data)){
    data$p[i] <- exp(-data$distance[i] * data$distance[i]/(2 * (Sigma[data$site[i], data$year[i]]^2)))
    data$obs[i] <- rbinom(n = 1, size = 1, prob = data$p[i])
  }
  
  # Retain only data from observations if "discard0"
  if(discard0){
    data <- subset(data, obs == 1)
  }
  
  # Summarise number of individuals observed per year-site combination
  #DS.count <- array(NA, nrow = Jmax, ncol = Tmax)
  DS.count <- array(NA, c(2, Jmax, Tmax))
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      data.sub <- subset(data, year == t & site == j)
      DS.count[1, j, t] <- sum(data.sub$no_juv)
      DS.count[2, j, t] <- sum(data.sub$no_ad)
    }
  }
  
  # Collate and return data
  DS.data <- list(d = data$distance, d_year = data$year, d_site = data$site,
                  DS.count = DS.count, L = L)
  return(DS.data)
}

## Simulate hierarchical distance sampling data
DS.data <- simulate.HDSdata(Jmax, Tmax, G.age = G.age, 
                            Mu.dd, sigmaT.dd, sigmaJ.dd,
                            W, discard0 = TRUE, 
                            min.Tlength, max.Tlength)

## Function to extract reproductive data from distance sampling data
extract.RepData <- function(Jmax, Tmax, DS.count){
  
  # Extract necessary data from DS counts
  sumR_obs <- c(DS.count[1,,])
  sumAd_obs <- c(DS.count[2,,])
  sumR_obs_year <- rep(1:Tmax, each = Jmax)
  
  # Discard entries with 0 adults present
  drop.idx <- which(sumAd_obs == 0)
  sumR_obs <- sumR_obs[-drop.idx]
  sumAd_obs <- sumAd_obs[-drop.idx]
  sumR_obs_year <- sumR_obs_year[-drop.idx]
  
  # Count observations
  N_sumR_obs <- length(sumR_obs)
  
  # Collate and return data
  Rep.data <- list(sumR_obs = sumR_obs, sumAd_obs = sumAd_obs,
                   sumR_obs_year = sumR_obs_year, N_sumR_obs = N_sumR_obs)
  return(Rep.data)
  
}

## Extract reproductive data
Rep.data <- extract.RepData(Jmax, Tmax, DS.count = DS.data$DS.count)


########################################
# KNOWN-FATE TELEMETRY DATA SIMULATION #
########################################

# NOTE: As of now, only one (= annual) survival interval is included in the data
#       simulation. However, extending this to conditional survival over two
#       seasonal periods (as in the real data) is straightforward. 

## Function for simulating known-fate telemetry data
simulate.RTdata <- function(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs){
  
  # Make vectors for storing data
  n.rel.S1 <- n.surv.S1 <- rep(NA, Tmax)
  n.rel.S2 <- n.surv.S2 <- rep(NA, Tmax)
  
  for(t in Tmin.RT:Tmax.RT){
    
    # Sample number of individuals fitted with transmitters each year
    n.rel.S1[t] <- rpois(1, nind.avg.RT)
    n.rel.S2[t] <- rpois(1, nind.avg.RT)
    
    # Simulate survivors to the next year
    n.surv.S1[t] <- rbinom(1, size = n.rel.S1[t], prob = sqrt(SurvProbs[t]))
    n.surv.S2[t] <- rbinom(1, size = n.rel.S2[t], prob = sqrt(SurvProbs[t]))
  }
  
  # Assemble data in a list and return
  RT.data <- list(Survs1 = cbind(n.rel.S1, n.surv.S1),
                  Survs2 = cbind(n.rel.S2, n.surv.S2))
  return(RT.data)
}

## Simulate known-fate telemetry data
RT.data <- simulate.RTdata(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs = VR.list$S[1,])


###############################
# NEST SURVEY DATA SIMULATION #
###############################

# ## Function for simulating nest survey data
# simulate.NSdata <- function(nind.avg.NS, Tmin.NS, Tmax.NS, Tmax, RepRates){
#   
#   # Determine number of nests surveyed in each year
#   n.nests <- rep(NA, Tmax)
#   for(t in Tmin.NS:Tmax.NS){
#     n.nests[t] <- rpois(1, nind.avg.NS)
#   }
#   
#   # Set up empty vectors for storing data
#   nChicks.obs <- year.obs <- vector()
#   
#   # Simulate nest survey data collection
#   for(t in Tmin.NS:Tmax.NS){
#     
#     # - Simulate number of chicks in each nest surveyed in year t
#     nChicks.obs.t <- rpois(n.nests[t], RepRates[t])
#     
#     # - Make a vector with year indices
#     year.obs.t <- rep(t, n.nests[t])
#     
#     # - Merge data of year t with data from previous years
#     nChicks.obs <- c(nChicks.obs, nChicks.obs.t)
#     year.obs <- c(year.obs, year.obs.t)
#   }
#   
#   # Assemble data in a list and return
#   NS.data <- list(nChicks.obs = nChicks.obs, year.obs = year.obs)
#   return(NS.data)
# }
# 
# ## Simulate nest survey data
# NS.data <- simulate.NSdata(nind.avg.NS, Tmin.NS, Tmax.NS, RepRates = VR.list$R[1,])


#########################
# FULL DATASET ASSEMBLY #
#########################

## Collect all simulated parameter values into a list
SimParams <- list(
  Seed = mySeed,
  Amax = Amax, Tmax = Tmax, Jmax = Jmax,
  avg_Gsize = avg_Gsize,
  N1 = N1,
  Mu.S = Mu.S, Mu.R = Mu.R, 
  sigmaT.S = sigmaT.S, sigmaT.R = sigmaT.R, 
  sigmaJ.S = sigmaJ.S, sigmaJ.R = sigmaJ.R, 
  
  min.Tlength = min.Tlength,
  max.Tlength = max.Tlength,
  W = W,
  Mu.dd = Mu.dd, sigmaT.dd = sigmaT.dd, sigmaJ.dd = sigmaJ.dd,
  
  Tmin.RT = Tmin.RT, Tmax.RT = Tmax.RT, nind.avg.RT = nind.avg.RT
  #Tmin.NS = Tmin.NS, Tmax.NS = Tmax.NS, nind.avg.NS = nind.avg.NS
)

## Collate and save all data
AllSimData <- list(
  SimParams = SimParams,
  VR.list = VR.list,
  N.data = SimData$N,
  group.data = G.age,
  DS.data = DS.data,
  Rep.data = Rep.data,
  RT.data = RT.data
)

saveRDS(AllSimData, "SimData_Full.rds")
