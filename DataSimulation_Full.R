
mySeed <- 0
set.seed(mySeed)


#########################
# SIMULATION PARAMETERS #
#########################


# General simulation parameters #
#-------------------------------#

Amax <- 2 # Number of age classes
Tmax <- 15 # Number of years
Jmax <- 50 # Number of sites/transect lines

# Population parameters #
#-----------------------#

# Initial population numbers per site
N1 <- matrix(NA, ncol = Amax, nrow = Jmax)
for(j in 1:Jmax){
  N1[j,1] <- round(runif(1, 100, 300)) # Number of juveniles
  N1[j,2] <- rpois(1, 0.5*N1[j,1]) # Number of adults
}


# Vital rate parameters #
#-----------------------#

## Annual survival
Mu.S <- 0.4 # Average annual survival probability
sigmaT.S <- 0 # SD of random year variation in survival
sigmaJ.S <- 0 # SD of random site variation in survival

## Reproduction
Mu.R <- 3 # Average number of chicks in August
sigmaT.R <- 0.5 # SD of random year variation in number of chicks
sigmaJ.R <- 0 # SD of random site variation in number of chicks

## Juvenile summer survival
#Mu.sJ <- 0.2 # Average summer survival of chicks
#sigmaT.sJ <- 0 # SD of random year variation in chick survival
#sigmaJ.sJ <- 0 # SD of random site variation in survival


# Data & observation parameters #
#-------------------------------#

## Line-transect distance sampling
min.Tlength <- 1000 # Minimum transect length
max.Tlength <- 1000  # Maximum transect length

W <- 200 # Truncation distance (max. distance at which observation is possible)
  
#pi <- 3.141593 # Approximation of pi

Mu.dd <- 75 # Average width parameter for half-normal detection function
sigmaT.dd <- 0 # SD of random year variation in detection probability
sigmaJ.dd <- 0 # SD of random year variation in detection probability

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
        N[j,1,t+1] <- rpois(1, lambda = N[j,2,t+1]*(VR.list$R[j,t+1])/2)
      
        # - Young-of-the-year recruiting in the end of summer
       # N[j,1,t+1] <- rbinom(1, size = Chicks[j,t+1], prob = VR.list$sJ[j,t+1])
    
    
      }else{
    
        # - Survivors
        N[j,2,t+1] <- sum(N[j,1:Amax,t])*VR.list$S[j,t]
      
        # - Reproduction of survivors
        N[j,1,t+1] <- N[j,2,t+1]*(VR.list$R[j,t+1]/2)
      
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


#################################
# LINE TRANSECT DATA SIMULATION #
#################################

# NOTE: As of now, the simulated DS data assumes independent observation of 
#       single animals, i.e. each observation = 1 individual. 
#       In the real ptarmigan data, an observation may pertain to single animals
#       or groups of multiple animals. Extension of the data simulation to 
#       accommodate group sampling will have to be made later. 

## Function for simulating hierarchical distance sampling data
simulate.HDSdata <- function(Jmax, Tmax, N.age,
                             Mu.dd, sigmaT.dd, sigmaJ.dd,
                             w, discard0 = TRUE, min.Tlength, max.Tlength){
  
  # Set site-specific transect lengths
  # NOTE: These are not currently influencing the data simulation in any way
  L <- matrix(runif(Jmax*Tmax, min.Tlength, max.Tlength), nrow = Jmax, ncol = Tmax)
  
  # Sample random year and site effects on sigma parameter
  epsilonT <- rnorm(Tmax, mean = 0, sd = sigmaT.dd)
  epsilonJ <- rnorm(Jmax, mean = 0, sd = sigmaJ.dd)
  
  # Calculate year- and site-specific sigma parameter
  Sigma <- exp(log(Mu.dd) + outer(epsilonJ, epsilonT, FUN = '+'))
  
  # Sum up true population sizes (adults + YOY) per year-site combination
  N.tot <- apply(N.age, c(1,3), sum)
  
  
  # Simulate line transect data (observation distances)
  data <- data.frame()
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      
      if(N.tot[j,t] == 0){
        data <- rbind(data, c(j, NA, NA, NA))
        next
      }
      
      # - Sample distance from transect line for each individual (assuming uniform distribution)
      #d <- runif(N.tot[j,t], 0, W) 
      d_juv <- runif(N.age[j, 1, t], 0, W)
      d_ad <- runif(N.age[j, 2, t], 0, W) 
      d <- as.matrix(cbind(c(d_ad, d_juv), c(rep(1, length(d_juv)), rep(2, length(d_ad)))))
      
      
      # - Calculate individual detection probabilities p based on distances d
      p <- exp(-d * d/(2 * (Sigma[j,t]^2)))
      
      # - Simulate observation process (whether each individual was observed or not)
      #y <- rbinom(n = N.tot[j,t], size = 1, prob = p)
      y_juv <- rbinom(n = N.age[j, 1, t], size = 1, prob = p)
      y_ad <- rbinom(n = N.age[j, 2, t], size = 1, prob = p)
      y <- as.matrix(cbind(c(y_ad, y_juv), c(rep(1, length(y_juv)), rep(2, length(y_ad)))))
      
      
      # - Retain only data from observations if "discard0"
      if(discard0){
        d <- d[y[,1] == 1,]
        # d_juv <- d_juv[y == 1]
        # d_ad <- d_ad[y == 1]
        y <- y[y[,1] == 1,]
        # y_juv <- y_juv[y == 1]
        # y_ad <- y_ad[y == 1]
        
      }
      
      
      # - Arrange data
      if(sum(y[,1]) > 0){
        data <- rbind(data, cbind(rep(j, sum(y[,1])), rep(t, sum(y[,1])), y[,1], d[,1], d[,2]))
      }else{
        data <- rbind(data, c(j, NA, NA, NA))
      } 
    }
  }
  colnames(data) <- c("site", "year", "obs", "distance", "age")
  
  # Summarise number of individuals observed per year-site combination
  #DS.count <- array(NA, nrow = Jmax, ncol = Tmax)
  DS.count <- array(NA, c(Jmax, 2, Tmax))
  for(t in 1:Tmax){
    for(j in 1:Jmax){
      for(a in 1:2){
        data.sub <- subset(data, year == t & site == j & age == a)
        DS.count[j, a, t] <- nrow(data.sub)
      }
    }
  }
  
  # Collate and return data
  DS.data <- list(d = data$d, d_year = data$year, d_site = data$site, d_age=data$age, 
                  DS.count = DS.count, L = L)
  return(DS.data)
}

## Simulate hierarchical distance sampling data
DS.data <- simulate.HDSdata(Jmax, Tmax, N.age = SimData$N, 
                            Mu.dd, sigmaT.dd, sigmaJ.dd,
                            W, discard0 = TRUE, 
                            min.Tlength, max.Tlength)



########################################
# KNOWN-FATE TELEMETRY DATA SIMULATION #
########################################

# NOTE: As of now, only one (= annual) survival interval is included in the data
#       simulation. However, extending this to conditional survival over two
#       seasonal periods (as in the real data) is straightforward. 

## Function for simulating known-fate telemetry data
simulate.RTdata <- function(nind.avg.RT, Tmin.RT, Tmax.RT, Tmax, SurvProbs){
  
  # Make vectors for storing data
  n.rel <- n.surv <- rep(NA, Tmax)
  
  for(t in Tmin.RT:Tmax.RT){
    
    # Sample number of individuals fitted with transmitters each year
    n.rel[t] <- rpois(1, nind.avg.RT)
    
    # Simulate survivors to the next year
    n.surv[t] <- rbinom(1, size = n.rel[t], prob = SurvProbs[t])
  }
  
  # Assemble data in a list and return
  RT.data <- list(n.rel = n.rel, n.surv = n.surv)
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
  DS.data = DS.data,
  RT.data = RT.data
)

