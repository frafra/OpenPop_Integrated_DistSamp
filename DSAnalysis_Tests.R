set.seed(0)

## Function for sourcing lines in a file (from https://gist.github.com/christophergandrud/1eb4e095974204b12af9)
source_lines <- function(file, lines){
  source(textConnection(readLines(file)[lines]))
}

## Source workflow without model fitting
source_lines(file = "Analysis_RealData.R", lines = c(1:70))

## Extract relevant data
y <- input_data$nim.data$y
nind <- length(y)
W <- input_data$nim.constants$W
zeros.dist <- input_data$nim.data$zeros.dist

## Augment data
nz <- 1000 # Augment observed data with nz = 1000 zeroes
aug <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y = 0 by definition
y.aug <- c(y, rep(NA, nz)) # Value of distance are missing for the augmented

## Bundle data
data.orig <-  list(y = y, nind = nind, W = W, zeros.dist = zeros.dist, pi = 3.141593)
data.aug <-  list(y = y.aug, nind = nind, W = W, nz = nz, aug = aug)

## Set initial values
inits.orig <- list(list(sigma = runif(1, 0, 50)), list(sigma = runif(1, 0, 50)), list(sigma = runif(1, 0, 50)))
inits.aug <- inits.orig
for(i in 1:3){
  inits.aug[[i]]$psi <- runif(1, 0, 1) 
  inits.aug[[i]]$z <- aug
}

## Set monitors
params <- c("sigma")


# Frequentist analysis using Distance package #
#---------------------------------------------#
library(Distance)

Distance.fit <- ds(data = y, key = "hn", adjustment = NULL)


# Bayesian analysis using zeros-trick (original model) #
#------------------------------------------------------#
library(nimble)

## Nimble code
Code_orig <- nimbleCode({
  
  for (i in 1:nind){ 
    #y[i] ~ dunif(0, W) # NOTE: this line is not actually necessary.
    L.f0[i] <- exp(-y[i]*y[i] / (2*sigma2)) * 1/esw
    nlogL.f0[i] <- -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
  }

  sigma2 <- sigma * sigma
  esw <- sqrt(pi * sigma2 / 2) 
  sigma ~ dunif(0, 200)
})

## Fit model
nim.fit <- nimbleMCMC(code = Code_orig, 
                      constants = data.orig, 
                      inits = inits.orig,
                      monitors = params,
                      niter = 11000,
                      nburnin = 1000,
                      nchains = 3,
                      thin = 2,
                      setSeed = 0,
                      samplesAsCodaMCMC = TRUE)


# Bayesian analysis using nimbleDistance::dHN #
#---------------------------------------------#
library(nimbleDistance)

## Nimble code
Code_nimDist <- nimbleCode({
  
  #y[1:nind] ~ dHN_V(sigma = sigma, Xmax = W, point = 0)
  
  for(i in 1:nind){
    y[i] ~ dHN(sigma = sigma, Xmax = W, point = 0)
  }
  
  sigma ~ dunif(0, 200)
  
})

## Fit model
nimDist.fit <- nimbleMCMC(code = Code_nimDist, 
                          constants = data.orig, 
                          inits = inits.orig,
                          monitors = params,
                          niter = 11000,
                          nburnin = 1000,
                          nchains = 3,
                          thin = 2,
                          setSeed = 0,
                          samplesAsCodaMCMC = TRUE)


# Bayesian analysis using data augmentation #
#-------------------------------------------#

## Nimble code (from AHM 1 chapter 8.3.1)
Code_nimDA <- nimbleCode({
  
  # Priors
  sigma ~ dunif(0, 200)  # Half-normal scale
  psi ~ dunif(0, 1)       # DA parameter
  
  # Likelihood
  for(i in 1:(nind + nz)){
    
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    y[i] ~ dunif(0, W)  # Distribution of distances
    
    # Observation model
    logp[i] <- -((y[i]*y[i])/(2*sigma*sigma)) # Half-normal detection fct.
    p[i] <- exp(logp[i])
    mu[i] <- z[i] * p[i]
    aug[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
  }
})

## Fit model
nimDA.fit <- nimbleMCMC(code = Code_nimDA, 
                        constants = data.aug, 
                        inits = inits.aug,
                        monitors = params,
                        niter = 11000,
                        nburnin = 1000,
                        nchains = 3,
                        thin = 2,
                        setSeed = 0,
                        samplesAsCodaMCMC = TRUE)


# Compare fits #
#--------------#
library(ggplot2)
library(viridis)

## Collate data
nim.out <- data.frame(
  sigma = c(as.matrix(nim.fit), as.matrix(nimDist.fit), as.matrix(nimDA.fit)),
  Fit = rep(c("Zeros trick", "nimbleDistance::dHN", "Data augmentation"), each = 15000)
)

## Plot results
ggplot(nim.out) + 
  geom_density(aes(x = sigma, color = Fit, fill = Fit), alpha = 0.5) + 
  geom_vline(aes(xintercept = exp(Distance.fit$ddf$par)), linetype = "dashed") + # Note: Distance uses log scale for parameters
  scale_fill_viridis(discrete = T) + 
  scale_color_viridis(discrete = T) + 
  theme_classic()
