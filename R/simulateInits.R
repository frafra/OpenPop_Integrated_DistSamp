simulateInits <- function(nim.data, nim.constants){
  
  mu.D1 <- runif(1, 3, 4)
  ratio.JA1 <- runif(1, 0.2, 0.6)
  
  N_exp <- array(NA, dim = c(nim.constants$N_ageC, nim.constants$N_sites, nim.constants$N_years))
  
  for(j in 1:nim.constants$N_sites){
    N_exp1 <- rpois(1, mu.D1*nim.data$L[j, 1]*(nim.constants$W/nim.constants$scale1)*2) ## Expected number of birds
    
    N_exp[1, j, 1] <- round(N_exp1*ratio.JA1)     
    N_exp[2, j, 1] <- N_exp1 - N_exp[1, j, 1]
  }
  
  list(
    mu.dd = runif(1, 4, 5), 
    sigma.dd = runif(1, 0.05, 2),
    #b = runif(1, 1, 50), 
    mu.D1 = mu.D1, 
    sigma.D = runif(1, 0.05, 2),
    Mu.R = runif(1, 1, 3), 
    sigma.R = runif(1, 0.05, 2),
    eps.dd = rep(0, nim.constants$N_years), 
    eps.R = rep(0, nim.constants$N_years), 
    eps.D1 = rep(0, nim.constants$N_sites),
    Mu.S1 = runif(1, 0.6, 0.7), 
    Mu.S2 = runif(1, 0.6, 0.7),
    N_exp = N_exp,
    ratio.JA1 = ratio.JA1
  )
}
