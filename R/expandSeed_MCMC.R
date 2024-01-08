expandSeed_MCMC <- function(seed, nchains){
  
  if(nchains > 5){
    stop("Seed expansion is currently only supported for up to 5 chains.")
  }
  
  seed.list <- c(
    seed,
    round((seed*70)/4),
    (seed + 10)*7,
    round(seed*5/3),
    seed + 37*2
  )
  
  return(seed.list[1:nchains])
}