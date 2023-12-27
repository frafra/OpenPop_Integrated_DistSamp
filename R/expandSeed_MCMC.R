expandSeed_MCMC <- function(seed, nchains){
  
  if(nchains > 6){
    stop("Seed expansion is currently only supported for up to 6 chains.")
  }
  
  seed.list <- c(
    seed,
    round((seed*70)/4),
    round(seed/3)+1,
    seed + 88,
    round(seed*5/3),
    seed + 37
  )
  
  return(seed.list[1:nchains])
}