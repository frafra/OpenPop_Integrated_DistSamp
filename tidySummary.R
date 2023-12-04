
library(tidybayes)
library(tidyverse)



#' Calculate summary stats for selected variables
#'
#' @param out.mcmc 
#'
#' @return a tibble containg summary stats for selected variables
#' @export
#'
#' @examples

tidySummary <- function(out.mcmc){


################################################################################
### Summarize muR
muR_tidy <- out.mcmc %>% tidybayes::spread_draws(Mu.R[area], 
                                                 sep = "[,]", regex = FALSE) %>% 
  group_by() %>%
  summarise(Mean = mean(Mu.R), 
            std = sd(Mu.R), 
            Median = median(Mu.R), 
            lower = quantile(Mu.R, probs = 0.025), 
            upper = quantile(Mu.R, probs = 0.975)) %>%
            mutate(model_parameter = "Mu_R")

################################################################################
### Summarize mu.S

muS_tidy <- out.mcmc %>% tidybayes::spread_draws(Mu.S[area], 
                                                 sep = "[,]", regex = FALSE) %>% 
  group_by() %>%
  summarise(Mean = mean(Mu.S), 
            std = sd(Mu.S), 
            Median = median(Mu.S), 
            lower = quantile(Mu.S, probs = 0.025), 
            upper = quantile(Mu.S, probs = 0.975)) %>%
            mutate(model_parameter = "Mu_S")

################################################################################
### Summarize mu.S1

muS1_tidy <- out.mcmc %>% tidybayes::spread_draws(Mu.S1, 
                                                 sep = "[,]", regex = FALSE) %>% 
  group_by() %>%
  summarise(Mean = mean(Mu.S1), 
            std = sd(Mu.S1), 
            Median = median(Mu.S1), 
            lower = quantile(Mu.S1, probs = 0.025), 
            upper = quantile(Mu.S1, probs = 0.975)) %>%
  mutate(model_parameter = "Mu_S1")

################################################################################
### Summarize mu.S2

muS2_tidy <- out.mcmc %>% tidybayes::spread_draws(Mu.S[], Mu.S1, 
                                                 sep = "[,]", regex = FALSE) %>% 
  mutate(Mu.S2 = Mu.S/Mu.S1) %>%
  group_by() %>%
  summarise(Mean = mean(Mu.S2), 
            std = sd(Mu.S2), 
            Median = median(Mu.S2), 
            lower = quantile(Mu.S2, probs = 0.025), 
            upper = quantile(Mu.S2, probs = 0.975))%>%
    mutate(model_parameter = "Mu_S2")

################################################################################
### Summarize density 
dens_tidy <- out.mcmc  %>% tidybayes::spread_draws(Density[Area, age, lineNR, year], 
                                                   sep = "[,]", regex = FALSE) %>% 
  group_by(Area, lineNR, year, .chain, .iteration, .draw) %>%
  summarise(sum_dens = sum(Density)) %>%
   ungroup() %>%
   group_by(year, .draw) %>%
   summarise(mean_it_density = mean(sum_dens)) %>%
   ungroup() %>%
  group_by(year) %>%
  summarise(Mean = mean(mean_it_density), 
            std = sd(mean_it_density), 
            Median = median(mean_it_density), 
            lower = quantile(mean_it_density, probs = 0.025), 
            upper = quantile(mean_it_density, probs = 0.975))%>%
   ungroup() %>%
    mutate(model_parameter = "Density") 
 
  

################################################################################
### Summarize R_year

R_year_tidy <- out.mcmc %>% tidybayes::spread_draws(R_year[area, year], sep = "[,]", regex = FALSE) %>%
  group_by(year) %>%
  summarise(Mean = mean(R_year), 
            std = sd(R_year), 
            Median = median(R_year), 
            lower = quantile(R_year, probs = 0.025), 
            upper = quantile(R_year, probs = 0.975)) %>%
  ungroup() %>%
  mutate(model_parameter = "R_year")

################################################################################
##### Effect of rodent covariate betaR.R

betaR.R_tidy <- out.mcmc %>% tidybayes::spread_draws(betaR.R[area], 
                                                     sep = "[,]", regex = FALSE) %>% 
  group_by() %>%
  summarise(Mean = mean(betaR.R), 
            std = sd(betaR.R), 
            Median = median(betaR.R), 
            lower = quantile(betaR.R, probs = 0.025), 
            upper = quantile(betaR.R, probs = 0.975)) %>%
  mutate(model_parameter = "betaR.R")


################################################################################
#### Putting together output

pooled_results <- bind_rows(dens_tidy, R_year_tidy, muR_tidy, 
                            muS_tidy, muS1_tidy, muS2_tidy, betaR.R_tidy)


pooled_results

}

