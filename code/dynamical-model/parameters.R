#### outputs ####

# notation:
# parameter value lowercase
# species group uppercase (A = annual, P = perennial adult, S = perennial first-year)
# disease condition underscore (_d = disease/ambient, _h = healthy/suppressed)


#### set-up ####

# load packages
library(tidyverse)


#### values from literature ####

sA <- 0.15 # annual seed survival (Redwood et al. 2018)
sP <- 0.05 # perennial seed survival (Garrison and Stier 2010)
d_d <- 0.50 # litter decomposition infected sites (Lane et al. in review)
d_h <- 0.47 # litter decomposition non-infected sites (Lane et al. in review)
bT <- 0 # set other sources of litter to zero
delta <- 1 # set fraction of perennial biomass that affects establishment as 1


#### import posterior draws ####

pS_draws <- read_csv("intermediate-data/pS_draws.csv")
pP_draws <- read_csv("intermediate-data/pP_draws.csv")
yA_draws <- read_csv("intermediate-data/yA_draws.csv")
yP_draws <- read_csv("intermediate-data/yP_draws.csv")
f_draws <- read_csv("intermediate-data/f_draws.csv")
alphaA_draws <- read_csv("intermediate-data/alphaA_draws.csv")
alphaP_draws <- read_csv("intermediate-data/alphaP_draws.csv")
gamma_draws <- read_csv("intermediate-data/gamma_draws.csv")
gA_fung_draws <- read_csv("intermediate-data/gA_fung_draws.csv")
gA_inf_draws <- read_csv("intermediate-data/gA_inf_draws.csv")
gP_draws <- read_csv("intermediate-data/gP_draws.csv")
eA_draws <- read_csv("intermediate-data/eA_draws.csv")
eP_draws <- read_csv("intermediate-data/eP_draws.csv")
betaA_draws <- read_csv("intermediate-data/betaA_draws.csv")
betaP_draws <- read_csv("intermediate-data/betaP_draws.csv")
bA_draws <- read_csv("intermediate-data/bA_draws.csv")
bP_draws <- read_csv("intermediate-data/bP_draws.csv")


#### parameter function ####

# pull parameter values based on number of iterations
# gA_type can be fungicide or infection
params_fun <- function(iters, gA_type = "fungicide"){
  
  # iterations to take from each posterior
  i <- sample(1:15000, iters, replace = FALSE)
  
  # perennial seedling interannual survival
  pS_h <- pS_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  pS_d <- pS_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial adult interannual survival
  pP_h <- pP_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  pP_d <- pP_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # combine survival
  s_h <- tibble(sA, sP, pS = pS_h, pP = pP_h)
  s_d <- tibble(sA, sP, pS = pS_d, pP = pP_d)
  
  # annual seed yield without competition
  yA_h <- yA_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  yA_d <- yA_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial seed yield without competition
  yP_h <- yP_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  yP_d <- yP_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial seed yield ratio
  f_h <- f_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  f_d <- f_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # combine seed yield
  y_h <- tibble(yA = yA_h, yP = yP_h, f = f_h)
  y_d <- tibble(yA = yA_d, yP = yP_d, f = f_d)
  
  # annual germination
  if(gA_type == "fungicide"){
    
    gA_h <- gA_fung_draws %>%
      filter(fungicide == 1 & .draw %in% i) %>%
      pull(value)
    gA_d <- gA_fung_draws %>%
      filter(fungicide == 0 & .draw %in% i) %>%
      pull(value)
    
  }else if(gA_type == "infection"){
    
    gA_h <- gA_inf_draws %>%
      filter(infection == 0 & .draw %in% i) %>%
      pull(value)
    gA_d <- gA_inf_draws %>%
      filter(infection == 1 & .draw %in% i) %>%
      pull(value)
    
  }else{
    
    stop("gA_type must be 'fungicide' or 'infection'")
    
  }
  
  # perennial germination
  gP_h <- gP_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  gP_d <- gP_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # combine germination
  g_h <- tibble(gA = gA_h, gP = gP_h)
  g_d <- tibble(gA = gA_d, gP = gP_d)
  
  # annual establishment without litter
  eA <- eA_draws %>%
    filter(.draw %in% i) %>%
    mutate(value = if_else(value > 1, 1, value)) %>%  # cap at 1
    pull(value)
  
  # perennial establishment without litter
  eP <- eP_draws %>%
    filter(.draw %in% i) %>%
    pull(value)
  
  # combine establishment
  e <- tibble(eA, eP)
  
  # annual biomass
  bA_h <- bA_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  bA_d <- bA_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial biomass
  bP_h <- bP_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  bP_d <- bP_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # combine decay
  decay_h <- tibble(bA = bA_h, bP = bP_h, d = d_h, bT, delta)
  decay_d <- tibble(bA = bA_d, bP = bP_d, d = d_d, bT, delta)
  
  # annual effect on seed yield 
  alphaA_h <- alphaA_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  alphaA_d <- alphaA_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)

  # perennial adult effect on seed yield
  alphaP_h <- alphaP_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  alphaP_d <- alphaP_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial effect ratio
  gamma_h <- gamma_draws %>%
    filter(fungicide == 1 & .draw %in% i) %>%
    pull(value)
  gamma_d <- gamma_draws %>%
    filter(fungicide == 0 & .draw %in% i) %>%
    pull(value)
  
  # combine alpha
  alpha_h <- tibble(alphaA = alphaA_h, alphaP = alphaP_h, gamma = gamma_h)
  alpha_d <- tibble(alphaA = alphaA_d, alphaP = alphaP_d, gamma = gamma_d)
  
  # annual response to litter 
  betaA_h <- betaA_draws %>%
    filter(sterilized == 1 & .draw %in% i) %>%
    pull(value)
  betaA_d <- betaA_draws %>%
    filter(sterilized == 0 & .draw %in% i) %>%
    pull(value)
  
  # perennial response to litter
  betaP <- betaP_draws %>%
    filter(.draw %in% i) %>%
    pull(value)

  # combine beta
  beta_h <- tibble(betaA = betaA_h, betaP)
  beta_d <- tibble(betaA = betaA_d, betaP)
    
  # combine all parameters
  parameters_h <- list(s = s_h, y = y_h, g = g_h, e = e, decay = decay_h, 
                       alpha = alpha_h, beta = beta_h)
  parameters_d <- list(s = s_d, y = y_d, g = g_d, e = e, decay = decay_d, 
                       alpha = alpha_d, beta = beta_d)
  
  # return list of parameters
  return(list(parameters_h, parameters_d))
}

test <- params_fun(100, gA_type = "fungicide")
parameters <- test[[1]]  # healthy parameters
