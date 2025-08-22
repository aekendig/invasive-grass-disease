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
alphaAA_draws <- read_csv("intermediate-data/alphaAA_draws.csv")
alphaAP_draws <- read_csv("intermediate-data/alphaAP_draws.csv")
alphaAS_draws <- read_csv("intermediate-data/alphaAS_draws.csv")
alphaPA_draws <- read_csv("intermediate-data/alphaPA_draws.csv")
alphaPP_draws <- read_csv("intermediate-data/alphaPP_draws.csv")
alphaPS_draws <- read_csv("intermediate-data/alphaPS_draws.csv")
gA_fung_draws <- read_csv("intermediate-data/gA_fung_draws.csv")
gA_inf_draws <- read_csv("intermediate-data/gA_inf_draws.csv")
gP_draws <- read_csv("intermediate-data/gP_draws.csv")
eA_draws <- read_csv("intermediate-data/eA_draws.csv")
eP_draws <- read_csv("intermediate-data/eP_draws.csv")
betaA_draws <- read_csv("intermediate-data/betaA_draws.csv")
betaP_draws <- read_csv("intermediate-data/betaP_draws.csv")
bA_draws <- read_csv("intermediate-data/bA_draws.csv")
bP_draws <- read_csv("intermediate-data/bP_draws.csv")
yA_drawsP <- read_csv("intermediate-data/yA_draws_Stricker_priors.csv")
alphaAA_drawsP <- read_csv("intermediate-data/alphaAA_draws_Stricker_priors.csv")
alphaAP_drawsP <- read_csv("intermediate-data/alphaAP_draws_Stricker_priors.csv")
alphaAS_drawsP <- read_csv("intermediate-data/alphaAS_draws_Stricker_priors.csv")
bA_drawsP <- read_csv("intermediate-data/bA_draws_Stricker_priors.csv")


#### parameter function ####

# pull parameter values based on number of iterations
# gA_type can be fungicide or infection
params_fun <- function(iters, gA_type = "fungicide", draws = NULL,
                       priors = "baseline"){
  
  # iterations to take from each posterior
  if(is.null(draws)){
    
    draws <- tibble(.draw = sample(1:15000, iters, replace = FALSE))
    
  }
  
  # parameters affected by priors
  if(priors == "baseline"){
    
    # annual seed yield without competition
    yA_h <- yA_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    yA_d <- yA_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # annual biomass
    bA_h <- bA_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    bA_d <- bA_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAA
    alphaAA_h <- alphaAA_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAA_d <- alphaAA_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAP
    alphaAP_h <- alphaAP_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAP_d <- alphaAP_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAS
    alphaAS_h <- alphaAS_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAS_d <- alphaAS_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
  }else if(priors == "Stricker"){
    
    # annual seed yield without competition
    yA_h <- yA_drawsP %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    yA_d <- yA_drawsP %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # annual biomass
    bA_h <- bA_drawsP %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    bA_d <- bA_drawsP %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAA
    alphaAA_h <- alphaAA_drawsP %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAA_d <- alphaAA_drawsP %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAP
    alphaAP_h <- alphaAP_drawsP %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAP_d <- alphaAP_drawsP %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
    # alphaAS
    alphaAS_h <- alphaAS_drawsP %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    alphaAS_d <- alphaAS_drawsP %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
  }else{
    
    stop("priors must be 'baseline' or 'Stricker'")
    
  }
  
  # perennial seedling interannual survival
  pS_h <- pS_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  pS_d <- pS_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # perennial adult interannual survival
  pP_h <- pP_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  pP_d <- pP_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine survival
  s_h <- tibble(sA, sP, pS = pS_h, pP = pP_h)
  s_d <- tibble(sA, sP, pS = pS_d, pP = pP_d)
  
  # perennial seed yield without competition
  yP_h <- yP_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  yP_d <- yP_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # perennial seed yield ratio
  f_h <- f_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  f_d <- f_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine seed yield
  y_h <- tibble(yA = yA_h, yP = yP_h, f = f_h)
  y_d <- tibble(yA = yA_d, yP = yP_d, f = f_d)
  
  # annual germination
  if(gA_type == "fungicide"){
    
    gA_h <- gA_fung_draws %>%
      filter(fungicide == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    gA_d <- gA_fung_draws %>%
      filter(fungicide == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
  }else if(gA_type == "infection"){
    
    gA_h <- gA_inf_draws %>%
      filter(infection == 0) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    gA_d <- gA_inf_draws %>%
      filter(infection == 1) %>%
      inner_join(draws) %>%
      arrange(match(.draw, draws$.draw)) %>%
      pull(value)
    
  }else{
    
    stop("gA_type must be 'fungicide' or 'infection'")
    
  }
  
  # perennial germination
  gP_h <- gP_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  gP_d <- gP_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine germination
  g_h <- tibble(gA = gA_h, gP = gP_h)
  g_d <- tibble(gA = gA_d, gP = gP_d)
  
  # annual establishment without litter
  eA <- eA_draws %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    mutate(value = if_else(value > 1, 1, value),
           value = if_else(value < 0, 0, value)) %>% 
    pull(value)
  
  # perennial establishment without litter
  eP <- eP_draws %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    mutate(value = if_else(value > 1, 1, value),
           value = if_else(value < 0, 0, value)) %>% 
    pull(value)
  
  # combine establishment
  e <- tibble(eA, eP)
  
  # perennial biomass
  bP_h <- bP_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  bP_d <- bP_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine decay
  decay_h <- tibble(bA = bA_h, bP = bP_h, d = d_h, bT, delta)
  decay_d <- tibble(bA = bA_d, bP = bP_d, d = d_d, bT, delta)
  
  # alphaPA
  alphaPA_h <- alphaPA_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  alphaPA_d <- alphaPA_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # alphaPP
  alphaPP_h <- alphaPP_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  alphaPP_d <- alphaPP_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # alphaPS
  alphaPS_h <- alphaPS_draws %>%
    filter(fungicide == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  alphaPS_d <- alphaPS_draws %>%
    filter(fungicide == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine alpha
  alpha_h <- tibble(alphaAA = alphaAA_h, alphaAP = alphaAP_h, 
                    alphaAS = alphaAS_h, alphaPA = alphaPA_h, 
                    alphaPP = alphaPP_h, alphaPS = alphaPS_h,)
  alpha_d <- tibble(alphaAA = alphaAA_d, alphaAP = alphaAP_d, 
                    alphaAS = alphaAS_d, alphaPA = alphaPA_d, 
                    alphaPP = alphaPP_d, alphaPS = alphaPS_d,)
  
  # annual response to litter 
  betaA_h <- betaA_draws %>%
    filter(sterilized == 1) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  betaA_d <- betaA_draws %>%
    filter(sterilized == 0) %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # perennial response to litter
  betaP <- betaP_draws %>%
    inner_join(draws) %>%
    arrange(match(.draw, draws$.draw)) %>%
    pull(value)
  
  # combine beta
  beta_h <- tibble(betaA = betaA_h, betaP)
  beta_d <- tibble(betaA = betaA_d, betaP)
    
  # combine all parameters
  parameters_h <- list(s = s_h, y = y_h, g = g_h, e = e, decay = decay_h, 
                       alpha = alpha_h, beta = beta_h, draws = draws)
  parameters_d <- list(s = s_d, y = y_d, g = g_d, e = e, decay = decay_d, 
                       alpha = alpha_d, beta = beta_d, draws = draws)
  
  # return list of parameters
  return(list(healthy = parameters_h, disease = parameters_d))
}
