#### set-up ####

# load packages
library(tidyverse)
library(brms)
library(tidybayes)

# notation
# parameter value lowercase
# species group uppercase (A = annual, P = perennial adult, S = perennial first-year)
# disease condition underscore (_d = disease/ambient, _h = healthy/suppressed)


#### values from literature ####

sA <- 0.15 # annual seed survival (Redwood et al. 2018)
sP <- 0.05 # perennial seed survival (Garrison and Stier 2010)
d_d <- 0.50 # litter decomposition infected sites (Lane et al. in review)
d_h <- 0.47 # litter decomposition non-infected sites (Lane et al. in review)
b_T <- 0 # set other sources of litter to zero
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





#### parameter functions ####

# parameters with control conditions
### to do: gA_type uses either fungicide or infection ####
params_fun <- function(iters, gA_type){
  
  # iterations to take from each posterior
  i <- sample(1:15000, iters, replace = FALSE)
  
  # perennial seedling interannual survival
  pS_ctrl <- pred_pS %>%
    filter(fungicide == 0 & .draw %in% i)
  pS_fung <- pred_pS %>%
    filter(fungicide == 1 & .draw %in% i)
  
  # perennial adult interannual survival
  pA_ctrl <- pred_pA %>%
    filter(fungicide == 0 & .draw %in% i)
  pA_fung <- pred_pA %>%
    filter(fungicide == 1 & .draw %in% i)
  
 
  

  
  # seed yield
  yA <- mvSeedDraws[iter, ] %>% pull(b_Intercept)
  yP <- evASeedDraws[iter, ] %>% pull(b_Intercept)
  y <- c(yA, yP, f)
  
  g <- parameters[[3]]  # [gA, gP]
  gA <- g[1]; gP <- g[2]
  
  
  
  
  
  
  e <- parameters[[4]]  # [eA, eP]
  
  decay <- parameters[[5]]  # [bA, bP, d, bT, delta]
  bA <- decay[1]; bP <- decay[2]; d <- decay[3]; bT <- decay[4]; delta <- decay[5]
  
  alpha <- parameters[[6]]  # [alphaA, alphaP, gamma]
  alphaA <- alpha[1]; alphaP <- alpha[2]; gamma <- alpha[3]
  
  beta <- parameters[[7]] 
  
}

# survival
s <- c(sA, sP, pS, pP)