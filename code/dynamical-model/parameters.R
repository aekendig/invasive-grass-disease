#### set-up ####

# load packages
library(tidyverse)
library(brms)
library(tidybayes)


#### values from literature ####

sA <- 0.15 # annual seed survival (Redwood et al. 2018)
sP <- 0.05 # perennial seed survival (Garrison and Stier 2010)
# d <- 0.61 # litter decomposition (DeMeester and Richter 2010) - look at Brett's paper for this
# h <- 0.29 # seedling survival from germination to establishment (Emery et al. 2013) - probably don't need, use germination/establishment from litter field experiment


#### posterior draws ####

# import models
load("output/evS_survival_model_2018_2019_density_exp.rda")
load("output/evA_survival_model_2018_2019_density_exp.rda")
load("output/evS_background_seed_model_2019_density_exp.rda")
load("output/evA_background_seed_model_2019_density_exp.rda")
load("output/mv_background_seed_model_2019_density_exp.rda")
load("output/ev_seed_model_2019_density_exp.rda") # for f
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
load("output/mv_germination_infection_model_2018_density_exp.rda")
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")
load("output/evS_establishment_model_2018_2019_density_exp.rda") # not from seed, combine with below
load("output/mv_establishment_model_2018_2019_density_exp.rda") # not from seed, combine with below
load("output/ev_establishment_model_2019_litter_exp.rda") # field germination/establishment without litter
load("output/mv_establishment_model_2018_litter_exp.rda") # field germination/establishment without litter, live litter effect
load("output/mv_establishment_model_2018_greenhouse_exp.rda") # for litter sensitivity
load("output/ev_establishment_model_2018_greenhouse_exp.rda") # for litter sensitivity

# no background plot data
# pred_dat_no_bg <- tibble(fungicide = 0, background_density = 0,
#                          background = "Mv seedling") %>%
#   add_row(fungicide = 1, background_density = 0,
#           background = "Mv seedling")

pred_dat_trt <- tibble(fungicide = c(0, 1))

# posterior draws on response scale
pred_pS <- pred_dat_trt %>% 
  add_epred_draws(evSSurvMod, re_formula = NA) %>% ungroup()
pred_pA <- pred_dat_trt %>% 
  add_epred_draws(evASurvMod, re_formula = NA) %>% ungroup()




evSSeedDraws <- as_draws_df(evSSeedMod)
evASeedDraws <- as_draws_df(evASeedMod)
mvSeedDraws <- as_draws_df(mvSeedMod)
mvGermD1Draws <- as_draws_df(mvGermD1Mod)
mvGermInfD1Draws <- as_draws_df(mvGermInfD1Mod)
evGermDraws <- as_draws_df(evGermMod)
evSEstDraws <- as_draws_df(evSEstMod)
mvEstDraws <- as_draws_df(mvEstMod)
evEstL2Draws <- as_draws_df(evEstL2Mod)
mvEstL1Draws <- as_draws_df(mvEstL1Mod)
mvEstGhDraws <- as_draws_df(mvEstGhMod)
evEstGhDraws <- as_draws_df(evEstGhMod)


#### parameter functions ####

# parameters with control conditions
params_fun <- function(iters){
  
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