#### status ####

# revised the response function
# tried it with gfung h and looks okay
# revise all the parameter sets and disease conditions in the response section
# make figures

# perhaps we put the effects and responses in one figure with the points and errors and forget the distributions?
# their tails are very long, so it's hard to see the majority of the distribution

# make a new figure with how the growth rates and outcomes of competition relate to effects and responses
# I think Nick was working on this at the end of the pdf and we may have to do something to the 
# competition-related factors for this figure/analysis, but okay to draft

#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidybayes)
library(ggtext)
library(patchwork)
library(janitor)
# scales packaged used within plotting functions below


#### parameters ####

# import parameters (loads tidyverse)
source("code/dynamical-model/parameters.R")

# parameter combinations
params_iters <- 15000

# parameters with fungicide effect on germination
params_gfung <- params_fun(iters = params_iters, gA_type = "fungicide")

# parameters with infection effect on germination
params_ginf <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws)

# parameters with infection effect and priors
params_infP <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws,
                          priors = "Stricker")

# remove raw parameter objects
rm(list = setdiff(ls(), c("params_iters", "params_gfung", "params_ginf",
                         "params_infP")))
gc()

# generations
gens <- 1000

# threshold for population present
pres_thresh <- 1e-2


#### load other scripts ####

# import simulation function
source("code/dynamical-model/kortessis_etal_2022_revised_model.R")

# figure settings
source("code/figure-prep/figure_settings.R")


#### functions ####

# combine healthy and disease simulations
sims_comb_fun <- function(sims_h, sims_d, widen = T){
  
  # combine simulations
  out_long <- bind_rows(sims_h, .id = "iteration") %>%
    mutate(disease = "h") %>%
    full_join(bind_rows(sims_d, .id = "iteration") %>%
                mutate(disease = "d")) %>%
    mutate(iteration = as.numeric(iteration),
           disease = fct_relevel(disease, "h"))
  
  # make wide by disease
  out_wide <- out_long %>%
    pivot_wider(names_from = disease, 
                values_from = -c(iteration, .draw, generation, disease))
  
  # output both formats
  return(list(long = out_long, wide = out_wide))

}

# calculate responses
resp_fun <- function(iter, sim_outcome, parameters, gens = 2,
                     L_add, CA_add, CP_add){
  
  # initial conditions (annual seeds, litter, perennial seeds, perennial adults)
  inits <- c(sim_outcome$annual_seeds,
             sim_outcome$litter_star,
             sim_outcome$perennial_seeds,
             sim_outcome$perennial_adults)
  
  # initial competition values (experienced by annual, perennial)
  inits_comp <- c(sim_outcome$C_A_star,
                  sim_outcome$C_P_star)
  
  # simulate one year under baseline conditions
  sim_base <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
                          parameters = parameters, init_C = inits_comp,
                          return_M = T)
  
  # re-save competition factors
  inits_L2 <- inits
  inits_compA2 <- inits_comp
  inits_compP2 <- inits_comp
  
  # inflate competition factors
  inits_L2[2] <- inits[2] + L_add
  inits_compA2[1] <- inits_comp[1] + CA_add
  inits_compP2[2] <- inits_comp[2] + CP_add
  
  # simulate one year with inflation for each species
  sim_L <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_L2,
                       parameters = parameters, init_C = inits_comp,
                       return_M = T)
  sim_compA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
                           parameters = parameters, init_C = inits_compA2)
  sim_compP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
                           parameters = parameters, init_C = inits_compP2,
                           return_M = T)
  
  # output dataframe
  out_base <- sim_base[[1]]
  out_L <- sim_L[[1]]
  
  # competition factors
  L_base <- inits[2]
  L_adj <- inits_L2[2]
  CA_base <- inits_comp[1]
  CA_adj <- inits_compA2[1]
  CP_base <- inits_comp[2]
  CP_adj <- inits_compP2[2]
  
  # check for presence of annual
  if(inits[1] >= pres_thresh){
    
    # calculate growth rates r = (ln[N(t+1)] -  ln[N(t)])/delta_t
    gr_baseA <- log(out_base$annual_seeds[2] / out_base$annual_seeds[1])
    gr_LA <- log(out_L$annual_seeds[2] / out_L$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    # calculate responses
    resp_LA <- -1 * (gr_LA - gr_baseA) / (L_adj - L_base)
    resp_CA <- -1 * (gr_CA - gr_baseA) / (CA_adj - CA_base)
    
  } else {
    
    gr_baseA <- gr_LA <- gr_CA <- resp_LA <- resp_CA <- NA
    
  }
  
  # check for presence of perennial
  if(inits[3] >= pres_thresh & inits[4] >= pres_thresh){
    
    # calculate growth rates dominant eigenvalue
    gr_baseP <- log(eigen(sim_base[[2]], only.values = T)$values[1])
    gr_LP <- log(eigen(sim_L[[2]], only.values = T)$values[1])
    gr_CP <- log(eigen(sim_compP[[2]], only.values = T)$values[1])
    
    # calculate responses
    resp_LP <- -1 * (gr_LP - gr_baseP) / (L_adj - L_base)
    resp_CP <- -1 * (gr_CP - gr_baseP) / (CP_adj - CP_base)
    
  } else {
    
    gr_baseP <- gr_LP <- gr_CP <- resp_LP <- resp_CP <- NA
    
  }
  
  # combine values
  out <- tibble(.draw = parameters[["draws"]][iter, ]$.draw,
                gr_baseA = gr_baseA, # check that these are zero
                gr_baseP = gr_baseP,
                gr_LA = gr_LA,
                gr_LP = gr_LP,
                gr_CA = gr_CA,
                gr_CP = gr_CP,
                L_base = L_base,
                L_adj = L_adj,
                CA_base = CA_base,
                CA_adj = CA_adj,
                CP_base = CP_base,
                CP_adj = CP_adj,
                resp_LA = resp_LA,
                resp_LP = resp_LP,
                resp_CA = resp_CA,
                resp_CP = resp_CP)
  
  return(out)
}


#### single species simulations ####

# initial conditions (annual seeds, litter, perennial seeds, perennial adults)
initsA <- c(1, 0, 0, 0)
initsP <- c(0, 0, 0, 1)

# output lits
sims_gfung_A_h <- list()
sims_gfung_A_d <- list()

sims_gfung_P_h <- list()
sims_gfung_P_d <- list()

sims_ginf_A_h <- list()
sims_ginf_A_d <- list()

sims_infP_A_h <- list()
sims_infP_A_d <- list()

# cycle through parameters
for(i in 1:params_iters){
  
  # fungicide effect, healthy
  sims_gfung_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["healthy"]])
  sims_gfung_P_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["healthy"]])

  # fungicide effect, disease
  sims_gfung_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["disease"]])
  sims_gfung_P_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["disease"]])

  # infection effect, healthy
  sims_ginf_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_ginf[["healthy"]])

  # infection effect, disease
  sims_ginf_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                    parameters = params_ginf[["disease"]])
  
  # infection effect and priors, healthy
  sims_infP_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                    parameters = params_infP[["healthy"]])
  
  # infection effect and priors, disease
  sims_infP_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                    parameters = params_infP[["disease"]])
  
}

# combine healthy and disease
sims_gfung_A <- sims_comb_fun(sims_gfung_A_h, sims_gfung_A_d)
sims_ginf_A <- sims_comb_fun(sims_ginf_A_h, sims_ginf_A_d)
sims_infP_A <- sims_comb_fun(sims_infP_A_h, sims_infP_A_d)
sims_gfung_P <- sims_comb_fun(sims_gfung_P_h, sims_gfung_P_d)

# save
save(sims_gfung_A, file = "output/time_series_gfung_parms_annual_only_20250719.rda")
save(sims_ginf_A, file = "output/time_series_ginf_parms_annual_only_20250719.rda")
save(sims_infP_A, file = "output/time_series_infP_parms_annual_only_20250719.rda")
save(sims_gfung_P, file = "output/time_series_gfung_parms_perennial_only_20250719.rda")

# reload if needed
load("output/time_series_gfung_parms_annual_only_20250719.rda")
load("output/time_series_ginf_parms_annual_only_20250719.rda")
load("output/time_series_infP_parms_annual_only_20250719.rda")
load("output/time_series_gfung_parms_perennial_only_20250719.rda")


#### single species time series figures ####

# select random rows for visualization
viz_draw_sample <- sample(1:15000, 1000)

# combine datasets
sims_gfung_single <- sims_gfung_A[["long"]] %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds, litter) %>%
  rename(annual_litter = litter,
         annual = annual_seeds) %>%
  full_join(sims_gfung_P[["long"]] %>% 
              filter(.draw %in% viz_draw_sample) %>%
              mutate(perennial = perennial_seeds + perennial_adults) %>%
              select(disease, generation, .draw, perennial, litter) %>%
              rename(perennial_litter = litter)) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, annual_litter, perennial, perennial_litter)) %>%
  mutate(species = str_replace(name, "_", " ") %>%
           str_replace("annual", "*M. vimineum*") %>%
           str_replace("perennial", "*E. virginicus*"),
         species = if_else(str_detect(species, "litter"), 
                           paste(species, "(g)"), 
                           paste(species, "density")),
         species = fct_reorder(species, as.numeric(as.factor(name))))

sims_ginf_single <- sims_ginf_A[["long"]] %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds, litter) %>%
  rename(annual_litter = litter,
         annual = annual_seeds) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, annual_litter)) %>%
  mutate(species = str_replace(name, "_", " ") %>%
           str_replace("annual", "*M. vimineum*"),
         species = if_else(str_detect(species, "litter"), 
                           paste(species, "(g)"), 
                           paste(species, "density")),
         species = fct_reorder(species, as.numeric(as.factor(name))))

sims_infP_single <- sims_infP_A[["long"]] %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds, litter) %>%
  rename(annual_litter = litter,
         annual = annual_seeds) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, annual_litter)) %>%
  mutate(species = str_replace(name, "_", " ") %>%
           str_replace("annual", "*M. vimineum*"),
         species = if_else(str_detect(species, "litter"), 
                           paste(species, "(g)"), 
                           paste(species, "density")),
         species = fct_reorder(species, as.numeric(as.factor(name))))

# figure
sims_gfung_single_fig <- ggplot(sims_gfung_single, 
       aes(x = generation, y = value, color = .draw, group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
    theme(axis.title.y = element_blank(),
          strip.placement.y = "outside",
          strip.text = element_markdown())

sims_ginf_single_fig <- ggplot(sims_ginf_single, 
       aes(x = generation, y = value, color = .draw, group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_markdown())

sims_infP_single_fig <- ggplot(sims_infP_single, 
                               aes(x = generation, y = value, color = .draw, group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_markdown())

# save
ggsave("output/time_series_gfung_parms_single_species.png", 
       sims_gfung_single_fig, width = 6, height = 8.5)
ggsave("output/time_series_ginf_parms_single_species.png", 
       sims_ginf_single_fig, width = 6, height = 4.5)
ggsave("output/time_series_infP_parms_single_species.png", 
       sims_infP_single_fig, width = 6, height = 4.5)


#### single species equilibria ####

# save last time point
# make sure order matches parameters
eq_gfung_A <- sims_gfung_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(C_AA = annual_competition,
         C_PA = perennial_competition,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(C_AA = annual_competition,
         C_PA = perennial_competition,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_infP_A <- sims_infP_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(C_AA = annual_competition,
         C_PA = perennial_competition,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_gfung_P <- sims_gfung_P[["long"]] %>%
  filter(generation == gens) %>%
  mutate(C_AP = annual_competition,
         C_PP = perennial_competition,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

# remove raw simulations
rm(sims_gfung_A, sims_gfung_P, sims_ginf_A, sims_infP_A)
gc()

# identify cases where species couldn't establish
non_gfung_A <- eq_gfung_A %>% 
  filter(annual_seeds < pres_thresh) %>% 
  select(.draw, disease)

non_ginf_A <- eq_ginf_A %>% 
  filter(annual_seeds < pres_thresh) %>% 
  select(.draw, disease)

non_infP_A <- eq_infP_A %>% 
  filter(annual_seeds < pres_thresh) %>% 
  select(.draw, disease)

non_gfung_P <- eq_gfung_P %>% 
  filter(perennial_seeds < pres_thresh | perennial_adults < pres_thresh) %>% 
  select(.draw, disease)

# collapse parameters into tables
params_gfung2 <- params_gfung %>% 
  imap_dfr(~ bind_cols(.x) %>%
             mutate(disease = str_sub(.y, 1, 1)))

params_ginf2 <- params_ginf %>% 
  imap_dfr(~ bind_cols(.x) %>%
             mutate(disease = str_sub(.y, 1, 1)))

params_infP2 <- params_infP %>% 
  imap_dfr(~ bind_cols(.x) %>%
             mutate(disease = str_sub(.y, 1, 1)))

# reference equilibria
eq_gfung <- eq_gfung_A %>% 
  select(.draw, disease, disease_name, litter_A = litter, C_AA, C_PA,
         annual_seeds) %>% 
  full_join(eq_gfung_P %>% 
              select(.draw, disease, disease_name, litter_P = litter,
                     C_AP, C_PP, perennial_seeds, perennial_adults)) %>% 
  mutate(litter_star = (litter_A + litter_P) / 2) %>% 
  left_join(params_gfung2) %>% 
  mutate(C_A_star = ((gA * eA * yA) / (1 - sA * (1 - gA))) * (1 / (1 + betaA * litter_star)),
         C_P_star = ((gP * eP * yP) / (1 - sP * (1 - gP))) * (1 / (1 + betaP * litter_star)) * (f + (pS / (1 - pP))),
         litter_A_effect = litter_A - litter_star,
         litter_P_effect = litter_P - litter_star,
         C_AA_effect = C_AA - C_A_star,
         C_AP_effect = C_AP - C_A_star,
         C_PA_effect = C_PA - C_P_star,
         C_PP_effect = C_PP - C_P_star)

eq_ginf <- eq_ginf_A %>% 
  select(.draw, disease, disease_name, litter_A = litter, C_AA, C_PA,
         annual_seeds) %>% 
  full_join(eq_gfung_P %>% 
              select(.draw, disease, disease_name, litter_P = litter,
                     C_AP, C_PP, perennial_seeds, perennial_adults)) %>% 
  mutate(litter_star = (litter_A + litter_P) / 2) %>% 
  left_join(params_ginf2) %>% 
  mutate(C_A_star = ((gA * eA * yA) / (1 - sA * (1 - gA))) * (1 / (1 + betaA * litter_star)),
         C_P_star = ((gP * eP * yP) / (1 - sP * (1 - gP))) * (1 / (1 + betaP * litter_star)) * (f + (pS / (1 - pP))),
         litter_A_effect = litter_A - litter_star,
         litter_P_effect = litter_P - litter_star,
         C_AA_effect = C_AA - C_A_star,
         C_AP_effect = C_AP - C_A_star,
         C_PA_effect = C_PA - C_P_star,
         C_PP_effect = C_PP - C_P_star)

eq_infP <- eq_infP_A %>% 
  select(.draw, disease, disease_name, litter_A = litter, C_AA, C_PA,
         annual_seeds) %>% 
  full_join(eq_gfung_P %>% 
              select(.draw, disease, disease_name, litter_P = litter,
                     C_AP, C_PP, perennial_seeds, perennial_adults)) %>% 
  mutate(litter_star = (litter_A + litter_P) / 2) %>% 
  left_join(params_infP2) %>% 
  mutate(C_A_star = ((gA * eA * yA) / (1 - sA * (1 - gA))) * (1 / (1 + betaA * litter_star)),
         C_P_star = ((gP * eP * yP) / (1 - sP * (1 - gP))) * (1 / (1 + betaP * litter_star)) * (f + (pS / (1 - pP))),
         litter_A_effect = litter_A - litter_star,
         litter_P_effect = litter_P - litter_star,
         C_AA_effect = C_AA - C_A_star,
         C_AP_effect = C_AP - C_A_star,
         C_PA_effect = C_PA - C_P_star,
         C_PP_effect = C_PP - C_P_star)


#### invasion simulations ####

# output lits (first letter is invader, second is resident)
sims_gfung_AP_h <- list()
sims_gfung_PA_h <- list()

sims_gfung_AP_d <- list()
sims_gfung_PA_d <- list()

sims_ginf_AP_h <- list()
sims_ginf_PA_h <- list()

sims_ginf_AP_d <- list()
sims_ginf_PA_d <- list()

sims_infP_AP_h <- list()
sims_infP_PA_h <- list()

sims_infP_AP_d <- list()
sims_infP_PA_d <- list()

grwr_gfung_PA_h <- list()
grwr_gfung_PA_d <- list()

grwr_ginf_PA_h <- list()
grwr_ginf_PA_d <- list()

grwr_infP_PA_h <- list()
grwr_infP_PA_d <- list()

# cycle through parameters
for(i in 1:params_iters){
  
  # get draw from parameters
  draw_gfung_h <- params_gfung[["healthy"]][["draws"]]$.draw[i]
  draw_gfung_d <- params_gfung[["disease"]][["draws"]]$.draw[i]
  draw_ginf_h <- params_ginf[["healthy"]][["draws"]]$.draw[i]
  draw_ginf_d <- params_ginf[["disease"]][["draws"]]$.draw[i]
  draw_infP_h <- params_infP[["healthy"]][["draws"]]$.draw[i] 
  draw_infP_d <- params_infP[["disease"]][["draws"]]$.draw[i] 
  
  # select equilibrium values
  inits_gfung_P_h <- eq_gfung_P %>%
    filter(disease == "h" & .draw == draw_gfung_h)
  inits_gfung_P_d <- eq_gfung_P %>%
    filter(disease == "d" & .draw == draw_gfung_h)

  inits_ginf_P_h <- eq_gfung_P %>%
    filter(disease == "h" & .draw == draw_ginf_h)
  inits_ginf_P_d <- eq_gfung_P %>%
    filter(disease == "d" & .draw == draw_ginf_h)
  
  inits_infP_P_h <- eq_gfung_P %>%
    filter(disease == "h" & .draw == draw_infP_h)
  inits_infP_P_d <- eq_gfung_P %>%
    filter(disease == "d" & .draw == draw_infP_h)
  
  inits_gfung_A_h <- eq_gfung_A %>%
    filter(disease == "h" & .draw == draw_gfung_h)
  inits_gfung_A_d <- eq_gfung_A %>%
    filter(disease == "d" & .draw == draw_gfung_h)

  inits_ginf_A_h <- eq_ginf_A %>%
    filter(disease == "h" & .draw == draw_ginf_h)
  inits_ginf_A_d <- eq_ginf_A %>%
    filter(disease == "d" & .draw == draw_ginf_h)
  
  inits_infP_A_h <- eq_infP_A %>% 
    filter(disease == "h" & .draw == draw_infP_h)
  inits_infP_A_d <- eq_infP_A %>% 
    filter(disease == "d" & .draw == draw_infP_h)
  
  # initial conditions based on resident's simulation 
  # (annual seeds, litter, perennial seeds, perennial adults)
  inits_gfung_resP_h <- c(1,
                          inits_gfung_P_h$litter,
                          inits_gfung_P_h$perennial_seeds,
                          inits_gfung_P_h$perennial_adults)
  inits_gfung_resP_d <- c(1,
                          inits_gfung_P_d$litter,
                          inits_gfung_P_d$perennial_seeds,
                          inits_gfung_P_d$perennial_adults)

  inits_ginf_resP_h <- c(1,
                          inits_ginf_P_h$litter,
                          inits_ginf_P_h$perennial_seeds,
                          inits_ginf_P_h$perennial_adults)
  inits_ginf_resP_d <- c(1,
                          inits_ginf_P_d$litter,
                          inits_ginf_P_d$perennial_seeds,
                          inits_ginf_P_d$perennial_adults)
  
  inits_infP_resP_h <- c(1,
                         inits_infP_P_h$litter,
                         inits_infP_P_h$perennial_seeds,
                         inits_infP_P_h$perennial_adults)
  inits_infP_resP_d <- c(1,
                         inits_infP_P_d$litter,
                         inits_infP_P_d$perennial_seeds,
                         inits_infP_P_d$perennial_adults)
  
  inits_gfung_resA_h <- c(inits_gfung_A_h$annual_seeds,
                          inits_gfung_A_h$litter,
                          0,
                          1)
  inits_gfung_resA_d <- c(inits_gfung_A_d$annual_seeds,
                          inits_gfung_A_d$litter,
                          0,
                          1)

  inits_ginf_resA_h <- c(inits_ginf_A_h$annual_seeds,
                          inits_ginf_A_h$litter,
                          0,
                          1)
  inits_ginf_resA_d <- c(inits_ginf_A_d$annual_seeds,
                          inits_ginf_A_d$litter,
                          0,
                          1)
  
  inits_infP_resA_h <- c(inits_infP_A_h$annual_seeds,
                         inits_infP_A_h$litter,
                         0,
                         1)
  inits_infP_resA_d <- c(inits_infP_A_d$annual_seeds,
                         inits_infP_A_d$litter,
                         0,
                         1)
  
  # fungicide effect, healthy
  sims_gfung_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resP_h,
                                      parameters = params_gfung[["healthy"]])
  sims_gfung_PA_h_out <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_gfung_resA_h,
                                     parameters = params_gfung[["healthy"]],
                                     return_M = T)
  sims_gfung_PA_h[[i]] <- sims_gfung_PA_h_out[[1]]
  grwr_gfung_PA_h[[i]] <- tibble(grwr = eigen(sims_gfung_PA_h_out[[2]],
                                              only.values = T)$values[1],
                                 .draw = draw_gfung_h)

  # fungicide effect, disease
  sims_gfung_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resP_d,
                                      parameters = params_gfung[["disease"]])
  sims_gfung_PA_d_out <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resA_d,
                                      parameters = params_gfung[["disease"]],
                                      return_M = T)
  sims_gfung_PA_d[[i]] <- sims_gfung_PA_d_out[[1]]
  grwr_gfung_PA_d[[i]] <- tibble(grwr = eigen(sims_gfung_PA_d_out[[2]],
                                              only.values = T)$values[1],
                                 .draw = draw_gfung_d)

  # infection effect, healthy
  sims_ginf_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_ginf_resP_h,
                                     parameters = params_ginf[["healthy"]])
  sims_ginf_PA_h_out <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_ginf_resA_h,
                                     parameters = params_ginf[["healthy"]],
                                     return_M = T)
  sims_ginf_PA_h[[i]] <- sims_ginf_PA_h_out[[1]]
  grwr_ginf_PA_h[[i]] <- tibble(grwr = eigen(sims_ginf_PA_h_out[[2]],
                                             only.values = T)$values[1],
                                .draw = draw_ginf_h)

  # infection effect, disease
  sims_ginf_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_ginf_resP_d,
                                     parameters = params_ginf[["disease"]])
  sims_ginf_PA_d_out <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_ginf_resA_d,
                                     parameters = params_ginf[["disease"]],
                                     return_M = T)
  sims_ginf_PA_d[[i]] <- sims_ginf_PA_d_out[[1]]
  grwr_ginf_PA_d[[i]] <- tibble(grwr = eigen(sims_ginf_PA_d_out[[2]],
                                             only.values = T)$values[1],
                                .draw = draw_ginf_h)
  
  # infection effect + priors, healthy
  sims_infP_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_infP_resP_h,
                                     parameters = params_infP[["healthy"]])
  sims_infP_PA_h_out <- dyn_mod_fun(iter = i, gen = gens, 
                                    init_cond = inits_infP_resA_h,
                                    parameters = params_infP[["healthy"]],
                                    return_M = T)
  sims_infP_PA_h[[i]] <- sims_infP_PA_h_out[[1]]
  grwr_infP_PA_h[[i]] <- tibble(grwr = eigen(sims_infP_PA_h_out[[2]], 
                                             only.values = T)$values[1],
                                .draw = draw_infP_h)
  
  # infection effect + priors, disease
  sims_infP_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_infP_resP_d,
                                     parameters = params_infP[["disease"]])
  sims_infP_PA_d_out <- dyn_mod_fun(iter = i, gen = gens, 
                                    init_cond = inits_infP_resA_d,
                                    parameters = params_infP[["disease"]],
                                    return_M = T)
  sims_infP_PA_d[[i]] <- sims_infP_PA_d_out[[1]]
  grwr_infP_PA_d[[i]] <- tibble(grwr = eigen(sims_infP_PA_d_out[[2]], 
                                             only.values = T)$values[1],
                                .draw = draw_infP_h)
  
}

# combine healthy and disease
sims_gfung_AP <- sims_comb_fun(sims_gfung_AP_h, sims_gfung_AP_d)
sims_gfung_PA <- sims_comb_fun(sims_gfung_PA_h, sims_gfung_PA_d)
sims_ginf_AP <- sims_comb_fun(sims_ginf_AP_h, sims_ginf_AP_d)
sims_ginf_PA <- sims_comb_fun(sims_ginf_PA_h, sims_ginf_PA_d)
sims_infP_AP <- sims_comb_fun(sims_infP_AP_h, sims_infP_AP_d)
sims_infP_PA <- sims_comb_fun(sims_infP_PA_h, sims_infP_PA_d)

grwr_gfung_PA <- bind_rows(grwr_gfung_PA_h, .id = "iteration") %>%
  mutate(disease = "h") %>%
  full_join(bind_rows(grwr_gfung_PA_d, .id = "iteration") %>%
              mutate(disease = "d")) %>%
  mutate(iteration = as.numeric(iteration),
         disease = fct_relevel(disease, "h"))
grwr_ginf_PA <- bind_rows(grwr_ginf_PA_h, .id = "iteration") %>%
  mutate(disease = "h") %>%
  full_join(bind_rows(grwr_ginf_PA_d, .id = "iteration") %>%
              mutate(disease = "d")) %>%
  mutate(iteration = as.numeric(iteration),
         disease = fct_relevel(disease, "h"))
grwr_infP_PA <- bind_rows(grwr_infP_PA_h, .id = "iteration") %>%
  mutate(disease = "h") %>%
  full_join(bind_rows(grwr_infP_PA_d, .id = "iteration") %>%
              mutate(disease = "d")) %>%
  mutate(iteration = as.numeric(iteration),
         disease = fct_relevel(disease, "h"))

# divide for smaller files
sims_gfung_AP_long <- sims_gfung_AP[["long"]]
sims_gfung_AP_wide <- sims_gfung_AP[["wide"]]
sims_gfung_PA_long <- sims_gfung_PA[["long"]]
sims_gfung_PA_wide <- sims_gfung_PA[["wide"]]
sims_ginf_AP_long <- sims_ginf_AP[["long"]]
sims_ginf_AP_wide <- sims_ginf_AP[["wide"]]
sims_ginf_PA_long <- sims_ginf_PA[["long"]]
sims_ginf_PA_wide <- sims_ginf_PA[["wide"]]
sims_infP_AP_long <- sims_infP_AP[["long"]]
sims_infP_AP_wide <- sims_infP_AP[["wide"]]
sims_infP_PA_long <- sims_infP_PA[["long"]]
sims_infP_PA_wide <- sims_infP_PA[["wide"]]

# save
save(sims_gfung_AP_long,
     file = "output/time_series_gfung_parms_annual_invasion_20250719.rda")
save(sims_gfung_AP_wide,
     file = "output/time_series_wide_gfung_parms_annual_invasion_20250719.rda")
save(sims_gfung_PA_long,
     file = "output/time_series_gfung_parms_perennial_invasion_20250719.rda")
save(sims_gfung_PA_wide,
     file = "output/time_series_wide_gfung_parms_perennial_invasion_20250719.rda")
save(sims_ginf_AP_long,
     file = "output/time_series_ginf_parms_annual_invasion_20250719.rda")
save(sims_ginf_AP_wide,
     file = "output/time_series_wide_ginf_parms_annual_invasion_20250719.rda")
save(sims_ginf_PA_long,
     file = "output/time_series_ginf_parms_perennial_invasion_20250719.rda")
save(sims_ginf_PA_wide,
     file = "output/time_series_wide_ginf_parms_perennial_invasion_20250719.rda")
save(sims_infP_AP_long, 
     file = "output/time_series_infP_parms_annual_invasion_20250719.rda")
save(sims_infP_AP_wide, 
     file = "output/time_series_wide_infP_parms_annual_invasion_20250719.rda")
save(sims_infP_PA_long, 
     file = "output/time_series_infP_parms_perennial_invasion_20250719.rda")
save(sims_infP_PA_wide, 
     file = "output/time_series_wide_infP_parms_perennial_invasion_20250719.rda")
save(grwr_gfung_PA,
     file = "output/grwr_gfung_parms_perennial_invasion_20250719.rda")
save(grwr_ginf_PA,
     file = "output/grwr_ginf_parms_perennial_invasion_20250719.rda")
save(grwr_infP_PA, 
     file = "output/grwr_infP_parms_perennial_invasion_20250719.rda")

# reload if needed
load("output/time_series_gfung_parms_annual_invasion_20250719.rda")
load("output/time_series_gfung_parms_perennial_invasion_20250719.rda")
load("output/time_series_ginf_parms_annual_invasion_20250719.rda")
load("output/time_series_ginf_parms_perennial_invasion_20250719.rda")
load("output/time_series_infP_parms_annual_invasion_20250719.rda")
load("output/time_series_infP_parms_perennial_invasion_20250719.rda")
load("output/grwr_gfung_parms_perennial_invasion_20250719.rda")
load("output/grwr_ginf_parms_perennial_invasion_20250719.rda")
load("output/grwr_infP_parms_perennial_invasion_20250719.rda")


#### invasion time series figures ####

# combine datasets
sims_gfung_inv <- sims_gfung_AP_long %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  rename(annual = annual_seeds) %>%
  full_join(sims_gfung_PA_long %>% 
              filter(.draw %in% viz_draw_sample) %>%
              mutate(perennial = perennial_seeds + perennial_adults) %>%
              select(disease, generation, .draw, perennial)) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, perennial)) %>%
  mutate(species = fct_recode(name, 
                              "*M. vimineum* density" = "annual",
                              "*E. virginicus* density" = "perennial"),
         species = fct_reorder(species, as.numeric(as.factor(name))))

sims_ginf_inv <- sims_ginf_AP_long %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  rename(annual = annual_seeds) %>%
  full_join(sims_ginf_PA_long %>% 
              filter(.draw %in% viz_draw_sample) %>%
              mutate(perennial = perennial_seeds + perennial_adults) %>%
              select(disease, generation, .draw, perennial)) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, perennial)) %>%
  mutate(species = fct_recode(name, 
                              "*M. vimineum* density" = "annual",
                              "*E. virginicus* density" = "perennial"),
         species = fct_reorder(species, as.numeric(as.factor(name))))

sims_infP_inv <- sims_infP_AP_long %>% 
  filter(.draw %in% viz_draw_sample) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  rename(annual = annual_seeds) %>%
  full_join(sims_infP_PA_long %>% 
              filter(.draw %in% viz_draw_sample) %>%
              mutate(perennial = perennial_seeds + perennial_adults) %>%
              select(disease, generation, .draw, perennial)) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease")) %>%
  pivot_longer(cols = c(annual, perennial)) %>%
  mutate(species = fct_recode(name, 
                              "*M. vimineum* density" = "annual",
                              "*E. virginicus* density" = "perennial"),
         species = fct_reorder(species, as.numeric(as.factor(name))))

# figure
sims_gfung_inv_fig <- ggplot(sims_gfung_inv, 
                                aes(x = generation, y = value, color = .draw, 
                                    group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_markdown())

sims_ginf_inv_fig <- ggplot(sims_ginf_inv, 
                               aes(x = generation, y = value, color = .draw, 
                                   group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_markdown())

sims_infP_inv_fig <- ggplot(sims_infP_inv, 
                            aes(x = generation, y = value, color = .draw, 
                                group = .draw)) +
  geom_line(linewidth = 0.5) +
  facet_grid(species ~ disease_name, scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", color = "Parameter draw") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.placement.y = "outside",
        strip.text = element_markdown())

# save
ggsave("output/time_series_gfung_parms_invasions.png", 
       sims_gfung_inv_fig, width = 6, height = 4.5)
ggsave("output/time_series_ginf_parms_invasions.png", 
       sims_ginf_inv_fig, width = 6, height = 4.5)
ggsave("output/time_series_infP_parms_invasions.png", 
       sims_infP_inv_fig, width = 6, height = 4.5)


#### invasion outcomes ####

# calculate growth rates r = (ln[N(t+1)] -  ln[N(t)])/delta_t (or log(lambda))
inv_gfung_AP <- sims_gfung_AP_long %>%
  filter(generation == 2) %>%
  mutate(annual_gr = log(annual_seeds),
         annual_invasion = if_else(annual_gr > 0, "yes", "no"))

inv_gfung_PA <- grwr_gfung_PA %>%
  mutate(perennial_gr = log(grwr),
         perennial_invasion = if_else(perennial_gr > 0, 
                                      "yes", "no"))

inv_ginf_AP <- sims_ginf_AP_long %>%
  filter(generation == 2) %>%
  mutate(annual_gr = log(annual_seeds),
         annual_invasion = if_else(annual_gr > 0, "yes", "no")) 

inv_ginf_PA <- grwr_ginf_PA %>%
  mutate(perennial_gr = log(grwr),
         perennial_invasion = if_else(perennial_gr > 0, 
                                      "yes", "no"))

inv_infP_AP <- sims_infP_AP_long %>%
  filter(generation == 2) %>%
  mutate(annual_gr = log(annual_seeds),
         annual_invasion = if_else(annual_gr > 0, "yes", "no"))

inv_infP_PA <- grwr_infP_PA %>%
  mutate(perennial_gr = log(grwr),
         perennial_invasion = if_else(perennial_gr > 0, 
                                      "yes", "no")) 

# combine invasion outcomes
# add final values from annual invasion
# add final values from each species alone
inv_gfung <- inv_gfung_AP %>%
  select(.draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_gfung_PA %>%
              select(.draw, disease, perennial_invasion, perennial_gr)) %>%
  left_join(non_gfung_A %>% mutate(annual_est = "no")) %>% 
  left_join(non_gfung_P %>% mutate(perennial_est = "no")) %>% 
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no"),
         annual_est = replace_na(annual_est, "yes"),
         perennial_est = replace_na(perennial_est, "yes"),
         outcome = case_when(annual_est == "no" | perennial_est == "no" ~
                               "non-establishment",
                             coexist == "yes" ~ "coexist",
                             annual_invasion == "yes" & 
                               perennial_invasion == "no" ~ "*M. vimineum* only",
                             perennial_invasion == "yes" &
                               annual_invasion == "no" ~ "*E. virginicus* only",
                             TRUE ~ "priority effect") %>%
           fct_relevel("coexist", "*M. vimineum* only")) %>%
  left_join(sims_gfung_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_gfung %>%
              select(.draw, disease, disease_name, 
                     annual_seeds_A = annual_seeds, 
                     perennial_seeds_P = perennial_seeds, 
                     perennial_adults_P = perennial_adults))

inv_ginf <- inv_ginf_AP %>%
  select(.draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_ginf_PA %>%
              select(.draw, disease, perennial_invasion, perennial_gr)) %>%
  left_join(non_ginf_A %>% mutate(annual_est = "no")) %>% 
  left_join(non_gfung_P %>% mutate(perennial_est = "no")) %>% 
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no"),
         annual_est = replace_na(annual_est, "yes"),
         perennial_est = replace_na(perennial_est, "yes"),
         outcome = case_when(annual_est == "no" | perennial_est == "no" ~
                               "non-establishment",
                             coexist == "yes" ~ "coexist",
                             annual_invasion == "yes" & 
                               perennial_invasion == "no" ~ "*M. vimineum* only",
                             perennial_invasion == "yes" &
                               annual_invasion == "no" ~ "*E. virginicus* only",
                             TRUE ~ "priority effect") %>%
           fct_relevel("coexist", "*M. vimineum* only")) %>%
  left_join(sims_ginf_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_ginf %>%
              select(.draw, disease, disease_name, 
                     annual_seeds_A = annual_seeds, 
                     perennial_seeds_P = perennial_seeds, 
                     perennial_adults_P = perennial_adults))

inv_infP <- inv_infP_AP %>%
  select(.draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_infP_PA %>%
              select(.draw, disease, perennial_invasion, perennial_gr)) %>%
  left_join(non_infP_A %>% mutate(annual_est = "no")) %>% 
  left_join(non_gfung_P %>% mutate(perennial_est = "no")) %>% 
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no"),
         annual_est = replace_na(annual_est, "yes"),
         perennial_est = replace_na(perennial_est, "yes"),
         outcome = case_when(annual_est == "no" | perennial_est == "no" ~
                               "non-establishment",
                             coexist == "yes" ~ "coexist",
                             annual_invasion == "yes" & 
                               perennial_invasion == "no" ~ "*M. vimineum* only",
                             perennial_invasion == "yes" &
                               annual_invasion == "no" ~ "*E. virginicus* only",
                             TRUE ~ "priority effect") %>%
           fct_relevel("coexist", "*M. vimineum* only")) %>%
  left_join(sims_infP_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_infP %>%
              select(.draw, disease, disease_name, 
                     annual_seeds_A = annual_seeds, 
                     perennial_seeds_P = perennial_seeds, 
                     perennial_adults_P = perennial_adults))

# remove sims
rm(sims_gfung_AP_long, sims_gfung_PA_long, sims_ginf_AP_long, 
   sims_ginf_PA_long, sims_infP_AP_long, sims_infP_PA_long)
gc()

# priority effect simulations
inv_gfung %>%
  filter(outcome == "priority effect") %>%
  select(.draw, disease, annual_seeds_A, perennial_seeds_P,
         perennial_adults_P) %>%
  mutate(perennial_P = perennial_seeds_P + perennial_adults_P) %>%
  filter(annual_seeds_A < 1 | perennial_P < 1)
# both can establish

inv_ginf %>%
  filter(outcome == "priority effect") %>%
  select(.draw, disease, annual_seeds_A, perennial_seeds_P, 
         perennial_adults_P) %>%
  mutate(perennial_P = perennial_seeds_P + perennial_adults_P) %>%
  filter(annual_seeds_A < 1 | perennial_P < 1)
# both can establish

inv_infP %>%
  filter(outcome == "priority effect") %>%
  select(.draw, disease, annual_seeds_A, perennial_seeds_P, 
         perennial_adults_P) %>%
  mutate(perennial_P = perennial_seeds_P + perennial_adults_P) %>%
  filter(annual_seeds_A < 1 | perennial_P < 1)
# both can establish

# summarize for figure
inv_gfung_sum <- inv_gfung %>%
  count(disease_name, outcome) %>%
  group_by(disease_name) %>%
  mutate(prop = n / (sum(n))) %>%
  ungroup() %>% 
  mutate(outcome = str_replace(outcome, " only", "<br>only") %>% 
           str_replace(" effect", "<br>effect") %>% 
           fct_relevel("coexist", "*M. vimineum*<br>only"))

inv_ginf_sum <- inv_ginf %>%
  count(disease_name, outcome) %>%
  group_by(disease_name) %>%
  mutate(prop = n / (sum(n))) %>%
  ungroup() %>% 
  mutate(outcome = str_replace(outcome, " only", "<br>only") %>% 
           str_replace(" effect", "<br>effect") %>% 
           fct_relevel("coexist", "*M. vimineum*<br>only"))

inv_infP_sum <- inv_infP %>%
  count(disease_name, outcome) %>%
  group_by(disease_name) %>%
  mutate(prop = n / (sum(n))) %>%
  ungroup() %>% 
  mutate(outcome = str_replace(outcome, " only", "<br>only") %>% 
           str_replace(" effect", "<br>effect") %>% 
           fct_relevel("coexist", "*M. vimineum*<br>only"))


#### invasion outcome figures ####

# summary figure
coex_gfung_sum_fig <- inv_gfung_sum %>% 
  filter(outcome != "non-establishment") %>% 
  ggplot(aes(x = outcome, y = prop, 
             fill = disease_name, color = disease_name)) +
  geom_col(position = "dodge", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = dis_pal) +
  scale_color_manual(values = dis_pal) +
  labs(y = "Set 1 draws") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 8),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.direction = "vertical",
        legend.position.inside = c(0.76, 0.75))

coex_ginf_sum_fig <- coex_gfung_sum_fig +
  filter(inv_ginf_sum, outcome != "non-establishment") +
  labs(y = "Set 2 draws") +
  theme(legend.position = "none")

coex_infP_sum_fig <- coex_gfung_sum_fig +
  filter(inv_infP_sum, outcome != "non-establishment") +
  labs(y = "Set 3 draws") +
  theme(legend.position = "none")

# growth rate figures
coex_infP_grd_fig <- inv_infP %>%
  filter(disease == "d" & outcome != "non-establishment") %>%
  ggplot(aes(y = annual_gr, x = perennial_gr, color = outcome, 
             shape = outcome)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 0.5) +
  scale_color_manual(values = col_pal4) +
  scale_shape_manual(values = shape_pal4) +
  labs(x = "*E. virginicus* GRWR", y = "*M. vimineum* GRWR",
       title = "Parameter Set 3: Ambient Disease") +
  fig_theme +
  theme(axis.title.x = element_markdown(size = 9, color = "black"),
        axis.title.y = element_markdown(size = 9, color = "black"),
        plot.title = element_text(size = 9, face = "plain"),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "inside",
        legend.position.inside = c(0.73, 0.85)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# combine
comb_sum_fig <- coex_gfung_sum_fig + coex_ginf_sum_fig + coex_infP_sum_fig +
  plot_layout(ncol = 1, axes = "collect_x") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

coex_fig <- coex_infP_grd_fig +
  theme(plot.tag = element_text(size = 9, face = "bold")) + 
  comb_sum_fig +
  plot_layout(nrow = 1, widths = c(0.8, 1)) +
  plot_annotation(tag_levels = "A")

# save
ggsave("output/coexistence_outcomes.png", coex_fig,
       width = 6.5, height = 6.5)

# table
inv_gfung_sum %>%
  mutate(parameter_set = 1) %>%
  full_join(inv_ginf_sum %>%
              mutate(parameter_set = 2)) %>%
  full_join(inv_infP_sum %>%
              mutate(parameter_set = 3)) %>%
  mutate(perc = janitor::round_half_up(prop * 100, 1),
         outcome = str_replace(outcome, "<br>", " ") %>% 
           str_remove_all("\\*") %>% 
           fct_relevel("coexist", "M. vimineum only", "E. virginicus only",
                       "priority effect", "non-establishment")) %>%
  select(-c(n, prop)) %>%
  pivot_wider(names_from = parameter_set, values_from = perc) %>%
  arrange(disease_name, outcome) %>%
  write_csv("output/coexistence_outcomes_all_parms.csv")


#### responses ####

# output lists
resp_gfung_h <- list()
resp_gfung_d <- list()

resp_ginf_h <- list()
resp_ginf_d <- list()

resp_infP_h <- list()
resp_infP_d <- list()

# calculate inflation amounts
L_add_gfung <- 0.001 * mean(eq_gfung$litter_star)
L_add_ginf <- 0.001 * mean(eq_ginf$litter_star)
L_add_infP <- 0.001 * mean(eq_infP$litter_star)

CA_add_gfung <- 0.001 * mean(eq_gfung$C_A_star)
CA_add_ginf <- 0.001 * mean(eq_ginf$C_A_star)
CA_add_infP <- 0.001 * mean(eq_infP$C_A_star)

CP_add_gfung <- 0.001 * mean(eq_gfung$C_P_star)
CP_add_ginf <- 0.001 * mean(eq_ginf$C_P_star)
CP_add_infP <- 0.001 * mean(eq_infP$C_P_star)

for(i in 1:params_iters){
  
  # get draw from parameters
  draw_gfung_h <- params_gfung[["healthy"]][["draws"]]$.draw[i]
  draw_gfung_d <- params_gfung[["disease"]][["draws"]]$.draw[i]
  draw_ginf_h <- params_ginf[["healthy"]][["draws"]]$.draw[i]
  draw_ginf_d <- params_ginf[["disease"]][["draws"]]$.draw[i]
  draw_infP_h <- params_infP[["healthy"]][["draws"]]$.draw[i] 
  draw_infP_d <- params_infP[["disease"]][["draws"]]$.draw[i] 
  
  # select equilibrium values
  eq_gfung_draw_h <- eq_gfung %>%
    filter(.draw == draw_gfung_h & disease == "h")
  eq_gfung_draw_d <- eq_gfung %>%
    filter(.draw == draw_gfung_d & disease == "d")

  eq_ginf_draw_h <- eq_ginf %>%
    filter(.draw == draw_ginf_h & disease == "h")
  eq_ginf_draw_d <- eq_ginf %>%
    filter(.draw == draw_ginf_d & disease == "d")
  
  eq_infP_draw_h <- eq_infP %>% 
    filter(.draw == draw_infP_h & disease == "h")
  eq_infP_draw_d <- eq_infP %>% 
    filter(.draw == draw_infP_d & disease == "d")
  
  # calculate responses for each parameter set and disease condition
  resp_gfung_h[[i]] <- resp_fun(iter = i, sim_outcome = eq_gfung_draw_h,
                                parameters = params_gfung[["healthy"]],
                                L_add = L_add_gfung, CA_add = CA_add_gfung,
                                CP_add = CP_add_gfung)
  resp_gfung_d[[i]] <- resp_fun(iter = i, sim_outcome = eq_gfung_draw_d,
                                parameters = params_gfung[["disease"]],
                                L_add = L_add_gfung, CA_add = CA_add_gfung,
                                CP_add = CP_add_gfung)

  resp_ginf_h[[i]] <- resp_fun(iter = i, sim_outcome = eq_ginf_draw_h,
                                parameters = params_ginf[["healthy"]],
                               L_add = L_add_ginf, CA_add = CA_add_ginf,
                               CP_add = CP_add_ginf)
  resp_ginf_d[[i]] <- resp_fun(iter = i, sim_outcome = eq_ginf_draw_d,
                                parameters = params_ginf[["disease"]],
                               L_add = L_add_ginf, CA_add = CA_add_ginf,
                               CP_add = CP_add_ginf)
  
  resp_infP_h[[i]] <- resp_fun(iter = i, sim_outcome = eq_infP_draw_h, 
                               parameters = params_infP[["healthy"]],
                               L_add = L_add_infP, CA_add = CA_add_infP,
                               CP_add = CP_add_infP)
  resp_infP_d[[i]] <- resp_fun(iter = i, sim_outcome = eq_infP_draw_d, 
                               parameters = params_infP[["disease"]],
                               L_add = L_add_infP, CA_add = CA_add_infP,
                               CP_add = CP_add_infP)
  
}

# combine healthy and disease
resp_gfung <- resp_gfung_h %>% 
  bind_rows() %>% 
  mutate(disease = "h") %>% 
  rbind(bind_rows(resp_gfung_d) %>% 
          mutate(disease = "d")) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

resp_ginf <- resp_ginf_h %>% 
  bind_rows() %>% 
  mutate(disease = "h") %>% 
  rbind(bind_rows(resp_ginf_d) %>% 
          mutate(disease = "d")) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

resp_infP <- resp_infP_h %>% 
  bind_rows() %>% 
  mutate(disease = "h") %>% 
  rbind(bind_rows(resp_infP_d) %>% 
          mutate(disease = "d")) %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

# save
save(resp_gfung, file = "output/responses_gfung_parms_20260722.rda")
save(resp_ginf, file = "output/responses_ginf_parms_20260722.rda")
save(resp_infP, file = "output/responses_infP_parms_20260722.rda")


#### response figures ####

# reload if needed
load("output/responses_gfung_parms_20260722.rda")
load("output/responses_ginf_parms_20260722.rda")
load("output/responses_infP_parms_20260722.rda")

# remove non-establishment parameter combinations
# make long by species
resp_gfung2 <- resp_gfung %>% 
  anti_join(non_gfung_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, starts_with("resp")) %>% 
  pivot_longer(cols = starts_with("resp"),
               names_to = c("competition", "species"),
               names_pattern = "resp_(.)(.)") %>% 
  mutate(competition = fct_recode(competition,
                                  "litter" = "L",
                                  "density" = "C"),
         species = fct_recode(species,
                              "M. vimineum" = "A",
                              "E. virginicus" = "P")) %>% 
  group_by(disease, competition, species) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            .groups = "drop")

resp_ginf2 <- resp_ginf %>% 
  anti_join(non_ginf_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, starts_with("resp")) %>% 
  pivot_longer(cols = starts_with("resp"),
               names_to = c("competition", "species"),
               names_pattern = "resp_(.)(.)") %>% 
  mutate(competition = fct_recode(competition,
                                  "litter" = "L",
                                  "density" = "C"),
         species = fct_recode(species,
                              "M. vimineum" = "A",
                              "E. virginicus" = "P")) %>% 
  group_by(disease, competition, species) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            .groups = "drop")

resp_infP2 <- resp_infP %>% 
  anti_join(non_infP_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, starts_with("resp")) %>% 
  pivot_longer(cols = starts_with("resp"),
               names_to = c("competition", "species"),
               names_pattern = "resp_(.)(.)") %>% 
  mutate(competition = fct_recode(competition,
                                  "litter" = "L",
                                  "density" = "C"),
         species = fct_recode(species,
                              "M. vimineum" = "A",
                              "E. virginicus" = "P")) %>% 
  group_by(disease, competition, species) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            .groups = "drop")

# gfung figures
resp_gfung_h_fig <- resp_gfung2 %>% 
  filter(disease == "h") %>% 
  ggplot(aes(x = competition, y = median, color = species)) +
  geom_errorbar(aes(ymin = L95, ymax = U95),
                linewidth = 0.3, width = 0,
                position = position_dodge(dodge_width)) +
  geom_point(aes(shape = species), size = 2,
             position = position_dodge(dodge_width)) +
  scale_y_log10(limits = c(min(resp_gfung2$L95),
                           max(resp_gfung2$U95))) +
  scale_color_manual(values = spp_pal) +
  scale_shape_manual(values = spp_shape_pal) +
  labs(x = "Competition type", 
       y = expression("Response ("*log[10]*")"),
       title = "Disease suppressed") +
  fig_theme +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

resp_gfung_d_fig <- resp_gfung_h_fig +
  filter(resp_gfung2, disease == "d") +
  labs(title = "Ambient disease")

resp_gfung_fig <- resp_gfung_h_fig + resp_gfung_d_fig + 
  plot_layout(nrow = 1, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# ginf figures
resp_ginf_h_fig <- resp_gfung_h_fig +
  filter(resp_ginf2, disease == "h") +
  scale_y_log10(limits = c(min(resp_ginf2$L95),
                           max(resp_ginf2$U95)))

resp_ginf_d_fig <- resp_gfung_d_fig +
  filter(resp_ginf2, disease == "d") +
  scale_y_log10(limits = c(min(resp_ginf2$L95),
                           max(resp_ginf2$U95)))

resp_ginf_fig <- resp_ginf_h_fig + resp_ginf_d_fig + 
  plot_layout(nrow = 1, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# infP figures
resp_infP_h_fig <- resp_gfung_h_fig +
  filter(resp_infP2, disease == "h") +
  scale_y_log10(limits = c(min(resp_infP2$L95),
                           max(resp_infP2$U95)))

resp_infP_d_fig <- resp_gfung_d_fig +
  filter(resp_infP2, disease == "d") +
  scale_y_log10(limits = c(min(resp_infP2$L95),
                           max(resp_infP2$U95)))

resp_infP_fig <- resp_infP_h_fig + resp_infP_d_fig + 
  plot_layout(nrow = 1, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# save
ggsave("output/responses_gfung_parms.png", resp_gfung_fig,
       width = 6, height = 2.7)
ggsave("output/responses_ginf_parms.png", resp_ginf_fig,
       width = 6, height = 2.7)
ggsave("output/responses_infP_parms.png", resp_infP_fig,
       width = 6, height = 2.7)


##### litter response distributions - delete or appendix #####

mv_litter_resp_dist <- ggplot(resp_gfung2, aes(x = resp_LA)) +
  stat_slab(aes(color = disease_name, fill = disease_name), alpha = 0.5) +
  scale_fill_manual(values = c("goldenrod", "#7C9B5B")) +
  scale_color_manual(values = c("goldenrod", "#7C9B5B")) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "*M. vimineum* response to litter-mediated competition") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

ev_litter_resp_dist <- mv_litter_resp_dist +
  aes(x = resp_LP) +
  labs(x = "*E. virginicus* response to litter-mediated competition")

mv_litter_resp_dist2 <- mv_litter_resp_dist +
  resp_ginf2

ev_litter_resp_dist2 <- ev_litter_resp_dist +
  resp_ginf2

mv_litter_resp_dist3 <- mv_litter_resp_dist +
  resp_infP2

ev_litter_resp_dist3 <- ev_litter_resp_dist +
  resp_infP2


##### competition response distributions - delete or appendix #####

mv_comp_resp_dist <- mv_litter_resp_dist +
  aes(x = resp_CA) +
  labs(x = "*M. vimineum* response to density-mediated competition")

ev_comp_resp_dist <- ev_litter_resp_dist +
  aes(x = resp_CP) +
  labs(x = "*E. virginicus* response to density-mediated competition")

mv_comp_resp_dist2 <- mv_comp_resp_dist +
  resp_ginf2

ev_comp_resp_dist2 <- ev_comp_resp_dist +
  resp_ginf2

mv_comp_resp_dist3 <- mv_comp_resp_dist +
  resp_infP2

ev_comp_resp_dist3 <- ev_comp_resp_dist +
  resp_infP2


##### litter response points - not updated, delete or appendix #####
mv_litter_resp_pt <- ggplot(resp_gfung2, 
                            aes(x = disease_name, y = resp_LA)) +
  stat_pointinterval(aes(color = disease_name), point_size = 2,
                     point_interval = median_hdci, .width = c(.66, .95)) +
  scale_color_manual(values = c("goldenrod", "#7C9B5B"), guide = "none") +
  # scale_y_continuous(n.breaks = 3) +
  fig_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))

ev_litter_resp_pt <- mv_litter_resp_pt %+%
  aes(y = resp_LP)

mv_comp_resp_pt <- mv_litter_resp_pt %+%
  aes(y = resp_CA)

ev_comp_resp_pt <- ev_litter_resp_pt %+%
  aes(y = resp_CP)

mv_litter_resp_pt2 <- mv_litter_resp_pt %+%
  resp_ginf_long

ev_litter_resp_pt2 <- mv_litter_resp_pt2 %+%
  aes(y = resp_LP)

mv_comp_resp_pt2 <- mv_litter_resp_pt2 %+%
  aes(y = resp_CA)

ev_comp_resp_pt2 <- mv_litter_resp_pt2 %+%
  aes(y = resp_CP)

mv_litter_resp_pt3 <- mv_litter_resp_pt %+%
  resp_infP_long

ev_litter_resp_pt3 <- mv_litter_resp_pt3 %+%
  aes(y = resp_LP)

mv_comp_resp_pt3 <- mv_litter_resp_pt3 %+%
  aes(y = resp_CA)

ev_comp_resp_pt3 <- mv_litter_resp_pt3 %+%
  aes(y = resp_CP)

# combine
resp_inset_coord <- c(0.25, 0.5, 0.75, 0.99)

mv_litter_resp_fig <- mv_litter_resp_dist + 
  inset_element(mv_litter_resp_pt, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_litter_resp_fig <- ev_litter_resp_dist + 
  inset_element(ev_litter_resp_pt, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
mv_comp_resp_fig <- mv_comp_resp_dist + 
  inset_element(mv_comp_resp_pt, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_comp_resp_fig <- ev_comp_resp_dist + 
  inset_element(ev_comp_resp_pt, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])

mv_litter_resp_fig2 <- mv_litter_resp_dist2 + 
  inset_element(mv_litter_resp_pt2, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_litter_resp_fig2 <- ev_litter_resp_dist2 + 
  inset_element(ev_litter_resp_pt2, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
mv_comp_resp_fig2 <- mv_comp_resp_dist2 + 
  inset_element(mv_comp_resp_pt2, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_comp_resp_fig2 <- ev_comp_resp_dist2 + 
  inset_element(ev_comp_resp_pt2, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])

mv_litter_resp_fig3 <- mv_litter_resp_dist3 + 
  inset_element(mv_litter_resp_pt3, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_litter_resp_fig3 <- ev_litter_resp_dist3 + 
  inset_element(ev_litter_resp_pt3, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
mv_comp_resp_fig3 <- mv_comp_resp_dist3 + 
  inset_element(mv_comp_resp_pt3, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])
ev_comp_resp_fig3 <- ev_comp_resp_dist3 + 
  inset_element(ev_comp_resp_pt3, resp_inset_coord[1], resp_inset_coord[2],
                resp_inset_coord[3], resp_inset_coord[4])


#### effect figures ####

# remove non-establishment parameter combinations
# make long by species
eq_gfung2 <- eq_gfung %>% 
  anti_join(non_gfung_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, ends_with("effect")) %>% 
  pivot_longer(cols = ends_with("effect"),
               names_to = "effect") %>% 
  mutate(competition = case_when(
    str_detect(effect, "litter") ~ "litter",
    str_detect(effect, "AA|PP") ~ "intrasp_density",
    str_detect(effect, "AP|PA") ~ "intersp_density"),
    species = case_when(
      str_detect(effect, "litter_A|C_AA|C_PA") ~ "M. vimineum",
      str_detect(effect, "litter_P|C_PP|C_AP") ~ "E. virginicus")) %>% 
  select(-effect) %>% 
  group_by(disease, species, competition) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            min = min(value),
            max = max(value),
            .groups = "drop") %>% 
  pivot_wider(names_from = competition, 
              values_from = c(median, L95, U95, min, max),
              names_glue = "{competition}_{.value}")

eq_ginf2 <- eq_ginf %>% 
  anti_join(non_ginf_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, ends_with("effect")) %>% 
  pivot_longer(cols = ends_with("effect"),
               names_to = "effect") %>% 
  mutate(competition = case_when(
    str_detect(effect, "litter") ~ "litter",
    str_detect(effect, "AA|PP") ~ "intrasp_density",
    str_detect(effect, "AP|PA") ~ "intersp_density"),
    species = case_when(
      str_detect(effect, "litter_A|C_AA|C_PA") ~ "M. vimineum",
      str_detect(effect, "litter_P|C_PP|C_AP") ~ "E. virginicus")) %>% 
  select(-effect) %>% 
  group_by(disease, species, competition) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            min = min(value),
            max = max(value),
            .groups = "drop") %>% 
  pivot_wider(names_from = competition, 
              values_from = c(median, L95, U95, min, max),
              names_glue = "{competition}_{.value}")

eq_infP2 <- eq_infP %>% 
  anti_join(non_infP_A) %>% 
  anti_join(non_gfung_P) %>% 
  select(.draw, disease, ends_with("effect")) %>% 
  pivot_longer(cols = ends_with("effect"),
               names_to = "effect") %>% 
  mutate(competition = case_when(
    str_detect(effect, "litter") ~ "litter",
    str_detect(effect, "AA|PP") ~ "intrasp_density",
    str_detect(effect, "AP|PA") ~ "intersp_density"),
    species = case_when(
      str_detect(effect, "litter_A|C_AA|C_PA") ~ "M. vimineum",
      str_detect(effect, "litter_P|C_PP|C_AP") ~ "E. virginicus")) %>% 
  select(-effect) %>% 
  group_by(disease, species, competition) %>% 
  summarize(median = median(value),
            L95 = median_hdci(value)$ymin,
            U95 = median_hdci(value)$ymax,
            min = min(value),
            max = max(value),
            .groups = "drop") %>% 
  pivot_wider(names_from = competition, 
              values_from = c(median, L95, U95, min, max),
              names_glue = "{competition}_{.value}")

# gfung figures
eff_gfung_hle_fig <- eq_gfung2 %>% 
  filter(disease == "h") %>% 
  ggplot(aes(x = litter_median, y = intersp_density_median, color = species)) +
  geom_errorbar(aes(xmin = litter_L95, xmax = litter_U95), 
                linewidth = 0.3, width = 0) +
  geom_errorbar(aes(ymin = intersp_density_L95, ymax = intersp_density_U95),
                linewidth = 0.3, width = 0) +
  geom_point(aes(shape = species), size = 2) +
  scale_color_manual(values = spp_pal) +
  scale_shape_manual(values = spp_shape_pal) +
  scale_x_continuous(limits = c(min(eq_gfung2$litter_L95),
                                max(eq_gfung2$litter_U95))) +
  scale_y_continuous(limits = c(min(eq_gfung2$intersp_density_L95),
                                max(eq_gfung2$intersp_density_U95))) +
  labs(x = "Litter-mediated competition", 
       y = "Interspecific density-mediated\ncompetition",
       title = "Disease suppressed") +
  fig_theme +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

eff_gfung_dle_fig <- eff_gfung_hle_fig +
  filter(eq_gfung2, disease == "d") +
  labs(title = "Ambient disease")

eff_gfung_hea_fig <- eq_gfung2 %>% 
  filter(disease == "h") %>% 
  ggplot(aes(x = intersp_density_median, y = intrasp_density_median, color = species)) +
  geom_errorbar(aes(xmin = intersp_density_L95, xmax = intersp_density_U95),
                linewidth = 0.3, width = 0) +
  geom_errorbar(aes(ymin = intrasp_density_L95, ymax = intrasp_density_U95),
                linewidth = 0.3, width = 0) +
  geom_point(aes(shape = species), size = 2) +
  scale_color_manual(values = spp_pal) +
  scale_shape_manual(values = spp_shape_pal) +
  scale_x_continuous(limits = c(min(eq_gfung2$intersp_density_L95),
                                max(eq_gfung2$intersp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_gfung2$intrasp_density_L95),
                                max(eq_gfung2$intrasp_density_U95))) +
  labs(x = "Interspecific density-mediated competition",
       y = "Intraspecific density-mediated\ncompetition") +
  fig_theme +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

eff_gfung_dea_fig <- eff_gfung_hea_fig +
  filter(eq_gfung2, disease == "d")

eff_gfung_hal_fig <- eq_gfung2 %>% 
  filter(disease == "h") %>% 
  ggplot(aes(x = intrasp_density_median, y = litter_median, color = species)) +
  geom_errorbar(aes(xmin = intrasp_density_L95, xmax = intrasp_density_U95),
                linewidth = 0.3, width = 0) +
  geom_errorbar(aes(ymin = litter_L95, ymax = litter_U95), 
                linewidth = 0.3, width = 0) +
  geom_point(aes(shape = species), size = 2) +
  scale_color_manual(values = spp_pal) +
  scale_shape_manual(values = spp_shape_pal) +
  scale_x_continuous(limits = c(min(eq_gfung2$intrasp_density_L95),
                                max(eq_gfung2$intrasp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_gfung2$litter_L95),
                                max(eq_gfung2$litter_U95))) +
  labs(x = "Intraspecific density-mediated competition",
       y = "Litter-mediated competition") +
  fig_theme +
  theme(legend.title = element_blank(),
        legend.text = element_text(face = "italic"))

eff_gfung_dal_fig <- eff_gfung_hal_fig +
  filter(eq_gfung2, disease == "d")

eff_gfung_fig <- eff_gfung_hle_fig + eff_gfung_dle_fig + 
  eff_gfung_hea_fig + eff_gfung_dea_fig + 
  eff_gfung_hal_fig + eff_gfung_dal_fig + 
  plot_layout(nrow = 3, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# ginf figures
eff_ginf_hle_fig <- eff_gfung_hle_fig +
  filter(eq_ginf2, disease == "h") +
  scale_x_continuous(limits = c(min(eq_ginf2$litter_L95),
                                max(eq_ginf2$litter_U95))) +
  scale_y_continuous(limits = c(min(eq_ginf2$intersp_density_L95),
                                max(eq_ginf2$intersp_density_U95)))

eff_ginf_dle_fig <- eff_ginf_hle_fig +
  filter(eq_ginf2, disease == "d") +
  labs(title = "Ambient disease")

eff_ginf_hea_fig <- eff_gfung_hea_fig +
  filter(eq_ginf2, disease == "h") + 
  scale_x_continuous(limits = c(min(eq_ginf2$intersp_density_L95),
                                max(eq_ginf2$intersp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_ginf2$intrasp_density_L95),
                                max(eq_ginf2$intrasp_density_U95)))

eff_ginf_dea_fig <- eff_ginf_hea_fig +
  filter(eq_ginf2, disease == "d")

eff_ginf_hal_fig <- eff_gfung_hal_fig +
  filter(eq_ginf2, disease == "h") +
  scale_x_continuous(limits = c(min(eq_ginf2$intrasp_density_L95),
                                max(eq_ginf2$intrasp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_ginf2$litter_L95),
                                max(eq_ginf2$litter_U95)))

eff_ginf_dal_fig <- eff_ginf_hal_fig +
  filter(eq_ginf2, disease == "d")

eff_ginf_fig <- eff_ginf_hle_fig + eff_ginf_dle_fig + 
  eff_ginf_hea_fig + eff_ginf_dea_fig + 
  eff_ginf_hal_fig + eff_ginf_dal_fig + 
  plot_layout(nrow = 3, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# infP figures
eff_infP_hle_fig <- eff_gfung_hle_fig +
  filter(eq_infP2, disease == "h") +
  scale_x_continuous(limits = c(min(eq_infP2$litter_L95),
                                max(eq_infP2$litter_U95))) +
  scale_y_continuous(limits = c(min(eq_infP2$intersp_density_L95),
                                max(eq_infP2$intersp_density_U95)))

eff_infP_dle_fig <- eff_infP_hle_fig +
  filter(eq_infP2, disease == "d") +
  labs(title = "Ambient disease")

eff_infP_hea_fig <- eff_gfung_hea_fig +
  filter(eq_infP2, disease == "h") + 
  scale_x_continuous(limits = c(min(eq_infP2$intersp_density_L95),
                                max(eq_infP2$intersp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_infP2$intrasp_density_L95),
                                max(eq_infP2$intrasp_density_U95)))

eff_infP_dea_fig <- eff_infP_hea_fig +
  filter(eq_infP2, disease == "d")

eff_infP_hal_fig <- eff_gfung_hal_fig +
  filter(eq_infP2, disease == "h") +
  scale_x_continuous(limits = c(min(eq_infP2$intrasp_density_L95),
                                max(eq_infP2$intrasp_density_U95))) +
  scale_y_continuous(limits = c(min(eq_infP2$litter_L95),
                                max(eq_infP2$litter_U95)))

eff_infP_dal_fig <- eff_infP_hal_fig +
  filter(eq_infP2, disease == "d")

eff_infP_fig <- eff_infP_hle_fig + eff_infP_dle_fig + 
  eff_infP_hea_fig + eff_infP_dea_fig + 
  eff_infP_hal_fig + eff_infP_dal_fig + 
  plot_layout(nrow = 3, guides = "collect", axes = "collect") +
  plot_annotation(tag_levels = "A")  &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal") 

# save
ggsave("output/effects_gfung_parms.png", eff_gfung_fig,
       width = 6, height = 8)
ggsave("output/effects_ginf_parms.png", eff_ginf_fig,
       width = 6, height = 8)
ggsave("output/effects_infP_parms.png", eff_infP_fig,
       width = 6, height = 8)


##### litter effect distributions - delete or appendix ######
mv_litter_eff_dist <- ggplot(eq_gfung, aes(x = litter_A_effect)) +
  stat_slab(aes(color = disease_name, fill = disease_name), alpha = 0.5) +
  scale_fill_manual(values = dis_pal) +
  scale_color_manual(values = dis_pal) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "*M. vimineum* effect on litter-mediated competition") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank())

ev_litter_eff_dist <- mv_litter_eff_dist +
  aes(x = litter_P_effect) +
  labs(x = "*E. virginicus* effect on litter-mediated competition")

mv_litter_eff_dist2 <- mv_litter_eff_dist +
  eq_ginf

ev_litter_eff_dist2 <- ev_litter_eff_dist +
  eq_ginf

mv_litter_eff_dist3 <- mv_litter_eff_dist +
  eq_infP

ev_litter_eff_dist3 <- ev_litter_eff_dist +
  eq_infP

##### competition effect distributions - delete or appendix #####
mv_mv_comp_eff_dist <- mv_litter_eff_dist +
  aes(x = C_AA_effect) +
  labs(x = "*M. vimineum* effect on density-mediated intraspecific competition")

ev_mv_comp_eff_dist <- mv_litter_eff_dist +
  aes(x = C_PA_effect) +
  labs(x = "*M. vimineum* effect on density-mediated interspecific competition")

ev_ev_comp_eff_dist <- ev_litter_eff_dist +
  aes(x = C_PP_effect) +
  labs(x = "*E. virginicus* effect on density-mediated intraspecific competition")

mv_ev_comp_eff_dist <- ev_litter_eff_dist +
  aes(x = C_AP_effect) +
  labs(x = "*E. virginicus* effect on density-mediated interspecific competition")

mv_mv_comp_eff_dist2 <- mv_mv_comp_eff_dist +
  eq_ginf

ev_mv_comp_eff_dist2 <- ev_mv_comp_eff_dist +
  eq_ginf

ev_ev_comp_eff_dist2 <- ev_ev_comp_eff_dist +
  eq_ginf

mv_ev_comp_eff_dist2 <- mv_ev_comp_eff_dist +
  eq_ginf

mv_mv_comp_eff_dist3 <- mv_mv_comp_eff_dist +
  eq_infP

ev_mv_comp_eff_dist3 <- ev_mv_comp_eff_dist +
  eq_infP

ev_ev_comp_eff_dist3 <- ev_ev_comp_eff_dist +
  eq_infP

mv_ev_comp_eff_dist3 <- mv_ev_comp_eff_dist +
  eq_infP

##### litter effect points - delete or appendix #####
mv_litter_eff_pt <- ggplot(eq_gfung, 
                           aes(x = disease_name, y = litter_A_effect)) +
  stat_pointinterval(aes(color = disease_name), point_size = 2,
                     point_interval = median_hdci, .width = c(.66, .95)) +
  scale_color_manual(values = c("goldenrod", "#7C9B5B"), guide = "none") +
  scale_y_continuous(labels = scales::comma, n.breaks = 3) +
  labs(y = "*M. vimineum* effect on litter-mediated competition") +
  fig_theme +
  theme(axis.title.y = element_markdown(),
        axis.title.x = element_blank())

ev_litter_eff_pt <- mv_litter_eff_pt +
  aes(y = litter_P_effect) +
  labs(y = "*E. virginicus* effect on litter-mediated competition")

mv_litter_eff_pt2 <- mv_litter_eff_pt +
  eq_ginf

ev_litter_eff_pt2 <- ev_litter_eff_pt  +
  eq_ginf

mv_litter_eff_pt3 <- mv_litter_eff_pt +
  eq_infP

ev_litter_eff_pt3 <- ev_litter_eff_pt +
  eq_infP

##### competition effect points - delete or appendix #####
mv_mv_comp_eff_pt <- mv_litter_eff_pt +
  aes(y = C_AA_effect) +
  labs(y = "*M. vimineum* effect on density-mediated intraspecific competition")

ev_mv_comp_eff_pt <- mv_litter_eff_pt +
  aes(y = C_PA_effect) +
  labs(y = "*M. vimineum* effect on density-mediated interspecific competition")

ev_ev_comp_eff_pt <- ev_litter_eff_pt +
  aes(y = C_PP_effect) +
  labs(y = "*E. virginicus* effect on density-mediated intraspecific competition")

mv_ev_comp_eff_pt <- ev_litter_eff_pt +
  aes(y = C_AP_effect) +
  labs(y = "*E. virginicus* effect on density-mediated interspecific competition")

mv_mv_comp_eff_pt2 <- mv_mv_comp_eff_pt +
  eq_ginf

ev_mv_comp_eff_pt2 <- ev_mv_comp_eff_pt +
  eq_ginf

ev_ev_comp_eff_pt2 <- ev_ev_comp_eff_pt +
  eq_ginf

mv_ev_comp_eff_pt2 <- mv_ev_comp_eff_pt +
  eq_ginf

mv_mv_comp_eff_pt3 <- mv_mv_comp_eff_pt +
  eq_infP

ev_mv_comp_eff_pt3 <- ev_mv_comp_eff_pt +
  eq_infP

ev_ev_comp_eff_pt3 <- ev_ev_comp_eff_pt +
  eq_infP

mv_ev_comp_eff_pt3 <- mv_ev_comp_eff_pt +
  eq_infP


#### combine litter and density figuresn - delete or appendix ####

eff_resp_litter_fig <- mv_litter_eff_fig + ev_litter_eff_fig + mv_litter_resp_fig + 
  ev_litter_resp_fig +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

eff_resp_comp_fig <- mv_comp_eff_fig + ev_comp_eff_fig + mv_comp_resp_fig + 
  ev_comp_resp_fig +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

eff_resp_litter_fig2 <- mv_litter_eff_fig2 + ev_litter_eff_fig + 
  mv_litter_resp_fig2 + ev_litter_resp_fig2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

eff_resp_comp_fig2 <- mv_comp_eff_fig2 + ev_comp_eff_fig + 
  mv_comp_resp_fig2 + ev_comp_resp_fig2 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

eff_resp_litter_fig3 <- mv_litter_eff_fig3 + ev_litter_eff_fig + 
  mv_litter_resp_fig3 + ev_litter_resp_fig3 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

eff_resp_comp_fig3 <- mv_comp_eff_fig3 + ev_comp_eff_fig + 
  mv_comp_resp_fig3 + ev_comp_resp_fig3 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave("output/litter_effects_responses_gfung_parms.png", eff_resp_litter_fig,
       width = 6.5, height = 6.5)
ggsave("output/competition_effects_responses_gfung_parms.png", eff_resp_comp_fig,
       width = 6.5, height = 6.5)
ggsave("output/litter_effects_responses_ginf_parms.png", eff_resp_litter_fig2,
       width = 6.5, height = 6.5)
ggsave("output/competition_effects_responses_ginf_parms.png", eff_resp_comp_fig2,
       width = 6.5, height = 6.5)
ggsave("output/litter_effects_responses_infP_parms.png", eff_resp_litter_fig3,
       width = 6.5, height = 6.5)
ggsave("output/competition_effects_responses_infP_parms.png", eff_resp_comp_fig3,
       width = 6.5, height = 6.5)


#### effects/responses values for text ####

# effect on litter
eq_infP_A %>%
  group_by(disease_name) %>%
  median_hdci(litter)

eq_gfung_P %>%
  group_by(disease_name) %>%
  median_hdci(litter)

# treatment difference in effect on litter
eq_infP_A %>%
  select(.draw, disease_name, litter) %>%
  pivot_wider(names_from = disease_name,
              values_from = litter) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

eq_gfung_P %>%
  select(.draw, disease_name, litter) %>%
  pivot_wider(names_from = disease_name,
              values_from = litter) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

# response to litter
resp_infP_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_LA)

resp_infP_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_LP)

# treatment difference in response to litter
resp_infP_long %>%
  select(.draw, disease_name, resp_LA) %>%
  pivot_wider(names_from = disease_name,
              values_from = resp_LA) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

resp_infP_long %>%
  select(.draw, disease_name, resp_LP) %>%
  pivot_wider(names_from = disease_name,
              values_from = resp_LP) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

# effect on competition
eq_infP_A %>%
  group_by(disease_name) %>%
  median_hdci(comp_eff)

eq_gfung_P %>%
  group_by(disease_name) %>%
  median_hdci(comp_eff)

# treatment difference in effect on competition
eq_infP_A %>%
  select(.draw, disease_name, comp_eff) %>%
  pivot_wider(names_from = disease_name,
              values_from = comp_eff) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

eq_gfung_P %>%
  select(.draw, disease_name, comp_eff) %>%
  pivot_wider(names_from = disease_name,
              values_from = comp_eff) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

# response to competition
resp_infP_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_CA)

resp_infP_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_CP)

# treatment difference in response to competition
resp_infP_long %>%
  select(.draw, disease_name, resp_CA) %>%
  pivot_wider(names_from = disease_name,
              values_from = resp_CA) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)

resp_infP_long %>%
  select(.draw, disease_name, resp_CP) %>%
  pivot_wider(names_from = disease_name,
              values_from = resp_CP) %>%
  mutate(diff = `Disease suppressed` - `Ambient disease`) %>%
  median_hdci(diff)


#### effect/response table ####

# effect on litter
litter_eff_tab <- eq_gfung_A %>%
  group_by(disease_name) %>%
  median_hdci(litter) %>% 
  mutate(parameter_set = 1) %>% 
  full_join(eq_ginf_A %>%
              group_by(disease_name) %>%
              median_hdci(litter) %>% 
              mutate(parameter_set = 2)) %>% 
  full_join(eq_infP_A %>%
              group_by(disease_name) %>%
              median_hdci(litter) %>% 
              mutate(parameter_set = 3)) %>% 
  mutate(species = "M. vimineum") %>% 
  full_join(eq_gfung_P %>%
              group_by(disease_name) %>%
              median_hdci(litter) %>% 
              cross_join(tibble(parameter_set = c(1, 2, 3))) %>% 
              mutate(species = "E. virginicus")) %>% 
  select(disease_name, species, litter, .lower, .upper, parameter_set) %>% 
  rename(var = litter) %>% 
  mutate(type = "effect",
         competition = "litter")

# response to litter
litter_resp_tab <- resp_gfung_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_LA) %>% 
  mutate(parameter_set = 1) %>% 
  full_join(resp_ginf_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_LA)  %>% 
              mutate(parameter_set = 2)) %>% 
  full_join(resp_infP_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_LA)  %>% 
              mutate(parameter_set = 3)) %>% 
  mutate(species = "M. vimineum") %>% 
  rename(var = resp_LA) %>% 
  full_join(resp_gfung_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_LP) %>% 
              mutate(parameter_set = 1) %>% 
              full_join(resp_ginf_long %>%
                          group_by(disease_name) %>%
                          median_hdci(resp_LP)  %>% 
                          mutate(parameter_set = 2)) %>% 
              full_join(resp_infP_long %>%
                          group_by(disease_name) %>%
                          median_hdci(resp_LP)  %>% 
                          mutate(parameter_set = 3)) %>% 
              mutate(species = "E. virginicus") %>% 
              rename(var = resp_LP)) %>% 
  select(disease_name, species, var, .lower, .upper, parameter_set) %>% 
  mutate(type = "response",
         competition = "litter")

# effect on density
density_eff_tab <- eq_gfung_A %>%
  group_by(disease_name) %>%
  median_hdci(comp_eff) %>% 
  mutate(parameter_set = 1) %>% 
  full_join(eq_ginf_A %>%
              group_by(disease_name) %>%
              median_hdci(comp_eff) %>% 
              mutate(parameter_set = 2)) %>% 
  full_join(eq_infP_A %>%
              group_by(disease_name) %>%
              median_hdci(comp_eff) %>% 
              mutate(parameter_set = 3)) %>% 
  mutate(species = "M. vimineum") %>% 
  full_join(eq_gfung_P %>%
              group_by(disease_name) %>%
              median_hdci(comp_eff) %>% 
              cross_join(tibble(parameter_set = c(1, 2, 3))) %>% 
              mutate(species = "E. virginicus")) %>% 
  select(disease_name, species, comp_eff, .lower, .upper, parameter_set) %>% 
  rename(var = comp_eff) %>% 
  mutate(type = "effect",
         competition = "density")

# response to density
density_resp_tab <- resp_gfung_long %>%
  group_by(disease_name) %>%
  median_hdci(resp_CA) %>% 
  mutate(parameter_set = 1) %>% 
  full_join(resp_ginf_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_CA)  %>% 
              mutate(parameter_set = 2)) %>% 
  full_join(resp_infP_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_CA)  %>% 
              mutate(parameter_set = 3)) %>% 
  mutate(species = "M. vimineum") %>% 
  rename(var = resp_CA) %>% 
  full_join(resp_gfung_long %>%
              group_by(disease_name) %>%
              median_hdci(resp_CP) %>% 
              mutate(parameter_set = 1) %>% 
              full_join(resp_ginf_long %>%
                          group_by(disease_name) %>%
                          median_hdci(resp_CP)  %>% 
                          mutate(parameter_set = 2)) %>% 
              full_join(resp_infP_long %>%
                          group_by(disease_name) %>%
                          median_hdci(resp_CP)  %>% 
                          mutate(parameter_set = 3)) %>% 
              mutate(species = "E. virginicus") %>% 
              rename(var = resp_CP)) %>% 
  select(disease_name, species, var, .lower, .upper, parameter_set) %>% 
  mutate(type = "response",
         competition = "density")

# combine
eff_resp_tab <- litter_eff_tab %>% 
  full_join(litter_resp_tab) %>% 
  full_join(density_eff_tab) %>% 
  full_join(density_resp_tab) %>% 
  relocate(competition, type, species, .before = 1) %>% 
  mutate(across(.cols = c(var, .lower, .upper),
                .fns = ~case_when(abs(.x) > 1 ~ as.character(round_half_up(.x, 0)),
                                  abs(round_half_up(.x, 2)) >= 0.01 ~ 
                                    as.character(round_half_up(.x, 2)),
                                  .x == 0 ~ "0",
                                  TRUE ~
                                    format(.x, scientific = T, digits = 2)))) %>% 
  mutate(var_all = paste0(var, "\n(", .lower, " to ", .upper, ")")) %>% 
  select(-c(var, .lower, .upper)) %>% 
  pivot_wider(names_from = parameter_set,
              values_from = var_all)
  
write_csv(eff_resp_tab, "output/effects_responses_all_parms.csv")
