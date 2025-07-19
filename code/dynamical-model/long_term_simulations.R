#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidybayes)
library(ggtext)

# import parameters (loads tidyvers)
source("code/dynamical-model/parameters.R")

# import simulation function
source("code/dynamical-model/APL_Sim_Tree.R")

# figure settings
source("code/figure-prep/figure_settings.R")


#### functions ####

# combine healthy and disease simulations
sims_comb_fun <- function(sims_h, sims_d){
  
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
                values_from = c(annual_seeds, litter, 
                                perennial_seeds, perennial_adults,
                                annual_competition, perennial_competition))
  
  # output both formats
  return(list(long = out_long, wide = out_wide))

}

# calculate responses
resp_fun <- function(iter, sim_outcome, parameters){
  
  # when species coexist
  if(sim_outcome$coexist[iter] == "yes"){
    
    # initial conditions (annual seeds, litter, perennial seeds, perennial adults)
    inits <- c(sim_outcome$annual_seeds[iter],
               sim_outcome$litter[iter],
               sim_outcome$perennial_seeds[iter],
               sim_outcome$perennial_adults[iter])
    
    # initial competition values (experienced by annual, perennial)
    inits_comp <- c(sim_outcome$annual_competition[iter],
                    sim_outcome$perennial_competition[iter])
    
    # simulate one year under baseline conditions
    sim_base <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits,
                            parameters = parameters)
    
    # inflate competition factors
    inits_L2 <- inits
    inits_L2[2] <- sim_outcome$litter[iter] * 1.01
    
    inits_compA2 <- inits_comp
    inits_compA2[1] <- sim_outcome$annual_competition[iter] * 1.01
    
    inits_compP2 <- inits_comp
    inits_compP2[2] <- sim_outcome$perennial_competition[iter] * 1.01
    
    # simulate one year with inflation for each species
    sim_L <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_L2,
                         parameters = parameters)
    sim_compA <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits,
                             parameters = parameters, init_C = inits_compA2)
    sim_compP <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits,
                             parameters = parameters, init_C = inits_compP2)
    
    # calculate growth rates r = lnN(t+1) -  lnN(t)
    gr_baseA <- log(sim_base$annual_seeds[2] / sim_base$annual_seeds[1])
    gr_LA <- log(sim_L$annual_seeds[2] / sim_L$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    gr_baseP <- log((sim_base$perennial_seeds[2] + sim_base$perennial_adults[2]) / 
                      (sim_base$perennial_seeds[1] + sim_base$perennial_adults[1]))
    gr_LP <- log((sim_L$perennial_seeds[2] + sim_L$perennial_adults[2]) / 
                   (sim_L$perennial_seeds[1] + sim_L$perennial_adults[1]))
    gr_CP <- log((sim_compP$perennial_seeds[2] + sim_compP$perennial_adults[2]) / 
                   (sim_compP$perennial_seeds[1] + sim_compP$perennial_adults[1]))
    
    # calculate responses
    resp_LA <- (gr_LA - gr_baseA) / (inits2[2] - inits[2])
    resp_LP <- (gr_LP - gr_baseP) / (inits2[2] - inits[2])
    resp_CA <- (gr_CA - gr_baseA) / (inits_compA2[1] - inits_comp[1])
    resp_CP <- (gr_CP - gr_baseP) / (inits_compP2[2] - inits_comp[2])
    
  } else {
    
    # first letter is invader, second is resident
    # initial conditions (annual seeds, litter, perennial seeds, perennial adults)
    inits_AP <- c(1,
                  sim_outcome$litter_P[iter],
                  sim_outcome$perennial_seeds_P[iter],
                  sim_outcome$perennial_adults_P[iter])
    
    inits_PA <- c(sim_outcome$annual_seeds_A[iter],
                  sim_outcome$litter_A[iter],
                  0,
                  1)
    
    # initial competition values (experienced by annual, perennial)
    inits_compA <- c(sim_outcome$annual_competition_A[iter],
                     0)
    inits_compP <- c(0,
                     sim_outcome$perennial_competition_P[iter])
    
    # simulate one year under baseline conditions
    sim_base_AP <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_AP,
                               parameters = parameters)
    sim_base_PA <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_PA,
                               parameters = parameters)
    
    # inflate competition factors
    inits_LP2 <- inits_AP
    inits_LP2[2] <- sim_outcome$litter_P[iter] * 1.01
    
    inits_LA2 <- inits_PA
    inits_LA2[2] <- sim_outcome$litter_A[iter] * 1.01
    
    inits_compA2 <- inits_compA
    inits_compA2[1] <- sim_outcome$annual_competition_A[iter] * 1.01
    
    inits_compP2 <- inits_compP
    inits_compP2[2] <- sim_outcome$perennial_competition_P[iter] * 1.01

    # simulate one year with inflation for each species
    sim_LP <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_LP2,
                          parameters = parameters)
    sim_LA <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_LA2,
                          parameters = parameters)
    sim_compA <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_AP,
                             parameters = parameters, init_C = inits_compA2)
    sim_compP <- dyn_mod_fun(iter = iter, gen = 2, init_cond = inits_PA,
                             parameters = parameters, init_C = inits_compP2)
    
    # calculate growth rates r = lnN(t+1) -  lnN(t)
    gr_baseA <- log(sim_base_AP$annual_seeds[2] / sim_base_AP$annual_seeds[1])
    gr_LA <- log(sim_LAP$annual_seeds[2] / sim_LAP$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    gr_baseP <- log((sim_base_PA$perennial_seeds[2] + 
                       sim_base_PA$perennial_adults[2]) / 
                      (sim_base_PA$perennial_seeds[1] + 
                         sim_base_PA$perennial_adults[1]))
    gr_LP <- log((sim_LPA$perennial_seeds[2] + sim_LPA$perennial_adults[2]) / 
                   (sim_LPA$perennial_seeds[1] + sim_LPA$perennial_adults[1]))
    gr_CP <- log((sim_compP$perennial_seeds[2] + sim_compP$perennial_adults[2]) / 
                   (sim_compP$perennial_seeds[1] + sim_compP$perennial_adults[1]))
    
    # calculate responses
    resp_LA <- (gr_LA - gr_baseA) / (inits_LP2[2] - inits_AP[2])
    resp_LP <- (gr_LP - gr_baseP) / (inits_LA2[2] - inits_PA[2])
    resp_CA <- (gr_CA - gr_baseA) / (inits_compA2[1] - inits_compA[1])
    resp_CP <- (gr_CP - gr_baseP) / (inits_compP2[2] - inits_compP[2])
    
  }
  
  out <- tibble(
    .draw = parameters[["draws"]][iter, ]$.draw,
    resp_LA = resp_LA,
    resp_LP = resp_LP,
    resp_CA = resp_CA,
    resp_CP = resp_CP)
}


#### parameters ####

# parameter combinations
params_iters <- 10

# parameters with fungicide effect on germination
params_gfung <- params_fun(iters = params_iters, gA_type = "fungicide")

# parameters with infection effect on germination
params_ginf <- params_fun(iters = params_iters, gA_type = "infection")

# generations
gens <- 100


#### single species simulations ####

# initial conditions (annual seeds, litter, perennial seeds, perennial adults)
initsA <- c(1, 0, 0, 0)
initsP <- c(0, 0, 0, 1)

# output lits
sims_gfung_A_h <- list()
sims_gfung_P_h <- list()

sims_gfung_A_d <- list()
sims_gfung_P_d <- list()

sims_ginf_A_h <- list()

sims_ginf_A_d <- list()

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
  
}

# combine healthy and disease
sims_gfung_A <- sims_comb_fun(sims_gfung_A_h, sims_gfung_A_d)
sims_ginf_A <- sims_comb_fun(sims_ginf_A_h, sims_ginf_A_d)
sims_gfung_P <- sims_comb_fun(sims_gfung_P_h, sims_gfung_P_d)

# have the population stabilized? 
# (memory too limited to show all -> subsample)
sims_gfung_A[["long"]] %>% 
  # filter(iteration %in% sample(1:15000, 1000)) %>%
  ggplot(aes(x = generation, y = annual_seeds,
           group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

sims_ginf_A[["long"]] %>%
  # filter(iteration %in% sample(1:15000, 1000)) %>%
  ggplot(aes(x = generation, y = annual_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

sims_gfung_P[["long"]] %>%
  # filter(iteration %in% sample(1:15000, 1000)) %>% 
  ggplot(aes(x = generation, y = perennial_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

sims_gfung_P[["long"]] %>%
  # filter(iteration %in% sample(1:15000, 1000)) %>% 
  ggplot(aes(x = generation, y = perennial_adults,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)


#### single species figures ####

# mean hdci time series figures
(fig_dens_gfung_A <- ggplot(sims_gfung_A[["long"]], aes(x = generation, y = annual_seeds)) +
    stat_lineribbon(aes(fill = disease, color = disease), 
                    point_interval = mean_hdci, # multiple simulations were multi-modal 
                    .width = 0.95, alpha = 0.5) +
    scale_fill_manual(values = c(coral_pal[2], grey_pal[2]), 
                      name = "Disease status") +
    scale_color_manual(values = c(coral_pal[3], grey_pal[3]),
                       name = "Disease status") +
    labs(x = "Year", 
         y = "*M. vimineum* seed density") +
    fig_theme +
    theme(axis.title = element_markdown()))

(fig_dens_ginf_A <- fig_dens_gfung_A %+%
    sims_ginf_A[["long"]])

(fig_dens_gfung_S <- fig_dens_gfung_A %+%
    sims_gfung_P[["long"]] +
    aes(y = perennial_seeds) +
    labs(y = "*E. virginicus* seed density"))

(fig_dens_gfung_P <- fig_dens_gfung_A %+%
    sims_gfung_P[["long"]] +
    aes(y = perennial_adults) +
    labs(y = "*E. virginicus* adult density"))

# healthy - disease figures
(fig_diff_gfung_A <- ggplot(sims_gfung_A[["wide"]], 
                            aes(x = generation, y = annual_seeds_d -
                                  annual_seeds_h)) +
    stat_lineribbon(point_interval = mean_hdci, # multiple simulations were multi-modal 
                    .width = 0.95, alpha = 0.5,
                    fill = coral_pal[2], color = coral_pal[3]) +
    labs(x = "Year", 
         y = "Effect of disease on *M. vimineum* seed density") +
    fig_theme +
    theme(axis.title = element_markdown()))

(fig_diff_ginf_A <- fig_diff_gfung_A %+%
    sims_ginf_A[["wide"]])

(fig_diff_gfung_S <- fig_diff_gfung_A %+%
    sims_gfung_P[["wide"]] +
    aes(y = perennial_seeds_d - perennial_seeds_h) +
    labs(y = "Effect of disease on *E. virginicus* seed density"))

(fig_diff_gfung_P <- fig_diff_gfung_A %+%
    sims_gfung_P[["wide"]] +
    aes(y = perennial_adults_d - perennial_adults_h) +
    labs(y = "Effect of disease on *E. virginicus* seed density"))


#### effects ####

# alpha values
alpha_gfung_A <- params_gfung[["healthy"]][["alpha"]] %>%
  select(alphaPA) %>%
  cbind(params_g)
  mutate(disease = "h",
         .draw = params_gfung[["healthy"]][["draws"]]$.draw) %>%
  full_join(params_gfung[["disease"]][["alpha"]] %>%
              select(alphaPA) %>%
              mutate(disease = "d",
                     .draw = params_gfung[["disease"]][["draws"]]$.draw))

# save last time point
eq_gfung_A <- sims_gfung_A[["long"]] %>%
  filter(generation == gens) %>%
  left_join(alpha_gfung_A) %>%
  mutate(comp_eff = annual_seeds * alphaPA)
  
eq_gfung_A_wide <- sims_gfung_A[["wide"]] %>%
  filter(generation == gens) %>%
  left_join(alpha_gfung_A) %>%
  mutate(comp_eff_diff = alphaPA * (annual_seeds_d - 
                                      annual_seeds_h))

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens)

eq_ginf_A_wide <- sims_ginf_A[["wide"]] %>%
  filter(generation == gens)

eq_gfung_P <- sims_gfung_P[["long"]] %>%
  filter(generation == gens)

eq_gfung_P_wide <- sims_gfung_P[["wide"]] %>%
  filter(generation == gens)

# litter
eq_gfung_A %>%
  group_by(disease) %>%
  mean_hdci(litter)

eq_ginf_A %>%
  group_by(disease) %>%
  mean_hdci(litter)

eq_gfung_P %>%
  group_by(disease) %>%
  mean_hdci(litter)

# competition
eq_gfung_A %>%
  group_by(disease) %>%
  mean_hdci(comp_eff)

# figure
ggplot(eq_gfung_A, aes(x = litter, y = disease)) +
  # stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
  #           .width = c(.66, .95, 1)) + # use limits argument to cut-off tail
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  labs(x = "*M. vimineum* effect on litter", y = "Disease treatment") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ggplot(eq_gfung_A, aes(x = comp_eff, y = disease)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi,
            .width = c(.66, .95, 1)) + # use limits argument to cut-off tail
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  labs(x = "*M. vimineum* effect on litter", y = "Disease treatment") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ggplot(eq_gfung_A_wide, aes(x = litter_d - litter_h)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi,
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  labs(x = "Disease mediation of *M. vimineum* effect on litter") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ggplot(eq_gfung_A_wide, aes(x = comp_eff_diff)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdci,
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  labs(x = "Disease mediation of *M. vimineum* competitive effect") +
  fig_theme +
  theme(axis.title.x = element_markdown())


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

# cycle through parameters
for(i in 1:params_iters){
  
  # initial conditions based on resident's simulation 
  # (annual seeds, litter, perennial seeds, perennial adults)
  inits_gfung_resP_h <- c(1,
                          eq_gfung_P_wide$litter_h[i],
                          eq_gfung_P_wide$perennial_seeds_h[i],
                          eq_gfung_P_wide$perennial_adults_h[i])
  inits_gfung_resP_d <- c(1,
                          eq_gfung_P_wide$litter_d[i],
                          eq_gfung_P_wide$perennial_seeds_d[i],
                          eq_gfung_P_wide$perennial_adults_d[i])
  
  inits_gfung_resA_h <- c(eq_gfung_A_wide$annual_seeds_h[i],
                          eq_gfung_A_wide$litter_h[i],
                          0,
                          1)
  inits_gfung_resA_d <- c(eq_gfung_A_wide$annual_seeds_d[i],
                          eq_gfung_A_wide$litter_d[i],
                          0,
                          1)
  
  inits_ginf_resA_h <- c(eq_ginf_A_wide$annual_seeds_h[i],
                          eq_ginf_A_wide$litter_h[i],
                          0,
                          1)
  inits_ginf_resA_d <- c(eq_ginf_A_wide$annual_seeds_d[i],
                          eq_ginf_A_wide$litter_d[i],
                          0,
                          1)
  
  # fungicide effect, healthy
  sims_gfung_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, 
                                      init_cond = inits_gfung_resP_h,
                                      parameters = params_gfung[["healthy"]])
  sims_gfung_PA_h[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resA_h,
                                      parameters = params_gfung[["healthy"]])
  
  # fungicide effect, disease
  sims_gfung_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resP_d,
                                      parameters = params_gfung[["disease"]])
  sims_gfung_PA_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                      init_cond = inits_gfung_resA_d,
                                      parameters = params_gfung[["disease"]])
  
  # infection effect, healthy
  sims_ginf_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, 
                                     init_cond = inits_gfung_resP_h,
                                     parameters = params_ginf[["healthy"]])
  sims_ginf_PA_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, 
                                     init_cond = inits_ginf_resA_h,
                                     parameters = params_ginf[["healthy"]])
  
  # infection effect, disease
  sims_ginf_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_gfung_resP_d,
                                     parameters = params_ginf[["disease"]])
  sims_ginf_PA_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, 
                                     init_cond = inits_ginf_resA_d,
                                     parameters = params_ginf[["disease"]])
  
}

# combine healthy and disease
sims_gfung_AP <- sims_comb_fun(sims_gfung_AP_h, sims_gfung_AP_d)
sims_gfung_PA <- sims_comb_fun(sims_gfung_PA_h, sims_gfung_PA_d)
sims_ginf_AP <- sims_comb_fun(sims_ginf_AP_h, sims_ginf_AP_d)
sims_ginf_PA <- sims_comb_fun(sims_ginf_PA_h, sims_ginf_PA_d)

sims_gfung_AP[["long"]] %>% 
  ggplot(aes(x = generation, y = annual_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

sims_gfung_PA[["long"]] %>% 
  ggplot(aes(x = generation, y = perennial_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)


#### competition outcomes ####

# select last year
# are both species present?
inv_gfung_AP <- sims_gfung_AP[["long"]] %>%
  filter(generation == gens) %>%
  mutate(final_spp = case_when(annual_seeds >= 1 & 
                                 (perennial_seeds >= 1 | 
                                    perennial_adults >= 1) ~ "both",
                               annual_seeds >= 1 ~ "annual",
                               perennial_seeds >= 1 | perennial_adults >= 1 ~
                                 "perennial",
                               TRUE ~ "neither"),
         invasion = if_else(annual_seeds >= 1, "yes", "no"))

inv_gfung_PA <- sims_gfung_PA[["long"]] %>%
  filter(generation == gens) %>%
  mutate(final_spp = case_when(annual_seeds >= 1 & 
                                 (perennial_seeds >= 1 | 
                                    perennial_adults >= 1) ~ "both",
                               annual_seeds >= 1 ~ "annual",
                               perennial_seeds >= 1 | perennial_adults >= 1 ~
                                 "perennial",
                               TRUE ~ "neither"),
         invasion = if_else(perennial_seeds >= 1 | perennial_adults >= 1 , 
                            "yes", "no"))

# combine invasion outcomes
inv_gfung <- eq_gfung_AP %>%
  mutate(invader = "ann") %>%
  full_join(eq_gfung_PA %>%
              mutate(invader = "per")) %>%
  pivot_wider(values_from = -c(iteration, .draw, generation, disease),
              names_from = "invader",
              names_sep = "_") 
# check for any NA's (neither invade)

# select equilibrium values based on invasion outcome
# litter of both species with coexistence or just one without
eq_gfung <- inv_gfung %>%
  left_join(eq_gfung_A %>%
              rename(annual_seeds_A = annual_seeds,
                     litter_A = litter) %>%
              select(disease, .draw, annual_seeds_A, litter_A)) %>%
  left_join(eq_gfung_P %>%
              rename(perennial_seeds_P = perennial_seeds,
                     perennial_adults_P = perennial_adults,
                     litter_P = litter) %>%
              select(disease, .draw, perennial_seeds_P, perennial_adults_P,
                     litter_P)) %>%
  mutate(litter_A_fin = if_else(invasion_ann == "yes" & invasion_per == "yes",
                                  litter_ann, litter_A), # litter_ann should equal litter_per -- check
         litter_P_fin = if_else(invasion_ann == "yes" & invasion_per == "yes",
                                litter_per, litter_P)) # litter_per should equal litter_ann -- check

# when the invader excludes the resident, is that reflected in the alternative invasion?
eq_gfung_AP %>%
  filter(final_spp == "annual") %>%
  select(.draw, disease) %>%
  inner_join(eq_gfung_PA) %>%
  select(final_spp) # yes

eq_gfung_PA %>%
  filter(final_spp == "perennial") %>%
  select(.draw, disease) %>%
  inner_join(eq_gfung_AP) %>%
  select(final_spp)

# when they coexist, are equilibrium values equal?
eq_gfung %>%
  filter(final_spp_ann == "both" & final_spp_per == "both") %>%
  filter(annual_seeds_ann != annual_seeds_per) %>%
  select(starts_with("annual_seeds"))

eq_gfung %>%
  filter(final_spp_ann == "both" & final_spp_per == "both") %>%
  filter(perennial_seeds_ann != perennial_seeds_per) %>%
  select(starts_with("perennial_seeds"))
# no, but check again with longer simulations

# figure
ggplot(eq_gfung_AP, aes(x = final_spp, fill = disease)) +
  geom_bar(position = "dodge")

ggplot(eq_gfung_AP, aes(x = invasion, fill = disease)) +
  geom_bar(position = "dodge")

ggplot(eq_gfung_PA, aes(x = final_spp, fill = disease)) +
  geom_bar(position = "dodge")

ggplot(eq_gfung_PA, aes(x = invasion, fill = disease)) +
  geom_bar(position = "dodge")


#### responses ####

# use the simulations where A invades P as "sim outcome" so that two-species equilibrium values come from this

# from Nick:
# To summarize, we need to calculate the sensitivities just like effects, at particular points. Here is how you could do it in simulation. 1. Simulate the appropriate scenario (resident equilibrium for each species or the two species equilibrium). 
# 
# 2. Save the equilibrium densities of CA, CP, and L. Calculate the growth rates r = ln N(t+1) -  lnN(t) under those equilibrium densities. Call them r̂.
# Amy note: I think use the equilibirum density of the focal species as N(t) if they coexist, use 1 as N(t) if they don't
# 
# 3. Increase the amount of a single factor by a small amount epsilon (e.g., L* = L̂ + ϵ) while leaving the others the same. 
# 
# 4. Put those equilibrium densities and the slightly modified density into the equation N(t+1)/N(t) to calculate the growth rate. Call it r*.
# Amy note: same as above, use the real equilibrium density of focal species as N(t), or 1 if they don't coexist (invading)
# 
# 5. The sensitivity to the factor you changed is then (r* - r̂)/ϵ. 

# change in growth rate due to change in litter amount

# if species coexist, litter starts at two-species equilibrium value
# the growth rate of both species is likely zero (at equilibrium, straight lines)
# initialize with higher litter amount and everything else at equilibrium
# measure change in plant densities (divide and log-transform for growth rate)

# if species don't coexist, litter starts as equilibrium value of alternative species
# the growth rate 

#### response to 