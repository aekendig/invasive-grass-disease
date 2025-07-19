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
resp_fun <- function(iter, sim_outcome, parameters, gens = 2){
  
  # when species coexist
  if(sim_outcome$coexist == "yes"){
    
    # initial conditions (annual seeds, litter, perennial seeds, perennial adults)
    inits <- c(sim_outcome$annual_seeds,
               sim_outcome$litter,
               sim_outcome$perennial_seeds,
               sim_outcome$perennial_adults)
    
    # initial competition values (experienced by annual, perennial)
    inits_comp <- c(sim_outcome$annual_competition,
                    sim_outcome$perennial_competition)
    
    # simulate one year under baseline conditions
    sim_base <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
                            parameters = parameters)
    
    # inflate competition factors
    inits_L2 <- inits
    inits_L2[2] <- inits[2] * 1.01
    
    inits_compA2 <- inits_comp
    inits_compA2[1] <- inits_comp[1] * 1.01
    
    inits_compP2 <- inits_comp
    inits_compP2[2] <- inits_comp[2] * 1.01
    
    # simulate one year with inflation for each species
    sim_L <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_L2,
                         parameters = parameters, init_C = inits_comp)
    sim_compA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
                             parameters = parameters, init_C = inits_compA2)
    sim_compP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits,
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
    resp_LA <- (gr_LA - gr_baseA) / (inits_L2[2] - inits[2])
    resp_LP <- (gr_LP - gr_baseP) / (inits_L2[2] - inits[2])
    resp_CA <- (gr_CA - gr_baseA) / (inits_compA2[1] - inits_comp[1])
    resp_CP <- (gr_CP - gr_baseP) / (inits_compP2[2] - inits_comp[2])
    
  } else {
    
    # first letter is invader, second is resident
    # initial conditions (annual seeds, litter, perennial seeds, perennial adults)
    inits_AP <- c(1,
                  sim_outcome$litter_P,
                  sim_outcome$perennial_seeds_P,
                  sim_outcome$perennial_adults_P)
    
    inits_PA <- c(sim_outcome$annual_seeds_A,
                  sim_outcome$litter_A,
                  0,
                  1)
    
    # simulate one year under baseline conditions
    sim_base_AP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_AP,
                               parameters = parameters)
    sim_base_PA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_PA,
                               parameters = parameters)
    
    # initial competition values (experienced by annual, perennial)
    inits_compA <- c(sim_base_AP$annual_competition[1],
                     sim_base_AP$perennial_competition[1])
    inits_compP <- c(sim_base_PA$annual_competition[1],
                     sim_base_PA$perennial_competition[1])
    
    # inflate competition factors
    inits_LA2 <- inits_AP
    inits_LA2[2] <- inits_AP[2] * 1.01
    
    inits_LP2 <- inits_PA
    inits_LP2[2] <- inits_PA[2] * 1.01
    
    inits_compA2 <- inits_compA
    inits_compA2[1] <- inits_compA[1] * 1.01
    
    inits_compP2 <- inits_compP
    inits_compP2[2] <- inits_compP[2] * 1.01

    # simulate one year with inflation for each species
    sim_LA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_LA2,
                          parameters = parameters, init_C = inits_compA)
    sim_LP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_LP2,
                          parameters = parameters, init_C = inits_compP)
    sim_compA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_AP,
                             parameters = parameters, init_C = inits_compA2)
    sim_compP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_PA,
                             parameters = parameters, init_C = inits_compP2)
    
    # calculate growth rates r = lnN(t+1) -  lnN(t)
    gr_baseA <- log(sim_base_AP$annual_seeds[2] / sim_base_AP$annual_seeds[1])
    gr_LA <- log(sim_LA$annual_seeds[2] / sim_LA$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    gr_baseP <- log((sim_base_PA$perennial_seeds[2] + 
                       sim_base_PA$perennial_adults[2]) / 
                      (sim_base_PA$perennial_seeds[1] + 
                         sim_base_PA$perennial_adults[1]))
    gr_LP <- log((sim_LP$perennial_seeds[2] + sim_LP$perennial_adults[2]) / 
                   (sim_LP$perennial_seeds[1] + sim_LP$perennial_adults[1]))
    gr_CP <- log((sim_compP$perennial_seeds[2] + sim_compP$perennial_adults[2]) / 
                   (sim_compP$perennial_seeds[1] + sim_compP$perennial_adults[1]))
    
    # calculate responses
    resp_LA <- (gr_LA - gr_baseA) / (inits_LA2[2] - inits_AP[2])
    resp_LP <- (gr_LP - gr_baseP) / (inits_LP2[2] - inits_PA[2])
    resp_CA <- (gr_CA - gr_baseA) / (inits_compA2[1] - inits_compA[1])
    resp_CP <- (gr_CP - gr_baseP) / (inits_compP2[2] - inits_compP[2])
    
  }
  
  # combine values
  # indicate cases where population always goes to zero
  out <- tibble(
    .draw = parameters[["draws"]][iter, ]$.draw,
    generation = gens,
    coexist = sim_outcome$coexist,
    resp_LA = resp_LA,
    resp_LP = resp_LP,
    resp_CA = resp_CA,
    resp_CP = resp_CP)
}


#### parameters ####

# parameter combinations
params_iters <- 100

# parameters with fungicide effect on germination
params_gfung <- params_fun(iters = params_iters, gA_type = "fungicide")

# parameters with infection effect on germination
params_ginf <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws)

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


#### single-species equilibrium and effects ####

# save last time point
# make sure order matches parameters
eq_gfung_A <- sims_gfung_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1)
  
eq_gfung_A_wide <- sims_gfung_A[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = perennial_competition_d - 1,
         comp_eff_h = perennial_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1)

eq_ginf_A_wide <- sims_ginf_A[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = perennial_competition_d - 1,
         comp_eff_h = perennial_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

eq_gfung_P <- sims_gfung_P[["long"]] %>%
  filter(generation == gens)%>%
  mutate(comp_eff = annual_competition - 1)

eq_gfung_P_wide <- sims_gfung_P[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = annual_competition_d - 1,
         comp_eff_h = annual_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

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

eq_ginf_A %>%
  group_by(disease) %>%
  mean_hdci(comp_eff)

eq_gfung_P %>%
  group_by(disease) %>%
  mean_hdci(comp_eff)

# figure
ggplot(eq_gfung_A, aes(x = litter, y = disease)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi,
            .width = c(.66, .95, 1)) + # use limits argument to cut-off tail
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
  labs(x = "*M. vimineum* effect on competition", y = "Disease treatment") +
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
  
  # get draw from parameters
  draw_gfung_h <- params_gfung[["healthy"]][["draws"]]$.draw[i] 
  draw_gfung_d <- params_gfung[["disease"]][["draws"]]$.draw[i] 
  draw_ginf_h <- params_ginf[["healthy"]][["draws"]]$.draw[i] 
  draw_ginf_d <- params_ginf[["disease"]][["draws"]]$.draw[i] 
  
  # select equilibrium values
  inits_gfung_P_h <- eq_gfung_P %>% 
    filter(disease == "h" & .draw == draw_gfung_h)
  inits_gfung_P_d <- eq_gfung_P %>% 
    filter(disease == "d" & .draw == draw_gfung_h)
  
  inits_ginf_P_h <- eq_gfung_P %>% 
    filter(disease == "h" & .draw == draw_ginf_h)
  inits_ginf_P_d <- eq_gfung_P %>% 
    filter(disease == "d" & .draw == draw_ginf_h)
  
  inits_gfung_A_h <- eq_gfung_A %>% 
    filter(disease == "h" & .draw == draw_gfung_h)
  inits_gfung_A_d <- eq_gfung_A %>% 
    filter(disease == "d" & .draw == draw_gfung_h)
  
  inits_ginf_A_h <- eq_ginf_A %>% 
    filter(disease == "h" & .draw == draw_ginf_h)
  inits_ginf_A_d <- eq_ginf_A %>% 
    filter(disease == "d" & .draw == draw_ginf_h)
  
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
                                     init_cond = inits_ginf_resP_h,
                                     parameters = params_ginf[["healthy"]])
  sims_ginf_PA_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, 
                                     init_cond = inits_ginf_resA_h,
                                     parameters = params_ginf[["healthy"]])
  
  # infection effect, disease
  sims_ginf_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens,
                                     init_cond = inits_ginf_resP_d, 
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

sims_gfung_PA[["long"]] %>% 
  ggplot(aes(x = generation, y = perennial_adults,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)


#### invasion outcomes ####

# generations to assess invasion
gens_inv <- 10

# select year
# was invasion successful?
inv_gfung_AP <- sims_gfung_AP[["long"]] %>%
  filter(generation == gens_inv) %>%
  mutate(annual_invasion = if_else(annual_seeds > 1, "yes", "no"))

inv_gfung_PA <- sims_gfung_PA[["long"]] %>%
  filter(generation == gens_inv) %>%
  mutate(perennial_invasion = if_else(perennial_adults > 1 , "yes", "no"))

inv_ginf_AP <- sims_ginf_AP[["long"]] %>%
  filter(generation == gens_inv) %>%
  mutate(annual_invasion = if_else(annual_seeds > 1, "yes", "no"))

inv_ginf_PA <- sims_ginf_PA[["long"]] %>%
  filter(generation == gens_inv) %>%
  mutate(perennial_invasion = if_else(perennial_adults > 1 , "yes", "no"))

# combine invasion outcomes
# add final values from annual invasion
# add final values from each species alone
inv_gfung <- inv_gfung_AP %>%
  select(iteration, .draw, disease, annual_invasion) %>%
  full_join(inv_gfung_PA %>%
              select(iteration, .draw, disease, perennial_invasion)) %>%
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no")) %>%
  left_join(sims_gfung_AP[["long"]] %>%
              filter(generation == gens)) %>%
  left_join(eq_gfung_P %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(iteration, .draw, disease))) %>%
  left_join(eq_gfung_A %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(iteration, .draw, disease)))

inv_ginf <- inv_ginf_AP %>%
  select(iteration, .draw, disease, annual_invasion) %>%
  full_join(inv_ginf_PA %>%
              select(iteration, .draw, disease, perennial_invasion)) %>%
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no")) %>%
  left_join(sims_ginf_AP[["long"]] %>%
              filter(generation == gens)) %>%
  left_join(eq_gfung_P %>% # don't need ginf parameters for P alone
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(iteration, .draw, disease))) %>%
  left_join(eq_ginf_A %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(iteration, .draw, disease)))

# split by disease (for response simulations)
# make sure order matches parameters
inv_gfung_h <- inv_gfung %>%
  filter(disease == "h")
inv_gfung_d <- inv_gfung %>%
  filter(disease == "d")

inv_ginf_h <- inv_ginf %>%
  filter(disease == "h")
inv_ginf_d <- inv_ginf %>%
  filter(disease == "d")

# figure
ggplot(inv_gfung, aes(x = annual_invasion, fill = disease)) +
  geom_bar(position = "dodge")

ggplot(inv_gfung, aes(x = perennial_invasion, fill = disease)) +
  geom_bar(position = "dodge")

ggplot(inv_gfung, aes(x = coexist, fill = disease)) +
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

# output lits
resp_gfung_h <- list()
resp_gfung_d <- list()

resp_ginf_h <- list()
resp_ginf_d <- list()

# cycle through parameters
for(i in 1:params_iters){
  
  # get draw from parameters
  draw_gfung_h <- params_gfung[["healthy"]][["draws"]]$.draw[i] 
  draw_gfung_d <- params_gfung[["disease"]][["draws"]]$.draw[i] 
  draw_ginf_h <- params_ginf[["healthy"]][["draws"]]$.draw[i] 
  draw_ginf_d <- params_ginf[["disease"]][["draws"]]$.draw[i] 
  
  # select equilibrium values
  inv_gfung_draw_h <- inv_gfung_h %>% 
    filter(.draw == draw_gfung_h)
  inv_gfung_draw_d <- inv_gfung_d %>% 
    filter(.draw == draw_gfung_d)
  
  inv_ginf_draw_h <- inv_ginf_h %>% 
    filter(.draw == draw_ginf_h)
  inv_ginf_draw_d <- inv_ginf_d %>% 
    filter(.draw == draw_ginf_d)
  
  # calculate responses for each parameter set and disease condition
  resp_gfung_h[[i]] <- resp_fun(iter = i, sim_outcome = inv_gfung_draw_h, 
                                parameters = params_gfung[["healthy"]])
  resp_gfung_d[[i]] <- resp_fun(iter = i, sim_outcome = inv_gfung_draw_d, 
                                parameters = params_gfung[["disease"]])
  
  resp_ginf_h[[i]] <- resp_fun(iter = i, sim_outcome = inv_ginf_draw_h, 
                                parameters = params_ginf[["healthy"]])
  resp_ginf_d[[i]] <- resp_fun(iter = i, sim_outcome = inv_ginf_draw_d, 
                                parameters = params_ginf[["disease"]])
  
}

# combine healthy and disease
resp_gfung <- sims_comb_fun(resp_gfung_h, resp_gfung_d)
resp_ginf <- sims_comb_fun(resp_ginf_h, resp_ginf_d)

# values
resp_gfung[[1]] %>%
  group_by(disease, coexist) %>%
  mean_hdci(resp_LA)

resp_gfung[[1]] %>%
  group_by(disease, coexist) %>%
  mean_hdci(resp_CA)

resp_gfung[[1]] %>%
  group_by(disease, coexist) %>%
  mean_hdci(resp_LP)

resp_gfung[[1]] %>%
  group_by(disease, coexist) %>%
  mean_hdci(resp_CA)
