#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidybayes)
library(ggtext)
library(patchwork)
# scales packaged used within plotting functions below

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
params_iters <- 15000

# parameters with fungicide effect on germination
params_gfung <- params_fun(iters = params_iters, gA_type = "fungicide")

# parameters with infection effect on germination
params_ginf <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws)

# generations
gens <- 1000


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

# save
save(sims_gfung_A, file = "output/time_series_gfung_parms_annual_only_20250719.rda")
save(sims_ginf_A, file = "output/time_series_ginf_parms_annual_only_20250719.rda")
save(sims_gfung_P, file = "output/time_series_gfung_parms_perennial_only_20250719.rda")

# reload if needed
load("output/time_series_gfung_parms_annual_only_20250719.rda")
load("output/time_series_ginf_parms_annual_only_20250719.rda")
load("output/time_series_gfung_parms_perennial_only_20250719.rda")


#### single species time series figures ####

# select random rows for visualization
viz_row_sample <- sample(1:15000, 1000)

# combine datasets
sims_gfung_single <- sims_gfung_A[["long"]] %>% 
  filter(iteration %in% viz_row_sample) %>%
  select(disease, generation, .draw, annual_seeds, litter) %>%
  rename(annual_litter = litter,
         annual = annual_seeds) %>%
  full_join(sims_gfung_P[["long"]] %>% 
              filter(iteration %in% viz_row_sample) %>%
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
  filter(iteration %in% viz_row_sample) %>%
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

# save
ggsave("output/time_series_gfung_parms_single_species.png", 
       sims_gfung_single_fig, width = 6, height = 8.5)
ggsave("output/time_series_ginf_parms_single_species.png", 
       sims_ginf_single_fig, width = 6, height = 4.5)


#### evaluate single species dynamics ####

# years to get 2 adults
perennial_2_adults <- sims_gfung_P[["long"]] %>%
  filter(perennial_adults >= 2) %>%
  group_by(disease, .draw) %>%
  mutate(min_gen = min(generation)) %>%
  filter(generation == min_gen)
# some never do (nrows < 30,000)

# visualize
ggplot(perennial_2_adults, aes(x = generation)) +
  geom_histogram(binwidth = 1) +
  geom_vline(xintercept = gens_inv_P, color = "red") +
  facet_wrap(~disease, scales = "free")

# quantiles
perennial_2_adults %>%
  group_by(disease) %>%
  reframe(quantile(generation))

# cut-off for invasion growth rate calculation
gens_inv_P <- 2
gens_inv_A <- 2

# density at gen years
perennial_gen_years <- sims_gfung_P[["long"]] %>%
  filter(generation == gens_inv_P) %>%
  mutate(perennial = perennial_adults + perennial_seeds)

# visualize
ggplot(perennial_gen_years, aes(x = perennial)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~disease, scales = "free")
# too many first-years, which can impact the annual

# collapse (very few seeds)
perennial_10_seeds <- sims_gfung_P[["long"]] %>%
  filter(generation > 1 & perennial_seeds <= 10) %>%
  distinct(disease, .draw) %>%
  inner_join(sims_gfung_P[["long"]]) %>%
  mutate(perennial = perennial_adults + perennial_seeds,
         iter = paste(disease, .draw))

# dynamics
ggplot(perennial_10_seeds, aes(x = generation, y = perennial,
                               color = iter)) +
  geom_line()
# yes, these all collapse

perennial_10_seeds %>% filter(generation == 31)
perennial_10_seeds %>% filter(generation == gens)

# cut-off for establishment
gens_est_P <- 1000
gens_est_A <- 31

# density at year two related to longer-term establishment?
perennial_inv_est_years <- sims_gfung_P[["long"]] %>%
  filter(generation %in% c(gens_inv_P, gens_est_P)) %>%
  mutate(perennial = perennial_adults + perennial_seeds,
         generation = if_else(generation == gens_inv_P, "inv", "est")) %>%
  select(disease, generation, .draw, starts_with("perennial")) %>%
  pivot_wider(names_from = generation,
              values_from = starts_with("perennial")) %>%
  # mutate(established = if_else(perennial_adults_est > 1,
  #                              "established",
  #                              "not established"))
  mutate(established = if_else(perennial_est > 1, 
                               "established",
                               "not established"))

annual_inv_est_years <- sims_gfung_A[["long"]] %>%
  filter(generation %in% c(gens_inv_A, gens_est_A)) %>%
  mutate(generation = if_else(generation == gens_inv_A, "inv", "est")) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  pivot_wider(names_from = generation,
              values_from = annual_seeds,
              names_glue = "annual_seeds_{generation}") %>%
  mutate(established = if_else(annual_seeds_est > 1, 
                               "established",
                               "not established"))

# visualize perennial
ggplot(perennial_inv_est_years, aes(x = log(perennial_inv) / (gens_inv_P - 1), 
                                 y = log(perennial_est) /  (gens_est_P - 1))) +
  geom_point(aes(color = established)) +
  facet_wrap(~ disease, scales = "free")
# disease has a distinct effect on perennial density at 2 years, but not 
# adult long-term. the two variables are not related

ggplot(perennial_inv_est_years, aes(x = log(perennial_adults_inv) / (gens_inv_P - 1), 
                                   y = log(perennial_est) /  (gens_est_P - 1))) +
  geom_point(aes(color = established)) +
  facet_wrap(~ disease, scales = "free")

ggplot(perennial_inv_est_years, aes(x = log(perennial_seeds_inv) / (gens_inv_P - 1), 
                                   y = log(perennial_est) /  (gens_est_P - 1))) +
  geom_point(aes(color = established)) +
  facet_wrap(~ disease, scales = "free")
# need to separately evaluate growth after one year and establishment

# visualize annual
ggplot(annual_inv_est_years, aes(x = log(annual_seeds_inv) / (gens_inv_A - 1), 
                                    y = log(annual_seeds_est) /  (gens_est_A - 1))) +
  geom_point(aes(color = established)) +
  facet_wrap(~ disease, scales = "free")
# look completely aligned

# how well is the establishment cut-off related to 1,000 year value?
perennial_inv_est_years %>%
  filter(established == "not established") %>%
  select(disease, .draw, established) %>%
  full_join(sims_gfung_P[["long"]] %>%
               filter(generation == gens & perennial_adults <= 1)) %>%
  select(disease, .draw, established, perennial_adults) %>%
  data.frame()
# completely aligned

annual_inv_est_years %>%
  filter(established == "not established") %>%
  select(disease, .draw, established) %>%
  full_join(sims_gfung_A[["long"]] %>%
              filter(generation == gens & annual_seeds <= 1)) %>%
  select(disease, .draw, established, annual_seeds) %>%
  data.frame()
# completely aligned

# long-term dynamics
perennial_non_est <- perennial_inv_est_years %>%
  filter(established == "not established") %>%
  select(disease, .draw) %>%
  inner_join(sims_gfung_P[["long"]]) %>%
  mutate(perennial = perennial_adults + perennial_seeds,
         iter = paste(disease, .draw))

ggplot(perennial_non_est, aes(x = generation, y = perennial,
                               color = iter)) +
  geom_line()
# only a subset collapse, the others persist at reasonable densities


#### single-species equilibrium and effects ####

# save last time point
# make sure order matches parameters
eq_gfung_A <- sims_gfung_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))
  
eq_gfung_A_wide <- sims_gfung_A[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = perennial_competition_d - 1,
         comp_eff_h = perennial_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_ginf_A_wide <- sims_ginf_A[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = perennial_competition_d - 1,
         comp_eff_h = perennial_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

eq_gfung_P <- sims_gfung_P[["long"]] %>%
  filter(generation == gens)%>%
  mutate(comp_eff = annual_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_gfung_P_wide <- sims_gfung_P[["wide"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff_d = annual_competition_d - 1,
         comp_eff_h = annual_competition_h - 1,
         comp_eff_diff = comp_eff_d - comp_eff_h)

# figure
# patchwork bug requires letters built into main figures instead of added later
mv_litter_eff_dist <- ggplot(eq_gfung_A, aes(x = litter)) +
  stat_slab(aes(color = disease_name, fill = disease_name), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 15000)) +
  scale_fill_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_color_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "*M. vimineum* effect on litter", title = "A") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot")

mv_litter_eff_dist2 <- mv_litter_eff_dist %+%
  eq_ginf_A

ev_litter_eff_dist <- mv_litter_eff_dist %+%
  eq_gfung_P +
  coord_cartesian(xlim = c(0, 15000)) +
  labs(x = "*E. virginicus* effect on litter", title = "B")

mv_comp_eff_dist <- mv_litter_eff_dist %+%
  aes(x = comp_eff) +
  coord_cartesian(xlim = range(eq_gfung_A$comp_eff)) +
  labs(x = "*M. vimineum* effect on competition")

mv_comp_eff_dist2 <- mv_comp_eff_dist %+%
  eq_ginf_A

ev_comp_eff_dist <- ev_litter_eff_dist %+%
  aes(x = comp_eff) +
  coord_cartesian(xlim = c(0, 600)) +
  labs(x = "*E. virginicus* effect on competition")

mv_litter_eff_pt <- ggplot(eq_gfung_A, aes(x = disease_name, y = litter)) +
  stat_pointinterval(aes(color = disease_name), point_size = 2,
                     point_interval = median_hdci, .width = c(.66, .95)) +
  scale_color_manual(values = c(grey_pal[2], coral_pal[2]), guide = "none") +
  scale_y_continuous(labels = scales::comma, n.breaks = 3) +
  fig_theme +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

mv_litter_eff_pt2 <- mv_litter_eff_pt %+%
  eq_ginf_A

ev_litter_eff_pt <- mv_litter_eff_pt %+%
  eq_gfung_P

mv_comp_eff_pt <- mv_litter_eff_pt %+%
  aes(y = comp_eff)

mv_comp_eff_pt2 <- mv_litter_eff_pt2 %+%
  aes(y = comp_eff)

ev_comp_eff_pt <- ev_litter_eff_pt %+%
  aes(y = comp_eff)

# combine
eff_inset_coord <- c(0.25, 0.5, 0.75, 0.99)

mv_litter_eff_fig <- mv_litter_eff_dist + 
  inset_element(mv_litter_eff_pt, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])
mv_litter_eff_fig2 <- mv_litter_eff_dist2 + 
  inset_element(mv_litter_eff_pt2, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])
ev_litter_eff_fig <- ev_litter_eff_dist + 
  inset_element(ev_litter_eff_pt, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])
mv_comp_eff_fig <- mv_comp_eff_dist + 
  inset_element(mv_comp_eff_pt, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])
mv_comp_eff_fig2 <- mv_comp_eff_dist2 + 
  inset_element(mv_comp_eff_pt2, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])
ev_comp_eff_fig <- ev_comp_eff_dist + 
  inset_element(ev_comp_eff_pt, eff_inset_coord[1], eff_inset_coord[2],
                eff_inset_coord[3], eff_inset_coord[4])


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

# save
save(sims_gfung_AP, file = "output/time_series_gfung_parms_annual_invasion_20250719.rda")
save(sims_gfung_PA, file = "output/time_series_gfung_parms_perennial_invasion_20250719.rda")
save(sims_ginf_AP, file = "output/time_series_ginf_parms_annual_invasion_20250719.rda")
save(sims_ginf_PA, file = "output/time_series_ginf_parms_perennial_invasion_20250719.rda")

# reload if needed
load("output/time_series_gfung_parms_annual_invasion_20250719.rda")
load("output/time_series_gfung_parms_perennial_invasion_20250719.rda")
load("output/time_series_ginf_parms_annual_invasion_20250719.rda")
load("output/time_series_ginf_parms_perennial_invasion_20250719.rda")


#### invasion time series figures ####

# combine datasets
sims_gfung_inv <- sims_gfung_AP[["long"]] %>% 
  filter(iteration %in% viz_row_sample) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  rename(annual = annual_seeds) %>%
  full_join(sims_gfung_PA[["long"]] %>% 
              filter(iteration %in% viz_row_sample) %>%
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

sims_ginf_inv <- sims_ginf_AP[["long"]] %>% 
  filter(iteration %in% viz_row_sample) %>%
  select(disease, generation, .draw, annual_seeds) %>%
  rename(annual = annual_seeds) %>%
  full_join(sims_ginf_PA[["long"]] %>% 
              filter(iteration %in% viz_row_sample) %>%
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

# save
ggsave("output/time_series_gfung_parms_invasions.png", 
       sims_gfung_inv_fig, width = 6, height = 4.5)
ggsave("output/time_series_ginf_parms_invasions.png", 
       sims_ginf_inv_fig, width = 6, height = 4.5)


#### invasion outcomes ####

# generations to assess invasion
gens_inv_A <- 2
gens_inv_P <- 2

# perennial densities at end
inv_gfung_PA_fin <- sims_gfung_PA[["long"]] %>%
  filter(generation == gens) %>%
  mutate(perennial_fin = perennial_adults + perennial_seeds) %>%
  distinct(disease, .draw, perennial_fin)

inv_ginf_PA_fin <- sims_ginf_PA[["long"]] %>%
  filter(generation == gens) %>%
  mutate(perennial_fin = perennial_adults + perennial_seeds) %>%
  distinct(disease, .draw, perennial_fin)

# select year
# was invasion successful?
inv_gfung_AP <- sims_gfung_AP[["long"]] %>%
  filter(generation == gens_inv_A) %>%
  mutate(annual_gr = log(annual_seeds)/(gens_inv_A - 1),
         annual_invasion = if_else(annual_gr > 0, "yes", "no"))

inv_gfung_PA <- sims_gfung_PA[["long"]] %>%
  filter(generation == gens_inv_P) %>%
  left_join(inv_gfung_PA_fin) %>%
  mutate(perennial = perennial_adults + perennial_seeds,
         perennial_gr = log(perennial)/(gens_inv_P - 1),
         perennial_invasion = if_else(perennial_gr > 0 & perennial_fin > 1 , 
                                      "yes", "no"))

inv_ginf_AP <- sims_ginf_AP[["long"]] %>%
  filter(generation == gens_inv_A) %>%
  mutate(annual_gr = log(annual_seeds)/(gens_inv_A - 1),
         annual_invasion = if_else(annual_gr > 0, "yes", "no"))

inv_ginf_PA <- sims_ginf_PA[["long"]] %>%
  filter(generation == gens_inv_P) %>%
  left_join(inv_ginf_PA_fin) %>%
  mutate(perennial = perennial_adults + perennial_seeds,
         perennial_gr = log(perennial)/(gens_inv_P - 1),
         perennial_invasion = if_else(perennial_gr > 0 & perennial_fin > 1 , 
                                      "yes", "no"))

# combine invasion outcomes
# add final values from annual invasion
# add final values from each species alone
inv_gfung <- inv_gfung_AP %>%
  select(iteration, .draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_gfung_PA %>%
              select(iteration, .draw, disease, perennial_invasion, 
                     perennial_gr)) %>%
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no"),
         outcome = case_when(coexist == "yes" ~ "coexist",
                             annual_invasion == "yes" & 
                               perennial_invasion == "no" ~ "annual only",
                             perennial_invasion == "yes" &
                               annual_invasion == "no" ~ "perennial only",
                             TRUE ~ "priority effect") %>%
           fct_relevel("coexist")) %>%
  left_join(sims_gfung_AP[["long"]] %>%
              filter(generation == gens)) %>%
  left_join(eq_gfung_P %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(iteration, .draw, disease, disease_name))) %>%
  left_join(eq_gfung_A %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(iteration, .draw, disease, disease_name)))

inv_ginf <- inv_ginf_AP %>%
  select(iteration, .draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_ginf_PA %>%
              select(iteration, .draw, disease, perennial_invasion, 
                     perennial_gr)) %>%
  mutate(coexist = if_else(annual_invasion == "yes" & 
                             perennial_invasion == "yes",
                           "yes", "no"),
         outcome = case_when(coexist == "yes" ~ "coexist",
                             annual_invasion == "yes" & 
                               perennial_invasion == "no" ~ "annual only",
                             perennial_invasion == "yes" &
                               annual_invasion == "no" ~ "perennial only",
                             TRUE ~ "priority effect") %>%
           fct_relevel("coexist")) %>%
  left_join(sims_ginf_AP[["long"]] %>%
              filter(generation == gens)) %>%
  left_join(eq_gfung_P %>% # don't need ginf parameters for P alone
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(iteration, .draw, disease, disease_name))) %>%
  left_join(eq_ginf_A %>%
              select(-generation) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(iteration, .draw, disease, disease_name)))

# priority effect simulations
inv_gfung %>%
  filter(outcome == "priority effect") %>%
  select(.draw, disease, annual_seeds_A, perennial_seeds_P, 
         perennial_adults_P) %>%
  head(n = 50) %>%
  data.frame()
# can establish

inv_ginf %>%
  filter(outcome == "priority effect") %>%
  select(.draw, disease, annual_seeds_A, perennial_seeds_P, 
         perennial_adults_P) %>%
  head(n = 50) %>%
  data.frame()

# split by disease (for response simulations)
inv_gfung_h <- inv_gfung %>%
  filter(disease == "h")
inv_gfung_d <- inv_gfung %>%
  filter(disease == "d")

inv_ginf_h <- inv_ginf %>%
  filter(disease == "h")
inv_ginf_d <- inv_ginf %>%
  filter(disease == "d")

# summarize for figure
inv_gfung_sum <- inv_gfung %>%
  count(disease_name, outcome) %>%
  group_by(disease_name) %>%
  mutate(prop = n / (sum(n))) %>%
  ungroup()

inv_ginf_sum <- inv_ginf %>%
  count(disease_name, outcome) %>%
  group_by(disease_name) %>%
  mutate(prop = n / (sum(n))) %>%
  ungroup()

# summary figure
coex_gfung_sum_fig <- ggplot(inv_gfung_sum, aes(x = disease_name, y = prop, 
                          fill = outcome)) +
  geom_col(position = "dodge", show.legend = F) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = col_pal4, name = "Invasion outcome") +
  labs(y = "Parameter draws", title = "A") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot")

coex_ginf_sum_fig <- coex_gfung_sum_fig %+%
  inv_ginf_sum

# growth rate figures
coex_gfung_grd_fig <- inv_gfung %>%
  filter(disease == "d") %>%
  ggplot(aes(x = annual_gr, y = perennial_gr, color = outcome, 
             shape = outcome)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 0.5) +
  facet_wrap(~disease_name) +
  scale_color_manual(values = col_pal4, name = "Invasion outcome") +
  scale_shape_manual(values = shape_pal4, name = "Invasion outcome") +
  labs(y = "*E. virginicus* GRWR", x = "*M. vimineum* GRWR",
       title = "B") +
  fig_theme +
  theme(axis.title = element_markdown(),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot",
        strip.placement = "inside") +
  guides(color = guide_legend(override.aes = list(size = 3)))

coex_gfung_grh_fig <- coex_gfung_grd_fig %+%
  filter(inv_gfung, disease == "h") +
  labs(title = "C")

coex_ginf_grd_fig <- coex_gfung_grd_fig %+%
  filter(inv_ginf, disease == "d")

coex_ginf_grh_fig <- coex_ginf_grd_fig %+%
  filter(inv_ginf, disease == "h") +
  labs(title = "C")

# combine
coex_gfung_fig <- coex_gfung_sum_fig / (coex_gfung_grd_fig + coex_gfung_grh_fig) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

coex_ginf_fig <- coex_ginf_sum_fig / (coex_ginf_grd_fig + coex_ginf_grh_fig) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# save
ggsave("output/coexistence_outcomes_gfung_parms.png", coex_gfung_fig,
       width = 6.5, height = 6.5)
ggsave("output/coexistence_outcomes_ginf_parms.png", coex_ginf_fig,
       width = 6.5, height = 6.5)


#### responses ####

# output lists
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

# save
save(resp_gfung, file = "output/responses_gfung_parms_20250719.rda")
save(resp_ginf, file = "output/responses_ginf_parms_20250719.rda")

# reload if needed
load("output/responses_gfung_parms_20250719.rda")
load("output/responses_ginf_parms_20250719.rda")

# format for figures
resp_gfung_long <- resp_gfung[[1]] %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

resp_ginf_long <- resp_ginf[[1]] %>%
  mutate(disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))
  
# distribution figures
mv_litter_resp_dist <- ggplot(resp_gfung_long, aes(x = resp_LA)) +
  stat_slab(aes(color = disease_name, fill = disease_name), alpha = 0.5) +
  scale_fill_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_color_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "*M. vimineum* response to litter", title = "C") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot")

ev_litter_resp_dist <- mv_litter_resp_dist %+%
  aes(x = resp_LP) +
  labs(x = "*E. virginicus* response to litter", title = "D")
# lots of errors, but seems like most accurate representation from all the geoms
# tried histogram, freqpoly, and dots

mv_comp_resp_dist <- mv_litter_resp_dist %+%
  aes(x = resp_CA) +
  labs(x = "*M. vimineum* response to competition")

ev_comp_resp_dist <- ev_litter_resp_dist %+%
  aes(x = resp_CP) +
  labs(x = "*E. virginicus* response to competition")


mv_litter_resp_dist2 <- mv_litter_resp_dist %+%
  resp_ginf_long

ev_litter_resp_dist2 <- ev_litter_resp_dist %+%
  resp_ginf_long

mv_comp_resp_dist2 <- mv_comp_resp_dist %+%
  resp_ginf_long

ev_comp_resp_dist2 <- ev_comp_resp_dist %+%
  resp_ginf_long

# point figures
mv_litter_resp_pt <- ggplot(resp_gfung_long, 
                            aes(x = disease_name, y = resp_LA)) +
  stat_pointinterval(aes(color = disease_name), point_size = 2,
                     point_interval = median_hdci, .width = c(.66, .95)) +
  scale_color_manual(values = c(grey_pal[2], coral_pal[2]), guide = "none") +
  scale_y_continuous(n.breaks = 3) +
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

#### combine litter figures, combine response figures ####

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

ggsave("output/litter_effects_responses_gfung_parms.png", eff_resp_litter_fig,
       width = 6.5, height = 6.5)
ggsave("output/competition_effects_responses_gfung_parms.png", eff_resp_comp_fig,
       width = 6.5, height = 6.5)
ggsave("output/litter_effects_responses_ginf_parms.png", eff_resp_litter_fig2,
       width = 6.5, height = 6.5)
ggsave("output/competition_effects_responses_ginf_parms.png", eff_resp_comp_fig2,
       width = 6.5, height = 6.5)
