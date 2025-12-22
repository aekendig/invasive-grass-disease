#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidybayes)
library(ggtext)
library(patchwork)
library(janitor)
# scales packaged used within plotting functions below

# import parameters (loads tidyverse)
source("code/dynamical-model/parameters.R")

# import simulation function
source("code/dynamical-model/APL_SIm_Tree.R")

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
                            parameters = parameters,
                            return_M = T)
    
    # inflate competition factors
    inits_L2 <- inits
    inits_L2[2] <- inits[2] * 1.01
    
    inits_compA2 <- inits_comp
    inits_compA2[1] <- inits_comp[1] * 1.01
    
    inits_compP2 <- inits_comp
    inits_compP2[2] <- inits_comp[2] * 1.01
    
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

    # calculate growth rates r = (ln[N(t+1)] -  ln[N(t)])/delta_t
    gr_baseA <- log(out_base$annual_seeds[2] / out_base$annual_seeds[1])
    gr_LA <- log(out_L$annual_seeds[2] / out_L$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    # calculate growth rates dominant eigenvalue
    gr_baseP <- eigen(sim_base[[2]], only.values = T)$values[1]
    gr_LP <- eigen(sim_L[[2]], only.values = T)$values[1]
    gr_CP <- eigen(sim_compP[[2]], only.values = T)$values[1]
    
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
                               parameters = parameters, return_M = T)
    
    # initial competition values (experienced by annual, perennial)
    inits_compA <- c(sim_base_AP$annual_competition[1],
                     sim_base_AP$perennial_competition[1])
    inits_compP <- c(sim_base_PA[[1]]$annual_competition[1],
                     sim_base_PA[[1]]$perennial_competition[1])
    
    # inflate litter and competition factors for annual
    inits_LA2 <- inits_AP
    inits_LA2[2] <- inits_AP[2] * 1.01
    
    inits_compA2 <- inits_compA
    inits_compA2[1] <- inits_compA[1] * 1.01
    
    # inflate factors for perennial (if no annual, add 0.01)
    inits_LP2 <- inits_PA
    if(inits_PA[2] > 0){
      
      inits_LP2[2] <- inits_PA[2] * 1.01
    
    } else {
      
      inits_LP2[2] <- 0.01
        
      }
    
    inits_compP2 <- inits_compP
    if(inits_compP[2] > 0){
      
      inits_compP2[2] <- inits_compP[2] * 1.01
      
    } else {
      
      inits_compP2[2] <- 0.01
      
    }

    # simulate one year with inflation for each species
    sim_LA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_LA2,
                          parameters = parameters, init_C = inits_compA)
    sim_LP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_LP2,
                          parameters = parameters, init_C = inits_compP,
                          return_M = T)
    sim_compA <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_AP,
                             parameters = parameters, init_C = inits_compA2)
    sim_compP <- dyn_mod_fun(iter = iter, gen = gens, init_cond = inits_PA,
                             parameters = parameters, init_C = inits_compP2,
                             return_M = T)
    
    # calculate growth rates r = lnN(t+1) -  lnN(t)
    gr_baseA <- log(sim_base_AP$annual_seeds[2] / sim_base_AP$annual_seeds[1])
    gr_LA <- log(sim_LA$annual_seeds[2] / sim_LA$annual_seeds[1])
    gr_CA <- log(sim_compA$annual_seeds[2] / sim_compA$annual_seeds[1])
    
    # calculate growth rates dominant eigenvalue
    gr_baseP <- eigen(sim_base_PA[[2]], only.values = T)$values[1]
    gr_LP <- eigen(sim_LP[[2]], only.values = T)$values[1]
    gr_CP <- eigen(sim_compP[[2]], only.values = T)$values[1]
    
    # calculate responses
    resp_LA <- (gr_LA - gr_baseA) / (inits_LA2[2] - inits_AP[2])
    resp_LP <- (gr_LP - gr_baseP) / (inits_LP2[2] - inits_PA[2])
    resp_CA <- (gr_CA - gr_baseA) / (inits_compA2[1] - inits_compA[1])
    resp_CP <- (gr_CP - gr_baseP) / (inits_compP2[2] - inits_compP[2])
    
  }
  
  # combine values
  # indicate cases where population always goes to zero
  out <- tibble(.draw = parameters[["draws"]][iter, ]$.draw,
                generation = gens,
                coexist = sim_outcome$coexist,
                resp_LA = resp_LA,
                resp_LP = resp_LP,
                resp_CA = resp_CA,
                resp_CP = resp_CP)
  
  return(out)
}


#### parameters ####

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

# generations
gens <- 1000


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

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_infP_A <- sims_infP_A[["long"]] %>%
  filter(generation == gens) %>%
  mutate(comp_eff = perennial_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

eq_gfung_P <- sims_gfung_P[["long"]] %>%
  filter(generation == gens)%>%
  mutate(comp_eff = annual_competition - 1,
         disease_name = fct_recode(disease,
                                   "Ambient disease" = "d",
                                   "Disease suppressed" = "h") %>%
           fct_relevel("Ambient disease"))

# figure
# patchwork bug requires letters built into main figures instead of added later
mv_litter_eff_dist <- ggplot(eq_gfung_A, aes(x = litter)) +
  stat_slab(aes(color = disease_name, fill = disease_name), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 15000)) +
  scale_fill_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_color_manual(values = c(grey_pal[2], coral_pal[2])) +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "*M. vimineum* effect on litter-mediated competition", title = "A") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot")

mv_litter_eff_dist2 <- mv_litter_eff_dist %+%
  eq_ginf_A

mv_litter_eff_dist3 <- mv_litter_eff_dist %+%
  eq_infP_A +
  coord_cartesian(xlim = range(eq_infP_A$litter))

ev_litter_eff_dist <- mv_litter_eff_dist %+%
  eq_gfung_P +
  coord_cartesian(xlim = c(0, 15000)) +
  labs(x = "*E. virginicus* effect on litter-mediated competition", title = "B")

mv_comp_eff_dist <- mv_litter_eff_dist %+%
  aes(x = comp_eff) +
  coord_cartesian(xlim = range(eq_gfung_A$comp_eff)) +
  labs(x = "*M. vimineum* effect on density-mediated competition")

mv_comp_eff_dist2 <- mv_comp_eff_dist %+%
  eq_ginf_A +
  coord_cartesian(xlim = range(eq_ginf_A$comp_eff))

mv_comp_eff_dist3 <- mv_comp_eff_dist %+%
  eq_infP_A +
  coord_cartesian(xlim = range(eq_infP_A$comp_eff))

ev_comp_eff_dist <- ev_litter_eff_dist %+%
  aes(x = comp_eff) +
  coord_cartesian(xlim = c(0, 600)) +
  labs(x = "*E. virginicus* effect on density-mediated competition")

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

mv_litter_eff_pt3 <- mv_litter_eff_pt %+%
  eq_infP_A

ev_litter_eff_pt <- mv_litter_eff_pt %+%
  eq_gfung_P

mv_comp_eff_pt <- mv_litter_eff_pt %+%
  aes(y = comp_eff)

mv_comp_eff_pt2 <- mv_litter_eff_pt2 %+%
  aes(y = comp_eff)

mv_comp_eff_pt3 <- mv_litter_eff_pt3 %+%
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
mv_litter_eff_fig3 <- mv_litter_eff_dist3 + 
  inset_element(mv_litter_eff_pt3, eff_inset_coord[1], eff_inset_coord[2],
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
mv_comp_eff_fig3 <- mv_comp_eff_dist3 + 
  inset_element(mv_comp_eff_pt3, eff_inset_coord[1], eff_inset_coord[2],
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
  left_join(sims_gfung_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_gfung_P %>%
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(.draw, disease, disease_name))) %>%
  left_join(eq_gfung_A %>%
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(.draw, disease, disease_name)))

inv_ginf <- inv_ginf_AP %>%
  select(.draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_ginf_PA %>%
              select(.draw, disease, perennial_invasion, perennial_gr)) %>%
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
  left_join(sims_ginf_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_gfung_P %>% # don't need ginf parameters for P alone
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(.draw, disease, disease_name))) %>%
  left_join(eq_ginf_A %>%
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(.draw, disease, disease_name)))

inv_infP <- inv_infP_AP %>%
  select(.draw, disease, annual_invasion, annual_gr) %>%
  full_join(inv_infP_PA %>%
              select(.draw, disease, perennial_invasion, perennial_gr)) %>%
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
  left_join(sims_infP_AP_long %>%
              filter(generation == gens) %>%
              select(-c(generation, iteration))) %>%
  left_join(eq_gfung_P %>% # don't need infP parameters for P alone
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_P"),
                          .cols = -c(.draw, disease, disease_name))) %>%
  left_join(eq_infP_A %>%
              select(-c(generation, iteration)) %>%
              rename_with(.fn = ~paste0(.x, "_A"),
                          .cols = -c(.draw, disease, disease_name)))

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

# split by disease (for response simulations)
inv_gfung_h <- inv_gfung %>%
  filter(disease == "h")
inv_gfung_d <- inv_gfung %>%
  filter(disease == "d")

inv_ginf_h <- inv_ginf %>%
  filter(disease == "h")
inv_ginf_d <- inv_ginf %>%
  filter(disease == "d")

inv_infP_h <- inv_infP %>%
  filter(disease == "h")
inv_infP_d <- inv_infP %>%
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

inv_infP_sum <- inv_infP %>%
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

coex_infP_sum_fig <- coex_gfung_sum_fig %+%
  inv_infP_sum

# change outcome labels for figure
inv_gfung2 <- inv_gfung %>%
  mutate(outcome = str_replace(outcome, "annual", "*M. vimineum*") %>%
           str_replace("perennial", "*E. virginicus*") %>%
           fct_relevel("coexist", "*M. vimineum* only"))

inv_ginf2 <- inv_ginf %>%
  mutate(outcome = str_replace(outcome, "annual", "*M. vimineum*") %>%
           str_replace("perennial", "*E. virginicus*") %>%
           fct_relevel("coexist", "*M. vimineum* only"))

inv_infP2 <- inv_infP %>%
  mutate(outcome = str_replace(outcome, "annual", "*M. vimineum*") %>%
           str_replace("perennial", "*E. virginicus*") %>%
           fct_relevel("coexist", "*M. vimineum* only"))

# growth rate figures
coex_gfung_grd_fig <- inv_gfung2 %>%
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
        strip.placement = "inside",
        legend.text = element_markdown()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

coex_gfung_grh_fig <- coex_gfung_grd_fig %+%
  filter(inv_gfung2, disease == "h") +
  labs(title = "C")

coex_ginf_grd_fig <- coex_gfung_grd_fig %+%
  filter(inv_ginf2, disease == "d")

coex_ginf_grh_fig <- coex_ginf_grd_fig %+%
  filter(inv_ginf2, disease == "h") +
  labs(title = "C")

coex_infP_grd_fig <- coex_gfung_grd_fig %+%
  filter(inv_infP2, disease == "d")

coex_infP_grh_fig <- coex_infP_grd_fig %+%
  filter(inv_infP2, disease == "h") +
  labs(title = "C")

# combine
coex_gfung_fig <- coex_gfung_sum_fig / (coex_gfung_grd_fig + coex_gfung_grh_fig) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

coex_ginf_fig <- coex_ginf_sum_fig / (coex_ginf_grd_fig + coex_ginf_grh_fig) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

coex_infP_fig <- coex_infP_sum_fig / (coex_infP_grd_fig + coex_infP_grh_fig) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# save
ggsave("output/coexistence_outcomes_gfung_parms.png", coex_gfung_fig,
       width = 6.5, height = 6.5)
ggsave("output/coexistence_outcomes_ginf_parms.png", coex_ginf_fig,
       width = 6.5, height = 6.5)
ggsave("output/coexistence_outcomes_infP_parms.png", coex_infP_fig,
       width = 6.5, height = 6.5)

# table
inv_gfung_sum %>%
  mutate(parameter_set = 1) %>%
  full_join(inv_ginf_sum %>%
              mutate(parameter_set = 2)) %>%
  full_join(inv_infP_sum %>%
              mutate(parameter_set = 3)) %>%
  mutate(perc = janitor::round_half_up(prop * 100)) %>%
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
  inv_gfung_draw_h <- inv_gfung_h %>%
    filter(.draw == draw_gfung_h)
  inv_gfung_draw_d <- inv_gfung_d %>%
    filter(.draw == draw_gfung_d)

  inv_ginf_draw_h <- inv_ginf_h %>%
    filter(.draw == draw_ginf_h)
  inv_ginf_draw_d <- inv_ginf_d %>%
    filter(.draw == draw_ginf_d)
  
  inv_infP_draw_h <- inv_infP_h %>% 
    filter(.draw == draw_infP_h)
  inv_infP_draw_d <- inv_infP_d %>% 
    filter(.draw == draw_infP_d)
  
  # calculate responses for each parameter set and disease condition
  resp_gfung_h[[i]] <- resp_fun(iter = i, sim_outcome = inv_gfung_draw_h,
                                parameters = params_gfung[["healthy"]])
  resp_gfung_d[[i]] <- resp_fun(iter = i, sim_outcome = inv_gfung_draw_d,
                                parameters = params_gfung[["disease"]])

  resp_ginf_h[[i]] <- resp_fun(iter = i, sim_outcome = inv_ginf_draw_h,
                                parameters = params_ginf[["healthy"]])
  resp_ginf_d[[i]] <- resp_fun(iter = i, sim_outcome = inv_ginf_draw_d,
                                parameters = params_ginf[["disease"]])
  
  resp_infP_h[[i]] <- resp_fun(iter = i, sim_outcome = inv_infP_draw_h, 
                               parameters = params_infP[["healthy"]])
  resp_infP_d[[i]] <- resp_fun(iter = i, sim_outcome = inv_infP_draw_d, 
                               parameters = params_infP[["disease"]])
  
}

# combine healthy and disease
resp_gfung <- sims_comb_fun(resp_gfung_h, resp_gfung_d)
resp_ginf <- sims_comb_fun(resp_ginf_h, resp_ginf_d)
resp_infP <- sims_comb_fun(resp_infP_h, resp_infP_d)

# save
save(resp_gfung, file = "output/responses_gfung_parms_20250719.rda")
save(resp_ginf, file = "output/responses_ginf_parms_20250719.rda")
save(resp_infP, file = "output/responses_infP_parms_20250719.rda")

# reload if needed
load("output/responses_gfung_parms_20250719.rda")
load("output/responses_ginf_parms_20250719.rda")
load("output/responses_infP_parms_20250719.rda")

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

resp_infP_long <- resp_infP[[1]] %>%
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
  labs(x = "*M. vimineum* response to litter-mediated competition", 
       title = "C") +
  fig_theme +
  theme(axis.title.x = element_markdown(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot")

ev_litter_resp_dist <- mv_litter_resp_dist %+%
  aes(x = resp_LP) +
  labs(x = "*E. virginicus* response to litter-mediated competition", 
       title = "D")

mv_comp_resp_dist <- mv_litter_resp_dist %+%
  aes(x = resp_CA) +
  labs(x = "*M. vimineum* response to density-mediated competition") +
  theme(axis.title.x = element_markdown(hjust = 0.75))

ev_comp_resp_dist <- ev_litter_resp_dist %+%
  aes(x = resp_CP) +
  labs(x = "*E. virginicus* response to density-mediated competition")

mv_litter_resp_dist2 <- mv_litter_resp_dist %+%
  resp_ginf_long

ev_litter_resp_dist2 <- ev_litter_resp_dist %+%
  resp_ginf_long

mv_comp_resp_dist2 <- mv_comp_resp_dist %+%
  resp_ginf_long

ev_comp_resp_dist2 <- ev_comp_resp_dist %+%
  resp_ginf_long

mv_litter_resp_dist3 <- mv_litter_resp_dist %+%
  resp_infP_long

ev_litter_resp_dist3 <- ev_litter_resp_dist %+%
  resp_infP_long

mv_comp_resp_dist3 <- mv_comp_resp_dist %+%
  resp_infP_long

ev_comp_resp_dist3 <- ev_comp_resp_dist %+%
  resp_infP_long

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

#### combine litter and density figures ####

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
