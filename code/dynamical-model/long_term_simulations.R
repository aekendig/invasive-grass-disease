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
    mutate(disease = "suppressed") %>%
    full_join(bind_rows(sims_d, .id = "iteration") %>%
                mutate(disease = "ambient")) %>%
    mutate(iteration = as.numeric(iteration),
           disease = fct_relevel(disease, "suppressed"))
  
  # make wide by disease
  out_wide <- out_long %>%
    pivot_wider(names_from = disease, 
                values_from = c(annual_seeds, litter, 
                                perennial_seeds, perennial_adults))
  
  # output both formats
  return(list(long = out_long, wide = out_wide))

}


#### parameters ####

# parameter combinations
params_iters <- 100

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
# initsAP <- c(1, 0, 0, 1)

# output lits
sims_gfung_A_h <- list()
sims_gfung_P_h <- list()
# sims_gfung_AP_h <- list()

sims_gfung_A_d <- list()
sims_gfung_P_d <- list()
# sims_gfung_AP_d <- list()

sims_ginf_A_h <- list()
# sims_ginf_AP_h <- list()

sims_ginf_A_d <- list()
# sims_ginf_AP_d <- list()

# cycle through parameters
for(i in 1:params_iters){
  
  # fungicide effect, healthy
  sims_gfung_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["healthy"]])
  sims_gfung_P_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["healthy"]])
  # sims_gfung_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
  #                                     parameters = params_gfung[["healthy"]])
  
  # fungicide effect, disease
  sims_gfung_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["disease"]])
  sims_gfung_P_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["disease"]])
  # sims_gfung_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
  #                                     parameters = params_gfung[["disease"]])
  
  # infection effect, healthy
  sims_ginf_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_ginf[["healthy"]])
  # sims_ginf_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
  #                                    parameters = params_ginf[["healthy"]])
  
  # infection effect, disease
  sims_ginf_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                    parameters = params_ginf[["disease"]])
  # sims_ginf_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
  #                                    parameters = params_ginf[["disease"]])
  
}

# combine healthy and disease
sims_gfung_A <- sims_comb_fun(sims_gfung_A_h, sims_gfung_A_d)
sims_ginf_A <- sims_comb_fun(sims_ginf_A_h, sims_ginf_A_d)
sims_gfung_P <- sims_comb_fun(sims_gfung_P_h, sims_gfung_P_d)
# sims_gfung_AP <- sims_comb_fun(sims_gfung_AP_h, sims_gfung_AP_d)
# sims_ginf_AP <- sims_comb_fun(sims_ginf_AP_h, sims_ginf_AP_d)

# have the population stabilized? 
# (memory too limited to show all -> subsample)
sims_gfung_A[["long"]] %>% 
  filter(iteration %in% sample(1:15000, 1000)) %>%
  ggplot(aes(x = generation, y = annual_seeds,
           group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

# ggplot(sims_gfung_AP[["long"]], aes(x = generation, y = annual_seeds,
#                                    group = iteration, color = iteration)) +
#   geom_line() +
#   facet_wrap(~ disease)

sims_ginf_A[["long"]] %>%
  filter(iteration %in% sample(1:15000, 1000)) %>%
  ggplot(aes(x = generation, y = annual_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

# ggplot(sims_ginf_AP[["long"]], aes(x = generation, y = annual_seeds,
#                                    group = iteration, color = iteration)) +
#   geom_line() +
#   facet_wrap(~ disease) +
#   coord_cartesian(ylim = c(0, 70000))

sims_gfung_P[["long"]] %>%
  filter(iteration %in% sample(1:15000, 1000)) %>% 
  ggplot(aes(x = generation, y = perennial_seeds,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

sims_gfung_P[["long"]] %>%
  filter(iteration %in% sample(1:15000, 1000)) %>% 
  ggplot(aes(x = generation, y = perennial_adults,
             group = iteration, color = iteration)) +
  geom_line() +
  facet_wrap(~ disease)

# ggplot(sims_gfung_AP[["long"]], aes(x = generation, y = perennial_seeds,
#                                    group = iteration, color = iteration)) +
#   geom_line() +
#   facet_wrap(~ disease) # no
# 
# ggplot(sims_gfung_AP[["long"]], aes(x = generation, y = perennial_adults,
#                                    group = iteration, color = iteration)) +
#   geom_line() +
#   facet_wrap(~ disease)

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

# (fig_dens_gfung_AP_A <- fig_dens_gfung_A %+%
#     sims_gfung_AP[["long"]])
# 
# (fig_dens_ginf_AP_A <- fig_dens_gfung_A %+%
#     sims_ginf_AP[["long"]])
# 
# (fig_dens_gfung_AP_S <- fig_dens_gfung_S %+%
#     sims_gfung_AP[["long"]])
# 
# (fig_dens_gfung_AP_P <- fig_dens_gfung_P %+%
#     sims_gfung_AP[["long"]])
# 
# (fig_dens_ginf_AP_S <- fig_dens_gfung_S %+%
#     sims_ginf_AP[["long"]])
# 
# (fig_dens_ginf_AP_P <- fig_dens_gfung_P %+%
#     sims_ginf_AP[["long"]])

# healthy - disease figures
(fig_diff_gfung_A <- ggplot(sims_gfung_A[["wide"]], 
                            aes(x = generation, y = annual_seeds_infected -
                                  annual_seeds_suppressed)) +
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
    aes(y = perennial_seeds_infected - perennial_seeds_suppressed) +
    labs(y = "Effect of disease on *E. virginicus* seed density"))

(fig_diff_gfung_P <- fig_diff_gfung_A %+%
    sims_gfung_P[["wide"]] +
    aes(y = perennial_adults_infected - perennial_adults_suppressed) +
    labs(y = "Effect of disease on *E. virginicus* seed density"))


#### single species values ####

# alpha values
alpha_gfung_A <- params_gfung[["healthy"]][["alpha"]] %>%
  select(alphaPA) %>%
  mutate(disease = "suppressed",
         .draw = params_gfung[["healthy"]][["draws"]]$.draw) %>%
  full_join(params_gfung[["disease"]][["alpha"]] %>%
              select(alphaPA) %>%
              mutate(disease = "ambient",
                     .draw = params_gfung[["disease"]][["draws"]]$.draw))

# save last time point
eq_gfung_A <- sims_gfung_A[["long"]] %>%
  filter(generation == gens) %>%
  left_join(alpha_gfung_A) %>%
  mutate(comp_eff = annual_seeds * alphaPA)
  
eq_gfung_A_wide <- sims_gfung_A[["wide"]] %>%
  filter(generation == gens) %>%
  left_join(alpha_gfung_A) %>%
  mutate(comp_eff_diff = alphaPA * (annual_seeds_ambient - 
                                      annual_seeds_suppressed))

eq_ginf_A <- sims_ginf_A[["long"]] %>%
  filter(generation == gens)

eq_gfung_P <- sims_gfung_P[["long"]] %>%
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

ggplot(eq_gfung_A_wide, aes(x = litter_ambient - litter_suppressed)) +
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
