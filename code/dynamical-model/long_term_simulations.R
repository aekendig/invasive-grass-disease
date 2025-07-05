#### set-up ####

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
                mutate(disease = "infected")) %>%
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


#### simple simulations ####

# initial conditions (annual seeds, litter, perennial seeds, perennial adults)
initsA <- c(1, 0, 0, 0)
initsP <- c(0, 0, 0, 1)
initsAP <- c(1, 0, 0, 1)

# output lits
sims_gfung_A_h <- list()
sims_gfung_P_h <- list()
sims_gfung_AP_h <- list()

sims_gfung_A_d <- list()
sims_gfung_P_d <- list()
sims_gfung_AP_d <- list()

sims_ginf_A_h <- list()
sims_ginf_AP_h <- list()

sims_ginf_A_d <- list()
sims_ginf_AP_d <- list()

# cycle through parameters
for(i in 1:params_iters){
  
  # fungicide effect, healthy
  sims_gfung_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["healthy"]])
  sims_gfung_P_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["healthy"]])
  sims_gfung_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
                                      parameters = params_gfung[["healthy"]])
  
  # fungicide effect, disease
  sims_gfung_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_gfung[["disease"]])
  sims_gfung_P_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsP,
                                     parameters = params_gfung[["disease"]])
  sims_gfung_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
                                      parameters = params_gfung[["disease"]])
  
  # infection effect, healthy
  sims_ginf_A_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                     parameters = params_ginf[["healthy"]])
  sims_ginf_AP_h[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
                                     parameters = params_ginf[["healthy"]])
  
  # infection effect, disease
  sims_ginf_A_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsA,
                                    parameters = params_ginf[["disease"]])
  sims_ginf_AP_d[[i]] <- dyn_mod_fun(iter = i, gen = gens, init_cond = initsAP,
                                     parameters = params_ginf[["disease"]])
  
}

# combine healthy and disease
sims_gfung_A <- sims_comb_fun(sims_gfung_A_h, sims_gfung_A_d)
sims_ginf_A <- sims_comb_fun(sims_ginf_A_h, sims_ginf_A_d)
sims_gfung_P <- sims_comb_fun(sims_gfung_P_h, sims_gfung_P_d)
sims_gfung_AP <- sims_comb_fun(sims_gfung_AP_h, sims_gfung_AP_d)
sims_ginf_AP <- sims_comb_fun(sims_ginf_AP_h, sims_ginf_AP_d)

# mean hdci figures
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

(fig_dens_gfung_AP_A <- fig_dens_gfung_A %+%
    sims_gfung_AP[["long"]])

(fig_dens_ginf_AP_A <- fig_dens_gfung_A %+%
    sims_ginf_AP[["long"]])

(fig_dens_gfung_AP_S <- fig_dens_gfung_S %+%
    sims_gfung_AP[["long"]])

(fig_dens_gfung_AP_P <- fig_dens_gfung_P %+%
    sims_gfung_AP[["long"]])

(fig_dens_ginf_AP_S <- fig_dens_gfung_S %+%
    sims_ginf_AP[["long"]])

(fig_dens_ginf_AP_P <- fig_dens_gfung_P %+%
    sims_ginf_AP[["long"]])

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
