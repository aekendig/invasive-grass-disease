##### outputs ####

# mv_establishment_model_2018_litter_exp.rda
# mv_establishment_model_2018_litter_exp.csv

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)
library(ggtext)

# import data
estL1Dat <- read_csv("data/both_germination_disease_jul_2018_litter_exp.csv")
plots <- read_csv("data/plot_treatments_2018_litter_exp.csv")

# load germination model and data
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
mvGermDat <- read_csv("intermediate-data/mv_germination_disease_2018_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# Beverton-Holt function
bh_fun <- function(dat_in, b){
  
  # extract values
  xmin = min(dat_in$litter.g.m2)
  xmax = max(dat_in$litter.g.m2)
  V = filter(dat_in, litter.g.m2 == xmin) %>%
    pull(prop_est_adj2) %>%
    mean()
  print(V)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(v = V / (1 + b * x))
  
  # plot
  print(ggplot(dat_in, aes(x = litter.g.m2, y = prop_est_adj2)) +
          stat_summary(geom = "point", fun = "mean") +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = v)))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# edit variables
# remove unnecessary variables
plots2 <- plots %>%
  mutate(live = if_else(litter == "live", 1, 0),
         sterilized = ifelse(live == 0, "sterilized", "live") %>%
           fct_relevel("sterilized"),
         litter.g.m2 = litter_weight.g) %>%
  select(-c(flag_color, justification, litter_weight.g))

# average germination with fungicide
prop_germ <- tibble(fungicide = 1,
                  seeds = round(mean(mvGermDat$seeds))) %>%
  add_epred_draws(mvGermD1Mod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds) %>%
  pull(germ_frac) %>%
  mean()
  
# select plots with seeds added only
mvEstL1Dat <- estL1Dat %>%
  filter(seeds_added == "yes") %>%
  left_join(plots2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  mutate(mv_est_planted = mv_germ - mv_germ_ev, # (planted + background) - only background (none planted in Ev section)
         prop_est = mv_est_planted/200,
         prop_est_adj = if_else(prop_est < 0, 0, prop_est),
         prop_est_adj2 = prop_est_adj/prop_germ)

range(mvEstL1Dat$prop_est) # some negative - added in adjustment above


#### model ####

# initial visualization
ggplot(mvEstL1Dat, aes(x = litter.g.m2, y = prop_est_adj, color = sterilized)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0)

ggplot(mvEstL1Dat, aes(x = litter.g.m2, y = prop_est_adj2, color = sterilized)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0)

# beverton-holt
mvEstL1Dat %>% filter(sterilized == "sterilized") %>%
  bh_fun(b = 0.008)
mvEstL1Dat %>% filter(sterilized == "live" | litter.g.m2 == 0) %>%
  bh_fun(b = 0.01)

# check prior distribution
val <- seq(0, 1, length.out = 50)
dens <- dexp(val, 7)
plot(val, dens, type = "l")

# fit model
mvEstL1Mod <- brm(data = mvEstL1Dat, family = gaussian,
                      bf(prop_est_adj2 ~ e0/(1 + beta * litter.g.m2),
                         e0 ~ 1 + (1 | site), 
                         beta ~ sterilized + 0, 
                         nl = T),
                      prior <- c(prior(normal(0.8, 0.5), coef = 'Intercept', 
                                       nlpar = "e0"),
                                 prior(exponential(7), lb = 0, nlpar = "beta")),
                      iter = 6000, warmup = 1000, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.99))
mod_check_fun(mvEstL1Mod)

# save
save(mvEstL1Mod, file = "output/mv_establishment_model_2018_litter_exp.rda")


#### figures and tables ####

# load
load("output/mv_establishment_model_2018_litter_exp.rda")

# table
write_csv(tidy(mvEstL1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_establishment_model_2018_litter_exp.csv")

# prediction data
mvEstL1Draws <- mvEstL1Dat %>%
  distinct(live, sterilized) %>%
  expand_grid(litter.g.m2 = seq(0, max(mvEstL1Dat$litter.g.m2), 
                                length.out = 100)) %>%
  add_epred_draws(mvEstL1Mod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(sterilized, "ambient" = "live"))

# figure
mv_est_fig <- ggplot(mvEstL1Draws, aes(x = litter.g.m2, y = .epred)) +
  stat_lineribbon(aes(fill = trt, color = trt), point_interval = mean_hdci, 
                  .width = 0.95, alpha = 0.5) +
  scale_fill_manual(values = c(coral_pal[2], grey_pal[2]), 
                    name = "Disease treatment") +
  scale_color_manual(values = c(coral_pal[3], grey_pal[3]),
                     name = "Disease treatment") +
  labs(x = "Litter (g/m<sup>2</sup>)", 
       y = "*M. vimineum* establishment") +
  fig_theme +
  theme(axis.title = element_markdown())

# beta values
mvEstL1Beta <- gather_draws(mvEstL1Mod,
                            c(b_beta_sterilizedsterilized, 
                              b_beta_sterilizedlive)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "live"),
                       "ambient", "sterilized") %>%
           fct_relevel("sterilized"))

mv_beta_fig <- ggplot(mvEstL1Beta, aes(x = .value, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdci, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(x = "*M. vimineum* response to litter", y = "Disease treatment") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

# save to combine with EV establishment
save(mv_est_fig, file = "output/mv_establishment_figure_2018_litter_exp.rda")
save(mv_beta_fig, file = "output/mv_establishment_beta_figure_2018_litter_exp.rda")
