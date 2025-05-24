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
    pull(prop_germ_adj) %>%
    mean()
  print(V)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(v = V / (1 + b * x))
  
  # plot
  print(ggplot(dat_in, aes(x = litter.g.m2, y = prop_germ_adj)) +
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

# select plots with seeds added only
mvEstL1Dat <- estL1Dat %>%
  filter(seeds_added == "yes") %>%
  left_join(plots2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  group_by(site) %>%
  mutate(mv_germ_ev_avg = mean(mv_germ_ev)) %>%
  ungroup() %>%
  mutate(mv_germ_planted = mv_germ - mv_germ_ev, # (planted + background) - only background (none planted in Ev section)
         prop_germ = mv_germ_planted/200,
         prop_germ_adj = if_else(prop_germ < 0, 0, prop_germ))



#### model ####

# initial visualization
ggplot(mvEstL1Dat, aes(x = litter.g.m2, y = prop_germ, color = sterilized)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0)

range(mvEstL1Dat$prop_germ) # some negative

ggplot(mvEstL1Dat, aes(x = litter.g.m2, y = prop_germ_adj, color = sterilized)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0)

ggplot(mvEstL1Dat, aes(mv_germ_ev, mv_germ)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

# beverton-holt
mvEstL1Dat %>% filter(sterilized == "sterilized") %>%
  bh_fun(b = 0.008)
mvEstL1Dat %>% filter(sterilized == "live") %>%
  bh_fun(b = 0.005)

# check prior distribution
val <- seq(0, 100, length.out = 50)
dens <- dexp(val, 6)
plot(val, dens, type = "l")

# fit models
mvEstL1Mod <- brm(data = mvEstL1Dat, family = gaussian,
                      bf(prop_germ_adj ~ e0/(1 + beta * litter.g.m2),
                         e0 ~ 1 + (1 | site), 
                         beta ~ sterilized + 0, 
                         nl = T),
                      prior <- c(prior(normal(0.6, 0.5), coef = 'Intercept', 
                                       nlpar = "e0"),
                                 prior(exponential(6), lb = 0, nlpar = "beta")),
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

# save to combine with EV establishment
save(mv_est_fig, file = "output/mv_establishment_figure_2018_litter_exp.rda")
