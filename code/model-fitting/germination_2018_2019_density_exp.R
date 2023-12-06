##### outputs ####

# mv_germination_fungicide_model_2018_density_exp.rda
# mv_germination_fungicide_model_data_2018_density_exp.csv
# ev_germination_fungicide_model_2018_2019_density_exp.rda
# ev_germination_fungicide_model_data_2018_2019_density_exp.rda
# mv_germination_infection_model_2018_density_exp.csv
# mv_seed_infection_dark_fungicide_model_2018_density_exp.csv
# mv_germination_infection_figure_2018_density_exp.pdf


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)
library(car)
library(tidybayes)
library(broom.mixed)

# import data
mvGermDat <- read_csv("intermediate-data/mv_germination_disease_2018_density_exp.csv")
evGermDat <- read_csv("intermediate-data/ev_germination_2018_2019_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# Mv data
# calculate proportions and make study design columns
mvGermD1Dat <- mvGermDat %>%
  mutate(prop_germ = germination_final / seeds,
         prop_dark = seeds_dark / seeds,
         prop_light = seeds_light / seeds,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)))

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 ) %>%
  data.frame()

# trials per plot
mvGermD1Dat %>%
  count(plotf) %>%
  rename(trials = "n") %>%
  count(trials)

# Ev data
# calculate proportions and make study design columns
evGermDat2 <- evGermDat %>%
  mutate(prop_germ = germinants/seeds_planted,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         yearf = as.factor(year))
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()


#### Mv models ####

# initial visualization
mvGermD1Dat %>%
  select(prop_dark, prop_light, prop_germ) %>%
  ggpairs()

ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# light/dark infection correlation
cor.test(~ prop_dark + prop_light, data = mvGermD1Dat) # not correlated

# fungicide model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ fungicide + (1|plotf), # can't converge with site/plot
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermD1Mod)
save(mvGermD1Mod, file = "output/mv_germination_fungicide_model_2018_density_exp.rda")

# seed infection model
mvGermInfD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ prop_dark + prop_light + (1|plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermInfD1Mod)
save(mvGermInfD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")

# proportion of seeds with dark infection
mvPropDarkMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ fungicide + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod)
save(mvPropDarkMod, file = "output/mv_seed_infection_dark_model_2018_density_exp.rda")

# proportion of seeds with light infection
mvPropLightMod <- brm(data = mvGermD1Dat, family = binomial,
                      seeds_light | trials(seeds) ~ fungicide + (1|plotf),
                      prior <- c(prior(normal(0, 10), class = "Intercept"),
                                 prior(normal(0, 10), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropLightMod)
save(mvPropLightMod, file = "output/mv_seed_infection_light_model_2018_density_exp.rda")

# load
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
load("output/mv_germination_infection_model_2018_density_exp.rda")
load("output/mv_seed_infection_dark_model_2018_density_exp.rda")
load("output/mv_seed_infection_light_model_2018_density_exp.rda")

# tables
write_csv(tidy(mvGermD1Mod), "output/mv_germination_fungicide_model_2018_density_exp.csv")
write_csv(tidy(mvGermInfD1Mod), "output/mv_germination_infection_model_2018_density_exp.csv")
write_csv(tidy(mvPropDarkMod), "output/mv_seed_infection_dark_model_2018_density_exp.csv")


#### Ev models ####

# initial visualization
ggplot(evGermDat2, aes(year, prop_germ, color = age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(0.2))

ggplot(evGermDat2, aes(treatment, prop_germ)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                 germinants | trials(seeds_planted) ~ fungicide + (1|site) + (1|yearf),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, 
                 control = list(adapt_delta = 0.999, max_treedepth = 15))
# mod_check_fun(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# table
write_csv(tidy(evGermMod), "output/ev_germination_fungicide_model_2018_2019_density_exp.csv")


#### fungicide effect figure ####

# posterior draws
mvGermD1Draws <- as_draws_df(mvGermD1Mod)
evGermDraws <- as_draws_df(evGermMod)

# combine
germDraws <- tibble(sp = "M. vimineum",
                    fung_eff = mvGermD1Draws$b_fungicide) %>%
  full_join(tibble(sp = "E. virginicus",
                   fung_eff = evGermDraws$b_fungicide))

# figure
ggplot(germDraws, aes(x = sp, y = fung_eff)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = NA) +
  labs(y = "Fungicide effect on germination") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))


#### seed infection figure ####

# posterior draws
mvGermInfD1Draws <- as_draws_df(mvGermInfD1Mod)
mvPropDarkDraws <- as_draws_df(mvPropDarkMod)
mvPropLightDraws <- as_draws_df(mvPropLightMod)

# combine
infDraws <- tibble(fungi = "dark",
                   response = "Infection effect on germination",
                   eff = mvGermInfD1Draws$b_prop_dark) %>%
  full_join(tibble(fungi = "light",
                   response = "Infection effect on germination",
                   eff = mvGermInfD1Draws$b_prop_light)) %>%
  full_join(tibble(fungi = "dark",
                   response = "Fungicide effect on infection",
                   eff = mvPropDarkDraws$b_fungicide) %>%
              full_join(tibble(fungi = "light",
                               response = "Fungicide effect on infection",
                               eff = mvPropLightDraws$b_fungicide)))

# figure
ggplot(infDraws, aes(x = fungi, y = eff)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = NA) +
  facet_wrap(~ response, strip.position = "left",
             scales = "free") +
  labs(x = "Seed fungi color") +
  fig_theme +
  theme(axis.title.y = element_blank())


#### START HERE: values for text ####

# fungicide effects
mean_hdci(evGermDraws$b_fungicide)

# Mv fungicide effect on infection
hypothesis(mvPropDarkMod, "plogis(Intercept + fungicide) - plogis(Intercept) = 0")

# Mv infection effect on germination
filter(mvGermSim, prop_dark == min(prop_dark))
filter(mvGermSim, prop_dark == max(prop_dark))

# Mv infection effect on germination
set.seed(184)
posterior_predict(mvGermD1Mod,
                  newdata = filter(mvGermD1Dat, 
                                   prop_dark %in% c(max(mvGermD1Dat$prop_dark), min(mvGermD1Dat$prop_dark))) %>%
                    select(prop_dark) %>%
                    unique() %>%
                    mutate(prop_light = 0,
                           seeds = 30,
                           plotf = "A"),
                  allow_new_levels = T) %>%
  as_tibble(.name_repair = ~ c("min_inf", "max_inf")) %>% # min_inf is first row
  mutate(min_inf_prop = min_inf / 30,
         max_inf_prop = max_inf / 30,
         inf_eff = 100 * (max_inf_prop - min_inf_prop) / min_inf_prop) %>%
  select(min_inf_prop, max_inf_prop, inf_eff) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)
