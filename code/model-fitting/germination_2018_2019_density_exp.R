#### outputs ####
# models
# output/mv_germination_fungicide_model_2018_density_exp.rda
# output/mv_germination_infection_model_2018_density_exp.rda
# output/mv_seed_infection_dark_model_2018_density_exp.rda
# output/mv_seed_infection_light_model_2018_density_exp.rda
# output/ev_germination_fungicide_model_2018_2019_density_exp.rda
# tables
# output/mv_germination_fungicide_model_2018_density_exp.csv
# output/mv_germination_infection_model_2018_density_exp.csv
# output/mv_seed_infection_dark_model_2018_density_exp.csv
# output/ev_germination_fungicide_model_2018_2019_density_exp.csv
# figures
# output/germination_fungicide_figure_2018_2019_density_exp.rda
# output/germination_infection_figure_2018_density_exp.png


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(GGally)
library(broom.mixed)
library(ggtext)
library(patchwork)

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
         plotID = paste(site, plot, fungicide, sep = "_"))

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 ) %>%
  data.frame()

# trials per plot
mvGermD1Dat %>%
  count(plotID) %>%
  rename(trials = "n") %>%
  count(trials)

# Ev data
# calculate proportions and make study design columns
evGermDat2 <- evGermDat %>%
  mutate(prop_germ = germinants/seeds_planted,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotID = paste(site, plot, fungicide, sep = "_"),
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

ggplot(mvGermD1Dat, aes(treatment, prop_germ)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

ggplot(mvGermD1Dat, aes(site, prop_germ, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# light/dark infection correlation
cor.test(~ prop_dark + prop_light, data = mvGermD1Dat) # not correlated

# fungicide model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ fungicide + (1|plotID), # can't converge with site/plot
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermD1Mod)
save(mvGermD1Mod, file = "output/mv_germination_fungicide_model_2018_density_exp.rda")

# seed infection model
mvGermInfD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ prop_dark + prop_light + (1|plotID),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermInfD1Mod)
save(mvGermInfD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")

# proportion of seeds with dark infection
mvPropDarkMod <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ fungicide + (1|plotID),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvPropDarkMod)
save(mvPropDarkMod, file = "output/mv_seed_infection_dark_model_2018_density_exp.rda")

# proportion of seeds with light infection
mvPropLightMod <- brm(data = mvGermD1Dat, family = binomial,
                      seeds_light | trials(seeds) ~ fungicide + (1|plotID),
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
write_csv(tidy(mvGermD1Mod, conf.method = "HPDinterval"), 
               "output/mv_germination_fungicide_model_2018_density_exp.csv")
write_csv(tidy(mvGermInfD1Mod, conf.method = "HPDinterval"), 
          "output/mv_germination_infection_model_2018_density_exp.csv")
write_csv(tidy(mvPropDarkMod, conf.method = "HPDinterval"), 
          "output/mv_seed_infection_dark_model_2018_density_exp.csv")


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
                 germinants | trials(seeds_planted) ~ fungicide + yearf + age + (1|site/plotID),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.999, max_treedepth = 15))
mod_check_fun(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# table
write_csv(tidy(evGermMod, conf.method = "HPDinterval"), 
          "output/ev_germination_fungicide_model_2018_2019_density_exp.csv")


#### figures ####

# prediction data
pred_dat_trt <- mvGermD1Dat %>%
  distinct(fungicide, treatment) %>%
  mutate(treatment = fct_relevel(treatment, "water")) %>%
  expand_grid(tibble(seeds = round(mean(mvGermD1Dat$seeds)),
                     seeds_planted = round(mean(evGermDat2$seeds_planted)),
                     age = NA, # predictions are averaged over levels
                     yearf = NA))

pred_dat_inf_dark <- mvGermD1Dat %>%
  distinct(prop_dark) %>%
  expand_grid(tibble(seeds = round(mean(mvGermD1Dat$seeds)),
                     prop_light = 0))

pred_dat_inf_light <- mvGermD1Dat %>%
  distinct(prop_light) %>%
  expand_grid(tibble(seeds = round(mean(mvGermD1Dat$seeds)),
                     prop_dark = 0))

# posterior draws
mvGermD1Draws <- pred_dat_trt %>%
  add_epred_draws(mvGermD1Mod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds)
mvGermInfDarkD1Draws <- pred_dat_inf_dark %>%
  add_epred_draws(mvGermInfD1Mod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds)
mvGermInfLightD1Draws <- pred_dat_inf_light %>%
  add_epred_draws(mvGermInfD1Mod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds)
mvPropDarkDraws <- pred_dat_trt %>%
  add_epred_draws(mvPropDarkMod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds)
mvPropLightDraws <- pred_dat_trt %>%
  add_epred_draws(mvPropLightMod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds)
evGermDraws <- pred_dat_trt %>%
  add_epred_draws(evGermMod, re_formula = NA) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds_planted)

# mv fig
mv_germ_fig <- ggplot(mvGermD1Draws, aes(x = treatment, y = germ_frac)) +
  stat_pointinterval(fatten_point = 2,
                     point_interval = mean_hdi,
                     .width = c(0.80, 0.95)) +
  labs(x = "Disease treatment", y = "*M. vimineum* germination fraction") +
  fig_theme +
  theme(axis.title.y = element_markdown())

mv_dark_fig <- ggplot(mvGermInfDarkD1Draws, 
                      aes(x = prop_dark, y = germ_frac)) +
  stat_lineribbon(point_interval = mean_hdi,
                  .width = 0.95, alpha = 0.5, fill = "black") +
  labs(x = "Proportion infected with dark fungi", 
       y = "Germination fraction") +
  fig_theme +
  theme(axis.title.y = element_markdown())

mv_light_fig <- mv_dark_fig %+%
  mvGermInfLightD1Draws +
  aes(x = prop_light) +
  labs(x = "Proportion infected with light fungi")

mv_prop_dark_fig <- mv_germ_fig %+%
  mvPropDarkDraws +
  labs(y = "Proportion infected with dark fungi")

mv_prop_light_fig <- mv_germ_fig %+%
  mvPropLightDraws +
  labs(y = "Proportion infected with light fungi")

ev_germ_fig <- mv_germ_fig %+%
  evGermDraws +
  labs(y = "*E. virginicus* germination fraction")

# combine figures and save
germ_fung_fig <- mv_germ_fig + ev_germ_fig + plot_annotation(tag_levels = "A")
ggsave("output/germination_fungicide_figure_2018_2019_density_exp.png",
       germ_fung_fig, width = 6, height = 3)

germ_inf_fig <- mv_prop_dark_fig + mv_dark_fig + 
  mv_prop_light_fig + mv_light_fig +
  plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = "A")
ggsave("output/germination_infection_figure_2018_density_exp.png",
       germ_inf_fig, width = 6, height = 6)
