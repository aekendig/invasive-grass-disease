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


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(GGally)
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
write_csv(tidy(evGermMod, conf.method = "HPDinterval"), 
          "output/ev_germination_fungicide_model_2018_2019_density_exp.csv")


#### fungicide effect figure ####

# posterior draws
mvGermD1Draws <- as_draws_df(mvGermD1Mod)
evGermDraws <- as_draws_df(evGermMod)

# combine
germDraws <- tibble(sp = "M. vimineum",
                    int = mvGermD1Draws$b_Intercept,
                    beta = mvGermD1Draws$b_fungicide) %>%
  full_join(tibble(sp = "E. virginicus",
                   int = evGermDraws$b_Intercept,
                   beta = evGermDraws$b_fungicide)) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         prob_int = 100 * plogis(int),
         prob_fung = 100 * plogis(int + beta),
         prob_diff = prob_fung - prob_int)

# figure
germ_fung_fig <- ggplot(germDraws, aes(x = sp, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  labs(y = "Fungicide effect on germination (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))

save(germ_fung_fig, file = "output/germination_fungicide_figure_2018_2019_density_exp.rda")

# values for text
germDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_int)

germDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_fung)

germDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_diff)


#### seed infection figure ####

# posterior draws
mvGermInfD1Draws <- as_draws_df(mvGermInfD1Mod)
mvPropDarkDraws <- as_draws_df(mvPropDarkMod)
mvPropLightDraws <- as_draws_df(mvPropLightMod)

# combine
infDraws <- tibble(fungi = "dark",
                   response = "Infection effect on germination (log-odds)",
                   beta = mvGermInfD1Draws$b_prop_dark,
                   int = mvGermInfD1Draws$b_Intercept) %>%
  full_join(tibble(fungi = "light",
                   response = "Infection effect on germination (log-odds)",
                   beta = mvGermInfD1Draws$b_prop_light,
                   int = mvGermInfD1Draws$b_Intercept)) %>%
  full_join(tibble(fungi = "dark",
                   response = "Fungicide effect on infection (log-odds)",
                   beta = mvPropDarkDraws$b_fungicide,
                   int = mvPropDarkDraws$b_Intercept) %>%
              full_join(tibble(fungi = "light",
                               response = "Fungicide effect on infection (log-odds)",
                               beta = mvPropLightDraws$b_fungicide,
                               int = mvPropLightDraws$b_Intercept))) %>%
  mutate(prob_int = 100 * exp(int) / (1 + exp(int)),
         prob_beta = 100 * exp(int + beta) / (1 + exp(int + beta)),
         prob_diff = prob_beta - prob_int,
         response = fct_relevel(response, "Fungicide effect on infection (log-odds)"))

# figure
germ_inf_fig <- ggplot(infDraws, aes(x = fungi, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  facet_wrap(~ response, strip.position = "left",
             scales = "free_y", ncol = 1) +
  labs(x = "Seed fungi color") +
  fig_theme +
  theme(axis.title.y = element_blank())

ggsave("output/germination_infection_figure_2018_density_exp.png",
       germ_inf_fig,
       width = 3, height = 5.5)

# values for text
infDraws %>% 
  group_by(fungi, response) %>%
  mean_hdci(prob_int)

infDraws %>% 
  group_by(fungi, response) %>%
  mean_hdci(prob_beta)

infDraws %>% 
  group_by(fungi, response) %>%
  mean_hdci(prob_diff)
