##### outputs ####

# ev_establishment_model_2019_litter_exp.rda
# ev_establishment_model_2018_litter_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)

# import data
estL2Dat <- read_csv("data/both_germination_disease_jun_2019_litter_exp.csv")
plots <- read_csv("data/litter_weight_apr_2019_litter_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}


#### edit data ####

# edit plot data
# put the removed litter into the addition litter
# convert units
# remove unnecessary variables
plots2 <- plots %>%
  spread(key = treatment, value = litter_weight.lb) %>%
  mutate(addition = addition + removal,
         removal = 0) %>%
  gather(key = "treatment", value = "litter_weight.lb", -c(date, site, block)) %>%
  mutate(litter_weight.g = litter_weight.lb * 453.592,
         litter.g.m2 = litter_weight.g / 4, # 4 because the plots are 2m^2
         treatment = fct_relevel(treatment, "removal", "control"),
         plot = as.factor(paste(site, block, sep = "_"))) %>%
  select(-c(date))

# Ev planting data
plant <- tibble(treatment = c("removal", "control", "addition"),
                ev_tot = c(50, 26, 26)) 

# June germination data
# none of the germinants were infected
estL2Dat2 <- estL2Dat %>%
  select(-c(date, flag_color, mv_germ, mv_infec)) %>%
  full_join(plots2) %>%
  full_join(plant) %>%
  mutate(ev_prop_germ = ev_germ / ev_tot,
         treatment = fct_relevel(treatment, "removal", "control"))


#### fit model ####

# visualize
ggplot(estL2Dat2, aes(x = litter.g.m2, y = ev_prop_germ)) +
  geom_point()

# germination with litter
estL2Dat2 %>%
  filter(litter.g.m2 > 0 & ev_germ > 0)
# only one

# use data to estimate establishment without litter
estL2Dat3 <- filter(estL2Dat2, litter.g.m2 == 0)

# initial fit
evEstL2Mod <- brm(ev_germ ~ 1, data = estL2Dat3, family = poisson,
                  prior = c(prior(normal(200, 100), class = Intercept)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(evEstL2Mod)

# save model
save(evEstL2Mod, file = "output/ev_establishment_model_2019_litter_exp.rda")

# table
write_csv(tidy(evEstL2Mod), "output/ev_establishment_model_2019_litter_exp.csv")

# load
load("output/ev_establishment_model_2019_litter_exp.rda")


#### values for text ####

# posterior draws
evEstL2Draws <- as_draws_df(evEstL2Mod) %>% as_tibble()

# establishment
evEstL2Draws %>%
  mean_hdci(exp(b_Intercept) / 50)
