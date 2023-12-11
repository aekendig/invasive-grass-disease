##### outputs ####

# mv_litter_establishment_data_2018_litter_exp.csv
# mv_litter_establishment_model_2018_litter_exp.rda


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)

# import data
estL1Dat <- read_csv("data/both_germination_disease_jul_2018_litter_exp.csv")
plotsL <- read_csv("data/plot_treatments_2018_litter_exp.csv")

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

# edit variables
# remove unnecessary variables
plotsL2 <- plotsL %>%
  mutate(live = if_else(litter == "live", 1, 0),
         sterilized = ifelse(live == 0, "sterilized", "live") %>%
           fct_relevel("sterilized"),
         litter.g.m2 = litter_weight.g,
         litter.g.cm2 = litter.g.m2/10000) %>%
  select(-c(flag_color, justification, litter_weight.g))

# select plots with seeds added only
mvEstL1Dat <- estL1Dat %>%
  filter(seeds_added == "yes") %>%
  left_join(plotsL2) %>%
  select(-c(date, seeds_added, ev_germ, ev_infec)) %>%
  mutate(mv_germ_planted = mv_germ - mv_germ_ev, # (planted + background) - only background (none planted in Ev section)
         mv_seeds = 200 + mv_germ_ev,
         prop_germ = mv_germ_planted/mv_seeds)


#### model ####

# initial visualization
ggplot(mvEstL1Dat, aes(x = litter.g.cm2, y = prop_germ, color = sterilized)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0)

ggplot(mvEstL1Dat, aes(mv_germ_ev, mv_germ)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

# initial fit
mvEstL1Mod <- brm(mv_germ ~ mv_germ_ev + litter.g.cm2 +  litter.g.cm2:live + (1|site),
                  data = mvEstL1Dat, family = negbinomial,
                  prior = c(prior(normal(0, 10), class = Intercept),
                            prior(normal(0, 10), class = b),
                            prior(normal(1, 10), coef = "mv_germ_ev")),
                  iter = 6000, warmup = 1000, chains = 1)
summary(mvEstL1Mod)

# increase chains
mvEstL1Mod <- update(mvEstL1Mod, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.999))
mod_check_fun(mvEstL1Mod)

# save
save(mvEstL1Mod, file = "output/mv_litter_establishment_model_2018_litter_exp.rda")

# table
write_csv(tidy(mvEstL1Mod), "output/mv_litter_establishment_model_2018_litter_exp.csv")


#### figure ####

# posterior draws
mvEstL1Draws <- as_draws_df(mvEstL1Mod)

# edit
estDraws <- tibble(beta = mvEstL1Draws$b_litter.g.cm2,
                   treatment = "sterilized") %>%
  full_join(mvEstL1Draws %>%
              transmute(beta = b_litter.g.cm2 + `b_litter.g.cm2:live`,
                        treatment = "live"))

# figure
ggplot(estDraws, aes(x = treatment, y = `beta`)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  labs(y = "Litter sensitivity") +
  fig_theme +
  theme(axis.title.x = element_blank())

#### start here ####
# values for text
# update code tracking sheet
