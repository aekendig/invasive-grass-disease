##### outputs ####


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(broom.mixed)
library(ggtext)
library(patchwork)

# import plot information
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
  print(plot(mod))
  
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### format data ####

# missing data
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling

# add columns
# remove missing
mvBioDat <- mvBioD2Dat %>%
  filter(!is.na(biomass_weight.g)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotID = paste(site, plot, fungicide, sep = "_"),
         log_bio = log(biomass_weight.g))

evBioDat <- evBioD2Dat %>%
  filter(!is.na(weight)) %>%
  rename(biomass_weight.g = weight) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotID = paste(site, plot, fungicide, sep = "_"),
         log_bio = log(biomass_weight.g))


#### Mv model ####

# initial visualizations
ggplot(mvBioDat, aes(fungicide, biomass_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2)

ggplot(mvBioDat, aes(as.factor(plot), biomass_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ treatment)

ggplot(mvBioDat, aes(x = biomass_weight.g)) +
  geom_density() +
  facet_wrap(~ treatment)

ggplot(mvBioDat, aes(x = log_bio)) +
  geom_density()

mvBioDat %>%
  filter(fungicide == 0) %>%
  pull(log_bio) %>%
  mean() 

# model
mvBioD2Mod <- brm(data = mvBioDat, family = gaussian,
                   log_bio ~ fungicide + (1|site/plotID), 
                   prior <- c(prior(normal(2.5, 0.5), class = "Intercept"),
                              prior(normal(0, 1), class = "b"),
                              prior(exponential(1), class = "sd")), # use default for sigma
                  control = list(adapt_delta = 0.9999), 
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvBioD2Mod)
save(mvBioD2Mod, 
     file = "output/mv_biomass_fungicide_model_2019_density_exp.rda")


#### Ev model ####

# initial visualizations
ggplot(evBioDat, aes(fungicide, biomass_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ age)

ggplot(evBioDat, aes(x = biomass_weight.g)) +
  geom_density() +
  facet_wrap(~ treatment)

ggplot(evBioDat, aes(x = log_bio)) +
  geom_density()

evBioDat %>%
  filter(fungicide == 0 & age == "adult") %>%
  pull(log_bio) %>%
  mean() 

# model
# use default prior for sigma
evBioD2Mod <- brm(data = evBioDat, family = gaussian,
                  log_bio ~ fungicide + age + (1|site/plotID), 
                  prior <- c(prior(normal(1.4, 0.5), class = "Intercept"), # reduce extreme predictions
                             prior(normal(0, 2), class = "b"), # broad enough for age effect
                             prior(exponential(1), class = "sd")), # reduce extreme predictions
                  control = list(adapt_delta = 0.999), 
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(evBioD2Mod)
save(evBioD2Mod, 
     file = "output/ev_biomass_fungicide_model_2019_density_exp.rda")


#### tables and figures ####

# load
load("output/mv_biomass_fungicide_model_2019_density_exp.rda")
load("output/ev_biomass_fungicide_model_2019_density_exp.rda")

# tables
write_csv(tidy(mvBioD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_biomass_fungicide_model_2019_density_exp.csv")
write_csv(tidy(evBioD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/ev_biomass_fungicide_model_2019_density_exp.csv")

# prediction data
pred_dat_trt <- mvBioDat %>%
  distinct(fungicide, treatment) %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"),
         age = NA) # predictions are averaged over levels

# posterior draws
mvBioD2Draws <- pred_dat_trt %>%
  add_epred_draws(mvBioD2Mod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(biomass = exp(.epred))

evBioD2Draws <- pred_dat_trt %>%
  add_epred_draws(evBioD2Mod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(biomass = exp(.epred))

# figures
mv_bio_fig <- ggplot(mvBioD2Draws, aes(x = biomass, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(y = "Disease treatment", x = "*M. vimineum* biomass (g)") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ev_bio_fig <- ggplot(evBioD2Draws, aes(x = biomass, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(y = "Disease treatment", x = "*E. virginicus* biomass (g)") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

# combine figures
bio_fung_fig <- mv_bio_fig + ev_bio_fig + 
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 1, axes = "collect", guides = "collect") &
  theme(legend.position = "bottom")
  
ggsave("output/biomass_fungicide_figure_2019_density_exp.png",
       bio_fung_fig, width = 6, height = 3.2)
  
