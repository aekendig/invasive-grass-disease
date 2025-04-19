#### outputs #####

# models
# output/evS_background_seed_model_2019_density_exp.rda
# output/evA_background_seed_model_2019_density_exp.rda
# output/mv_background_seed_model_2019_density_exp.rda
# output/ev_seed_model_2019_density_exp.rda
# tables
# output/evS_background_seed_model_2019_density_exp.csv
# output/evA_background_seed_model_2019_density_exp.csv
# output/mv_background_seed_model_2019_density_exp.csv
# output/ev_seed_model_2019_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)

# import data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# Beverton-Holt function
bh_fun <- function(dat_in, a){
  
  # extract values
  xmin = min(dat_in$background_density)
  xmax = max(dat_in$background_density)
  yE = filter(dat_in, sp == "Ev" & background_density == 0) %>%
    summarise(mean_seeds = mean(seeds)) %>%
    as.numeric() %>%
    round()
  print(yE)
  
  yM = filter(dat_in, sp == "Mv" & background_density == 0) %>%
    summarise(mean_seeds = mean(seeds)) %>%
    as.numeric() %>%
    round()
  print(yM)
  
  # create data
  datE<- tibble(x = seq(xmin, xmax, length.out = 100),
                 sp = "Ev") %>%
    mutate(y = yE / (1 + a * x))
  datM <- tibble(x = seq(xmin, xmax, length.out = 100),
                 sp = "Mv") %>%
    mutate(y = yM / (1 + a * x))
  dat = full_join(datE, datM)
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean") +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = y)) +
          facet_wrap(~ sp, scales = "free_y"))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# 2019 list of all plants
# all dead plants were replaced
focD2Dat <- plots %>%
  select(plot, treatment) %>%
  expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
  expand_grid(ID = as.character(c(1, 2, 3))) %>%
  mutate(sp = "Mv",
         age = "seedling") %>%
  full_join(plots %>%
              select(plot, treatment) %>%
              expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
              expand_grid(tibble(ID = c("1", "2", "3", "A"),
                                 age = c(rep("seedling", 3), "adult"))) %>%
              mutate(sp = "Ev")) %>%
  mutate(fungicide = if_else(treatment == "fungicide", 1, 0),
         treatment = fct_relevel(treatment, "water"),
         plotID = paste(site, plot, fungicide, sep = "_"))

# Ev seeds 2019
evSeedD2Dat2 <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(focD2Dat %>%
              filter(sp == "Ev")) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = replace_na(seeds, 0),
         seeds1 = seeds + 1)

# Mv seeds 2019
mvSeedDat <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = replace_na(seeds, 0))

# combine Ev adults and Mv
evASeedDat <- evSeedD2Dat2 %>% filter(age == "adult")
combSeedDat <- evASeedDat %>% full_join(mvSeedDat)

# split by background
combEvASeedDat <- combSeedDat %>% filter(background %in% c("none", "Ev adult"))
combEvSSeedDat <- combSeedDat %>% filter(background %in% c("none", "Ev seedling"))
combMvSeedDat <- combSeedDat %>% filter(background %in% c("none", "Mv seedling"))


#### fit Beverton-Holt models ####

# initial visualizations
combEvASeedDat %>% filter(treatment == "water") %>%
  bh_fun(a = 0)
combEvASeedDat %>% filter(treatment == "fungicide") %>%
  bh_fun(a = 0)

combEvSSeedDat %>% filter(treatment == "water") %>%
  bh_fun(a = 0)
combEvSSeedDat %>% filter(treatment == "fungicide") %>%
  bh_fun(a = 0.01)

combMvSeedDat %>% filter(treatment == "water") %>%
  bh_fun(a = 0.03)
combMvSeedDat %>% filter(treatment == "fungicide") %>%
  bh_fun(a = 0.03)

# check prior distribution
val <- seq(0, 100, length.out = 50)
dens <- dexp(val, 5)
plot(val, dens, type = "l")

# fit models
combEvASeedMod <- brm(data = combEvASeedDat, family = gaussian,
                  bf(seeds ~ s0/(1 + alpha * background_density),
                     s0 ~ sp * fungicide + (1 | site/plotID), 
                     alpha ~ treatment + 0, 
                     nl = T),
                  prior <- c(prior(normal(81, 10), coef = 'Intercept', 
                                   nlpar = "s0"),
                             prior(normal(1119, 100), coef = 'spMv',
                                   nlpar = "s0"),
                             prior(normal(52, 10), coef = 'fungicide', 
                                   nlpar = "s0"),
                             prior(normal(-358, 10), coef = 'spMv:fungicide', 
                                   nlpar = "s0"),
                             prior(exponential(1), lb = 0, nlpar = "alpha"),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(combEvASeedMod)

combEvSSeedMod <- update(combEvASeedMod, newdata = combEvSSeedDat)
mod_check_fun(combEvSSeedMod)

combMvSeedMod <- update(combEvASeedMod, newdata = combMvSeedDat)
mod_check_fun(combMvSeedMod)

# save models
save(combEvASeedMod, file = "output/evA_background_seed_model_2019_density_exp.rda")
save(combEvSSeedMod, file = "output/evS_background_seed_model_2019_density_exp.rda")
save(combMvSeedMod, file = "output/mv_background_seed_model_2019_density_exp.rda")

# save tables
write_csv(tidy(combEvASeedMod, conf.method = "HPDinterval"), 
          "output/evA_background_seed_model_2019_density_exp.csv")
write_csv(tidy(combEvSSeedMod, conf.method = "HPDinterval"), 
          "output/evS_background_seed_model_2019_density_exp.csv")
write_csv(tidy(combMvSeedMod, conf.method = "HPDinterval"), 
          "output/mv_background_seed_model_2019_density_exp.csv")


#### fit Ev age model ####

# initial visualization
ggplot(evSeedD2Dat2, aes(x = age, y = seeds)) +
  geom_boxplot() +
  facet_wrap(~ treatment)

evSeedD2Dat2 %>%
  filter(age == "adult" & fungicide == 0) %>%
  pull(seeds1) %>%
  log() %>%
  mean()

ggplot(evSeedD2Dat2, aes(x = seeds1)) +
  geom_density()


# fit model
evSeedMod <- brm(data = evSeedD2Dat2, family = lognormal,
                  seeds1 ~ age * fungicide + (1|site/plotID),
                  prior <- c(prior(normal(4, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.99))
mod_check_fun(evSeedMod)

# save model
save(evSeedMod, file = "output/ev_seed_model_2019_density_exp.rda")

# save table
write_csv(tidy(evSeedMod, conf.method = "HPDinterval"), 
          "output/ev_seed_model_2019_density_exp.csv")
