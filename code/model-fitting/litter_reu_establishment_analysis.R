#### outputs ####

# litter_reu_mv_establishment_model.rda
# litter_reu_ev_establishment_model.rda

#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)
library(ggtext)
library(patchwork)

# import data
mvEstDat <- read_csv("data/litter_reu_mv_establishment_data.csv")
evEstDat <- read_csv("data/litter_reu_ev_establishment_data.csv")

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
    pull(prop_est) %>%
    mean()
  print(V)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(v = V / (1 + b * x))
  
  # plot
  print(ggplot(dat_in, aes(x = litter.g.m2, y = prop_est)) +
          stat_summary(geom = "point", fun = "mean") +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = v)))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

litterDat <- mvEstDat %>%
  select(Litter.g) %>%
  unique() %>%
  arrange() %>%
  mutate(litter.g.m2 = c(0, 50, 100, 200))

# remove competition treatments
# rename columns
# add litter column
mvEstDat2 <- mvEstDat %>%
  filter(Competition == 0) %>%
  rename(prop_est = PropEstMvDenCor) %>%
  left_join(litterDat)

evEstDat2 <- evEstDat %>%
  filter(Competition == 0) %>%
  rename(prop_est = PropEstEv) %>%
  left_join(litterDat)


#### fit models ####

# initial visualizations
bh_fun(mvEstDat2, b = 0.001)
bh_fun(evEstDat2, b = 0.0004)
# these don't work very well because the litter isn't high enough to get
# the establishment fraction to go down to zero


#### old code ####
ggplot(MvEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()

ggplot(EvEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()

# beta prior
x <- 0:200
y <- 0.85/(1 + 0.001 * x)
plot(x, y, type = "l")

val <- seq(0, 1, length.out = 50)
dens <- dexp(val, 1)
plot(val, dens, type = "l")

# Mv
mv_bh_mod <- brm(data = MvEstDat2, family = gaussian,
                 bf(PropEst ~ e0/(1 + beta * Litter.g.m2),
                    e0 ~ 1, 
                    beta ~ 1, 
                    nl = T),
                 prior <- c(prior(normal(0.85, 1), nlpar = "e0", lb = 0),
                            prior(exponential(1), nlpar = "beta", lb = 0),
                            prior(cauchy(0, 1), class = sigma)),
                 iter = 6000, warmup = 1000, chains = 1, cores = 1)

prior_summary(mv_bh_mod)
summary(mv_bh_mod, digits = 5)
mv_bh_mod <- update(mv_bh_mod, chains = 3)
summary(mv_bh_mod)
plot(mv_bh_mod)
fixef(mv_bh_mod)

# Ev
ev_bh_mod <- update(mv_bh_mod, newdata = EvEstDat2)

summary(ev_bh_mod)
plot(ev_bh_mod)
fixef(ev_bh_mod)

# save
save(mv_bh_mod, file = "output/litter_reu_mv_establishment_model.rda")
save(ev_bh_mod, file = "output/litter_reu_ev_establishment_model.rda")
