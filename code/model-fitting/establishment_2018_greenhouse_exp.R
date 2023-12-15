#### outputs ####

# mv_establishment_model_2018_greenhouse_exp.rda
# mv_establishment_model_2018_greenhouse_exp.csv
# ev_establishment_model_2018_greenhouse_exp.rda
# ev_establishment_model_2018_greenhouse_exp.csv


#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)

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


#### edit data ####

litterDat <- mvEstDat %>%
  select(Litter.g) %>%
  unique() %>%
  arrange() %>%
  mutate(Litter.g.m2 = c(0, 50, 100, 200))

# remove competition treatments
# rename columns
# add litter column
mvEstDat2 <- mvEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstMvDenCor) %>%
  left_join(litterDat)

evEstDat2 <- evEstDat %>%
  filter(Competition == 0) %>%
  rename(PropEst = PropEstEv) %>%
  left_join(litterDat)

# initial visualizations
ggplot(mvEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()

ggplot(evEstDat2, aes(x = Litter.g.m2, y = PropEst)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0, position = position_dodge(0.1)) +
  stat_summary(fun = "mean", geom = "point", position = position_dodge(0.1)) +
  theme_bw()


#### fit models ####

# beta prior
x <- 0:200
y <- 0.85/(1 + 0.001 * x)
plot(x, y, type = "l")

val <- seq(0, 1, length.out = 50)
dens <- dexp(val, 1)
plot(val, dens, type = "l")

# Mv
mvEstGhMod <- brm(data = mvEstDat2, family = gaussian,
                  bf(PropEst ~ e0/(1 + beta * Litter.g.m2),
                     e0 ~ 1, 
                     beta ~ 1, 
                     nl = T),
                  prior <- c(prior(normal(0.85, 1), nlpar = "e0", lb = 0),
                             prior(exponential(1), nlpar = "beta", lb = 0),
                             prior(cauchy(0, 1), class = sigma)),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvEstGhMod)

# Ev
evEstGhMod <- update(mvEstGhMod, newdata = evEstDat2)
mod_check_fun(evEstGhMod)

# save
save(mvEstGhMod, file = "output/mv_establishment_model_2018_greenhouse_exp.rda")
save(evEstGhMod, file = "output/ev_establishment_model_2018_greenhouse_exp.rda")

# tables
write_csv(tidy(mvEstGhMod), "output/mv_establishment_model_2018_greenhouse_exp.csv")
write_csv(tidy(evEstGhMod), "output/ev_establishment_model_2018_greenhouse_exp.csv")


#### values for text ####

# posterior draws
mvEstGhDraws <- as_draws_df(mvEstGhMod) %>% as_tibble()
evEstGhDraws <- as_draws_df(evEstGhMod) %>% as_tibble()

# litter sensitivity
mvEstGhDraws %>%
  mean_hdci(b_beta_Intercept)

evEstGhDraws %>%
  mean_hdci(b_beta_Intercept)
