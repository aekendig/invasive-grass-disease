##### outputs ####

# models
# "output/evS_severity_model_2018_2019_density_exp.rda"
# "output/evA_severity_model_2018_2019_density_exp.rda"
# "output/mv_severity_model_2018_2019_density_exp.rda"
# tables
# "output/evS_severity_model_2018_2019_density_exp.csv"
# "output/evA_severity_model_2018_2019_density_exp.csv"
# "output/mv_severity_model_2018_2019_density_exp.csv"


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)
library(tidybayes)
library(janitor)
library(ggtext)
library(patchwork)

# import data
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv",
                     guess_max = 3000)
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# severity
sevD1Dat2 <- sevD1Dat %>%
  filter(focal == 1 & bp_example == 0) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec, 
         leaf_area.pix, lesion_area.pix) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(leaves_infec = if_else(leaves_infec == 0 & lesion_area.pix > 0, 1, # add an infected leaf if scan found lesions
                                leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = if_else(leaves_infec == 0 & is.na(severity), 0, # make severity zero if no leaves infected and no severity info
                            severity),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         month = as.factor(month),
         plotID = paste(site, plot, fungicide, sep = "_"),
         plantID = paste(plotID, sp, age, ID, sep = "_")) %>%
  filter(!is.na(severity)) %>%
  mutate(severity_t = transform01(severity),
         log_severity_t = log(severity_t))

sevD2Dat2 <- sevD2Dat %>%
  filter(focal == 1) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec, 
         leaf_area.pix, lesion_area.pix) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  filter(!(month %in% c("may", "sep"))) %>% # too much data missing
  mutate(leaves_infec = if_else(leaves_infec == 0 & lesion_area.pix > 0, 1, # add an infected leaf if scan found lesions
                                leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = if_else(leaves_infec == 0 & is.na(severity), 0, # make severity zero if no leaves infected and no severity info
                            severity),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         month = as.factor(month) %>% fct_relevel("jun", "jul"),
         plotID = paste(site, plot, fungicide, sep = "_"),
         plantID = paste(plotID, sp, age, ID, sep = "_")) %>%
  filter(!is.na(severity)) %>%
  mutate(severity_t = transform01(severity),
         log_severity_t = log(severity_t))

# transformations
ggplot(sevD1Dat2, aes(x = severity, y = severity_t)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

ggplot(sevD2Dat2, aes(x = severity, y = severity_t)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

sevD1Dat2 %>% 
  mutate(diff = severity_t - severity)  %>%
  ggplot(aes(x = severity, y = diff)) +
  geom_point()

# sample sizes
sampsD1 <- sevD1Dat2 %>% 
  group_by(sp, age, background, month) %>% 
  summarize(n = n_distinct(plotID))

sampsD2 <- sevD2Dat2 %>% 
  group_by(sp, age, background, month) %>% 
  summarize(n = n_distinct(plotID))

# separate by focal and background
evSEvSSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Ev seedling"))
evSEvASevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Ev adult"))
evSMvSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Mv seedling"))
evAEvSSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Ev seedling"))
evAEvASevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Ev adult"))
evAMvSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Mv seedling"))
mvEvSSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Ev seedling"))
mvEvASevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Ev adult"))
mvMvSevD1Dat <- sevD1Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Mv seedling"))

evSEvSSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Ev seedling"))
evSEvASevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Ev adult"))
evSMvSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "seedling" & 
           background %in% c("none", "Mv seedling"))
evAEvSSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Ev seedling"))
evAEvASevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Ev adult"))
evAMvSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Ev" & age == "adult" & 
           background %in% c("none", "Mv seedling"))
mvEvSSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Ev seedling"))
mvEvASevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Ev adult"))
mvMvSevD2Dat <- sevD2Dat2 %>% 
  filter(sp == "Mv" & background %in% c("none", "Mv seedling"))


#### fit models ####

# tried to fit beta models with severity_t and there were a lot of issues
# with Mv as the background

# initial visualizations
ggplot(sevD1Dat2, aes(x = severity)) +
  geom_density() +
  facet_grid(sp + age ~ month, scales = "free")

ggplot(sevD1Dat2, aes(x = log(severity))) +
  geom_density() +
  facet_grid(sp + age ~ month, scales = "free")

(init_fig <- ggplot(evSEvSSevD1Dat, aes(x = background_density, 
                                        y = severity, 
                                        color = treatment)) +
   geom_point() +
   geom_smooth(method = "glm") +
   facet_wrap(~ month))
# late aug incomplete

init_fig %+% evSEvASevD1Dat # late aug incomplete
init_fig %+% evSMvSevD1Dat # late aug and sep incomplete
init_fig %+% evAEvSSevD1Dat
init_fig %+% evAEvASevD1Dat
init_fig %+% evAMvSevD1Dat # late aug and sep incomplete
init_fig %+% mvEvSSevD1Dat # severity increases over months
init_fig %+% mvEvASevD1Dat
init_fig %+% mvMvSevD1Dat

init_fig %+% evSEvSSevD2Dat
init_fig %+% evSEvASevD2Dat
init_fig %+% evSMvSevD2Dat # all background highest late august
init_fig %+% evAEvSSevD2Dat
init_fig %+% evAEvASevD2Dat # late august incomplete
init_fig %+% evAMvSevD2Dat # all background highest late august
init_fig %+% mvEvSSevD2Dat
init_fig %+% mvEvASevD2Dat
init_fig %+% mvMvSevD2Dat # all background highest late august

# intercepts
sevD1Dat2 %>% 
  filter(month == "jul" & background == "none" & fungicide == 0) %>%
  group_by(sp, age) %>%
  summarize(mean = mean(log_severity_t))

sevD2Dat2 %>% 
  filter(month == "jun" & background == "none" & fungicide == 0) %>%
  group_by(sp, age) %>%
  summarize(mean = mean(log_severity_t))

# fit models
# initially used site as a random effect, but it caused issues with
# fitting models, especially for Mv. Site is captured in plotID.

# year 1 EvS
evSEvSSevD1Mod <- brm(data = evSEvSSevD1Dat, family = gaussian,
                      log_severity_t ~ fungicide * background_density * month + 
                        (1|plotID/plantID),
                  prior <- c(prior(normal(-3.5, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evSEvSSevD1Mod)
save(evSEvSSevD1Mod, file = "output/evSEvS_severity_model_2018_density_exp.rda")

evSEvASevD1Mod <- update(evSEvSSevD1Mod, newdata = evSEvASevD1Dat, cores = 3)
mod_check_fun(evSEvASevD1Mod)
save(evSEvASevD1Mod, file = "output/evSEvA_severity_model_2018_density_exp.rda")

evSMvSevD1Mod <- update(evSEvSSevD1Mod, newdata = evSMvSevD1Dat, cores = 3)
mod_check_fun(evSMvSevD1Mod)
save(evSMvSevD1Mod, file = "output/evSMv_severity_model_2018_density_exp.rda")

# year 1 EvA
evAEvSSevD1Mod <- update(evSEvSSevD1Mod, newdata = evAEvSSevD1Dat, cores = 3,
                         prior = c(prior(normal(-3, 10), class = "Intercept"),
                                   prior(normal(0, 10), class = "b")))
mod_check_fun(evAEvSSevD1Mod)
save(evAEvSSevD1Mod, file = "output/evAEvS_severity_model_2018_density_exp.rda")

evAEvASevD1Mod <- update(evAEvSSevD1Mod, newdata = evAEvASevD1Dat, cores = 3)
mod_check_fun(evAEvASevD1Mod)
save(evAEvASevD1Mod, file = "output/evAEvA_severity_model_2018_density_exp.rda")

evAMvSevD1Mod <- update(evAEvSSevD1Mod, newdata = evAMvSevD1Dat, cores = 3)
mod_check_fun(evAMvSevD1Mod)
save(evAMvSevD1Mod, file = "output/evAMv_severity_model_2018_density_exp.rda")

# year 1 Mv
mvEvSSevD1Mod <- update(evSEvSSevD1Mod, newdata = mvEvSSevD1Dat, cores = 3,
                        prior = c(prior(normal(-6.2, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")))
mod_check_fun(mvEvSSevD1Mod)
save(mvEvSSevD1Mod, file = "output/mvEvS_severity_model_2018_density_exp.rda")

mvEvASevD1Mod <- update(mvEvSSevD1Mod, newdata = mvEvASevD1Dat, cores = 3)
mod_check_fun(mvEvASevD1Mod)
save(mvEvASevD1Mod, file = "output/mvEvA_severity_model_2018_density_exp.rda")

mvMvSevD1Mod <- update(mvEvSSevD1Mod, newdata = mvMvSevD1Dat, cores = 3)
mod_check_fun(mvMvSevD1Mod)
save(mvMvSevD1Mod, file = "output/mvMv_severity_model_2018_density_exp.rda")

# year 2 EvS
evSEvSSevD2Mod <- update(evSEvSSevD1Mod, newdata = evSEvSSevD2Dat, cores = 3,
                        prior = c(prior(normal(-3.2, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")))
mod_check_fun(evSEvSSevD2Mod)
save(evSEvSSevD2Mod, file = "output/evSEvS_severity_model_2019_density_exp.rda")

evSEvASevD2Mod <- update(evSEvSSevD2Mod, newdata = evSEvASevD2Dat, cores = 3)
mod_check_fun(evSEvASevD2Mod)
save(evSEvASevD2Mod, file = "output/evSEvA_severity_model_2019_density_exp.rda")

evSMvSevD2Mod <- update(evSEvSSevD2Mod, newdata = evSMvSevD2Dat, cores = 3)
mod_check_fun(evSMvSevD2Mod)
save(evSMvSevD2Mod, file = "output/evSMv_severity_model_2019_density_exp.rda")

# year 2 EvA
evAEvSSevD2Mod <- update(evSEvSSevD2Mod, newdata = evAEvSSevD2Dat, cores = 3,
                        prior = c(prior(normal(-5.4, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")))
mod_check_fun(evAEvSSevD2Mod)
save(evAEvSSevD2Mod, file = "output/evAEvS_severity_model_2019_density_exp.rda")

evAEvASevD2Mod <- update(evAEvSSevD2Mod, newdata = evAEvASevD2Dat, cores = 3)
mod_check_fun(evAEvASevD2Mod)
save(evAEvASevD2Mod, file = "output/evAEvA_severity_model_2019_density_exp.rda")

evAMvSevD2Mod <- update(evAEvSSevD2Mod, newdata = evAMvSevD2Dat, cores = 3)
mod_check_fun(evAMvSevD2Mod)
save(evAMvSevD2Mod, file = "output/evAMv_severity_model_2019_density_exp.rda")

# year 2 Mv
mvEvSSevD2Mod <- update(evSEvSSevD2Mod, newdata = mvEvSSevD2Dat, cores = 3,
                        prior = c(prior(normal(-5.7, 10), class = "Intercept"),
                                  prior(normal(0, 10), class = "b")))
mod_check_fun(mvEvSSevD2Mod)
save(mvEvSSevD2Mod, file = "output/mvEvS_severity_model_2019_density_exp.rda")

mvEvASevD2Mod <- update(mvEvSSevD2Mod, newdata = mvEvASevD2Dat, cores = 3)
mod_check_fun(mvEvASevD2Mod)
save(mvEvASevD2Mod, file = "output/mvEvA_severity_model_2019_density_exp.rda")

mvMvSevD2Mod <- update(mvEvSSevD2Mod, newdata = mvMvSevD2Dat, cores = 3)
mod_check_fun(mvMvSevD2Mod)
save(mvMvSevD2Mod, file = "output/mvMv_severity_model_2019_density_exp.rda")


#### tables and figures ####

# load
load("output/evSEvS_severity_model_2018_density_exp.rda")
load("output/evSEvA_severity_model_2018_density_exp.rda")
load("output/evSMv_severity_model_2018_density_exp.rda")

load("output/evAEvS_severity_model_2018_density_exp.rda")
load("output/evAEvA_severity_model_2018_density_exp.rda")
load("output/evAMv_severity_model_2018_density_exp.rda")

load("output/mvEvS_severity_model_2018_density_exp.rda")
load("output/mvEvA_severity_model_2018_density_exp.rda")
load("output/mvMv_severity_model_2018_density_exp.rda")

load("output/evSEvS_severity_model_2019_density_exp.rda")
load("output/evSEvA_severity_model_2019_density_exp.rda")
load("output/evSMv_severity_model_2019_density_exp.rda")

load("output/evAEvS_severity_model_2019_density_exp.rda")
load("output/evAEvA_severity_model_2019_density_exp.rda")
load("output/evAMv_severity_model_2019_density_exp.rda")

load("output/mvEvS_severity_model_2019_density_exp.rda")
load("output/mvEvA_severity_model_2019_density_exp.rda")
load("output/mvMv_severity_model_2019_density_exp.rda")

# tables
write_csv(tidy(evSEvSSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSEvS_severity_model_2018_density_exp.csv")
write_csv(tidy(evSEvASevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSEvA_severity_model_2018_density_exp.csv")
write_csv(tidy(evSMvSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSMv_severity_model_2018_density_exp.csv")

write_csv(tidy(evAEvSSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAEvS_severity_model_2018_density_exp.csv")
write_csv(tidy(evAEvASevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAEvA_severity_model_2018_density_exp.csv")
write_csv(tidy(evAMvSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAMv_severity_model_2018_density_exp.csv")

write_csv(tidy(mvEvSSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvEvS_severity_model_2018_density_exp.csv")
write_csv(tidy(mvEvASevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvEvA_severity_model_2018_density_exp.csv")
write_csv(tidy(mvMvSevD1Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvMv_severity_model_2018_density_exp.csv")

write_csv(tidy(evSEvSSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSEvS_severity_model_2019_density_exp.csv")
write_csv(tidy(evSEvASevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSEvA_severity_model_2019_density_exp.csv")
write_csv(tidy(evSMvSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evSMv_severity_model_2019_density_exp.csv")

write_csv(tidy(evAEvSSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAEvS_severity_model_2019_density_exp.csv")
write_csv(tidy(evAEvASevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAEvA_severity_model_2019_density_exp.csv")
write_csv(tidy(evAMvSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evAMv_severity_model_2019_density_exp.csv")

write_csv(tidy(mvEvSSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvEvS_severity_model_2019_density_exp.csv")
write_csv(tidy(mvEvASevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvEvA_severity_model_2019_density_exp.csv")
write_csv(tidy(mvMvSevD2Mod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mvMv_severity_model_2019_density_exp.csv")

# posterior draws
evSEvSSevD1Draws <- evSEvSSevD1Dat %>% 
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSEvSSevD1Mod, re_formula = ~0) %>% 
  ungroup()
evAEvSSevD1Draws <- evAEvSSevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAEvSSevD1Mod, re_formula = ~0) %>% 
  ungroup()
mvEvSSevD1Draws <- mvEvSSevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvEvSSevD1Mod, re_formula = ~0) %>% 
  ungroup()

evSEvASevD1Draws <- evSEvASevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSEvASevD1Mod, re_formula = ~0) %>% 
  ungroup()
evAEvASevD1Draws <- evAEvASevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAEvASevD1Mod, re_formula = ~0) %>% 
  ungroup()
mvEvASevD1Draws <- mvEvASevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvEvASevD1Mod, re_formula = ~0) %>% 
  ungroup()

evSMvSevD1Draws <- evSMvSevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSMvSevD1Mod, re_formula = ~0) %>% 
  ungroup()
evAMvSevD1Draws <- evAMvSevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAMvSevD1Mod, re_formula = ~0) %>% 
  ungroup()
mvMvSevD1Draws <- mvMvSevD1Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvMvSevD1Mod, re_formula = ~0) %>% 
  ungroup()

evSEvSSevD2Draws <- evSEvSSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSEvSSevD2Mod, re_formula = ~0) %>% 
  ungroup()
evAEvSSevD2Draws <- evAEvSSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAEvSSevD2Mod, re_formula = ~0) %>% 
  ungroup()
mvEvSSevD2Draws <- mvEvSSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvEvSSevD2Mod, re_formula = ~0) %>% 
  ungroup()

evSEvASevD2Draws <- evSEvASevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSEvASevD2Mod, re_formula = ~0) %>% 
  ungroup()
evAEvASevD2Draws <- evAEvASevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAEvASevD2Mod, re_formula = ~0) %>% 
  ungroup()
mvEvASevD2Draws <- mvEvASevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvEvASevD2Mod, re_formula = ~0) %>% 
  ungroup()

evSMvSevD2Draws <- evSMvSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evSMvSevD2Mod, re_formula = ~0) %>% 
  ungroup()
evAMvSevD2Draws <- evAMvSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(evAMvSevD2Mod, re_formula = ~0) %>% 
  ungroup()
mvMvSevD2Draws <- mvMvSevD2Dat %>%
  group_by(month, fungicide, treatment) %>%
  expand(background_density = full_seq(background_density, period = 1)) %>% 
  ungroup() %>%
  add_epred_draws(mvMvSevD2Mod, re_formula = ~0) %>% 
  ungroup()

# combine
sevD1Draws <- evSEvSSevD1Draws %>% 
  mutate(focal = "first-year *E. virginicus*") %>% 
  full_join(evAEvSSevD1Draws %>% 
              mutate(focal = "adult *E. virginicus*")) %>%
  full_join(mvEvSSevD1Draws %>%
              mutate(focal = "*M. vimineum*")) %>%
  mutate(background = "first-year *E. virginicus*") %>%
  full_join(evSEvASevD1Draws %>% 
              mutate(focal = "first-year *E. virginicus*") %>%
              full_join(evAEvASevD1Draws %>% 
                          mutate(focal = "adult *E. virginicus*")) %>%
              full_join(mvEvASevD1Draws %>%
                          mutate(focal = "*M. vimineum*")) %>%
              mutate(background = "adult *E. virginicus*")) %>%
  full_join(evSMvSevD1Draws %>%
              mutate(focal = "first-year *E. virginicus*") %>%
              full_join(evAMvSevD1Draws %>% 
                          mutate(focal = "adult *E. virginicus*")) %>%
              full_join(mvMvSevD1Draws %>%
                          mutate(focal = "*M. vimineum*")) %>%
              mutate(background = "*M. vimineum*")) %>% 
  mutate(severity = exp(.epred),
         trt = fct_recode(treatment, "ambient" = "water"))

sevD2Draws <- evSEvSSevD2Draws %>% 
  mutate(focal = "first-year *E. virginicus*") %>% 
  full_join(evAEvSSevD2Draws %>% 
              mutate(focal = "adult *E. virginicus*")) %>%
  full_join(mvEvSSevD2Draws %>%
              mutate(focal = "*M. vimineum*")) %>%
  mutate(background = "first-year *E. virginicus*") %>%
  full_join(evSEvASevD2Draws %>% 
              mutate(focal = "first-year *E. virginicus*") %>%
              full_join(evAEvASevD2Draws %>% 
                          mutate(focal = "adult *E. virginicus*")) %>%
              full_join(mvEvASevD2Draws %>%
                          mutate(focal = "*M. vimineum*")) %>%
              mutate(background = "adult *E. virginicus*")) %>%
  full_join(evSMvSevD2Draws %>%
              mutate(focal = "first-year *E. virginicus*") %>%
              full_join(evAMvSevD2Draws %>% 
                          mutate(focal = "adult *E. virginicus*")) %>%
              full_join(mvMvSevD2Draws %>%
                          mutate(focal = "*M. vimineum*")) %>%
              mutate(background = "*M. vimineum*")) %>% 
  mutate(severity = exp(.epred),
         trt = fct_recode(treatment, "ambient" = "water"))

# figure
julSevD1Fig <- ggplot(filter(sevD1Draws, month == "jul"), 
                      aes(x = background_density, y = severity)) +
  stat_lineribbon(aes(fill = trt, color = trt), point_interval = median_hdci, 
                  .width = 0.95, alpha = 0.5) +
  facet_grid(focal ~ background, scales = "free") +
  scale_fill_manual(values = c(coral_pal[2], grey_pal[2]), 
                    name = "Disease treatment") +
  scale_color_manual(values = c(coral_pal[3], grey_pal[3]),
                     name = "Disease treatment") +
  labs(x = "Background plant density", 
       y = "Focal plant disease severity") +
  fig_theme +
  theme(strip.text = element_markdown())

lateAugSevD1Fig <- julSevD1Fig %+%
  filter(sevD1Draws, month == "late_aug")

sepSevD1Fig <- julSevD1Fig %+%
  filter(sevD1Draws, month == "sep")

junSevD2Fig <- julSevD1Fig %+%
  filter(sevD2Draws, month == "jun")

julSevD2Fig <- julSevD1Fig %+%
  filter(sevD2Draws, month == "jul")

earlyAugSevD2Fig <- julSevD1Fig %+%
  filter(sevD2Draws, month == "early_aug")

lateAugSevD2Fig <- julSevD1Fig %+%
  filter(sevD2Draws, month == "late_aug")

# save
ggsave("output/severity_jul_figure_2018_density_exp.png",
       julSevD1Fig, width = 6.5, height = 6.5)
ggsave("output/severity_late_aug_figure_2018_density_exp.png",
       lateAugSevD1Fig, width = 6.5, height = 6.5)
ggsave("output/severity_sep_figure_2018_density_exp.png",
       sepSevD1Fig, width = 6.5, height = 6.5)
ggsave("output/severity_jun_figure_2019_density_exp.png",
       junSevD2Fig, width = 6.5, height = 6.5)
ggsave("output/severity_jul_figure_2019_density_exp.png",
       julSevD2Fig, width = 6.5, height = 6.5)
ggsave("output/severity_early_aug_figure_2019_density_exp.png",
       earlyAugSevD2Fig, width = 6.5, height = 6.5)
ggsave("output/severity_late_aug_figure_2019_density_exp.png",
       lateAugSevD2Fig, width = 6.5, height = 6.5)


#### values for text ####

