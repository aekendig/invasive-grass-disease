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
# parameter draws
# intermediate-data/yA_draws.csv
# intermediate-data/yP_draws.csv
# intermediate-data/f_draws.csv
# intermediate-data/alphaA_draws.csv
# intermediate-data/alphaP_draws.csv
# intermediate-data/gamma_draws.csv
# figures
# output/seed_fungicide_figure_2019_density_exp.png
# output/seed_alphas_ratio_figure_2018_2019_density_exp.png


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

# import data
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv")
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
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
  y0 = filter(dat_in, background_density == 0) %>%
    pull(seeds) %>%
    mean()
  print(y0)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(y = y0 / (1 + a * x))
  
  # plot
  print(ggplot(dat_in, aes(x = background_density, y = seeds)) +
          stat_summary(geom = "point", fun = "mean") +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = y)))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### format 2019 data ####

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
         seeds1 = seeds + 1,
         log_seeds = log(seeds1))

# Mv seeds 2019
mvSeedDat <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = replace_na(seeds, 0))

# isolate adult see data
evASeedDat <- evSeedD2Dat2 %>% filter(age == "adult")

# split by background
evAEvASeedDat <- evASeedDat %>% 
  filter(background %in% c("none", "Ev adult"))
evAEvSSeedDat <- evASeedDat %>% 
  filter(background %in% c("none", "Ev seedling"))
evAMvSeedDat <- evASeedDat %>% 
  filter(background %in% c("none", "Mv seedling"))

mvEvASeedDat <- mvSeedDat %>% 
  filter(background %in% c("none", "Ev adult"))
mvEvSSeedDat <- mvSeedDat %>% 
  filter(background %in% c("none", "Ev seedling"))
mvMvSeedDat <- mvSeedDat %>% 
  filter(background %in% c("none", "Mv seedling"))


#### format Ev 2018 data ####

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal)) %>%
  filter(!is.na(survival))

# Ev seeds 2018
evSeedD1Dat2 <- evSeedD1Dat %>%
  filter(focal == 1 & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev" & survival == 1) %>%
              select(-survival)) %>%
  mutate(seeds = replace_na(seeds, 0),
         seeds1 = seeds + 1,
         log_seeds = log(seeds1),
         fungicide = if_else(treatment == "fungicide", 1, 0),
         treatment = fct_relevel(treatment, "water"),
         plotID = paste(site, plot, fungicide, sep = "_"))

# combine Ev seeds
evSeedDat <- evSeedD2Dat2 %>%
  mutate(yearf = "2019") %>%
  full_join(evSeedD1Dat2 %>%
              mutate(yearf = "2018"))


#### fit Beverton-Holt models ####

# initial visualizations

ggplot(evAEvASeedDat, aes(x = seeds)) +
  geom_density() +
  facet_wrap(~ treatment, scales = "free")

ggplot(mvEvASeedDat, aes(x = seeds)) +
  geom_density() +
  facet_wrap(~ treatment, scales = "free")

# check prior distribution
val <- seq(0, 100, length.out = 50)
dens <- dexp(val, 5)
plot(val, dens, type = "l")

# fit EvA models
evAEvASeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0)
evAEvASeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0)

evAEvASeedMod <- brm(data = evAEvASeedDat, family = gaussian,
                  bf(seeds ~ s0/(1 + alpha * background_density),
                     s0 ~ fungicide + (1 | site), 
                     alpha ~ treatment + 0, 
                     nl = T),
                  prior <- c(prior(normal(81, 10), coef = 'Intercept', 
                                   nlpar = "s0"),
                             prior(normal(52, 10), coef = 'fungicide', 
                                   nlpar = "s0"),
                             prior(exponential(1), lb = 0, nlpar = "alpha")),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evAEvASeedMod)

evAEvSSeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0)
evAEvSSeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0.01)

evAEvSSeedMod <- update(evAEvASeedMod, newdata = evAEvSSeedDat, cores = 3)
mod_check_fun(evAEvSSeedMod)

evAMvSeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0.03)
evAMvSeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0.03)

evAMvSeedMod <- update(evAEvASeedMod, newdata = evAMvSeedDat, cores = 3)
mod_check_fun(evAMvSeedMod)

# fit Mv models
mvEvASeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0)
mvEvASeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0)

mvEvASeedMod <- brm(data = mvEvASeedDat, family = gaussian,
                     bf(seeds ~ s0/(1 + alpha * background_density),
                        s0 ~ fungicide + (1 | site/plotID), 
                        alpha ~ treatment + 0, 
                        nl = T),
                     prior <- c(prior(normal(1200, 100), coef = 'Intercept', 
                                      nlpar = "s0"),
                                prior(normal(-306, 100), coef = 'fungicide', 
                                      nlpar = "s0"),
                                prior(exponential(1), lb = 0, nlpar = "alpha")),
                     iter = 6000, warmup = 1000, chains = 3, cores = 3,
                     control = list(adapt_delta = 0.99))
mod_check_fun(mvEvASeedMod)

mvEvSSeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0)
mvEvSSeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0.01)

mvEvSSeedMod <- update(mvEvASeedMod, newdata = mvEvSSeedDat, cores = 3)
mod_check_fun(mvEvSSeedMod)

mvMvSeedDat %>% filter(treatment == "water") %>% bh_fun(a = 0.03)
mvMvSeedDat %>% filter(treatment == "fungicide") %>% bh_fun(a = 0.03)

mvMvSeedMod <- update(mvEvASeedMod, newdata = mvMvSeedDat, cores = 3)
mod_check_fun(mvMvSeedMod)

# save models
save(evAEvASeedMod, 
     file = "output/evA_focal_evA_background_seed_model_2019_density_exp.rda")
save(evAEvSSeedMod, 
     file = "output/evA_focal_evS_background_seed_model_2019_density_exp.rda")
save(evAMvSeedMod, 
     file = "output/evA_focal_mv_background_seed_model_2019_density_exp.rda")

save(mvEvASeedMod, 
     file = "output/mv_focal_evA_background_seed_model_2019_density_exp.rda")
save(mvEvSSeedMod, 
     file = "output/mv_focal_evS_background_seed_model_2019_density_exp.rda")
save(mvMvSeedMod, 
     file = "output/mv_focal_mv_background_seed_model_2019_density_exp.rda")


#### fit Ev age model ####

# initial visualization
ggplot(evSeedDat, aes(x = age, y = seeds)) +
  geom_boxplot() +
  facet_wrap(~ treatment + yearf)

ggplot(evSeedDat, aes(x = seeds1)) +
  geom_density() +
  facet_wrap(~ age + treatment)

ggplot(evSeedDat, aes(x = log_seeds)) +
  geom_density()  +
  facet_wrap(~ age + treatment)

evSeedDat %>%
  filter(age == "adult" & fungicide == 0 & yearf == "2018") %>%
  pull(log_seeds) %>%
  mean()

# fit model
evSeedMod <- brm(data = evSeedDat, family = Gamma(link = "log"),
                  seeds1 ~ age * fungicide + yearf + (1|site/plotID),
                  prior <- c(prior(normal(3.7, 1), class = "Intercept"),
                             prior(normal(0, 2), class = "b"),
                             prior(exponential(1), class = "sd")),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.99))
mod_check_fun(evSeedMod)

# save model
save(evSeedMod, file = "output/ev_seed_model_2019_density_exp.rda")


#### figures and tables ####

# load models
load("output/evA_focal_evA_background_seed_model_2019_density_exp.rda")
load("output/evA_focal_evS_background_seed_model_2019_density_exp.rda")
load("output/evA_focal_mv_background_seed_model_2019_density_exp.rda")
load("output/mv_focal_evA_background_seed_model_2019_density_exp.rda")
load("output/mv_focal_evS_background_seed_model_2019_density_exp.rda")
load("output/mv_focal_mv_background_seed_model_2019_density_exp.rda")
load("output/ev_seed_model_2019_density_exp.rda")

# save tables
write_csv(tidy(evAEvASeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evA_focal_evA_background_seed_model_2019_density_exp.csv")
write_csv(tidy(evAEvSSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evA_focal_evS_background_seed_model_2019_density_exp.csv")
write_csv(tidy(evAMvSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evA_focal_mv_background_seed_model_2019_density_exp.csv")
write_csv(tidy(mvEvASeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_focal_evA_background_seed_model_2019_density_exp.csv")
write_csv(tidy(mvEvSSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_focal_evS_background_seed_model_2019_density_exp.csv")
write_csv(tidy(mvMvSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_focal_mv_background_seed_model_2019_density_exp.csv")
write_csv(tidy(evSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/ev_seed_model_2019_density_exp.csv")

# predicted seed yield
evAEvASeedDraws <- evAEvASeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(evAEvASeedDat$background_density)) %>%
  add_epred_draws(evAEvASeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

evAEvSSeedDraws <- evAEvSSeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(evAEvSSeedDat$background_density)) %>%
  add_epred_draws(evAEvSSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

evAMvSeedDraws <- evAMvSeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 0:max(evAMvSeedDat$background_density)) %>%
  add_epred_draws(evAMvSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

mvEvASeedDraws <- mvEvASeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(mvEvASeedDat$background_density)) %>%
  add_epred_draws(mvEvASeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

mvEvSSeedDraws <- mvEvSSeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(mvEvSSeedDat$background_density)) %>%
  add_epred_draws(mvEvSSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

mvMvSeedDraws <- mvMvSeedDat %>%
  distinct(fungicide, treatment) %>%
  expand_grid(background_density = 0:max(mvMvSeedDat$background_density)) %>%
  add_epred_draws(mvMvSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

# perennial seed ratio
evSeedDraws <- evSeedD2Dat2 %>%
  distinct(age, fungicide, treatment) %>%
  mutate(yearf = NA) %>%
  add_epred_draws(evSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  select(.draw, fungicide, treatment, age, .epred) %>%
  pivot_wider(names_from = age, values_from = .epred) %>%
  mutate(ratio = seedling/adult,
         trt = fct_recode(treatment, "ambient" = "water"))

# competition coefficients
evAEvAAlphaDraws <- evAEvASeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

evAEvSAlphaDraws <- evAEvSSeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

evAMvAlphaDraws <- evAMvSeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

mvEvAAlphaDraws <- mvEvASeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

mvEvSAlphaDraws <- mvEvSSeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

mvMvAlphaDraws <- mvMvSeedMod %>%
  gather_draws(c(b_alpha_treatmentfungicide, 
                 b_alpha_treatmentwater)) %>%
  ungroup() %>%
  mutate(trt = if_else(str_detect(.variable, "water"), "ambient", 
                       "fungicide") %>%
           fct_relevel("fungicide"))

# check seed values
range(evSeedDraws$adult)
range(evSeedDraws$seedling)

# save draws
mvYDraws <- mvMvSeedDraws %>%
  filter(background_density == 0) %>%
  select(fungicide, .draw, .epred) %>%
  rename(value = .epred)

evYDraws <- evAMvSeedDraws %>%
  filter(background_density == 0) %>%
  select(fungicide, .draw, .epred) %>%
  rename(value = .epred)

evSeedDraws2 <- evSeedDraws %>%
  select(fungicide, .draw, ratio) %>%
  rename(value = ratio)

evAMvAlphaDraws2 <- evAMvAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

evAEvAAlphaDraws2 <- evAEvAAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

evAEvSAlphaDraws2 <- evAEvSAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

mvMvAlphaDraws2 <- mvMvAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

mvEvAAlphaDraws2 <- mvEvAAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

mvEvSAlphaDraws2 <- mvEvSAlphaDraws %>%
  mutate(fungicide = if_else(trt == "fungicide", 1, 0)) %>%
  select(fungicide, .draw, .value) %>%
  rename(value = .value)

write_csv(mvYDraws, "intermediate-data/yA_draws.csv")
write_csv(evYDraws, "intermediate-data/yP_draws.csv")
write_csv(evSeedDraws2, "intermediate-data/f_draws.csv")
write_csv(evAMvAlphaDraws2, "intermediate-data/alphaPA_draws.csv")
write_csv(evAEvAAlphaDraws2, "intermediate-data/alphaPP_draws.csv")
write_csv(evAEvSAlphaDraws2, "intermediate-data/alphaPS_draws.csv")
write_csv(mvMvAlphaDraws2, "intermediate-data/alphaAA_draws.csv")
write_csv(mvEvAAlphaDraws2, "intermediate-data/alphaAP_draws.csv")
write_csv(mvEvSAlphaDraws2, "intermediate-data/alphaAS_draws.csv")

# density figures
mvEvA_seed_fig <- mvEvASeedDraws %>%
  ggplot(aes(x = background_density, y = .epred)) +
  stat_lineribbon(aes(fill = trt, color = trt), point_interval = mean_hdi, 
                  .width = 0.95, alpha = 0.5) +
  scale_fill_manual(values = c(coral_pal[2], grey_pal[2]), 
                    name = "Disease treatment") +
  scale_color_manual(values = c(coral_pal[3], grey_pal[3]),
                     name = "Disease treatment") +
  labs(x = "Adult *E. virginicus* density", 
       y = "*M. vimineum* seed yield") +
  fig_theme +
  theme(axis.title = element_markdown())

evEvA_seed_fig <- mvEvA_seed_fig %+%
  evAEvASeedDraws +
  labs(y = "*E. virginicus* seed yield")

mvEvS_seed_fig <- mvEvA_seed_fig %+%
  mvEvSSeedDraws +
  labs(x = "First-year *E. virginicus* density")

evEvS_seed_fig <- mvEvS_seed_fig %+%
  evAEvSSeedDraws +
  labs(y = "*E. virginicus* seed yield")

mvMv_seed_fig <- mvEvA_seed_fig %+%
  mvMvSeedDraws +
  labs(x = "*M. vimineum* density")

evMv_seed_fig <- evAMvSeedDraws %>%
  ggplot(aes(x = background_density, y = .epred)) +
  stat_lineribbon(aes(fill = trt, color = trt), point_interval = mean_hdci, 
                  .width = 0.95, alpha = 0.5) +
  scale_fill_manual(values = c(coral_pal[2], grey_pal[2]), 
                    name = "Disease treatment") +
  scale_color_manual(values = c(coral_pal[3], grey_pal[3]),
                     name = "Disease treatment") +
  labs(x = "*M. vimineum* density", 
       y = "*E. virginicus* seed yield") +
  fig_theme +
  theme(axis.title = element_markdown())

# combine plots
seed_fung_fig <- (mvMv_seed_fig + evMv_seed_fig +
  plot_layout(axis_titles = "collect")) /
  (mvEvA_seed_fig + evEvA_seed_fig +
     plot_layout(axis_titles = "collect")) /
  (mvEvS_seed_fig + evEvS_seed_fig +
     plot_layout(axis_titles = "collect")) /
  plot_layout(guides = "collect")  + 
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") 

ggsave("output/seed_fungicide_figure_2019_density_exp.png",
       seed_fung_fig, width = 6, height = 8.2)

# alpha figures
evEvA_alpha_fig <- ggplot(evAEvAAlphaDraws, aes(x = .value, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(x = "*&alpha;~PP~*", 
       y = "Disease treatment") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

evEvS_alpha_fig <- evA_evA_alpha_fig %+% evAEvSAlphaDraws +
  labs(x = "*&alpha;~PS~*")

evMv_alpha_fig <- evA_evA_alpha_fig %+% evAMvAlphaDraws +
  labs(x = "*&alpha;~PA~*")

mvEvA_alpha_fig <- evA_evA_alpha_fig %+% mvEvAAlphaDraws +
  labs(x = "*&alpha;~AP~*")

mvEvS_alpha_fig <- evA_evA_alpha_fig %+% mvEvSAlphaDraws +
  labs(x = "*&alpha;~AS~*")

mvMv_alpha_fig <- evA_evA_alpha_fig %+% mvMvAlphaDraws +
  labs(x = "*&alpha;~AA~*")

# combine plots
alpha_fung_fig <- (mvMv_alpha_fig + evMv_alpha_fig +
    plot_layout(axes = "collect_y")) /
  (mvEvA_alpha_fig + evEvA_alpha_fig +
     plot_layout(axes = "collect_y")) /
  (mvEvS_alpha_fig + evEvS_alpha_fig +
     plot_layout(axes = "collect_y")) /
  plot_layout(guides = "collect")  + 
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom") 

ggsave("output/alpha_fungicide_figure_2019_density_exp.png",
       alpha_fung_fig, width = 6, height = 8.2)

# seed ratio figure
ratio_fung_fig <- ggplot(evSeedDraws, aes(x = ratio, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(y = "Disease treatment", x = "*E. virginicus* seed ratio") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ggsave("output/ev_seed_ratio_fungicide_figure_2018_2019_density_exp.png",
       ratio_fung_fig, width = 3, height = 3.2) 
