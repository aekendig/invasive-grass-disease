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
# figures
# output/seed_fungicide_figure_2019_density_exp.png


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
  yE = filter(dat_in, sp == "Ev" & background_density == 0) %>%
    pull(seeds) %>%
    mean()
  print(yE)
  
  yM = filter(dat_in, sp == "Mv" & background_density == 0) %>%
    pull(seeds) %>%
    mean()
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

# combine Ev adults and Mv
evASeedDat <- evSeedD2Dat2 %>% filter(age == "adult")
combSeedDat <- evASeedDat %>% full_join(mvSeedDat)

# split by background
combEvASeedDat <- combSeedDat %>% filter(background %in% c("none", "Ev adult"))
combEvSSeedDat <- combSeedDat %>% filter(background %in% c("none", "Ev seedling"))
combMvSeedDat <- combSeedDat %>% filter(background %in% c("none", "Mv seedling"))


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

ggplot(combEvASeedDat, aes(x = seeds)) +
  geom_density() +
  facet_wrap(sp ~ treatment, scales = "free")

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
                             prior(exponential(1), lb = 0, nlpar = "alpha")),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(combEvASeedMod)

combEvSSeedMod <- update(combEvASeedMod, newdata = combEvSSeedDat, cores = 3)
mod_check_fun(combEvSSeedMod)

combMvSeedMod <- update(combEvASeedMod, newdata = combMvSeedDat, cores = 3)
mod_check_fun(combMvSeedMod)

# save models
save(combEvASeedMod, file = "output/evA_background_seed_model_2019_density_exp.rda")
save(combEvSSeedMod, file = "output/evS_background_seed_model_2019_density_exp.rda")
save(combMvSeedMod, file = "output/mv_background_seed_model_2019_density_exp.rda")


#### fit Ev age model ####

# initial visualization
ggplot(evSeedDat, aes(x = age, y = seeds)) +
  geom_boxplot() +
  facet_wrap(~ treatment + yearf)

ggplot(evSeedDat, aes(x = seeds1)) +
  geom_density()

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
load("output/evA_background_seed_model_2019_density_exp.rda")
load("output/evS_background_seed_model_2019_density_exp.rda")
load("output/mv_background_seed_model_2019_density_exp.rda")
load("output/ev_seed_model_2019_density_exp.rda")

# save tables
write_csv(tidy(combEvASeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evA_background_seed_model_2019_density_exp.csv")
write_csv(tidy(combEvSSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evS_background_seed_model_2019_density_exp.csv")
write_csv(tidy(combMvSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/mv_background_seed_model_2019_density_exp.csv")
write_csv(tidy(evSeedMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/ev_seed_model_2019_density_exp.csv")

# prediction data
combEvASeedDraws <- combEvASeedDat %>%
  distinct(sp, fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(combEvASeedDat$background_density)) %>%
  add_epred_draws(combEvASeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

combEvSSeedDraws <- combEvSSeedDat %>%
  distinct(sp, fungicide, treatment) %>%
  expand_grid(background_density = 
                0:max(combEvSSeedDat$background_density)) %>%
  add_epred_draws(combEvSSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

combMvSeedDraws <- combMvSeedDat %>%
  distinct(sp, fungicide, treatment) %>%
  expand_grid(background_density = 0:max(combMvSeedDat$background_density)) %>%
  add_epred_draws(combMvSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

evSeedDraws <- evSeedD2Dat2 %>%
  distinct(age, fungicide, treatment) %>%
  mutate(yearf = NA) %>%
  add_epred_draws(evSeedMod, re_formula = ~0) %>% 
  ungroup() %>%
  select(.draw, fungicide, treatment, age, .epred) %>%
  pivot_wider(names_from = age, values_from = .epred) %>%
  mutate(ratio = seedling/adult,
         trt = fct_recode(treatment, "ambient" = "water"))

# check seed values
range(evSeedDraws$adult)
range(evSeedDraws$seedling)

# density figures
mvEvA_seed_fig <- filter(combEvASeedDraws, sp == "Mv") %>%
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
  filter(combEvASeedDraws, sp == "Ev") +
  labs(y = "*E. virginicus* seed yield")

mvEvS_seed_fig <- mvEvA_seed_fig %+%
  filter(combEvSSeedDraws, sp == "Mv") +
  labs(x = "First-year *E. virginicus* density")

evEvS_seed_fig <- mvEvS_seed_fig %+%
  filter(combEvSSeedDraws, sp == "Ev") +
  labs(y = "*E. virginicus* seed yield")

mvMv_seed_fig <- mvEvA_seed_fig %+%
  filter(combMvSeedDraws, sp == "Mv") +
  labs(x = "*M. vimineum* density")

evMv_seed_fig <- mvMv_seed_fig %+%
  filter(combMvSeedDraws, sp == "Ev") +
  labs(y = "*E. virginicus* seed yield")

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

# ratio figure
ratio_fung_fig <- ggplot(evSeedDraws, aes(x = ratio, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(y = "Disease treatment", x = "*E. virginicus* seed ratio") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

ggsave("output/seed_ratio_fungicide_figure_2018_2019_density_exp.png",
       ratio_fung_fig, width = 3, height = 3.2)  
