##### outputs ####

# Figure 2A (output/biomass_combined_figure_2019_density_exp.pdf)
# Table S1 (output/background_biomass_density_model_2019_dens_exp.csv)


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

# import data
bgBioD2Dat <- read_csv("data/bg_biomass_2019_density_exp.csv") 
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# figure settings
source("code/figure-prep/figure_settings.R")


##### edit data ####

# background biomass: combine bags (multiple if lots of biomass)
bgBioD2Dat2 <- bgBioD2Dat %>%
  group_by(site, plot, treatment) %>%
  summarise(biomass_bg = sum(biomass.g)) %>%
  ungroup()

# add 0 data to 1 plots
bgBioD2Dat3 <- tibble(site = rep(c("D1", "D2", "D3", "D4"), each = 2),
                      treatment = rep(c("water", "fungicide"), 4)) %>%
  mutate(plot = 1,
         biomass_bg = 0) %>%
  full_join(bgBioD2Dat2) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(background_species = case_when(background == "Mv seedling" ~ "*M. vimineum*",
                                        background == "Ev seedling" ~ "1st yr *E. virginicus*",
                                        background == "Ev adult" ~ "Adult *E. virginicus*") %>%
           as.factor(),
         treatment_fig = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"),
         fungicide = ifelse(treatment == "fungicide", 1, 0))

# divide by species
mvBgBioD2Dat <- filter(bgBioD2Dat3, background == "Mv seedling")
evSBgBioD2Dat <- filter(bgBioD2Dat3, background == "Ev seedling")
evABgBioD2Dat <- filter(bgBioD2Dat3, background == "Ev adult")


#### fit models ####

# initial visualization
ggplot(bgBioD2Dat3, aes(x = background_density, y = biomass_bg, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ background, scales = "free")

# distributions
ggplot(bgBioD2Dat3, aes(x = biomass_bg)) +
  geom_density() +
  facet_wrap(~ background, scales = "free")

bgBioD2Dat3 %>%
  filter(density_level != "none") %>%
  ggplot(aes(x = biomass_bg)) +
  geom_density() +
  geom_point(aes(x = biomass_bg, color = treatment), y = 0, shape = 108) +
  facet_wrap(~ density_level + background, scales = "free")

bgBioD2Dat3 %>%
  filter(density_level != "none") %>%
  ggplot(aes(x = log(biomass_bg))) +
  geom_density() +
  geom_point(aes(x = log(biomass_bg), color = treatment), y = 0, shape = 108) +
  facet_wrap(~ density_level + background, scales = "free")

# priors
bgBioD2Dat3 %>%
  filter(density_level == "low") %>%
  group_by(background) %>%
  summarize(avg_bio = mean(biomass_bg / background_density))

x <- seq(-1, 30, length.out = 100)
y <- dgamma(x, shape = 24, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")
y <- dgamma(x, shape = 10, scale = 1)
plot(x, y, type = "l")
y <- dgamma(x, shape = 2, scale = 1)
plot(x, y, type = "l")

# models
mvBgBioMod <- brm(data = mvBgBioD2Dat, family = gaussian,
                    bf(biomass_bg ~ (background_density * b0)/(1 + alpha * background_density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(24, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(mvBgBioMod)

evSBgBioMod <- brm(data = evSBgBioD2Dat, family = gaussian,
                   bf(biomass_bg ~ (background_density * b0)/(1 + alpha * background_density), 
                      b0 ~ 0 + treatment + (1|site), 
                      alpha ~ 0 + treatment, 
                      nl = T),
                   prior <- c(prior(gamma(2, 1), nlpar = "b0", lb = 0),
                              prior(exponential(0.5), nlpar = "alpha", lb = 0)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(evSBgBioMod)

evABgBioMod <- brm(data = evABgBioD2Dat, family = gaussian,
                   bf(biomass_bg ~ (background_density * b0)/(1 + alpha * background_density), 
                      b0 ~ 0 + treatment + (1|site), 
                      alpha ~ 0 + treatment, 
                      nl = T),
                   prior <- c(prior(gamma(10, 1), nlpar = "b0", lb = 0),
                              prior(exponential(0.5), nlpar = "alpha", lb = 0)),
                   iter = 4000, warmup = 2000, chains = 3, cores = 3,
                     control = list(adapt_delta = 0.99)) 
mod_check_fun(evABgBioMod)

# save models
save(mvBgBioMod, file = "output/mv_background_biomass_density_model_2019_density_exp.rda")
save(evSBgBioMod, file = "output/evS_background_biomass_density_model_2019_density_exp.rda")
save(evABgBioMod, file = "output/evA_background_biomass_density_model_2019_density_exp.rda")

# load models
load("output/mv_background_biomass_density_model_2019_density_exp.rda")
load("output/evS_background_biomass_density_model_2019_density_exp.rda")
load("output/evA_background_biomass_density_model_2019_density_exp.rda")


#### model table ####

bgBioTab <- tidy(mvBgBioMod) %>%
  mutate(plant_group = "M. vimineum") %>%
  full_join(tidy(evSBgBioMod) %>%
              mutate(plant_group = "1st yr E. virginicus")) %>%
  full_join(tidy(evABgBioMod) %>%
              mutate(plant_group = "adult E. virginicus")) %>%
  mutate(term = str_replace(term, "_treatment", " "),
         term = str_replace(term, "water", "control"),
         term = fct_recode(term, "b0 random intercept: site" = "sd__b0_(Intercept)"),
         term = fct_recode(term, "sigma" = "sd__Observation"),
         plant_group = fct_relevel(plant_group, "M. vimineum",
                                   "1st yr E. virginicus",
                                   "adult E. virginicus"),
         term = fct_relevel(term, "b0 control",
                            "b0 fungicide",
                            "alpha control",
                            "alpha fungicide",
                            "b0 random intercept: site"),
         across(c(estimate, std.error, conf.low, conf.high),
                ~ round_half_up(.x, 2))) %>%
  relocate(plant_group) %>%
  arrange(plant_group, term) %>%
  select(-c(effect, component, group))

# output tables
write_csv(bgBioTab, "output/background_biomass_density_model_2019_dens_exp.csv")


#### predicted values ####

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}

# prediction dataset template
bgBioPredDatTemplate <- bgBioD2Dat3 %>%
  group_by(treatment, treatment_fig, background, background_species) %>%
  summarize(min_dens = min(background_density),
            max_dens = max(background_density)) %>%
  ungroup() %>%
  mutate(background_density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(background_density) %>%
  mutate(site = "A") 

# apply to each plant group
mvBgBioPred <- bgBioPredDatTemplate %>%
  filter(background == "Mv seedling") %>%
  mutate(value = fitted(mvBgBioMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(mvBgBioMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(mvBgBioMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evSBgBioPred <- bgBioPredDatTemplate %>%
  filter(background == "Ev seedling") %>%
  mutate(value = fitted(evSBgBioMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evSBgBioMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evSBgBioMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evABgBioPred <- bgBioPredDatTemplate %>%
  filter(background == "Ev adult") %>%
  mutate(value = fitted(evABgBioMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evABgBioMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evABgBioMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

# combine plant groups
bgBioPred <- full_join(mvBgBioPred, evSBgBioPred) %>%
  full_join(evABgBioPred)


#### figure ####

# figure
bgBioFig <- ggplot(bgBioPred, aes(x = background_density, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment_fig), alpha = 0.3) +
  geom_line(aes(color = treatment_fig)) +
  geom_point(data = bgBioD2Dat3, 
             aes(y = biomass_bg, color = treatment_fig), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width), 
             alpha = 0.5, size = 0.75) +
  facet_wrap(~ background_species, scales = "free") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  labs(x = expression(paste("Density (", m^-2, ")", sep = "")),
       y = expression(paste("Plot biomass (g ", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(strip.text.x = element_markdown(),
        legend.box = "vertical")


save(bgBioFig, file = "output/background_biomass_density_figure_2019_dens_exp.rda")
