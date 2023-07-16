##### outputs ####

# Figure 2 (output/biomass_combined_figure_2019_density_exp.pdf)
# Table S1 (output/plot_biomass_density_model_2019_dens_exp.csv)
# Table S2 (output/focal_growth_interaction_coefficients_2019_density_exp.csv)
# Table S3 (output/focal_growth_interaction_coefficients_no_plot_1_high_EvA_2019_density_exp.csv)
# Table S4 (output/focal_growth_biomass_model_2019_dens_exp.csv)
# Table SX (output/focal_growth_biomass_model_no_plot_1_high_EvA_2019_dens_exp.csv)

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)
library(tidybayes)
library(ggtext)
library(janitor)
library(patchwork)

# import data
bgBioD2Dat <- read_csv("data/bg_biomass_2019_density_exp.csv") 
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, nsamples = 100))
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
  full_join(bgBioD2Dat2)

# plant group densities
plotDens <- plots %>%
  mutate(density = case_when(background == "Mv seedling" ~ background_density + 3,
                             background == "Ev seedling" ~ background_density + 3,
                             background == "Ev adult" ~ background_density + 1,
                             TRUE ~ background_density),
         density_level = fct_relevel(density_level, "none", "low", "medium", "high")) %>%
  select(plot, treatment, background, density, density_level)

# add focal biomass to background
# use average of others in plot if plant is missing biomass
plotBioD2Dat <- bgBioD2Dat3 %>%
  left_join(mvBioD2Dat %>% # add focal biomass
              group_by(site, plot, treatment) %>%
              mutate(biomass_weight_adj.g = mean(biomass_weight.g, na.rm = T)) %>%
              ungroup() %>%
              mutate(biomass_weight.g = case_when(is.na(biomass_weight.g) ~ biomass_weight_adj.g,
                                                  TRUE ~ biomass_weight.g)) %>%
              group_by(site, plot, treatment) %>%
              summarise(biomass_foc_mv = sum(biomass_weight.g)) %>%
              ungroup() %>%
              full_join(evBioD2Dat %>%
                          filter(ID %in% c("1", "2", "3")) %>%
                          group_by(site, plot, treatment) %>%
                          mutate(weight_adj = mean(weight, na.rm = T)) %>%
                          ungroup() %>%
                          mutate(weight = case_when(is.na(weight) ~ weight_adj,
                                                    TRUE ~ weight)) %>%
                          group_by(site, plot, treatment) %>%
                          summarise(biomass_foc_evS = sum(weight)) %>%
                          ungroup()) %>%
              full_join(evBioD2Dat %>%
                          filter(ID == "A") %>%
                          select(site, plot, treatment, weight) %>%
                          rename(biomass_foc_evA = weight))) %>%
  full_join(plotDens, relationship = "many-to-many") %>%
  mutate(biomass = case_when(background == "Mv seedling" ~ biomass_bg + biomass_foc_mv,
                             background == "Ev seedling" ~ biomass_bg + biomass_foc_evS,
                             background == "Ev adult" ~ biomass_bg + biomass_foc_evA),
         background_species = case_when(background == "Mv seedling" ~ "*M. vimineum*",
                                        background == "Ev seedling" ~ "1st yr *E. virginicus*",
                                        background == "Ev adult" ~ "Adult *E. virginicus*") %>%
           as.factor(),
         treatment_fig = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# divide by species
mvPlotBioD2Dat <- filter(plotBioD2Dat, background == "Mv seedling")
evSPlotBioD2Dat <- filter(plotBioD2Dat, background == "Ev seedling")
evAPlotBioD2Dat <- filter(plotBioD2Dat, background == "Ev adult")

# individual biomass
growthD2Dat <- mvBioD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  left_join(plotDens, relationship = "many-to-many") %>% 
  mutate(plant_growth = log(biomass_weight.g),
         age = ifelse(ID == "A", "adult", "seedling"),
         focal = paste(sp, age, sep = " "),
         foc = fct_recode(focal, m = "Mv seedling", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         background = str_replace(background, "_", " "),
         bg = fct_recode(background, m = "Mv seedling", a = "Ev adult", s = "Ev seedling") %>%
           fct_relevel("m"),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, str_sub(treatment, 1, 1)),
         focal_species = case_when(focal == "Mv seedling" ~ "*M. vimineum*",
                                   focal == "Ev seedling" ~ "1st yr *E. virginicus*",
                                   focal == "Ev adult" ~ "Adult *E. virginicus*"),
         intra = case_when(foc ==  bg ~ "yes",
                           str_detect(focal, "Ev") == T &  str_detect(background, "Ev") == T ~ "yes",
                           TRUE ~ "no")) %>%
  filter(!is.na(biomass_weight.g)) %>%
  left_join(plotBioD2Dat %>%
              select(site, treatment, treatment_fig, plot, background, background_species, biomass) %>%
              rename(plot_biomass = biomass)) %>%
  mutate(plot_biomass = if_else(as.character(focal) == as.character(background), 
                                plot_biomass - biomass_weight.g, plot_biomass))


#### fit biomass-density models ####

# initial visualization
ggplot(plotBioD2Dat, aes(x = density, y = biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ background, scales = "free")

# priors
plotBioD2Dat %>%
  filter(plot == 1) %>%
  group_by(treatment, background) %>%
  summarise(b0 = mean(biomass/density))

x <- seq(-1, 20, length.out = 100)
y <- dgamma(x, shape = 5, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

y <- dgamma(x, shape = 14, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# models
mvBioDensMod <- brm(data = mvPlotBioD2Dat, family = gaussian,
                    bf(biomass ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(14, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(mvBioDensMod)

evSBioDensMod <- brm(data = evSPlotBioD2Dat, family = gaussian,
                    bf(biomass ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(1, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 
mod_check_fun(evSBioDensMod)

evABioDensMod <- brm(data = evAPlotBioD2Dat, family = gaussian,
                     bf(biomass ~ (density * b0)/(1 + alpha * density), 
                        b0 ~ 0 + treatment + (1|site), 
                        alpha ~ 0 + treatment, 
                        nl = T),
                     prior <- c(prior(gamma(5, 1), nlpar = "b0", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 
mod_check_fun(evABioDensMod)

# save models
save(mvBioDensMod, file = "output/mv_plot_biomass_density_model_2019_density_exp.rda")
save(evSBioDensMod, file = "output/evS_plot_biomass_density_model_2019_density_exp.rda")
save(evABioDensMod, file = "output/evA_plot_biomass_density_model_2019_density_exp.rda")

# load models
load("output/mv_plot_biomass_density_model_2019_density_exp.rda")
load("output/evS_plot_biomass_density_model_2019_density_exp.rda")
load("output/evA_plot_biomass_density_model_2019_density_exp.rda")


#### biomass-density model table ####

bioDensTab <- tidy(mvBioDensMod) %>%
  mutate(plant_group = "M. vimineum") %>%
  full_join(tidy(evSBioDensMod) %>%
              mutate(plant_group = "1st yr E. virginicus")) %>%
  full_join(tidy(evABioDensMod) %>%
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
write_csv(bioDensTab, "output/plot_biomass_density_model_2019_dens_exp.csv")


#### biomass-density predicted values ####

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}

# prediction dataset template
plotPredDatTemplate <- plotBioD2Dat %>%
  group_by(treatment, treatment_fig, background, background_species) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(site = "A") 

# apply to each plant group
mvBioPredD2Dat <- plotPredDatTemplate %>%
  filter(background == "Mv seedling") %>%
  mutate(value = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(mvBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evSBioPredD2Dat <- plotPredDatTemplate %>%
  filter(background == "Ev seedling") %>%
  mutate(value = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evABioPredD2Dat <- plotPredDatTemplate %>%
  filter(background == "Ev adult") %>%
  mutate(value = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

# combine plant groups
bioPredD2Dat <- full_join(mvBioPredD2Dat, evSBioPredD2Dat) %>%
  full_join(evABioPredD2Dat)


#### fit individual growth models ####

# initial visualization
ggplot(growthD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_grid(focal ~ background, scales = "free")

ggplot(growthD2Dat, aes(plot_biomass, plant_growth, color = treatment)) +
  geom_point() +
  facet_grid(focal ~ background, scales = "free")

# remove plot 1 (no competitors)
growthD2Dat2 <- growthD2Dat %>%
  filter(plot > 1)

ggplot(growthD2Dat, aes(plot_biomass, plant_growth, color = treatment)) +
  geom_point() +
  facet_grid(foc ~ bg, scales = "free")

# fit models
growthD2Mod <- brm(data = growthD2Dat, family = gaussian,
                    plant_growth ~ foc * fungicide * (plot_biomass + plot_biomass:bg) + (1|plotf),
                    prior <- c(prior(normal(3, 1), class = "Intercept"),
                               prior(normal(0, 1), class = "b")), # use default for sigma
                    iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                    control = list(adapt_delta = 0.99999, max_treedepth = 15)) 
mod_check_fun(growthD2Mod)

# remove plot one (very low growth for some)
growthD2Mod2 <- update(growthD2Mod, newdata = growthD2Dat2)
mod_check_fun(growthD2Mod2)

# subset data to remove high EvA biomass
growthD2Dat2b <- growthD2Dat2 %>%
  filter(!(bg == "a" & plot_biomass > 150))

growthD2Mod2b <- update(growthD2Mod2, newdata = growthD2Dat2b)
mod_check_fun(growthD2Mod2b)

# save models and data
save(growthD2Mod, file = "output/focal_growth_biomass_model_2019_density_exp.rda")
save(growthD2Mod2, file = "output/focal_growth_biomass_model_no_plot_1_2019_density_exp.rda")
save(growthD2Mod2b, file = "output/focal_growth_biomass_model_no_plot_1_high_EvA_2019_density_exp.rda")

# load models
load("output/focal_growth_biomass_model_2019_density_exp.rda")
load("output/focal_growth_biomass_model_no_plot_1_2019_density_exp.rda")
load("output/focal_growth_biomass_model_no_plot_1_high_EvA_2019_density_exp.rda")


#### individual growth model table ####

# format tables
growthD2Tab <- tidy(growthD2Mod) %>%
  mutate(term = str_replace(term, "\\(Intercept\\)", "intercept"),
         term = str_replace(term, "foca", "adult Ev focal"),
         term = str_replace(term, "focs", "1st yr Ev focal"),
         term = str_replace(term, "plot_biomass", "biomass"),
         term = str_replace(term, "bga", "adult Ev bckgrd"),
         term = str_replace(term, "bgs", "1st yr Ev bckgrd"),
         term = str_replace(term, "sd__intercept", "random intercept: plot"),
         term = str_replace(term, "sd__Observation", "sigma")) %>%
  select(-c(effect, component, group))

growthD2Tab2b <- tidy(growthD2Mod2b) %>%
  mutate(term = str_replace(term, "\\(Intercept\\)", "intercept"),
         term = str_replace(term, "foca", "adult Ev focal"),
         term = str_replace(term, "focs", "1st yr Ev focal"),
         term = str_replace(term, "plot_biomass", "biomass"),
         term = str_replace(term, "bga", "adult Ev bckgrd"),
         term = str_replace(term, "bgs", "1st yr Ev bckgrd"),
         term = str_replace(term, "sd__intercept", "random intercept: plot"),
         term = str_replace(term, "sd__Observation", "sigma")) %>%
  select(-c(effect, component, group))

# output tables
write_csv(growthD2Tab, "output/focal_growth_biomass_model_2019_dens_exp.csv")
write_csv(growthD2Tab2b, "output/focal_growth_biomass_model_no_plot_1_high_EvA_2019_dens_exp.csv")


#### interaction coefficients (alphas) ####

# Mv background
mv_mv_ctrl_alpha = "plot_biomass = 0"
mv_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass = 0"
evS_mv_ctrl_alpha = "plot_biomass + focs:plot_biomass = 0"
evS_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + focs:plot_biomass + focs:fungicide:plot_biomass = 0"
evA_mv_ctrl_alpha = "plot_biomass + foca:plot_biomass = 0"
evA_mv_fung_alpha = "plot_biomass + fungicide:plot_biomass + foca:plot_biomass + foca:fungicide:plot_biomass = 0"

# EvS background
evS_evS_ctrl_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs = 0"
evS_evS_fung_alpha = "plot_biomass + focs:plot_biomass + plot_biomass:bgs + focs:plot_biomass:bgs + fungicide:plot_biomass + focs:fungicide:plot_biomass + fungicide:plot_biomass:bgs + focs:fungicide:plot_biomass:bgs = 0"
mv_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs = 0"
mv_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + fungicide:plot_biomass + fungicide:plot_biomass:bgs = 0"
evA_evS_ctrl_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs = 0"
evA_evS_fung_alpha = "plot_biomass +  plot_biomass:bgs + foca:plot_biomass + foca:plot_biomass:bgs + fungicide:plot_biomass +  fungicide:plot_biomass:bgs + foca:fungicide:plot_biomass + foca:fungicide:plot_biomass:bgs = 0"

# EvA background
evA_evA_ctrl_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga = 0"
evA_evA_fung_alpha = "plot_biomass + foca:plot_biomass + plot_biomass:bga + foca:plot_biomass:bga + fungicide:plot_biomass + foca:fungicide:plot_biomass + fungicide:plot_biomass:bga + foca:fungicide:plot_biomass:bga = 0"
mv_evA_ctrl_alpha = "plot_biomass +  plot_biomass:bga = 0"
mv_evA_fung_alpha = "plot_biomass +  plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga = 0"
evS_evA_ctrl_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga = 0"
evS_evA_fung_alpha = "plot_biomass + plot_biomass:bga + focs:plot_biomass + focs:plot_biomass:bga + fungicide:plot_biomass + fungicide:plot_biomass:bga + focs:fungicide:plot_biomass + focs:fungicide:plot_biomass:bga = 0"

# apply to each model
biomassD2alphas <- hypothesis(growthD2Mod, 
                              c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                                evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                                evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                                evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                                mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                                evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                                evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                                mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                                evS_evA_ctrl_alpha, evS_evA_fung_alpha))

biomassD2alphas2 <- hypothesis(growthD2Mod2, 
                              c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                                evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                                evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                                evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                                mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                                evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                                evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                                mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                                evS_evA_ctrl_alpha, evS_evA_fung_alpha))

biomassD2alphas2b <- hypothesis(growthD2Mod2b, 
                              c(mv_mv_ctrl_alpha, mv_mv_fung_alpha, 
                                evS_mv_ctrl_alpha, evS_mv_fung_alpha, 
                                evA_mv_ctrl_alpha, evA_mv_fung_alpha,
                                evS_evS_ctrl_alpha, evS_evS_fung_alpha,
                                mv_evS_ctrl_alpha, mv_evS_fung_alpha, 
                                evA_evS_ctrl_alpha, evA_evS_fung_alpha, 
                                evA_evA_ctrl_alpha, evA_evA_fung_alpha,
                                mv_evA_ctrl_alpha, mv_evA_fung_alpha, 
                                evS_evA_ctrl_alpha, evS_evA_fung_alpha))

# alpha processing
alpha_fun <- function(alphas) {
  
  out <- alphas %>%
    mutate(foc_bg_trt = c("m_m_ctrl", "m_m_fung", "s_m_ctrl", "s_m_fung", 
                          "a_m_ctrl", "a_m_fung", "s_s_ctrl", "s_s_fung",
                          "m_s_ctrl", "m_s_fung", "a_s_ctrl", "a_s_fung", 
                          "a_a_ctrl", "a_a_fung","m_a_ctrl", "m_a_fung", 
                          "s_a_ctrl", "s_a_fung")) %>%
    select(-Hypothesis) %>%
    rowwise() %>%
    mutate(foc = str_split(foc_bg_trt, "_")[[1]][1],
           bg = str_split(foc_bg_trt, "_")[[1]][2],
           trt = str_split(foc_bg_trt, "_")[[1]][3]) %>%
    ungroup() %>%
    left_join(growthD2Dat %>%
                distinct(foc, bg, focal_species, background_species)) %>%
    mutate(treatment_fig = fct_recode(trt, "control (water)" = "ctrl",
                                      "fungicide" = "fung"),
           sig = case_when((CI.Lower < 0 & CI.Upper < 0) | (CI.Lower > 0 & CI.Upper > 0) ~ "omits 0",
                           TRUE ~ "includes 0"),
           comp = as.character(sprintf("%.3f",round(Estimate, 3))))
  
  return(out)
  
}

alphaDat <- alpha_fun(biomassD2alphas[[1]])
alphaDat2 <- alpha_fun(biomassD2alphas2[[1]])
alphaDat2b <- alpha_fun(biomassD2alphas2b[[1]])

# save
write_csv(alphaDat, "output/focal_growth_interaction_coefficients_2019_density_exp.csv")
write_csv(alphaDat2b, "output/focal_growth_interaction_coefficients_no_plot_1_high_EvA_2019_density_exp.csv")


#### individual growth predicted values ####

# biomass gradient function
bio_fun <- function(f, b, trt){
  
  dat <- growthD2Dat %>% filter(foc == f & bg == b & treatment == trt)
  dat_out <- seq(min(dat$plot_biomass), max(dat$plot_biomass), length.out = 100)

  return(dat_out)
}

bio_fun2 <- function(f, b, trt){
  
  dat <- growthD2Dat2 %>% filter(foc == f & bg == b & treatment == trt)
  dat_out <- seq(min(dat$plot_biomass), max(dat$plot_biomass), length.out = 100)
  
  return(dat_out)
}

bio_fun2b <- function(f, b, trt){
  
  dat <- growthD2Dat2b %>% filter(foc == f & bg == b & treatment == trt)
  dat_out <- seq(min(dat$plot_biomass), max(dat$plot_biomass), length.out = 100)
  
  return(dat_out)
}

# predicted data
predD2Dat <- growthD2Dat %>%
  distinct(foc, bg, treatment, treatment_fig, fungicide, focal_species, background_species) %>%
  mutate(plotf = "A") %>%
  mutate(plot_biomass = pmap(list(foc, bg, treatment), bio_fun)) %>%
  unnest(plot_biomass) %>%
  mutate(plant_growth = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
  left_join(alphaDat %>%
              select(focal_species, background_species, treatment_fig, sig))

predD2Dat2 <- growthD2Dat2 %>%
  distinct(foc, bg, treatment, treatment_fig, fungicide, focal_species, background_species) %>%
  mutate(plotf = "A") %>%
  mutate(plot_biomass = pmap(list(foc, bg, treatment), bio_fun2)) %>%
  unnest(plot_biomass) %>%
  mutate(plant_growth = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod2, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
  left_join(alphaDat2 %>%
              select(focal_species, background_species, treatment_fig, sig))

predD2Dat2b <- growthD2Dat2b %>%
  distinct(foc, bg, treatment, treatment_fig, fungicide, focal_species, background_species) %>%
  mutate(plotf = "A") %>%
  mutate(plot_biomass = pmap(list(foc, bg, treatment), bio_fun2b)) %>%
  unnest(plot_biomass) %>%
  mutate(plant_growth = fitted(growthD2Mod2b, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(growthD2Mod2b, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(growthD2Mod2b, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
  left_join(alphaDat2b %>%
              select(focal_species, background_species, treatment_fig, sig))


#### figures ####

# plot biomass vs. density
bio_dens_fig <- ggplot(bioPredD2Dat, aes(x = density, y = value)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment_fig), alpha = 0.3) +
  geom_line(aes(color = treatment_fig)) +
  geom_point(data = plotBioD2Dat, 
             aes(y = biomass , color = treatment_fig, shape = density_level), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width), 
             alpha = 0.5, size = 0.75) +
  facet_wrap(~ background_species, scales = "free") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density level") +
  labs(x = expression(paste("Density (", m^-2, ")", sep = "")),
       y = expression(paste("Plot biomass (g ", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(strip.text.x = element_markdown(),
        legend.box = "vertical")

# plot individual vs. plot biomass
bio_bio_fig <- ggplot(predD2Dat, aes(x = plot_biomass, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment_fig), alpha = 0.3) +
  geom_line(aes(color = treatment_fig)) + 
  geom_point(data = growthD2Dat, 
             aes(color = treatment_fig, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat, sig == "omits 0"),
            x = -Inf, y = Inf, hjust = 0, vjust = 1,
            aes(label = paste("alpha", "==", comp, sep = ""),
                color = treatment_fig),
            parse = T, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal_species),
             cols = vars(background_species),
             scales = "free",
             switch = "both") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density level") +
  xlab(expression(paste("Plot biomass - focal biomass (g/", m^2, ")", sep = ""))) +
  ylab("Focal biomass (ln[g])") +
  fig_theme +
  theme(strip.text.x = element_markdown(),
        strip.text.y.left = element_markdown(),
        legend.box = "vertical")

bio_bio_fig2 <- ggplot(predD2Dat2, aes(x = plot_biomass, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment_fig), alpha = 0.3) +
  geom_line(aes(color = treatment_fig)) +
  geom_point(data = growthD2Dat2, 
             aes(color = treatment_fig, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat2, sig == "omits 0"),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
            aes(label = paste("alpha", "==", comp, sep = ""),
                color = treatment_fig),
            parse = T, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal_species),
             cols = vars(background_species),
             scales = "free",
             switch = "both") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density level") +
  xlab(expression(paste("Plot biomass - focal biomass (g/", m^2, ")", sep = ""))) +
  ylab("Focal biomass (ln[g])") +
  fig_theme +
  theme(strip.text.x = element_markdown(),
        strip.text.y.left = element_markdown(),
        legend.box = "vertical")

bio_bio_fig2b <- ggplot(predD2Dat2b, aes(x = plot_biomass, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment_fig), alpha = 0.3) +
  geom_line(aes(color = treatment_fig)) +
  geom_point(data = growthD2Dat2b, 
             aes(color = treatment_fig, shape = density_level), alpha = 0.5, size = 0.75) +
  geom_text(data = filter(alphaDat2b, sig == "omits 0"),
            x = -Inf, y = Inf, hjust = 0, vjust = 1,
            aes(label = paste("alpha", "==", comp, sep = ""),
                color = treatment_fig),
            parse = T, show.legend = F, size = textSize) +
  facet_grid(rows = vars(focal_species),
             cols = vars(background_species),
             scales = "free",
             switch = "both") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  scale_fill_manual(values = col_pal, name = "Disease treatment") +
  scale_shape(name = "Density level") +
  xlab(expression(paste("Plot biomass - focal biomass (g/", m^2, ")", sep = ""))) +
  ylab("Focal biomass (ln[g])") +
  fig_theme +
  theme(strip.text.x = element_markdown(),
        strip.text.y.left = element_markdown(),
        legend.box = "vertical")

# combine
bio_comb_fig <- bio_dens_fig + bio_bio_fig +
  plot_layout(ncol = 1, heights = c(0.3, 1), 
              guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  scale_shape_discrete(drop = F, name = "Density level") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        plot.tag = element_text(size = 10, face = "bold"))

pdf("output/biomass_combined_figure_2019_density_exp.pdf", width = 7.09, height = 8.5)
bio_comb_fig
dev.off()

pdf("output/focal_biomass_growth_no_plot_1_high_EvA_figure_2019_density_exp.pdf", width = 6.5, height = 6)
bio_bio_fig2b
dev.off()

pdf("output/focal_biomass_growth_no_plot_1_figure_2019_density_exp.pdf", width = 6.5, height = 6)
bio_bio_fig2
dev.off()


#### treatment effects ####

# treatment effects
b0_eff = "b0_treatmentfungicide - b0_treatmentwater = 0"
alpha_eff = "alpha_treatmentfungicide - alpha_treatmentwater = 0"

mvBioEff <- hypothesis(mvBioDensMod, c(b0_eff, alpha_eff))
evSBioEff <- hypothesis(evSBioDensMod, c(b0_eff, alpha_eff))
evABioEff <- hypothesis(evABioDensMod, c(b0_eff, alpha_eff))



