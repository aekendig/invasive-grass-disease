##### outputs ####

# Figure 2A

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(broom.mixed)
library(tidybayes)
library(ggtext)

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
                             TRUE ~ background_density)) %>%
  select(plot, treatment, background, density)

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
  mutate(biomass.g_m2 = case_when(background == "Mv seedling" ~ biomass_bg + biomass_foc_mv,
                                  background == "Ev seedling" ~ biomass_bg + biomass_foc_evS,
                                  background == "Ev adult" ~ biomass_bg + biomass_foc_evA),
         background_species = case_when(background == "Mv seedling" ~ "*M. vimineum*",
                                        background == "Ev seedling" ~ "1st yr *E. virginicus*",
                                        background == "Ev adult" ~ "Adult *E. virginicus*") %>%
           as.factor())

#### start here ####
# above dataset can be used for biomass vs. density
# add in data to do focal biomass vs. biomass


#### fit models ####

# initial visualization
ggplot(plotD2Dat, aes(x = density, y = biomass, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ age, scales = "free")

ggplot(plotD2Dat, aes(x = density, y = seeds, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, position = position_dodge(2)) +
  stat_summary(geom = "point", fun = "mean", size = 1.5, position = position_dodge(2)) +
  facet_wrap(~ age, scales = "free")

# priors
plotD2Dat %>%
  filter(plot == 1) %>%
  group_by(treatment, age) %>%
  summarise(b0 = mean(biomass/density),
            s0 = mean(seeds/density))

x <- seq(-1, 20, length.out = 100)
y <- dgamma(x, shape = 5, scale = 1) # note that this scale is 1/(stan scale)
plot(x, y, type = "l")

# models
evSBioDensMod <- brm(data = seedlingD2Dat, family = gaussian,
                    bf(biomass ~ (density * b0)/(1 + alpha * density), 
                       b0 ~ 0 + treatment + (1|site), 
                       alpha ~ 0 + treatment, 
                       nl = T),
                    prior <- c(prior(gamma(1, 1), nlpar = "b0", lb = 0),
                               prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                    iter = 6000, warmup = 1000, chains = 3,
                    control = list(adapt_delta = 0.99)) 

evABioDensMod <- brm(data = adultD2Dat, family = gaussian,
                     bf(biomass ~ (density * b0)/(1 + alpha * density), 
                        b0 ~ 0 + treatment + (1|site), 
                        alpha ~ 0 + treatment, 
                        nl = T),
                     prior <- c(prior(gamma(5, 1), nlpar = "b0", lb = 0),
                                prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma and sd
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999)) 

# mod_check_fun(evSBioDensMod)
# mod_check_fun(evABioDensMod)

# save models
save(evSBioDensMod, file = "output/evS_plot_biomass_density_model_2019_density_exp.rda")
save(evABioDensMod, file = "output/evA_plot_biomass_density_model_2019_density_exp.rda")

# load models
load("output/evS_plot_biomass_density_model_2019_density_exp.rda")
load("output/evA_plot_biomass_density_model_2019_density_exp.rda")

# output tables
write_csv(tidy(evSBioDensMod), "output/evS_plot_biomass_density_model_2019_dens_exp.csv")
write_csv(tidy(evABioDensMod), "output/evA_plot_biomass_density_model_2019_dens_exp.csv")


#### predicted values ####

# density gradient function
dens_fun <- function(min_dens, max_dens){
  
  density <- seq(min_dens, max_dens, length.out = 100)
  
  return(density)
}

# prediction dataset
plotPredDatTemplate <- plotD2Dat %>%
  group_by(treatment, fungicide, age) %>%
  summarize(min_dens = min(density),
            max_dens = max(density)) %>%
  ungroup() %>%
  mutate(density = pmap(list(min_dens, max_dens), dens_fun)) %>%
  unnest(density) %>%
  mutate(plotf = "A") 

evSBioPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "seedling") %>%
  mutate(value = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evSBioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evABioPredD2Dat <- plotPredDatTemplate %>%
  filter(age == "adult") %>%
  mutate(value = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Estimate"],
         lower = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q2.5"],
         upper = fitted(evABioDensMod, newdata = ., allow_new_levels = T)[, "Q97.5"])

evPredD2Dat <- evSSeedPredD2Dat %>%
  mutate(response = "seeds") %>%
  full_join(evASeedPredD2Dat %>%
              mutate(response = "seeds")) %>%
  full_join(evSBioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  full_join(evABioPredD2Dat %>%
              mutate(response = "biomass")) %>%
  mutate(plant_group = if_else(age == "adult", "Adult competitor (Ev)",
                               "1st yr competitor (Ev)") %>%
           fct_relevel("1st yr competitor (Ev)"))


#### figure ####


# labels
plot_labels <- c(biomass = "Biomass~(g~m^-2)",
                 seeds = "Seed~production~(m^-2)")

# figure
# pdf("output/ev_density_figure_2019_density_exp.pdf", width = 3.54, height = 3.54)
ggplot(plotBioD2Dat, aes(x = density, y = biomass.g_m2)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  # geom_line(aes(color = treatment)) +
  geom_point(aes(color = treatment), 
             position = position_jitterdodge(jitter.width = 0.1, 
                                             jitter.height = 0,
                                             dodge.width = dodge_width),
             size = 0.5) +
  facet_wrap(~ background_species, scales = "free") +
  scale_color_manual(values = col_pal, name = "Disease treatment") +
  # scale_fill_manual(values = col_pal, name = "Disease treatment") +
  labs(x = expression(paste("Planted density (", m^-2, ")", sep = "")),
       y = expression(paste("Plot biomass (g ", m^-2, ")", sep = ""))) +
  fig_theme +
  theme(strip.text.x = element_markdown())
# dev.off()


#### treatment effects ####

# treatment effects
b0_eff = "b0_treatmentfungicide - b0_treatmentcontrol = 0"
alpha_eff = "alpha_treatmentfungicide - alpha_treatmentcontrol = 0"

evSBioEff <- hypothesis(evSBioDensMod, c(b0_eff, alpha_eff))
evSSeedEff <- hypothesis(evSSeedDensMod, c(b0_eff, alpha_eff))
evABioEff <- hypothesis(evABioDensMod, c(b0_eff, alpha_eff))
evASeedEff <- hypothesis(evASeedDensMod, c(b0_eff, alpha_eff))

# combine treatment effects
evSBioEff[[1]] %>%
  mutate(response = "biomass",
         age = "seedling") %>%
  full_join(evSSeedEff[[1]] %>%
              mutate(response = "seeds",
                     age = "seedling")) %>%
  full_join(evABioEff[[1]] %>%
              mutate(response = "biomass",
                     age = "adult")) %>%
  full_join(evASeedEff[[1]] %>%
              mutate(response = "seeds",
                     age = "adult")) %>%
  mutate(parameter = rep(c("b0", "alpha"), 4)) %>%
  select(-c(Hypothesis, Evid.Ratio, Post.Prob)) %>%
  relocate(age, response, parameter)



