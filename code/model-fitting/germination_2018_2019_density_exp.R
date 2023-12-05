##### outputs ####

# mv_germination_fungicide_model_2018_density_exp.rda
# mv_germination_fungicide_model_data_2018_density_exp.csv
# ev_germination_fungicide_model_2018_2019_density_exp.rda
# ev_germination_fungicide_model_data_2018_2019_density_exp.rda
# mv_germination_infection_model_2018_density_exp.csv
# mv_seed_infection_dark_fungicide_model_2018_density_exp.csv
# mv_germination_infection_figure_2018_density_exp.pdf


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(GGally)
library(cowplot)
library(car)
library(tidybayes)
library(broom.mixed)

# import data
mvGermDat <- read_csv("intermediate-data/mv_germination_disease_2018_density_exp.csv")
evGermDat <- read_csv("intermediate-data/ev_germination_2018_2019_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}


#### edit data ####

# Mv data
# calculate proportions and make study design columns
mvGermD1Dat <- mvGermDat %>%
  mutate(prop_germ = germination_final / seeds,
         prop_dark = seeds_dark / seeds,
         prop_light = seeds_light / seeds,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)))

# check
filter(mvGermD1Dat, prop_germ > 1 | prop_dark > 1 | prop_light > 1 ) %>%
  data.frame()

# trials per plot
mvGermD1Dat %>%
  count(plotf) %>%
  rename(trials = "n") %>%
  count(trials)

# Ev data
# calculate proportions and make study design columns
evGermDat2 <- evGermDat %>%
  mutate(prop_germ = germinants/seeds_planted,
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         yearf = as.factor(year))
  
# sample size
evGermDat2 %>%
  group_by(year, age) %>%
  count()


#### Mv models ####

# initial visualization
mvGermD1Dat %>%
  select(prop_dark, prop_light, prop_germ) %>%
  ggpairs()

ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# light/dark infection correlation
cor.test(~ prop_dark + prop_light, data = mvGermD1Dat) # not correlated

# fungicide model
mvGermD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ fungicide + (1|plotf), # can't converge with site/plot
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
mod_check_fun(mvGermD1Mod)
save(mvGermD1Mod, file = "output/mv_germination_fungicide_model_2018_density_exp.rda")

#### start here ####
# can import models from microstegium-bipolaris folder and resave
# note name change on next one

# model
mvGermFungD1Mod <- brm(data = mvGermD1Dat, family = binomial,
                   germination_final | trials(seeds) ~ prop_dark + prop_light + (1|plotf),
                   prior <- c(prior(normal(0, 10), class = "Intercept"),
                              prior(normal(0, 10), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, cores = 3)
save(mvGermFungD1Mod, file = "output/mv_germination_infection_model_2018_density_exp.rda")
# mod_check_fun(mvGermD1Mod)
# prop dark decreases germination
# prop light increases germination



mvPropDarkMod2 <- brm(data = mvGermD1Dat, family = binomial,
                     seeds_dark | trials(seeds) ~ fungicide + (1|plotf),
                     prior <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3, cores = 3)
# mod_check_fun(mvPropDarkMod2)
# fungicide decreases prop dark

mvPropLightMod2 <- brm(data = mvGermD1Dat, family = binomial,
                      seeds_light | trials(seeds) ~ fungicide + (1|plotf),
                      prior <- c(prior(normal(0, 10), class = "Intercept"),
                                 prior(normal(0, 10), class = "b")), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3)
# mod_check_fun(mvPropLightMod2)
# fungicide doesn't affect prop light

# save


save(mvPropDarkMod2, file = "output/mv_seed_infection_dark_fungicide_model_2018_density_exp.rda")
save(mvPropLightMod2, file = "output/mv_seed_infection_light_fungicide_model_2018_density_exp.rda")

# load
load("output/mv_germination_infection_model_2018_density_exp.rda")
load("output/mv_germination_fungicide_model_2018_density_exp.rda")
load("output/mv_seed_infection_dark_fungicide_model_2018_density_exp.rda")
load("output/mv_seed_infection_light_fungicide_model_2018_density_exp.rda")

# tables
write_csv(tidy(mvGermD1Mod), "output/mv_germination_infection_model_2018_density_exp.csv")
write_csv(tidy(mvPropDarkMod2), "output/mv_seed_infection_dark_fungicide_model_2018_density_exp.csv")

# save corresponding data
write_csv(mvGermD1Dat, "output/mv_germination_fungicide_model_data_2018_density_exp.csv")


#### Ev models ####

# initial visualization
ggplot(evGermDat2, aes(year, prop_germ, color = age)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(0.2)) +
  stat_summary(geom = "point", fun = "mean", position = position_dodge(0.2))

ggplot(evGermDat2, aes(treatment, prop_germ)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# model
evGermMod <- brm(data = evGermDat2, family = binomial,
                 germinants | trials(seeds_planted) ~ fungicide + (1|site) + (1|yearf),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, 
                 control = list(adapt_delta = 0.999, max_treedepth = 15))
# mod_check_fun(evGermMod)

# save
save(evGermMod, file = "output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")

# save corresponding data
write_csv(evGermDat2, "output/ev_germination_fungicide_model_data_2018_2019_density_exp.rda")


#### figure ####

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.box = "vertical",
        plot.title = element_text(size = 7, hjust = 0.5))

# simulated data
mvFungSim <- tibble(fungicide = c(0, 1)) %>%
  mutate(plotf = "A",
         treatment = c("control", "fungicide"),
         seeds = round(mean(mvGermD1Dat$seeds))) %>%
  mutate(prop_dark = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvPropDarkMod2, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

mvGermSim <- tibble(prop_dark = seq(min(mvGermD1Dat$prop_dark), max(mvGermD1Dat$prop_dark), length.out = 100)) %>%
  mutate(plotf = "A",
         seeds = round(mean(mvGermD1Dat$seeds)),
         prop_light = mean(mvGermD1Dat$prop_light)) %>%
  mutate(prop_germ = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Estimate"]/seeds,
         lower = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Q2.5"]/seeds,
         upper = fitted(mvGermD1Mod, newdata = ., allow_new_levels = T)[, "Q97.5"]/seeds)

# figures
mvFungFig <- ggplot(mvGermD1Dat, aes(treatment, prop_dark)) +
    # geom_point(alpha = 0.3, size = 0.7, position = position_jitter(width = 0.35)) +
  # geom_violin() +
  ggdist::stat_dots(dotsize = .4, 
                    binwidth = 0.01,
                    color = "#238A8DFF",
                    fill = "#238A8DFF",
                    alpha = 0.7,
                    side = "both") +
    geom_errorbar(data = mvFungSim, width = 0, size = 0.25, aes(ymin = lower, ymax = upper)) +
    geom_point(data = mvFungSim, size = 1.25) +
  fig_theme +
  theme(axis.title.x = element_blank()) +
  ylab("Infected fraction")

mvGermFig <- ggplot(mvGermD1Dat, aes(prop_dark, prop_germ)) +
  geom_point(alpha = 0.3, size = 0.7, color = "#238A8DFF") +
  geom_line(data = mvGermSim) +
  geom_ribbon(data = mvGermSim, aes(ymin = lower, ymax = upper), alpha = 0.4) +
  fig_theme +
  xlab("Infected fraction") +
  ylab("Germination fraction")

pdf("output/mv_germination_infection_figure_2018_density_exp.pdf", width = 3.54, height = 1.57)
plot_grid(mvFungFig, mvGermFig,
          nrow = 1,
          labels = LETTERS[1:2],
          label_size = 10)
dev.off()


#### values for text ####

# Mv fungicide effect on infection
hypothesis(mvPropDarkMod2, "plogis(Intercept + fungicide) - plogis(Intercept) = 0")

# Mv infection effect on germination
filter(mvGermSim, prop_dark == min(prop_dark))
filter(mvGermSim, prop_dark == max(prop_dark))

# Mv infection effect on germination
set.seed(184)
posterior_predict(mvGermD1Mod,
                  newdata = filter(mvGermD1Dat, 
                                   prop_dark %in% c(max(mvGermD1Dat$prop_dark), min(mvGermD1Dat$prop_dark))) %>%
                    select(prop_dark) %>%
                    unique() %>%
                    mutate(prop_light = 0,
                           seeds = 30,
                           plotf = "A"),
                  allow_new_levels = T) %>%
  as_tibble(.name_repair = ~ c("min_inf", "max_inf")) %>% # min_inf is first row
  mutate(min_inf_prop = min_inf / 30,
         max_inf_prop = max_inf / 30,
         inf_eff = 100 * (max_inf_prop - min_inf_prop) / min_inf_prop) %>%
  select(min_inf_prop, max_inf_prop, inf_eff) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  median_hdi(value)
