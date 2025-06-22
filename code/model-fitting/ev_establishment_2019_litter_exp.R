##### outputs ####

# model
# output/ev_establishment_model_2019_litter_exp.rda
# table
# output/ev_establishment_model_2018_litter_exp.csv
# posterior draws
# intermediate-data/eP_draws.csv
# intermediate-data/betaP_draws.csv
# figure
# output/establishment_figure_2018_2019_litter_exp.png


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(brms)
library(broom.mixed)
library(ggtext)
library(patchwork)

# import data
estL2Dat <- read_csv("data/both_germination_disease_jun_2019_litter_exp.csv")
plots <- read_csv("data/litter_weight_apr_2019_litter_exp.csv")

# load germination model and data
load("output/ev_germination_fungicide_model_2018_2019_density_exp.rda")
evGermDat <- read_csv("intermediate-data/ev_germination_2018_2019_density_exp.csv")

# load Mv figures
load("output/mv_establishment_figure_2018_litter_exp.rda")
load("output/mv_establishment_beta_figure_2018_litter_exp.rda")

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
    pull(prop_est_adj) %>%
    mean()
  print(V)
  
  # create data
  dat <- tibble(x = seq(xmin, xmax, length.out = 100)) %>%
    mutate(v = V / (1 + b * x))
  
  # plot
  print(ggplot(dat_in, aes(x = litter.g.m2, y = prop_est_adj)) +
          stat_summary(geom = "point", fun = "mean") +
          stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.1) +
          geom_line(data = dat, aes(x = x, y = v)))
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# edit plot data
# put the removed litter into the addition litter
# convert units
# remove unnecessary variables
plots2 <- plots %>%
  pivot_wider(names_from = "treatment", values_from = "litter_weight.lb") %>%
  mutate(addition = addition + removal,
         removal = 0) %>%
  pivot_longer(cols = c(addition, control, removal),
               names_to = "treatment",
               values_to = "litter_weight.lb") %>%
  mutate(litter_weight.g = litter_weight.lb * 453.592,
         litter.g.m2 = litter_weight.g / 4, # 4 because the plots are 2m^2
         site_block = paste(site, block, sep = "_")) %>%
  select(-date)

# Ev planting data
plant <- tibble(treatment = c("removal", "control", "addition"),
                ev_tot = c(50, 26, 26)) 

# average germination with fungicide
prop_germ <- tibble(fungicide = 0,
                    yearf = "2018",
                    age = "adult",
                    seeds_planted = round(mean(evGermDat$seeds_planted))) %>%
  add_epred_draws(evGermMod, re_formula = ~0) %>% 
  ungroup() %>%
  mutate(germ_frac = .epred / seeds_planted) %>%
  pull(germ_frac) %>%
  mean()

# June germination data
# none of the germinants were infected
evEstL2Dat <- estL2Dat %>%
  select(-c(date, flag_color, mv_germ, mv_infec)) %>%
  full_join(plots2) %>%
  full_join(plant) %>%
  mutate(prop_est = ev_germ / ev_tot,
         prop_est_adj = prop_est / prop_germ,
         treatment = fct_relevel(treatment, "removal", "control"))


#### fit model ####

# visualize
ggplot(evEstL2Dat, aes(x = litter.g.m2, y = prop_est_adj)) +
  geom_point()

ggplot(evEstL2Dat, aes(x = litter.g.m2, y = prop_est_adj)) +
  geom_point(position = position_jitter(width = 1)) +
  facet_wrap(~ site_block)

bh_fun(evEstL2Dat, b = 0.1)

# germination with litter
evEstL2Dat %>%
  filter(litter.g.m2 > 0 & ev_germ > 0)
# only one

# check prior distribution
val <- seq(0, 1, length.out = 50)
dens <- dexp(val, 1)
plot(val, dens, type = "l")

# fit model
evEstL2Mod <- brm(data = evEstL2Dat, family = gaussian,
                  bf(prop_est_adj ~ e0/(1 + beta * litter.g.m2),
                     e0 ~ 1 + (1 | site_block),
                     beta ~ 1,
                     nl = T),
                  prior <- c(prior(normal(0.16, 0.05), coef = 'Intercept', 
                                   nlpar = "e0"),
                             prior(exponential(1), lb = 0, nlpar = "beta")),
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))

mod_check_fun(evEstL2Mod)

# save model
save(evEstL2Mod, file = "output/ev_establishment_model_2019_litter_exp.rda")


#### figures and tables ####

# load
load("output/ev_establishment_model_2019_litter_exp.rda")

# table
write_csv(tidy(evEstL2Mod, conf.method = "HPDinterval", rhat = T, ess = T),
          "output/ev_establishment_model_2019_litter_exp.csv")

# prediction data
evEstL2Draws <- tibble (litter.g.m2 = seq(0, max(evEstL2Dat$litter.g.m2), 
                                          length.out = 100)) %>%
  add_epred_draws(evEstL2Mod, re_formula = ~0) %>% 
  ungroup()

# beta draws
evEstL2Beta <- spread_draws(evEstL2Mod, b_beta_Intercept) %>%
  mutate(trt = "ambient")

# save draws
evEstL2Draws2 <- evEstL2Draws %>%
  filter(litter.g.m2 == 0) %>%
  select(.draw, .epred) %>%
  rename(value = .epred)

evEstL2Beta2 <- evEstL2Beta %>%
  rename(value = b_beta_Intercept) %>%
  select(.draw, value)

write_csv(evEstL2Draws2, "intermediate-data/eP_draws.csv")
write_csv(evEstL2Beta2, "intermediate-data/betaP_draws.csv")

# figure
ev_est_fig <- ggplot(evEstL2Draws, aes(x = litter.g.m2, y = .epred)) +
  stat_lineribbon(color = grey_pal[3], fill = grey_pal[2],
                  point_interval = mean_hdci, 
                  .width = 0.95, alpha = 0.5) +
  labs(x = "Litter (g/m<sup>2</sup>)", 
       y = "*E. virginicus* establishment") +
  fig_theme +
  theme(axis.title = element_markdown())

# beta values
ev_beta_fig <- ggplot(evEstL2Beta, aes(x = b_beta_Intercept, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdci, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdci, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(x = "*E. virginicus* response to litter", y = "Disease treatment") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.x = element_markdown())

# combine
est_litter_fig <- mv_est_fig + mv_beta_fig + ev_est_fig + ev_beta_fig +
  plot_layout(ncol = 2, guides = "collect")  + 
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.box = "vertical") 

ggsave("output/establishment_figure_2018_2019_litter_exp.png",
       est_litter_fig, width = 6, height = 6.2)
