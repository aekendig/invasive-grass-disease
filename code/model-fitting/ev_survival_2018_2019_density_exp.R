##### outputs ####

# models
# output/evS_survival_model_2018_2019_density_exp.rda
# output/evA_survival_model_2018_2019_density_exp.rda
# tables
# output/evS_survival_model_2018_2019_density_exp.csv
# output/evA_survival_model_2018_2019_density_exp.csv
# parameter draws
# intermediate-data/pS_draws.csv
# intermediate-data/pP_draws.csv
# figures
# output/survival_fungicide_figure_2018_2019_density_exp.png

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
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")

# model functions
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# annual perennial survival 2018-2019
# must be alive in July of 2018
evSurvDat <- survD1Dat %>%
  filter(month %in% c("July", "April") & sp == "Ev") %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(site, plot, treatment, sp, age, ID, focal, month, survival) %>%
  pivot_wider(names_from = month, values_from = survival) %>%
  filter(July == 1 & !is.na(April)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotID = paste(site, plot, fungicide, sep = "_")) %>%
  select(-July) %>%
  rename(survival = April)
# includes non-focal

# look at data
count(evSurvDat, age, survival) # pretty high for adult

# sample sizes for paper
count(evSurvDat, age)

# divide data by age
evSSurvDat <- filter(evSurvDat, age == "seedling")
evASurvDat <- filter(evSurvDat, age == "adult")

#### fit perennial survival models ####

# initial visualization
ggplot(evSSurvDat, aes(x = treatment, y = survival)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

ggplot(evASurvDat, aes(x = treatment, y = survival)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", fun = "mean")

# prior means
evSurvDat %>%
  group_by(age) %>%
  summarize(prop = sum(survival)/n(),
            .groups = "drop") %>%
  mutate(logodds = car::logit(prop))

# fit models
evSSurvMod <- brm(data = evSSurvDat, family = bernoulli,
                  survival ~ fungicide + (1|site/plotID),
                  prior <- c(prior(normal(0.9, 1), class = "Intercept"),
                             prior(normal(0, 1), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evSSurvMod)

evASurvMod <- brm(data = evASurvDat, family = bernoulli,
                  survival ~ fungicide + (1|site/plotID),
                  prior <- c(prior(normal(3.8, 1), class = "Intercept"),
                             prior(normal(0, 1), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evASurvMod)

# save models
save(evSSurvMod, file = "output/evS_survival_fungicide_model_2018_2019_density_exp.rda")
save(evASurvMod, file = "output/evA_survival_fungicide_model_2018_2019_density_exp.rda")


#### tables and figures ####

# load
load("output/evS_survival_fungicide_model_2018_2019_density_exp.rda")
load("output/evA_survival_fungicide_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSSurvMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evS_survival_fungicide_model_2018_2019_density_exp.csv")
write_csv(tidy(evASurvMod, conf.method = "HPDinterval", rhat = T, ess = T), 
          "output/evA_survival_fungicide_model_2018_2019_density_exp.csv")

# prediction data
pred_dat_trt <- evSurvDat %>%
  distinct(fungicide, treatment) %>%
  mutate(trt = fct_recode(treatment, "ambient" = "water"))

# posterior draws
evSSurvDraws <- pred_dat_trt %>%
  add_epred_draws(evSSurvMod, re_formula = ~0) %>% 
  ungroup()
evASurvDraws <- pred_dat_trt %>%
  add_epred_draws(evASurvMod, re_formula = ~0) %>% 
  ungroup()

# save draws
evSSurvDraws2 <- evSSurvDraws %>%
  select(fungicide, .draw, .epred) %>%
  rename(value = .epred)

evASurvDraws2 <- evASurvDraws %>%
  select(fungicide, .draw, .epred) %>%
  rename(value = .epred)

write_csv(evSSurvDraws2, "intermediate-data/pS_draws.csv")
write_csv(evASurvDraws2, "intermediate-data/pP_draws.csv")

# figure
evS_surv_fig <- ggplot(evSSurvDraws, aes(x = .epred, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(x = "Disease treatment", y = "First-year *E. virginicus* survival") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.y = element_markdown())

evA_surv_fig <- ggplot(evASurvDraws, aes(x = .epred, y = trt)) +
  stat_slab(aes(fill = after_stat(level)), point_interval = mean_hdi, 
            .width = c(.66, .95, 1)) +
  stat_pointinterval(point_interval = mean_hdi, .width = c(.66, .95),
                     shape = 21, fill = "white", point_size = 1.5) +
  labs(x = "Disease treatment", y = "Adult *E. virginicus* survival") +
  scale_fill_manual(values = coral_pal, name = "HDI") +
  fig_theme +
  theme(axis.title.y = element_markdown())

# combine
surv_fung_fig <- evS_surv_fig + evA_surv_fig + 
  plot_annotation(tag_levels = "A") +
  plot_layout(nrow = 1, axes = "collect", guides = "collect") &
  theme(legend.position = "bottom")

ggsave("output/survival_fungicide_figure_2018_2019_density_exp.png",
       surv_fung_fig, width = 6, height = 3.2)