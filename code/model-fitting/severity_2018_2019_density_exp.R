##### outputs ####




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
sevD1Dat <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sevD2Dat <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# model function
mod_check_fun <- function(mod){
  
  print(prior_summary(mod))
  print(summary(mod))
  print(pp_check(mod, ndraws = 100))
  print(plot(mod))
  
}

# figure settings
source("code/figure-prep/figure_settings.R")


#### edit data ####

# severity
sevD1Dat2 <- sevD1Dat %>%
  filter(focal == 1 & bp_example == 0) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec, leaf_area.pix, lesion_area.pix) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity),
         year = "1",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background)) %>%
  filter(!is.na(severity))

sevD2Dat2 <- sevD2Dat %>%
  filter(focal == 1) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec, leaf_area.pix, lesion_area.pix) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  filter(month != "sep") %>% # too much data missing in sep
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity),
         year = "2",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background)) %>%
  filter(!is.na(severity))

# combine by plant group
evSSevDat <- filter(sevD1Dat2, sp == "Ev" & age == "seedling") %>%
  full_join(filter(sevD2Dat2, sp == "Ev" & age == "seedling"))

evASevDat <- filter(sevD1Dat2, sp == "Ev" & age == "adult") %>%
  full_join(filter(sevD2Dat2, sp == "Ev" & age == "adult"))

mvSevDat <- filter(sevD1Dat2, sp == "Mv") %>%
  full_join(filter(sevD2Dat2, sp == "Mv"))


#### start here ####

#### initial visualizations ####

# Mv density
sevD1Dat3 %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% 1:4) %>%
  ggplot(aes(x = Mv_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# Mv continuous
sevD1Dat3 %>%
  ggplot(aes(x = Mv_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = Mv_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# EvS density
sevD1Dat3 %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = EvS_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% c(1, 5:7)) %>%
  ggplot(aes(x = EvS_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# EvS continuous
sevD1Dat3 %>%
  ggplot(aes(x = EvS_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = EvS_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# EvA density
sevD1Dat3 %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = EvA_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

sevD2Dat3 %>%
  filter(plot %in% c(1, 8:10)) %>%
  ggplot(aes(x = EvA_dens, y = next_severity, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = mean_cl_boot, 
               position = position_dodge(1), width = 0) +
  stat_summary(geom = "point", fun = "mean",
               position = position_dodge(1), size = 2) +
  stat_summary(geom = "line", fun = "mean",
               position = position_dodge(1)) +
  facet_grid(plant_group ~ month, scales = "free_y")

# EvA continuous
sevD1Dat3 %>%
  ggplot(aes(x = EvA_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = EvA_sev_dens, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# edge
ggplot(sevD1Dat3, aes(x = edge_severity, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

ggplot(sevD2Dat3, aes(x = edge_severity, y = next_severity, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# severity over time
sevD1Dat3 %>%
  ggplot(aes(x = month, y = next_severity, color = plant_group)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = mean_cl_boot,
               position = position_dodge(0.2)) +
  stat_summary(geom = "point", size = 2, fun = mean,
               position = position_dodge(0.2)) +
  facet_wrap(~ treatment)

sevD2Dat3 %>%
  mutate(month = fct_relevel(month, "jun", "jul")) %>%
  ggplot(aes(x = month, y = next_severity, color = plant_group)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = mean_cl_boot,
               position = position_dodge(0.2)) +
  stat_summary(geom = "point", size = 2, fun = mean,
               position = position_dodge(0.2)) +
  facet_wrap(~ treatment)

# edge severity and dens_sev
sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = Mv_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvS_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD1Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvA_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = Mv_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvS_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

sevD2Dat3 %>%
  ggplot(aes(x = edge_severity, y = EvA_sev_dens)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month, scales = "free")

# logit transformation
sevD2Dat3 %>%
  ggplot(aes(x = next_severity, y = next_severity_l)) +
  geom_point()
# really skews severity changes

sevD2Dat3 %>%
  ggplot(aes(x = next_severity, y = next_severity_t)) +
  geom_point()
# linear change

# environmental variables
sevD1Dat3 %>%
  ggplot(aes(x = soil_moisture_c, y = next_severity_t, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)

sevD1Dat3 %>%
  ggplot(aes(x = canopy_cover_c, y = next_severity_t, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)

sevD2Dat3 %>%
  ggplot(aes(x = dew_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_grid(plant_group ~ month)


#### fit models ####

# remove missing data
sevD1Dat3_aug2 <- sevD1Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD1Dat3_jul2 <- sevD1Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity)) %>%
  mutate(next_severity_t = transform01(next_severity))

sevD2Dat3_aug2 <- sevD2Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jul2 <- sevD2Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jun2 <- sevD2Dat3_jun %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens)) %>%
  mutate(next_severity_t = transform01(next_severity))

# fit models
sevD1Mod_sev_dens_jul <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity) + (1|plotf),
                             data = sevD1Dat3_jul2, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD1Mod_sev_dens_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens) + (1|plotf),
                             data = sevD1Dat3_aug2, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_sev_dens_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity) + (1|plotf),
                         data = sevD2Dat3_aug2, family = "beta",
                         prior <- c(prior(normal(0, 1), class = "Intercept"),
                                    prior(normal(0, 1), class = "b")), # use default for others
                         iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_sev_dens_jul <- update(sevD2Mod_sev_dens_aug, newdata = sevD2Dat3_jul2)

sevD2Mod_sev_dens_jun <- update(sevD2Mod_sev_dens_aug, newdata = sevD2Dat3_jun2)

# check models
# mod_check_fun(sevD1Mod_sev_dens_jul)
# mod_check_fun(sevD1Mod_sev_dens_aug)
# mod_check_fun(sevD2Mod_sev_dens_jun)
# mod_check_fun(sevD2Mod_sev_dens_jul)
# mod_check_fun(sevD2Mod_sev_dens_aug)

# save models
save(sevD1Mod_sev_dens_aug, file = "output/focal_severity_model_aug_2018_dens_exp.rda")
save(sevD1Mod_sev_dens_jul, file = "output/focal_severity_model_jul_2018_dens_exp.rda")

save(sevD2Mod_sev_dens_aug, file = "output/focal_severity_model_aug_2019_dens_exp.rda")
save(sevD2Mod_sev_dens_jul, file = "output/focal_severity_model_jul_2019_dens_exp.rda")
save(sevD2Mod_sev_dens_jun, file = "output/focal_severity_model_jun_2019_dens_exp.rda")

# save corresponding raw data
write_csv(sevD2Dat3_aug2, "output/focal_severity_model_data_aug_2019_dens_exp.csv")
write_csv(sevD2Dat3_jul2, "output/focal_severity_model_data_jul_2019_dens_exp.csv")

# load models
load("output/focal_severity_model_aug_2018_dens_exp.rda")
load("output/focal_severity_model_jul_2018_dens_exp.rda")

load("output/focal_severity_model_aug_2019_dens_exp.rda")
load("output/focal_severity_model_jul_2019_dens_exp.rda")
load("output/focal_severity_model_jun_2019_dens_exp.rda")

# model summary tables
sevD1_sev_dens_sum <- tidy(sevD1Mod_sev_dens_jul) %>%
  filter(effect == "fixed") %>%
  select(term:conf.high) %>%
  rename(Covariate = term,
         Estimate_16 = estimate,
         Est.Error_16 = std.error,
         "l-95% CI_16" = conf.low,
         "u-95% CI_16" = conf.high) %>%
  full_join(tidy(sevD1Mod_sev_dens_aug) %>%
              filter(effect == "fixed") %>%
              select(term:conf.high) %>%
              rename(Covariate = term,
                     Estimate_20 = estimate,
                     Est.Error_20 = std.error,
                     "l-95% CI_20" = conf.low,
                     "u-95% CI_20" = conf.high)) %>%
  mutate(Covariate = str_replace_all(Covariate, "plant_groupEvS", "1st yr comp"),
         Covariate = str_replace_all(Covariate, "plant_groupEvA", "adult comp"),
         Covariate = str_replace_all(Covariate, "Mv_sev_dens", "src invader"),
         Covariate = str_replace_all(Covariate, "EvS_sev_dens", "src 1st yr"),
         Covariate = str_replace_all(Covariate, "EvA_sev_dens", "src adult"),
         Covariate = str_replace_all(Covariate, "edge_severity", "src surr"),
         Covariate = str_replace_all(Covariate, "fungicide", "fung"))
write_csv(sevD1_sev_dens_sum, "output/focal_severity_model_2018_dens_exp.csv")

sevD2_sev_dens_sum <- tidy(sevD2Mod_sev_dens_jun) %>%
  filter(effect == "fixed") %>%
  select(term:conf.high) %>%
  rename(Covariate = term,
         Estimate_8 = estimate,
         Est.Error_8 = std.error,
         "l-95% CI_8" = conf.low,
         "u-95% CI_8" = conf.high) %>%
  full_join(tidy(sevD2Mod_sev_dens_jul) %>%
              filter(effect == "fixed") %>%
              select(term:conf.high) %>%
              rename(Covariate = term,
                     Estimate_12 = estimate,
                     Est.Error_12 = std.error,
                     "l-95% CI_12" = conf.low,
                     "u-95% CI_12" = conf.high)) %>%
  full_join(tidy(sevD2Mod_sev_dens_aug) %>%
              filter(effect == "fixed") %>%
              select(term:conf.high) %>%
              rename(Covariate = term,
                     Estimate_16 = estimate,
                     Est.Error_16 = std.error,
                     "l-95% CI_16" = conf.low,
                     "u-95% CI_16" = conf.high)) %>%
  mutate(Covariate = str_replace_all(Covariate, "plant_groupEvS", "1st yr comp"),
         Covariate = str_replace_all(Covariate, "plant_groupEvA", "adult comp"),
         Covariate = str_replace_all(Covariate, "Mv_sev_dens", "src invader"),
         Covariate = str_replace_all(Covariate, "EvS_sev_dens", "src 1st yr"),
         Covariate = str_replace_all(Covariate, "EvA_sev_dens", "src adult"),
         Covariate = str_replace_all(Covariate, "edge_severity", "src surr"),
         Covariate = str_replace_all(Covariate, "fungicide", "fung"))
write_csv(sevD2_sev_dens_sum, "output/focal_severity_model_2019_dens_exp.csv")


#### transmission coefficients ####

# marginal trends  
Mv_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "Mv_sev_dens", regrid = "response")
EvS_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvS_sev_dens", regrid = "response")
EvA_sev_dens_D1_aug <- emtrends(sevD1Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvA_sev_dens", regrid = "response")

edge_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "edge_severity", regrid = "response")
Mv_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "Mv_sev_dens", regrid = "response")
EvS_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvS_sev_dens", regrid = "response")
EvA_sev_dens_D1_jul <- emtrends(sevD1Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvA_sev_dens", regrid = "response")

edge_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "edge_severity", regrid = "response")
Mv_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "Mv_sev_dens", regrid = "response")
EvS_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvS_sev_dens", regrid = "response")
EvA_sev_dens_aug <- emtrends(sevD2Mod_sev_dens_aug, c("plant_group", "fungicide"), var = "EvA_sev_dens", regrid = "response")

edge_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "edge_severity", regrid = "response")
Mv_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "Mv_sev_dens", regrid = "response")
EvS_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvS_sev_dens", regrid = "response")
EvA_sev_dens_jul <- emtrends(sevD2Mod_sev_dens_jul, c("plant_group", "fungicide"), var = "EvA_sev_dens", regrid = "response")

edge_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "edge_severity", regrid = "response")
Mv_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "Mv_sev_dens", regrid = "response")
EvS_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "EvS_sev_dens", regrid = "response")
EvA_sev_dens_jun <- emtrends(sevD2Mod_sev_dens_jun, c("plant_group", "fungicide"), var = "EvA_sev_dens", regrid = "response")

sev_dens_D1_coef <- as_tibble(summary(edge_sev_dens_D1_jul, point.est = mean)) %>% 
  mutate(source = "edge", weeks = 16) %>%
  rename(trend = edge_severity.trend) %>%
  full_join(as_tibble(summary(Mv_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "Invader (Mv)", weeks = 20) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "1st yr comp. (Ev)", weeks = 20) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_D1_aug, point.est = mean)) %>% 
                        mutate(source = "Adult comp. (Ev)", weeks = 20) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 16) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 16) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_D1_jul, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 16) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  mutate(treatment = if_else(fungicide == 0, "control", "fungicide"),
         plant_group = fct_recode(plant_group, "Invader (Mv)" = "Mv",
                                  "1st yr comp. (Ev)" = "EvS",
                                  "Adult comp. (Ev)" = "EvA"),
         source = fct_relevel(source, "Invader (Mv)") %>%
           fct_recode("Surrounding invader" = "edge"),
         sig = case_when(lower.HPD > 0 & upper.HPD > 0 ~ "yes",
                         lower.HPD < 0 & upper.HPD < 0 ~ "yes",
                         TRUE ~ "no"))

sev_dens_D2_coef <- as_tibble(summary(edge_sev_dens_aug, point.est = mean)) %>% 
  mutate(source = "edge", weeks = 16) %>%
  rename(trend = edge_severity.trend) %>%
  full_join(as_tibble(summary(Mv_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 16) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 16) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_aug, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 16) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(edge_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "edge", weeks = 12) %>%
              rename(trend = edge_severity.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 12) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 12) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_jul, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 12) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  full_join(as_tibble(summary(edge_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "edge", weeks = 8) %>%
              rename(trend = edge_severity.trend)) %>%
  full_join(as_tibble(summary(Mv_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "Invader (Mv)", weeks = 8) %>%
              rename(trend = Mv_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvS_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "1st yr comp. (Ev)", weeks = 8) %>%
              rename(trend = EvS_sev_dens.trend)) %>%
  full_join(as_tibble(summary(EvA_sev_dens_jun, point.est = mean)) %>% 
              mutate(source = "Adult comp. (Ev)", weeks = 8) %>%
              rename(trend = EvA_sev_dens.trend)) %>%
  mutate(treatment = if_else(fungicide == 0, "control", "fungicide"),
         plant_group = fct_recode(plant_group, "Invader (Mv)" = "Mv",
                                  "1st yr comp. (Ev)" = "EvS",
                                  "Adult comp. (Ev)" = "EvA"),
         source = fct_relevel(source, "Invader (Mv)") %>%
           fct_recode("Surrounding invader" = "edge"),
         sig = case_when(lower.HPD > 0 & upper.HPD > 0 ~ "yes",
                         lower.HPD < 0 & upper.HPD < 0 ~ "yes",
                         TRUE ~ "no"))

# save compiled trends
write_csv(sev_dens_D1_coef, "output/focal_severity_model_coefficients_2018_dens_exp.csv")
write_csv(sev_dens_D2_coef, "output/focal_severity_model_coefficients_2019_dens_exp.csv")

# load compiled trends
sev_dens_D1_coef <- read_csv("output/focal_severity_model_coefficients_2018_dens_exp.csv") %>%
  mutate(source = fct_relevel(source, "Invader (Mv)"),
         plant_group = fct_relevel(plant_group, "Invader (Mv)"))
sev_dens_D2_coef <- read_csv("output/focal_severity_model_coefficients_2019_dens_exp.csv") %>%
  mutate(source = fct_relevel(source, "Invader (Mv)"),
         plant_group = fct_relevel(plant_group, "Invader (Mv)"))


#### fungicide effects ####

# absolute difference

# hypotheses
mv_fung_hyp <- "inv_logit_scaled(Intercept + fungicide)*100 - inv_logit_scaled(Intercept)*100 = 0"
evS_fung_hyp <- "inv_logit_scaled(Intercept + fungicide + plant_groupEvS + plant_groupEvS:fungicide)*100 - inv_logit_scaled(Intercept + plant_groupEvS)*100  = 0"
evA_fung_hyp <- "inv_logit_scaled(Intercept + fungicide + plant_groupEvA + plant_groupEvA:fungicide)*100 - inv_logit_scaled(Intercept + plant_groupEvA)*100  = 0"

hypothesis(sevD2Mod_sev_dens_aug, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp)) # Mv
hypothesis(sevD2Mod_sev_dens_jul, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp)) # EvS and A
hypothesis(sevD2Mod_sev_dens_jun, c(mv_fung_hyp, evS_fung_hyp, evA_fung_hyp))


#### coefficient figure ####

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

col_pal <- c("#66a61e", "#d95f02", "#7570b3", "#e6ab02")

# figure
pdf("output/focal_severity_figure_2018_density_exp.pdf", width = 5.12, height = 3.94)
ggplot(sev_dens_D1_coef, aes(x = weeks, y = trend*100, color = source)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = lower.HPD*100, ymax = upper.HPD*100), 
                position = position_dodge(1), width = 0, size =0.35) +
  geom_point(aes(fill = source, shape = sig), position = position_dodge(1), size = 2) +
  facet_grid(treatment ~ plant_group, scales = "free") +
  scale_shape_manual(values = c(1, 19), guide = "none") +
  scale_color_manual(values = col_pal, name = "Disease source") +
  scale_fill_manual(values = col_pal, name = "Disease source") +
  scale_x_continuous(breaks = c(16, 20)) +
  labs(x = "Weeks post planting", y = "Percentage point difference in disease severity") +
  fig_theme
dev.off()

pdf("output/focal_severity_figure_2019_density_exp.pdf", width = 5.12, height = 3.54)
ggplot(sev_dens_D2_coef, aes(x = weeks, y = trend * 100, color = source)) +
  geom_hline(yintercept = 0, size = 0.25) +
  geom_errorbar(aes(ymin = lower.HPD * 100, ymax = upper.HPD * 100), 
                position = position_dodge(2.25), width = 0, size = 0.35) +
  geom_point(aes(fill = source, shape = sig), position = position_dodge(2.25), size = 2) +
  facet_grid(treatment ~ plant_group, scales = "free") +
  scale_shape_manual(values = c(1, 19), guide = "none") +
  scale_color_manual(values = col_pal, name = "Disease source") +
  scale_fill_manual(values = col_pal, name = "Disease source") +
  scale_x_continuous(breaks = c(8, 12, 16)) +
  labs(x = "Weeks post planting", y = "Percentage point difference in disease severity") +
  fig_theme
dev.off()


#### raw data figure ####

# figure settings
col_pal2 = c("black", "#238A8DFF")

# plot names
plot_names <- tibble(plot = 1:10,
                     plot_type = c("None",
                                   "Low invader",
                                   "Med. invader",
                                   "High invader",
                                   "Low 1st yr comp.",
                                   "Med 1st yr comp.",
                                   "High 1st yr comp.",
                                   "Low adult comp.",
                                   "Med adult comp.",
                                   "High adult comp."))

# edit raw data
sevD1Dat4 <- sevD1Dat2 %>%
  mutate(severity = severity * 100) %>%
  mutate(weeks = case_when(month == "jul" ~ 8,
                           month == "late_aug" ~ 16,
                           month == "sep" ~ 20),
         plant_group2 = case_when(plant_group == "Mv" ~ "Invader (Mv)",
                                  plant_group == "EvS" ~ "1st yr comp. (Ev)",
                                  plant_group == "EvA" ~ "Adult comp. (Ev)")) %>%
  full_join(edgeSevD1Dat2 %>%
              mutate(weeks = 8) %>%
              full_join(envD1Dat %>%
                          mutate(edge_severity = 100 * mv_inf_sep.prop,
                                 edge_severity = ifelse(edge_severity > 100, 100, edge_severity),
                                 weeks = 20)) %>%
              mutate(severity = edge_severity,
                     plant_group2 = "Surrounding invader")) %>%
  left_join(plot_names) %>%
  mutate(plant_group2 = fct_relevel(plant_group2, "Invader (Mv)"),
         treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control"),
         plot_type = fct_relevel(plot_type, "None",
                                 "Low invader",
                                 "Med. invader",
                                 "High invader",
                                 "Low 1st yr comp.",
                                 "Med 1st yr comp.",
                                 "High 1st yr comp.",
                                 "Low adult comp.",
                                 "Med adult comp.",
                                 "High adult comp."))

sevD2Dat4 <- sevD2Dat2 %>%
  mutate(severity = severity * 100) %>%
  full_join(edgeSevD2Dat2 %>%
              mutate(severity = edge_severity,
                     plant_group = "surr",
                     plant_group2 = "Surrounding invader")) %>%
  filter(month %in% c("jun", "jul", "early_aug", "late_aug")) %>%
  mutate(weeks = case_when(month == "jun" ~ 4,
                           month == "jul" ~ 8,
                           month == "early_aug" ~ 12,
                           month == "late_aug" ~ 16),
         plant_group2 = case_when(plant_group == "Mv" ~ "Invader (Mv)",
                                  plant_group == "EvS" ~ "1st yr comp. (Ev)",
                                  plant_group == "EvA" ~ "Adult comp. (Ev)",
                                  TRUE ~ plant_group2)) %>%
  
  left_join(plot_names) %>%
  mutate(plant_group2 = fct_relevel(plant_group2, "Invader (Mv)"),
         treatment = fct_recode(treatment, "control" = "water") %>%
           fct_relevel("control"),
         plot_type = fct_relevel(plot_type, "None",
                                 "Low invader",
                                 "Med. invader",
                                 "High invader",
                                 "Low 1st yr comp.",
                                 "Med 1st yr comp.",
                                 "High 1st yr comp.",
                                 "Low adult comp.",
                                 "Med adult comp.",
                                 "High adult comp."))

# figures
pdf("output/focal_severity_raw_data_figure_2018_density_exp.pdf", width = 6, height = 8)
ggplot(sevD1Dat4, aes(x = weeks, y = severity, color = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.5), 
             size = 0.5) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5)) +
  facet_grid(plot_type ~ plant_group2) +
  scale_color_manual(values = col_pal2, "Disease treatment") +
  labs(x = "Weeks post planting", y = "Disease severity (%)") +
  fig_theme +
  theme(legend.margin = margin(-0.3, 0, -0.3, 0, unit = "cm"),
        strip.text.y = element_text(angle = 0))
dev.off()

pdf("output/focal_severity_raw_data_figure_2019_density_exp.pdf", width = 6, height = 8)
ggplot(sevD2Dat4, aes(x = weeks, y = severity, color = treatment)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1,
                                             dodge.width = 0.5), 
             size = 0.5) +
  stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5)) +
  facet_grid(plot_type ~ plant_group2) +
  scale_color_manual(values = col_pal2, name = "Disease treatment") +
  labs(x = "Weeks post planting", y = "Disease severity (%)") +
  fig_theme +
  theme(legend.margin = margin(-0.3, 0, -0.3, 0, unit = "cm"),
        strip.text.y = element_text(angle = 0))
dev.off()


#### environmental models ####

# remove missing data
sevD1Dat3_aug_env <- sevD1Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(soil_moisture_c) & !is.na(canopy_cover_c)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD1Dat3_jul_env <- sevD1Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(soil_moisture_c) & !is.na(canopy_cover_c)) %>%
  mutate(next_severity_t = transform01(next_severity))

sevD2Dat3_aug_env <- sevD2Dat3_aug %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(dew_c)) %>%
  mutate(next_severity_t = transform01(next_severity))
sevD2Dat3_jul_env <- sevD2Dat3_jul %>% 
  filter(!is.na(Mv_sev_dens) & !is.na(EvA_sev_dens) & !is.na(EvS_sev_dens) & !is.na(edge_severity) & !is.na(dew_c)) %>%
  mutate(next_severity_t = transform01(next_severity))

# fit models
sevD1Mod_env_jul <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + soil_moisture_c + canopy_cover_c) + (1|plotf),
                             data = sevD1Dat3_jul_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD1Mod_env_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + soil_moisture_c + canopy_cover_c) + (1|plotf),
                             data = sevD1Dat3_aug_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_env_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + dew_c) + (1|plotf),
                             data = sevD2Dat3_aug_env, family = "beta",
                             prior <- c(prior(normal(0, 1), class = "Intercept"),
                                        prior(normal(0, 1), class = "b")), # use default for others
                             iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)

sevD2Mod_env_jul <- update(sevD2Mod_env_aug, newdata = sevD2Dat3_jul_env)

# check models
# mod_check_fun(sevD1Mod_env_aug)
# mod_check_fun(sevD1Mod_env_jul)
# mod_check_fun(sevD2Mod_env_aug)
# mod_check_fun(sevD2Mod_env_jul)

# save models
save(sevD1Mod_env_aug, file = "output/focal_severity_model_env_aug_2018_dens_exp.rda")
save(sevD1Mod_env_jul, file = "output/focal_severity_model_env_jul_2018_dens_exp.rda")
save(sevD2Mod_env_aug, file = "output/focal_severity_model_env_aug_2019_dens_exp.rda")
save(sevD2Mod_env_jul, file = "output/focal_severity_model_env_jul_2019_dens_exp.rda")

# load models
load("output/focal_severity_model_env_aug_2018_dens_exp.rda")
load("output/focal_severity_model_env_jul_2018_dens_exp.rda")
load("output/focal_severity_model_env_aug_2019_dens_exp.rda")
load("output/focal_severity_model_env_jul_2019_dens_exp.rda")

# marginal trends  
soil_env_D1_aug <- emtrends(sevD1Mod_env_aug, c("plant_group", "fungicide"), var = "soil_moisture_c", regrid = "response")
cano_env_D1_aug <- emtrends(sevD1Mod_env_aug, c("plant_group", "fungicide"), var = "canopy_cover_c", regrid = "response")
soil_env_D1_jul <- emtrends(sevD1Mod_env_jul, c("plant_group", "fungicide"), var = "soil_moisture_c", regrid = "response")
cano_env_D1_jul <- emtrends(sevD1Mod_env_jul, c("plant_group", "fungicide"), var = "canopy_cover_c", regrid = "response")

dew_env_D2_aug <- emtrends(sevD2Mod_env_aug, "plant_group", var = "dew_c", regrid = "response")
dew_env_D2_jul <- emtrends(sevD2Mod_env_jul, "plant_group", var = "dew_c", regrid = "response")

# combine
env_D1_coef <- as_tibble(summary(soil_env_D1_aug, point.est = mean)) %>% 
  mutate(Variable = "soil moisture", Weeks = 20) %>%
  rename(trend = soil_moisture_c.trend) %>%
  full_join(as_tibble(summary(cano_env_D1_aug, point.est = mean)) %>% 
              mutate(Variable = "canopy cover", Weeks = 20) %>%
              rename(trend = canopy_cover_c.trend)) %>%
  full_join(as_tibble(summary(soil_env_D1_jul, point.est = mean)) %>% 
              mutate(Variable = "soil moisture", Weeks = 16) %>%
              rename(trend = soil_moisture_c.trend)) %>%
  full_join(as_tibble(summary(cano_env_D1_jul, point.est = mean)) %>% 
              mutate(Variable = "canopy cover", Weeks = 16) %>%
              rename(trend = canopy_cover_c.trend)) %>%
  mutate(Treatment = if_else(fungicide == 0, "control", "fungicide"),
         Focal = fct_recode(plant_group, "invader" = "Mv",
                                  "first-year comp." = "EvS",
                                  "adult comp." = "EvA"),
         Estimate = trend * 100,
         l95 = lower.HPD * 100,
         u95 = upper.HPD * 100) %>%
  select(Weeks, Focal, Variable, Treatment, Estimate, l95, u95) %>%
  arrange(Weeks, Variable, Focal, Treatment)

env_D2_coef <- as_tibble(summary(dew_env_D2_aug, point.est = mean)) %>% 
  mutate(Weeks = 16) %>%
  full_join(as_tibble(summary(dew_env_D2_jul, point.est = mean)) %>% 
              mutate(Weeks = 12)) %>%
  mutate(Focal = fct_recode(plant_group, "invader" = "Mv",
                            "first-year comp." = "EvS",
                            "adult comp." = "EvA"),
         Estimate = dew_c.trend * 100,
         l95 = lower.HPD * 100,
         u95 = upper.HPD * 100) %>%
  select(Weeks, Focal, Estimate, l95, u95) %>%
  arrange(Weeks, Focal)

# save
write_csv(env_D1_coef, "output/focal_severity_model_env_coefficients_2018_density_exp.csv")
write_csv(env_D2_coef, "output/focal_severity_model_env_coefficients_2019_density_exp.csv")

# check sig relationships
sevD1Dat3_aug_env %>%
  filter(plant_group == "Mv" & treatment == "water") %>%
  ggplot(aes(x = canopy_cover_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

sevD2Dat3_jul_env %>%
  filter(plant_group == "EvS") %>%
  ggplot(aes(x = dew_c, y = next_severity_t)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)
# neither are driven by outliers

# check models with temp and humidity separate
sevD2Mod_temp_hum_aug <- brm(next_severity_t ~ plant_group * fungicide * (Mv_sev_dens + EvS_sev_dens + EvA_sev_dens + edge_severity + temp_avg_s + hum_avg_s) + (1|plotf),
                        data = sevD2Dat3_aug_env, family = "beta",
                        prior <- c(prior(normal(0, 1), class = "Intercept"),
                                   prior(normal(0, 1), class = "b")), # use default for others
                        iter = 6000, warmup = 1000, chains = 3, cores = 3, init_r = 0.1)
sevD2Mod_temp_hum_jul <- update(sevD2Mod_temp_hum_aug, newdata = sevD2Dat3_jul_env)

summary(sevD2Mod_temp_hum_aug)
summary(sevD2Mod_temp_hum_jul)
# neither have statistically significant effects


#### water effect - compare to edge ####

# summarize disease severity by plot
sevD2Dat_sum <- sevD2Dat2 %>%
  group_by(month, site, plot, treatment, sp) %>%
  summarize(severity = mean(severity)) %>%
  ungroup() %>%
  filter(sp == "Mv") %>%
  inner_join(edgeSevD2Dat2 %>%
               mutate(edge_severity = edge_severity / 100)) %>%
  pivot_longer(cols = c(severity, edge_severity),
               names_to = "plant_type",
               values_to = "severity") %>%
  mutate(plant_type = if_else(plant_type == "edge_severity", "surrounding", "inside") %>%
           fct_relevel("inside"),
         severity_t = transform01(severity),
         month = fct_relevel(month, "jun", "jul", "early_aug"),
         plotf = paste0(site, plot, substr(treatment, 1, 1)),
         treatment = fct_relevel(treatment, "water"))

# visualize
ggplot(sevD2Dat_sum, aes(x = month, y = severity_t, color = plant_type)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot", position = position_dodge(1)) +
  stat_summary(geom = "point", size = 2, fun = "mean", position = position_dodge(1)) +
  facet_wrap(~ treatment)

# model
ctrlD2Mod <- brm(severity_t ~ plant_type*treatment*month + (1|plotf),
                 data = sevD2Dat_sum, family = Beta,
                 prior <- c(prior(normal(0, 1), class = "Intercept"),
                            prior(normal(0, 1), class = "b")),
                 iter = 6000, warmup = 1000, chains = 3, cores = 3)

# mod_check_fun(ctrlD2Mod)

# save model
save(ctrlD2Mod, file = "output/edge_experiment_mv_model_2019_density_exp.rda")

# load model
load("output/edge_experiment_mv_model_2019_density_exp.rda")

# do experimental plants have different severity than edge?
jun_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding)*100 = inv_logit_scaled(Intercept)*100"
jun_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + plant_typesurrounding:treatmentfungicide)*100 = inv_logit_scaled(Intercept + treatmentfungicide)*100"
jul_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthjul + plant_typesurrounding:monthjul)*100 = inv_logit_scaled(Intercept + monthjul)*100"
jul_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthjul + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthjul + treatmentfungicide:monthjul + plant_typesurrounding:treatmentfungicide:monthjul)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthjul + treatmentfungicide:monthjul)*100"
eau_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthearly_aug + plant_typesurrounding:monthearly_aug)*100 = inv_logit_scaled(Intercept + monthearly_aug)*100"
eau_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthearly_aug + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthearly_aug + treatmentfungicide:monthearly_aug + plant_typesurrounding:treatmentfungicide:monthearly_aug)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthearly_aug + treatmentfungicide:monthearly_aug)*100"
lau_ctrl_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + monthlate_aug + plant_typesurrounding:monthlate_aug)*100 = inv_logit_scaled(Intercept + monthlate_aug)*100"
lau_fung_comp = "inv_logit_scaled(Intercept + plant_typesurrounding + treatmentfungicide + monthlate_aug + plant_typesurrounding:treatmentfungicide + plant_typesurrounding:monthlate_aug + treatmentfungicide:monthlate_aug + plant_typesurrounding:treatmentfungicide:monthlate_aug)*100 = inv_logit_scaled(Intercept + treatmentfungicide + monthlate_aug + treatmentfungicide:monthlate_aug)*100"

# estimates
ctrlD2hyps <- hypothesis(ctrlD2Mod, c(jun_ctrl_comp, jun_fung_comp, jul_ctrl_comp, jul_fung_comp,
                                      eau_ctrl_comp, eau_fung_comp, lau_ctrl_comp, lau_fung_comp))

ctrlD2hypsOut <- ctrlD2hyps[[1]] %>%
  as_tibble() %>%
  mutate(Weeks = rep(c(4, 8, 12, 16), each = 2),
         Treatment = rep(c("control", "fungicide"), 4))

write_csv(ctrlD2hypsOut, "output/edge_experiment_mv_2019_density_exp.csv")
