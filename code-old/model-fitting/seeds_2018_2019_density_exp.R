#### outputs #####

#### don't need 2018: seeds depend on density and treatments were messed up ####

# models
# output/evS_seed_model_2018_2019_density_exp.rda
# output/evA_seed_model_2018_2019_density_exp.rda
# output/mv_seed_model_2018_2019_density_exp.rda
# tables
# output/evS_seed_model_2018_2019_density_exp.csv
# output/evA_seed_model_2018_2019_density_exp.csv
# output/mv_seed_model_2018_2019_density_exp.csv
# figures
# output/seed_fungicide_figure_2018_2019_density_exp.rda
# output/seed_competition_figure_2018_2019_density_exp.rda


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
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
evSeedD1Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2018_density_exp.csv") 
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
tillerD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
plots <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

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

# 2018 survival
# make survival 1 if the plant produced seeds in summer
# remove NA's 
survD1Dat2 <- survD1Dat %>%
  filter(month == "September" & focal == 1) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival)) %>%
  select(-c(month, field_notes, seeds_produced, focal)) %>%
  filter(!is.na(survival))

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
              mutate(sp = "Ev"))

# Ev seeds 2018
evSeedD1Dat2 <- evSeedD1Dat %>%
  filter(focal == 1 & ID_unclear == 0) %>%
  group_by(site, plot, treatment, sp, age, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(survD1Dat2 %>%
              filter(sp == "Ev" & survival == 1) %>%
              select(-survival)) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = replace_na(seeds, 0),
         year = "1",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         seeds1 = seeds + 1)

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
         year = "2",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         seeds1 = seeds + 1)

# Mv seeds 2018
# quadrat was 0.49 x 0.25 m
mvSeedD1Dat2 <- mvSeedD1Dat %>% # none missing 
  left_join(tillerD1Dat %>%
              filter(sp == "Mv") %>%
              select(site, plot, treatment, sp, ID, tillers_jul)) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = seeds_per_stem * tillers_jul,
         year = "1",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         seeds1 = seeds + 1)

# Mv seeds 2019
mvSeedD2Dat2 <- mvSeedD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  select(site, plot, treatment, sp, ID, seeds) %>%
  full_join(focD2Dat %>%
              filter(sp == "Mv")) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  mutate(seeds = replace_na(seeds, 0),
         year = "2",
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         seeds1 = seeds + 1)

# combine by plant group
evSSeedDat <- evSeedD1Dat2 %>%
  filter(age == "seedling") %>%
  full_join(evSeedD2Dat2 %>%
              filter(age == "seedling"))

evASeedDat <- evSeedD1Dat2 %>%
  filter(age == "adult") %>%
  full_join(evSeedD2Dat2 %>%
              filter(age == "adult"))

mvSeedDat <- full_join(mvSeedD1Dat2, mvSeedD2Dat2)


#### fit models ####

# initial visualizations
ggplot(evSSeedDat, aes(x = background_density, y = seeds, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

ggplot(evSSeedDat, aes(x = log(seeds1))) +
  geom_density()

ggplot(evASeedDat, aes(x = background_density, y = seeds, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

ggplot(evASeedDat, aes(x = log(seeds1))) +
  geom_density()

ggplot(mvSeedDat, aes(x = background_density, y = seeds, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

ggplot(mvSeedDat, aes(x = log(seeds1))) +
  geom_density()


# fit models
evSSeedMod <- brm(data = evSSeedDat, family = lognormal,
                  seeds1 ~ fungicide + background_density:background + background_density:background:fungicide + year + 
                    (1|site/plot),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evSSeedMod)

mvSeedMod <- update(evSSeedMod, newdata = mvSeedDat,
                    control = list(adapt_delta = 0.999, max_treedepth = 15))
mod_check_fun(mvSeedMod)

evASeedMod <- brm(data = evASeedDat, family = lognormal,
                  seeds1 ~ fungicide + background_density:background + background_density:background:fungicide + year + 
                    (1|site),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evASeedMod)


# save models
save(evSSeedMod, file = "output/evS_seed_model_2018_2019_density_exp.rda")
save(evASeedMod, file = "output/evA_seed_model_2018_2019_density_exp.rda")
save(mvSeedMod, file = "output/mv_seed_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSSeedMod, conf.method = "HPDinterval"), 
          "output/evS_seed_model_2018_2019_density_exp.csv")
write_csv(tidy(evASeedMod, conf.method = "HPDinterval"), 
          "output/evA_seed_model_2018_2019_density_exp.csv")
write_csv(tidy(mvSeedMod, conf.method = "HPDinterval"), 
          "output/mv_seed_model_2018_2019_density_exp.csv")

# load
load("output/evS_seed_model_2018_2019_density_exp.rda")
load("output/evA_seed_model_2018_2019_density_exp.rda")
load("output/mv_seed_model_2018_2019_density_exp.rda")


#### fungicide effect without comp. figure ####

# posterior draws
evSSeedDraws <- as_draws_df(evSSeedMod)
mvSeedDraws <- as_draws_df(mvSeedMod)
evASeedDraws <- as_draws_df(evASeedMod)

# combine for fungicide effect without competition
SeedDraws <- tibble(sp = "E. virginicus",
                    age = "first-year",
                   int = evSSeedDraws$b_Intercept,
                   beta = evSSeedDraws$b_fungicide) %>%
  full_join(tibble(sp = "M. vimineum",
                   age = "first-year",
                   int = mvSeedDraws$b_Intercept,
                   beta = mvSeedDraws$b_fungicide)) %>%
  full_join(tibble(sp = "E. virginicus",
                   age = "adult",
                   int = evASeedDraws$b_Intercept,
                   beta = evASeedDraws$b_fungicide)) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         age = fct_relevel(age, "first-year"),
         sp_age = paste(age, sp) %>%
           fct_relevel("first-year M. vimineum",
                       "first-year E. virginicus"),
         resp_int = exp(int) - 1,
         resp_fung = exp(int + beta) - 1,
         resp_diff = resp_fung - resp_int,
         resp_perc = 100 * resp_diff / resp_int)

# fungicide effect without competition
seed_fung_fig <- ggplot(SeedDraws, aes(x = sp_age, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  scale_x_discrete(labels = c(~ atop(NA, atop(NA, textstyle(italic("M. vimineum")))),
                              ~ atop(NA, atop(textstyle("first-year"), 
                                              textstyle(italic("E. virginicus")))),
                              ~ atop(NA, atop(textstyle("adult"), 
                                              textstyle(italic("E. virginicus")))))) +
  labs(y = "Fungicide effect on seeds (log)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 2))

save(seed_fung_fig, file = "output/seed_fungicide_figure_2018_2019_density_exp.rda")


#### fungicide effect on competition ####

# edit
compSeedDraws <- tibble(fung = mvSeedDraws$`b_fungicide:background_density:backgroundEvadult`,
                        ctrl = mvSeedDraws$`b_background_density:backgroundEvadult`,
                        background = "adult *E. virginicus*") %>%
  full_join(tibble(fung = mvSeedDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = mvSeedDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = mvSeedDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = mvSeedDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*")) %>%
  mutate(sp_age = "*M. vimineum*") %>%
  full_join(tibble(fung = evSSeedDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evSSeedDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSSeedDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSSeedDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSSeedDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSSeedDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evASeedDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evASeedDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   sp_age = "adult *E. virginicus*")) %>%
  full_join(tibble(fung = evASeedDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evASeedDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   sp_age = "adult *E. virginicus*")) %>%
  full_join(tibble(fung = evASeedDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evASeedDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   sp_age = "adult *E. virginicus*")) %>%
  mutate(background = fct_relevel(background, "*M. vimineum*",
                                  "first-year *E. virginicus*"),
         sp_age = fct_relevel(sp_age, "*M. vimineum*",
                              "first-year *E. virginicus*"))

# figure
seed_comp_fig <- ggplot(compSeedDraws, aes(x = background, y = fung, color = sp_age)) +
  geom_hline(yintercept = 0, color = "grey") +
  stat_pointinterval(fatten_point = 5,
                     shape = 95,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1),
                     position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Focal group") +
  labs(x = "Background group",
       y = "Fungicide effect on\nseed response (log)") +
  fig_theme +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown())

save(seed_comp_fig, file = "output/seed_competition_figure_2018_2019_density_exp.rda")


#### values for text ####

# fungicide effects without neighbors
SeedDraws %>%
  group_by(sp, age) %>%
  mean_hdci(resp_int)

SeedDraws %>%
  group_by(sp, age) %>%
  mean_hdci(resp_diff)

SeedDraws %>%
  group_by(sp, age) %>%
  mean_hdci(resp_perc)

compSeedDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(fung)

compSeedDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(ctrl)

mvSeedDraws %>%
  mutate(seed_none = exp(b_Intercept) - 1,
         seed_Mv = exp(b_Intercept +
                         `b_background_density:backgroundMvseedling` * 30) - 1) %>%
  transmute(Mveffect = seed_Mv - seed_none) %>%
  mean_hdci()
