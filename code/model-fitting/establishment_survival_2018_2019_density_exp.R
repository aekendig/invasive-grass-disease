##### outputs ####

# models
# output/evS_establishment_model_2018_2019_density_exp.rda
# output/mv_establishment_model_2018_2019_density_exp.rda
# output/evS_survival_model_2018_2019_density_exp.rda
# output/evA_survival_model_2018_2019_density_exp.rda
# tables
# output/evS_establishment_model_2018_2019_density_exp.csv
# output/mv_establishment_model_2018_2019_density_exp.csv
# output/evS_survival_model_2018_2019_density_exp.csv
# output/evA_survival_model_2018_2019_density_exp.csv


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(broom.mixed)
library(ggtext)

# import data
survD1Dat <- read_csv("intermediate-data/all_processed_survival_2018_density_exp.csv")
survD2Dat <- read_csv("data/all_replacement_2019_density_exp.csv")
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
  left_join(plots) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         year = "1") %>%
  select(-c(month, field_notes, seeds_produced)) %>%
  filter(!is.na(survival))

# 2019 survival data needs list of all plants (only replacements recorded)
# merge ID lists with plot
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

# 2019 focal survival
# 2018 survival starts in June because plants were replaced through May 24
survD2Dat2 <- survD2Dat %>%
  filter(focal == 1 & replace_date > 20190531) %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(plantings = length(unique(replace_date)) + 1) %>%
  ungroup() %>%
  full_join(focD2Dat) %>%
  left_join(plots) %>%
  mutate(plantings = replace_na(plantings, 1),
         survival = ifelse(plantings > 1, 0, 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background),
         year = "2")

# annual perennial survival 2018-2019
evSurvD1Dat <- survD1Dat %>%
  filter(month == "April" & sp == "Ev") %>%
  left_join(plots) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background)) %>%
  select(-c(month, field_notes, seeds_produced)) %>%
  filter(!is.na(survival))
# includes non-focal

# divide data by focal and combine years
evSEstDat <- filter(survD1Dat2, sp == "Ev" & age == "seedling") %>%
  full_join(filter(survD2Dat2, sp == "Ev" & age == "seedling"))

mvEstDat <- filter(survD1Dat2, sp == "Mv" & age == "seedling") %>%
  full_join(filter(survD2Dat2, sp == "Mv" & age == "seedling"))

evSSurvDat <- filter(evSurvD1Dat, age == "seedling" & focal == 1)
evASurvDat <- filter(evSurvD1Dat, age == "adult" & focal == 1)


#### fit establishment models ####

# initial visualization
ggplot(evSEstDat, aes(x = background_density, y = survival, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

ggplot(mvEstDat, aes(x = background_density, y = survival, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

# fit models
evSEstMod <- brm(data = evSEstDat, family = bernoulli,
                 survival ~ fungicide + background_density:background + background_density:background:fungicide + year + 
                   (1|site/plot),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.99))
mod_check_fun(evSEstMod)

mvEstMod <- update(evSEstMod, newdata = mvEstDat)
mod_check_fun(mvEstMod)

# save models
save(evSEstMod, file = "output/evS_establishment_model_2018_2019_density_exp.rda")
save(mvEstMod, file = "output/mv_establishment_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSEstMod), "output/evS_establishment_model_2018_2019_density_exp.csv")
write_csv(tidy(mvEstMod), "output/mv_establishment_model_2018_2019_density_exp.csv")

# load
load("output/evS_establishment_model_2018_2019_density_exp.rda")
load("output/mv_establishment_model_2018_2019_density_exp.rda")


#### fit perennial survival models ####

# initial visualization
ggplot(evSSurvDat, aes(x = background_density, y = survival, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~ treatment)

ggplot(evASurvDat, aes(x = background_density, y = survival, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~ treatment)

# fit models
evSSurvMod <- brm(data = evSSurvDat, family = bernoulli,
                  survival ~ fungicide + background_density:background + background_density:background:fungicide + 
                    (1|site/plot),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evSSurvMod)

evASurvMod <- brm(data = evASurvDat, family = bernoulli,
                  survival ~ fungicide + background_density:background + background_density:background:fungicide + 
                    (1|site),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evASurvMod)

# save models
save(evSSurvMod, file = "output/evS_survival_model_2018_2019_density_exp.rda")
save(evASurvMod, file = "output/evA_survival_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSSurvMod), "output/evS_survival_model_2018_2019_density_exp.csv")
write_csv(tidy(evASurvMod), "output/evA_survival_model_2018_2019_density_exp.csv")

# load
load("output/evS_survival_model_2018_2019_density_exp.rda")
load("output/evA_survival_model_2018_2019_density_exp.rda")


#### fungicide effect without comp. figure ####

# posterior draws
evSEstDraws <- as_draws_df(evSEstMod)
mvEstDraws <- as_draws_df(mvEstMod)
evSSurvDraws <- as_draws_df(evSSurvMod)
evASurvDraws <- as_draws_df(evASurvMod)

# combine for fungicide effect without competition
estDraws <- tibble(sp = "E. virginicus",
                   int = evSEstDraws$b_Intercept,
                   beta = evSEstDraws$b_fungicide) %>%
  full_join(tibble(sp = "M. vimineum",
                   int = mvEstDraws$b_Intercept,
                   beta = mvEstDraws$b_fungicide)) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         prob_int = exp(int) / (1 + exp(int)),
         prob_fung = exp(int + beta) / (1 + exp(int + beta)),
         prob_change = 100 * (prob_fung - prob_int) / prob_int)

survDraws <- tibble(age = "first-year",
                    int = evSSurvDraws$b_Intercept,
                    beta = evSSurvDraws$b_fungicide) %>%
  full_join(tibble(age = "adult",
                   int = evASurvDraws$b_Intercept,
                   beta = evASurvDraws$b_fungicide)) %>%
  mutate(age = fct_relevel(age, "first-year"),
         prob_int = exp(int) / (1 + exp(int)),
         prob_fung = exp(int + beta) / (1 + exp(int + beta)),
         prob_change = 100 * (prob_fung - prob_int) / prob_int)

# fungicide effect without competition
ggplot(estDraws, aes(x = sp, y = beta)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  labs(y = "Effect of fungicide on establishment (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))
# looks better with beta than prob_change

ggplot(survDraws, aes(x = age, y = beta)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  labs(y = "Effect of fungicide on perennial survival (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank())


#### fungicide effect on competition ####

# edit
compEstDraws <- tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundEvadult`,
                       ctrl = mvEstDraws$`b_background_density:backgroundEvadult`,
                       background = "*E. virginicus* adult effects") %>%
  full_join(tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = mvEstDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects")) %>%
  full_join(tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = mvEstDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects")) %>%
  mutate(sp = "M. vimineum") %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evSEstDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects",
                   sp = "E. virginicus")) %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSEstDraws$`b_background_density:backgroundEvadult`,
                   background = "*E. virginicus* adult effects",
                   sp = "E. virginicus")) %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSEstDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects",
                   sp = "E. virginicus")) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         perc_change = 100 * fung / abs(ctrl))

compSurvDraws <- tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundEvseedling`,
                        ctrl = evSSurvDraws$`b_background_density:backgroundEvseedling`,
                        background = "*E. virginicus* first-year effects",
                        age = "first-year") %>%
  full_join(tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSSurvDraws$`b_background_density:backgroundEvadult`,
                   background = "*E. virginicus* adult effects",
                   age = "first-year")) %>%
  full_join(tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSSurvDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects",
                   age = "first-year")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evASurvDraws$`b_background_density:backgroundEvadult`,
                   background = "*E. virginicus* adult effects",
                   age = "adult")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evASurvDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects",
                   age = "adult")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evASurvDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects",
                   age = "adult")) %>%
  mutate(age = fct_relevel(age, "first-year"),
         perc_change = 100 * fung / abs(ctrl))

# figure
ggplot(compEstDraws, aes(x = background, y = fung)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(~ sp, scales = "free_y", ncol = 1) +
  labs(y = "Effect of fungicide on establishment response (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 7, color = "black"),
        strip.text = element_text(face = "italic"))

ggplot(compSurvDraws, aes(x = background, y = fung)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(~ age, scales = "free_y", ncol = 1) +
  labs(y = "Effect of fungicide on perennial survival response (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 7, color = "black"))


#### values for text ####

estDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_change)

survDraws %>%
  group_by(age) %>%
  mean_hdci(prob_change)

compEstDraws %>%
  group_by(sp, background) %>%
  mean_hdci(perc_change)

compSurvDraws %>%
  group_by(age, background) %>%
  mean_hdci(perc_change)
