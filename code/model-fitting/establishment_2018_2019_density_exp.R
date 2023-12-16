##### outputs ####


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes)
library(brms)
library(broom.mixed)

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

# winter survival 2018-2019
winSurvD1Dat <- survD1Dat %>%
  filter(month == "April" | month == "September") %>%
  select(-field_notes) %>%
  spread(month, survival) %>%
  left_join(plots) %>%
  mutate(September = case_when(seeds_produced == 1 ~ 1, TRUE ~ September),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste0(site, plot, substr(treatment, 1, 1))) %>%
  filter(September == 1 & !is.na(April)) %>%
  select(-September) %>%
  rename(survival = April)

# annual adult survival 2018-2019
adultSurvD1Dat <- survD1Dat %>%
  filter(month == "April" & age == "adult" & !is.na(survival)) %>%
  select(-field_notes) %>%
  left_join(plots) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0))
# includes non-focal

# divide data by focal and combine years
evSEstDat <- filter(survD1Dat2, sp == "Ev" & age == "seedling") %>%
  full_join(filter(survD2Dat2, sp == "Ev" & age == "seedling"))

evAEstDat <- filter(survD1Dat2, sp == "Ev" & age == "adult") %>%
  full_join(filter(survD2Dat2, sp == "Ev" & age == "adult"))

mvEstDat <- filter(survD1Dat2, sp == "Mv" & age == "seedling") %>%
  full_join(filter(survD2Dat2, sp == "Mv" & age == "seedling"))


#### fit models ####

# initial visualization
ggplot(evSEstDat, aes(x = background_density, y = survival, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year ~ treatment)

ggplot(evAEstDat, aes(x = background_density, y = survival, color = background)) +
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

evAEstMod <- brm(data = evAEstDat, family = bernoulli,
                 survival ~ fungicide + background_density:background + background_density:background:fungicide + year + 
                   (1|site),
                 prior <- c(prior(normal(0, 10), class = "Intercept"),
                            prior(normal(0, 10), class = "b")), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.99))
mod_check_fun(evAEstMod)

mvEstMod <- update(evSEstMod, newdata = mvEstDat)
mod_check_fun(mvEstMod)

# save models
save(evSEstMod, file = "output/evS_establishment_model_2018_2019_density_exp.rda")
save(evAEstMod, file = "output/evA_establishment_model_2018_2019_density_exp.rda")
save(mvEstMod, file = "output/mv_establishment_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSEstMod), "output/evS_establishment_model_2018_2019_density_exp.csv")
write_csv(tidy(evAEstMod), "output/evA_establishment_model_2018_2019_density_exp.csv")
write_csv(tidy(mvEstMod), "output/mv_establishment_model_2018_2019_density_exp.csv")

# load
load("output/evS_establishment_model_2018_2019_density_exp.rda")
load("output/evA_establishment_model_2018_2019_density_exp.rda")
load("output/mv_establishment_model_2018_2019_density_exp.rda")


#### figures ####

# posterior draws
evSEstDraws <- as_draws_df(evSEstMod)
evAEstDraws <- as_draws_df(evAEstMod)
mvEstDraws <- as_draws_df(mvEstMod)

# combine
estDraws <- tibble(sp = "E. virginicus",
                   age = "first-year",
                   int = evSEstDraws$b_Intercept,
                   beta = evSEstDraws$b_fungicide) %>%
  full_join(tibble(sp = "E. virginicus",
                   age = "adult",
                   int = evAEstDraws$b_Intercept,
                   beta = evAEstDraws$b_fungicide)) %>%
  full_join(tibble(sp = "M. vimineum",
                   age = "first-year",
                   int = mvEstDraws$b_Intercept,
                   beta = mvEstDraws$b_fungicide)) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         age = fct_relevel(age, "first-year"),
         prob_int = exp(int) / (1 + exp(int)),
         prob_fung = exp(int + beta) / (1 + exp(int + beta)),
         prob_change = 100 * (prob_fung - prob_int) / prob_int)

ggplot(estDraws, aes(x = sp, y = prob_change)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(~ age, scales = "free_x", strip.position = "bottom") +
  labs(y = "Change in establishment with fungicide (%)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))
# looks better on logit scale (y = beta)
# may also want to denote age differently (different colors? put all on x-axis?)


#### start here: values for text ####

# also: effect of fungicide on density effects
# adult full-year survival


#### adult survival ####

mean(adultSurvD1Dat$survival)

adultSurvD1Mod <- brm(data = adultSurvD1Dat, family = bernoulli,
                      survival ~ 1 + (1|site),
                      prior <- prior(normal(0, 1), class = "Intercept"), # use default for sigma
                      iter = 6000, warmup = 1000, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.99)) 
# mod_check_fun(adultSurvD1Mod)

adultSurvD1Mod2 <- brm(data = adultSurvD1Dat, family = bernoulli,
                      survival ~ fungicide + (1|site),
                      prior <- c(prior(normal(0, 10), class = "Intercept"), # use default for sigma
                                 prior(normal(0, 10), class = "b")),
                      iter = 6000, warmup = 1000, chains = 3, cores = 3,
                      control = list(adapt_delta = 0.999)) 
# mod_check_fun(adultSurvD1Mod2)

# save
save(adultSurvD1Mod, file = "output/ev_adult_survival_model_2018_2019_density_exp.rda")
save(adultSurvD1Mod2, file = "output/ev_adult_survival_fungicide_model_2018_2019_density_exp.rda")

# load
load("output/ev_adult_survival_model_2018_2019_density_exp.rda")
load("output/ev_adult_survival_fungicide_model_2018_2019_density_exp.rda")


#### fungicide effects ####

mv_fung_eff = "fungicide = 0"
evS_fung_eff = "fungicide + fungicide:focs = 0"
evA_fung_eff = "fungicide + fungicide:foca = 0"

hypothesis(survFungD1Mod,
                          c(mv_fung_eff, evS_fung_eff, evA_fung_eff)) [[1]] %>%
  mutate(Year = "2018", Season = "growing season") %>%
  full_join(hypothesis(survFungD2Mod,
                       c(mv_fung_eff, evS_fung_eff, evA_fung_eff)) [[1]] %>%
              mutate(Year = "2019", Season = "growing season")) %>%
  full_join(hypothesis(winSurvFungD1Mod,
                       c(mv_fung_eff, evS_fung_eff)) [[1]] %>%
              mutate(Year = "2018-2019", Season = "winter")) %>%
  mutate(Focal = c(rep(c("Mv", "Ev first-year", "Ev adult"), 2), "Ev adult", "Ev first-year")) %>%
  select(Year, Season, Focal, Estimate:CI.Upper) %>%
  arrange(Season, Year, Focal)

