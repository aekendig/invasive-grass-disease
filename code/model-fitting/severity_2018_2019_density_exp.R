##### outputs ####

# models
# "output/evS_severity_model_2018_2019_density_exp.rda"
# "output/evA_severity_model_2018_2019_density_exp.rda"
# "output/mv_severity_model_2018_2019_density_exp.rda"
# tables
# "output/evS_severity_model_2018_2019_density_exp.csv"
# "output/evA_severity_model_2018_2019_density_exp.csv"
# "output/mv_severity_model_2018_2019_density_exp.csv"


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

# function to transform data to account for 0's and 1's
transform01 <- function(x) {
  (x * (length(x) - 1) + 0.5) / (length(x))
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
         year_month = paste0("year", year, "_", month),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background)) %>%
  filter(!is.na(severity))

sevD2Dat2 <- sevD2Dat %>%
  filter(focal == 1) %>%
  select(month, site, plot, treatment, sp, ID, age, leaves_tot, leaves_infec, leaf_area.pix, lesion_area.pix) %>%
  left_join(plots, by = c("plot", "treatment"),
            relationship = "many-to-many") %>%
  filter(!(month %in% c("may", "sep"))) %>% # too much data missing
  mutate(leaves_infec = case_when(leaves_infec == 0 & lesion_area.pix > 0 ~ 1, # add an infected leaf if scan found lesions
                                  TRUE ~ leaves_infec),
         severity = (lesion_area.pix * leaves_infec) / (leaf_area.pix * leaves_tot),
         severity = case_when(leaves_infec == 0 & is.na(severity) ~ 0, # make severity zero if no leaves infected and no severity info
                              TRUE ~ severity),
         year = "2",
         year_month = paste0("year", year, "_", month),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         background = if_else(background == "none", "Mv seedling", background)) %>%
  filter(!is.na(severity))

# combine by plant group
evSSevDat <- filter(sevD1Dat2, sp == "Ev" & age == "seedling") %>%
  full_join(filter(sevD2Dat2, sp == "Ev" & age == "seedling")) %>%
  mutate(severity_t = transform01(severity))

evASevDat <- filter(sevD1Dat2, sp == "Ev" & age == "adult") %>%
  full_join(filter(sevD2Dat2, sp == "Ev" & age == "adult")) %>%
  mutate(severity_t = transform01(severity))

mvSevDat <- filter(sevD1Dat2, sp == "Mv") %>%
  full_join(filter(sevD2Dat2, sp == "Mv")) %>%
  mutate(severity_t = transform01(severity))


#### fit models ####

# initial visualizations
ggplot(evSSevDat, aes(x = background_density, y = severity, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year_month ~ treatment)

ggplot(evSSevDat, aes(x = severity_t)) +
  geom_density()

ggplot(evASevDat, aes(x = background_density, y = severity, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year_month ~ treatment)

ggplot(evASevDat, aes(x = severity_t)) +
  geom_density()

ggplot(mvSevDat, aes(x = background_density, y = severity, color = background)) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_grid(year_month ~ treatment)

ggplot(mvSevDat, aes(x = severity_t)) +
  geom_density()

# fit models
evSSevMod <- brm(data = evSSevDat, family = "beta",
                  severity_t ~ fungicide + background_density:background + background_density:background:fungicide + year_month + 
                    (1|site/plot),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3, 
                 control = list(adapt_delta = 0.9999))
mod_check_fun(evSSevMod)

mvSevMod <- update(evSSevMod, newdata = mvSevDat,
                   control = list(adapt_delta = 0.9999, max_treedepth = 15))
mod_check_fun(mvSevMod)

evASevMod <- brm(data = evASevDat, family = lognormal,
                 severity_t ~ fungicide + background_density:background + background_density:background:fungicide + year_month + 
                    (1|site),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.999999))
mod_check_fun(evASevMod)

# save models
save(evSSevMod, file = "output/evS_severity_model_2018_2019_density_exp.rda")
save(evASevMod, file = "output/evA_severity_model_2018_2019_density_exp.rda")
save(mvSevMod, file = "output/mv_severity_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSSevMod), "output/evS_severity_model_2018_2019_density_exp.csv")
write_csv(tidy(evASevMod), "output/evA_severity_model_2018_2019_density_exp.csv")
write_csv(tidy(mvSevMod), "output/mv_severity_model_2018_2019_density_exp.csv")

# load
load("output/evS_severity_model_2018_2019_density_exp.rda")
load("output/evA_severity_model_2018_2019_density_exp.rda")
load("output/mv_severity_model_2018_2019_density_exp.rda")


#### fungicide effect without comp. figure ####

# posterior draws
evSSevDraws <- as_draws_df(evSSevMod)
mvSevDraws <- as_draws_df(mvSevMod)
evASevDraws <- as_draws_df(evASevMod)

# combine for fungicide effect without competition
SevDraws <- tibble(sp = "E. virginicus",
                    age = "first-year",
                    int = evSSevDraws$b_Intercept,
                    beta = evSSevDraws$b_fungicide) %>%
  full_join(tibble(sp = "M. vimineum",
                   age = "first-year",
                   int = mvSevDraws$b_Intercept,
                   beta = mvSevDraws$b_fungicide)) %>%
  full_join(tibble(sp = "E. virginicus",
                   age = "adult",
                   int = evASevDraws$b_Intercept,
                   beta = evASevDraws$b_fungicide)) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"),
         age = fct_relevel(age, "first-year"),
         resp_int = plogis(int),
         resp_fung = plogis(int + beta),
         resp_diff = resp_fung - resp_int,
         resp_change = 100 * (resp_fung - resp_int) / resp_int)

# fungicide effect without competition
ggplot(SevDraws, aes(x = sp, y = beta)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(~ age, scales = "free_x") +
  labs(y = "Effect of fungicide on severity (log odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))


#### fungicide effect on competition ####

# edit
compSevDraws <- tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundEvadult`,
                        ctrl = mvSevDraws$`b_background_density:backgroundEvadult`,
                        background = "*E. virginicus* adult effects") %>%
  full_join(tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = mvSevDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects")) %>%
  full_join(tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = mvSevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects")) %>%
  mutate(sp = "M. vimineum") %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evSSevDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects",
                   sp = "E. virginicus\nfirst-year")) %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSSevDraws$`b_background_density:backgroundEvadult`,
                   background = "*E. virginicus* adult effects",
                   sp = "E. virginicus\nfirst-year")) %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSSevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects",
                   sp = "E. virginicus\nfirst-year")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evASevDraws$`b_background_density:backgroundEvadult`,
                   background = "*E. virginicus* adult effects",
                   sp = "E. virginicus\nadult")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evASevDraws$`b_background_density:backgroundEvseedling`,
                   background = "*E. virginicus* first-year effects",
                   sp = "E. virginicus\nadult")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evASevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum* effects",
                   sp = "E. virginicus\nadult")) %>%
  mutate(sp = fct_relevel(sp, "M. vimineum"))

# figure
ggplot(compSevDraws, aes(x = background, y = fung)) +
  geom_hline(yintercept = 0) +
  geom_violin(fill = "paleturquoise", color = "paleturquoise4", 
              draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(~ sp, scales = "free_y", ncol = 1) +
  labs(y = "Effect of fungicide on severity response (log odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 7, color = "black"),
        strip.text = element_text(face = "italic"))


#### values for text ####

# fungicide effects without neighbors
SevDraws %>%
  group_by(sp, age) %>%
  mean_hdci(resp_change)

compSevDraws %>%
  group_by(sp, background) %>%
  mean_hdci(ctrl)

compSevDraws %>%
  group_by(sp, background) %>%
  mean_hdci(fung + ctrl)
