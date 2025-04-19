##### outputs ####

#### maybe don't need establishment from here ####

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
# figures
# output/establishment_fungicide_figure_2018_2019_density_exp.rda
# output/survival_fungicide_figure_2018_2019_density_exp.rda
# output/establishment_competition_figure_2018_2019_density_exp.rda
# output/survival_competition_figure_2018_2019_density_exp.rda

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
# library(tidybayes)
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

# # 2018 survival
# # make survival 1 if the plant produced seeds in summer
# # remove NA's 
# survD1Dat2 <- survD1Dat %>%
#   filter(month == "September" & focal == 1) %>%
#   left_join(plots) %>%
#   mutate(survival = case_when(seeds_produced == 1 ~ 1, 
#                               TRUE ~ survival),
#          fungicide = ifelse(treatment == "fungicide", 1, 0),
#          background = if_else(background == "none", "Mv seedling", background),
#          year = "1") %>%
#   select(-c(month, field_notes, seeds_produced)) %>%
#   filter(!is.na(survival))
# 
# # 2019 survival data needs list of all plants (only replacements recorded)
# # merge ID lists with plot
# focD2Dat <- plots %>%
#   select(plot, treatment) %>%
#   expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
#   expand_grid(ID = as.character(c(1, 2, 3))) %>%
#   mutate(sp = "Mv",
#          age = "seedling") %>%
#   full_join(plots %>%
#               select(plot, treatment) %>%
#               expand_grid(site = c("D1", "D2", "D3", "D4")) %>%
#               expand_grid(tibble(ID = c("1", "2", "3", "A"),
#                                  age = c(rep("seedling", 3), "adult"))) %>%
#               mutate(sp = "Ev"))
# 
# # 2019 focal survival
# # 2018 survival starts in June because plants were replaced through May 24
# survD2Dat2 <- survD2Dat %>%
#   filter(focal == 1 & replace_date > 20190531) %>%
#   group_by(site, plot, treatment, sp, ID) %>%
#   summarise(plantings = length(unique(replace_date)) + 1) %>%
#   ungroup() %>%
#   full_join(focD2Dat) %>%
#   left_join(plots) %>%
#   mutate(plantings = replace_na(plantings, 1),
#          survival = ifelse(plantings > 1, 0, 1),
#          fungicide = ifelse(treatment == "fungicide", 1, 0),
#          background = if_else(background == "none", "Mv seedling", background),
#          year = "2")

# annual perennial survival 2018-2019
evSurvD1Dat <- survD1Dat %>%
  filter(month == "April" & sp == "Ev") %>%
  left_join(plots) %>%
  mutate(survival = case_when(seeds_produced == 1 ~ 1, 
                              TRUE ~ survival),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotID = paste(site, plot, fungicide, sep = "_")) %>%
         # background = if_else(background == "none", "Mv seedling", background)) %>% 
  select(-c(month, field_notes, seeds_produced)) %>%
  filter(!is.na(survival))
# includes non-focal

# divide data by focal and combine years
# evSEstDat <- filter(survD1Dat2, sp == "Ev" & age == "seedling") %>%
#   full_join(filter(survD2Dat2, sp == "Ev" & age == "seedling"))
# 
# mvEstDat <- filter(survD1Dat2, sp == "Mv" & age == "seedling") %>%
#   full_join(filter(survD2Dat2, sp == "Mv" & age == "seedling"))

evSSurvDat <- filter(evSurvD1Dat, age == "seedling" & focal == 1)
evASurvDat <- filter(evSurvD1Dat, age == "adult" & focal == 1)


#### fit establishment models ####

# # initial visualization
# ggplot(evSEstDat, aes(x = background_density, y = survival, color = background)) +
#   geom_point() +
#   geom_smooth(method = "glm") +
#   facet_grid(year ~ treatment)
# 
# ggplot(mvEstDat, aes(x = background_density, y = survival, color = background)) +
#   geom_point() +
#   geom_smooth(method = "glm") +
#   facet_grid(year ~ treatment)
# 
# # fit models
# evSEstMod <- brm(data = evSEstDat, family = bernoulli,
#                  survival ~ fungicide + background_density:background + background_density:background:fungicide + year + 
#                    (1|site/plot),
#                  prior <- c(prior(normal(0, 10), class = "Intercept"),
#                             prior(normal(0, 10), class = "b")), # use default for sigma
#                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
#                  control = list(adapt_delta = 0.99))
# mod_check_fun(evSEstMod)
# 
# mvEstMod <- update(evSEstMod, newdata = mvEstDat)
# mod_check_fun(mvEstMod)
# 
# # save models
# save(evSEstMod, file = "output/evS_establishment_model_2018_2019_density_exp.rda")
# save(mvEstMod, file = "output/mv_establishment_model_2018_2019_density_exp.rda")
# 
# # save tables
# write_csv(tidy(evSEstMod, conf.method = "HPDinterval"), 
#           "output/evS_establishment_model_2018_2019_density_exp.csv")
# write_csv(tidy(mvEstMod, conf.method = "HPDinterval"), 
#           "output/mv_establishment_model_2018_2019_density_exp.csv")
# 
# # load
# load("output/evS_establishment_model_2018_2019_density_exp.rda")
# load("output/mv_establishment_model_2018_2019_density_exp.rda")


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
                  survival ~ fungicide + (1|site/plotID),
                  prior <- c(prior(normal(0, 10), class = "Intercept"),
                             prior(normal(0, 10), class = "b")), # use default for sigma
                  iter = 6000, warmup = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.99))
mod_check_fun(evSSurvMod)

evASurvMod <- update(evSSurvMod, newdata = evASurvDat)
mod_check_fun(evASurvMod)

# save models
save(evSSurvMod, file = "output/evS_survival_model_2018_2019_density_exp.rda")
save(evASurvMod, file = "output/evA_survival_model_2018_2019_density_exp.rda")

# save tables
write_csv(tidy(evSSurvMod, conf.method = "HPDinterval"), 
          "output/evS_survival_model_2018_2019_density_exp.csv")
write_csv(tidy(evASurvMod, conf.method = "HPDinterval"), 
          "output/evA_survival_model_2018_2019_density_exp.csv")

# # load
# load("output/evS_survival_model_2018_2019_density_exp.rda")
# load("output/evA_survival_model_2018_2019_density_exp.rda")


#### maybe don't need to keep below ####

##### fungicide effect without comp. figure ####

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
         prob_int = 100 * plogis(int),
         prob_fung = 100 * plogis(int + beta),
         prob_diff = prob_fung - prob_int)

survDraws <- tibble(age = "first-year",
                    int = evSSurvDraws$b_Intercept,
                    beta = evSSurvDraws$b_fungicide) %>%
  full_join(tibble(age = "adult",
                   int = evASurvDraws$b_Intercept,
                   beta = evASurvDraws$b_fungicide)) %>%
  mutate(age = fct_relevel(age, "first-year"),
         prob_int = 100 * plogis(int),
         prob_fung = 100 * plogis(int + beta),
         prob_diff = prob_fung - prob_int)

# fungicide effect without competition
est_fung_fig <- ggplot(estDraws, aes(x = sp, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  labs(y = "Fungicide effect on establishment (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "italic"))

surv_fung_fig <- ggplot(survDraws, aes(x = age, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  scale_x_discrete(labels = c(~ atop(NA, atop(textstyle("first-year"), 
                                              textstyle(italic("E. virginicus")))),
                              ~ atop(NA, atop(textstyle("adult"), 
                                              textstyle(italic("E. virginicus")))))) +
  labs(y = "Fungicide effect on survival (log-odds)") +
  fig_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(vjust = 2.7))

save(est_fung_fig, file = "output/establishment_fungicide_figure_2018_2019_density_exp.rda")
save(surv_fung_fig, file = "output/survival_fungicide_figure_2018_2019_density_exp.rda")


#### fungicide effect on competition ####

# edit
compEstDraws <- tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundEvadult`,
                       ctrl = mvEstDraws$`b_background_density:backgroundEvadult`,
                       background = "adult *E. virginicus*") %>%
  full_join(tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = mvEstDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = mvEstDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = mvEstDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*")) %>%
  mutate(sp_age = "*M. vimineum*") %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evSEstDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSEstDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSEstDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSEstDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   sp_age = "first-year *E. virginicus*")) %>%
  mutate(background = fct_relevel(background, "*M. vimineum*",
                                  "first-year *E. virginicus*"))

compSurvDraws <- tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundEvseedling`,
                        ctrl = evSSurvDraws$`b_background_density:backgroundEvseedling`,
                        background = "first-year *E. virginicus*",
                        age = "first-year") %>%
  full_join(tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSSurvDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   age = "first-year")) %>%
  full_join(tibble(fung = evSSurvDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSSurvDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   age = "first-year")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evASurvDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   age = "adult")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evASurvDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   age = "adult")) %>%
  full_join(tibble(fung = evASurvDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evASurvDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   age = "adult")) %>%
  mutate(background = fct_relevel(background, "*M. vimineum*",
                                  "first-year *E. virginicus*"),
         sp_age = paste(age, "*E. virginicus*") %>%
           fct_relevel("first-year *E. virginicus*"))

# figure
est_comp_fig <- ggplot(compEstDraws, aes(x = background, y = fung, color = sp_age)) +
  geom_hline(yintercept = 0, color = "grey") +
  stat_pointinterval(fatten_point = 5,
                     shape = 95,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1),
                     position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Focal group") +
  labs(x = "Background group",
       y = "Fungicide effect on\nestablishment\nresponse (log-odds)") +
  fig_theme +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown())

surv_comp_fig <- ggplot(compSurvDraws, aes(x = background, y = fung, color = sp_age)) +
  geom_hline(yintercept = 0, color = "grey") +
  stat_pointinterval(fatten_point = 5,
                     shape = 95,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1),
                     position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Focal group") +
  labs(x = "Background group",
       y = "Fungicide effect on\nsurvival\nresponse (log-odds)") +
  fig_theme +
  theme(axis.text.x = element_markdown(),
        legend.text = element_markdown())

save(est_comp_fig, file = "output/establishment_competition_figure_2018_2019_density_exp.rda")
save(surv_comp_fig, file = "output/survival_competition_figure_2018_2019_density_exp.rda")


#### values for text ####

estDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_int)

estDraws %>%
  group_by(sp) %>%
  mean_hdci(prob_diff)

survDraws %>%
  group_by(age) %>%
  mean_hdci(prob_int)

survDraws %>%
  group_by(age) %>%
  mean_hdci(prob_diff)

compEstDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(ctrl)

evSEstDraws %>%
  mutate(prob_none = 100 * plogis(b_Intercept),
            prob_Mv = 100 * plogis(b_Intercept + 
                                   `b_background_density:backgroundMvseedling` * 10),
            prob_EvA = 100 * plogis(b_Intercept + 
                                  `b_background_density:backgroundEvadult` * 10)) %>%
  transmute(Mveffect = prob_Mv - prob_none,
            EvAeffect = prob_EvA - prob_none) %>%
  mean_hdci()

compEstDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(fung)

compSurvDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(ctrl)

evSSurvDraws %>%
  mutate(prob_none = 100 * plogis(b_Intercept),
         prob_Mv = 100 * plogis(b_Intercept + 
                                  `b_background_density:backgroundMvseedling` * 10),
         prob_EvA = 100 * plogis(b_Intercept + 
                                   `b_background_density:backgroundEvadult` * 10)) %>%
  transmute(Mveffect = prob_Mv - prob_none,
            EvAeffect = prob_EvA - prob_none) %>%
  mean_hdci()

compSurvDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(fung)

evSSurvDraws %>%
  mutate(prob_none = 100 * plogis(b_Intercept + b_fungicide),
         prob_Mv = 100 * plogis(b_Intercept + b_fungicide + 
                                  `b_background_density:backgroundMvseedling` * 10 + 
                                  `b_fungicide:background_density:backgroundMvseedling` * 10),
         prob_EvA = 100 * plogis(b_Intercept + b_fungicide + 
                                   `b_background_density:backgroundEvadult` * 10 + 
                                   `b_fungicide:background_density:backgroundEvadult` * 10)) %>%
  transmute(Mveffect = prob_Mv - prob_none,
            EvAeffect = prob_EvA - prob_none) %>%
  mean_hdci()

