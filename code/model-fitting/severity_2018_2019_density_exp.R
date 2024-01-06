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
library(patchwork)

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
write_csv(tidy(evSSevMod, conf.method = "HPDinterval"), 
          "output/evS_severity_model_2018_2019_density_exp.csv")
write_csv(tidy(evASevMod, conf.method = "HPDinterval"), 
          "output/evA_severity_model_2018_2019_density_exp.csv")
write_csv(tidy(mvSevMod, conf.method = "HPDinterval"), 
          "output/mv_severity_model_2018_2019_density_exp.csv")

# load
load("output/evS_severity_model_2018_2019_density_exp.rda")
load("output/evA_severity_model_2018_2019_density_exp.rda")
load("output/mv_severity_model_2018_2019_density_exp.rda")


#### disease severity over time ####

# posterior draws
evSSevDraws <- as_draws_df(evSSevMod)
mvSevDraws <- as_draws_df(mvSevMod)
evASevDraws <- as_draws_df(evASevMod)

# select time values
timeSevDraws <- evSSevDraws %>%
  select(b_Intercept, starts_with("b_year")) %>%
  mutate(sp_age = "first-year *E. virginicus*") %>%
  full_join(evASevDraws %>%
              select(b_Intercept, starts_with("b_year")) %>%
              mutate(sp_age = "adult *E. virginicus*")) %>%
  full_join(mvSevDraws %>%
              select(b_Intercept, starts_with("b_year")) %>%
              mutate(sp_age = "*M. vimineum*")) %>%
  transmute(sp_age = sp_age,
            n = rep(1:15000, 3),
            year_1_Jul = b_Intercept,
            year_1_late_Aug = b_Intercept + b_year_monthyear1_late_aug,
            year_1_Sep = b_Intercept + b_year_monthyear1_sep,
            year_2_early_Aug = b_Intercept + b_year_monthyear2_early_aug,
            year_2_Jul = b_Intercept + b_year_monthyear2_jul,
            year_2_Jun = b_Intercept + b_year_monthyear2_jun,
            year_2_late_Aug = b_Intercept + b_year_monthyear2_late_aug) %>%
  pivot_longer(cols = -c(sp_age, n),
               names_to = "time",
               values_to = "sev") %>%
  mutate(sev_perc = 100 * plogis(sev),
         month = str_replace_all(time, "_", " ") %>%
           fct_relevel("year 2 Jul", after = 5) %>%
           fct_relevel("year 2 early Aug", after = 5),
         sp_age = fct_relevel(sp_age, "*M. vimineum*",
                              "first-year *E. virginicus*"))

sev_time_fig <- ggplot(timeSevDraws, aes(x = month, y = sev_perc, color = sp_age)) +
  stat_pointinterval(fatten_point = 5,
                     shape = 95,
                     point_interval = mean_hdci,
                     .width = 0.95,
                     position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Focal group") +
  labs(x = "Month",
       y = "Disease severity (%)") +
  fig_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.text = element_markdown())


#### fungicide effect without comp. figure ####

# combine for fungicide effect without competition
SevDraws <- tibble(sp_age = "first-year *E. virginicus*",
                    int = evSSevDraws$b_Intercept,
                    beta = evSSevDraws$b_fungicide) %>%
  full_join(tibble(sp_age = "*M. vimineum*",
                   int = mvSevDraws$b_Intercept,
                   beta = mvSevDraws$b_fungicide)) %>%
  full_join(tibble(sp_age = "adult *E. virginicus*",
                   int = evASevDraws$b_Intercept,
                   beta = evASevDraws$b_fungicide)) %>%
  mutate(sp_age = fct_relevel(sp_age, "*M. vimineum*",
                              "first-year *E. virginicus*"),
         resp_int = 100 * plogis(int),
         resp_fung = 100 * plogis(int + beta),
         resp_diff = resp_fung - resp_int)

# fungicide effect without competition
sev_fung_fig <- ggplot(SevDraws, aes(x = sp_age, y = beta)) +
  geom_hline(yintercept = 0) +
  stat_pointinterval(fatten_point = 3,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1)) +
  scale_x_discrete(labels = c(~ atop(NA, atop(NA, textstyle(italic("M. vimineum")))),
                              ~ atop(NA, atop(textstyle("first-year"), 
                                              textstyle(italic("E. virginicus")))),
                              ~ atop(NA, atop(textstyle("adult"), 
                                              textstyle(italic("E. virginicus")))))) +
  labs(x = "Focal group",
       y = "Fungicide effect on disease\nseverity (log-odds)") +
  fig_theme +
  theme(axis.text.x = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 5))


#### fungicide effect on competition ####

# edit
compSevDraws <- tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundEvadult`,
                        ctrl = mvSevDraws$`b_background_density:backgroundEvadult`,
                        background = "adult *E. virginicus*") %>%
  full_join(tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = mvSevDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = mvSevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = mvSevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*")) %>%
  mutate(sp_age = "*M. vimineum*") %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evSSevDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evSSevDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evSSevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evSSevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   sp_age = "first-year *E. virginicus*")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundEvadult`,
                   ctrl = evASevDraws$`b_background_density:backgroundEvadult`,
                   background = "adult *E. virginicus*",
                   sp_age = "adult *E. virginicus*")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundEvseedling`,
                   ctrl = evASevDraws$`b_background_density:backgroundEvseedling`,
                   background = "first-year *E. virginicus*",
                   sp_age = "adult *E. virginicus*")) %>%
  full_join(tibble(fung = evASevDraws$`b_fungicide:background_density:backgroundMvseedling`,
                   ctrl = evASevDraws$`b_background_density:backgroundMvseedling`,
                   background = "*M. vimineum*",
                   sp_age = "adult *E. virginicus*")) %>%
  mutate(background = fct_relevel(background, "*M. vimineum*",
                                  "first-year *E. virginicus*"),
         sp_age = fct_relevel(sp_age, "*M. vimineum*",
                              "first-year *E. virginicus*"))

# figure
sev_comp_fig <- ggplot(compSevDraws, aes(x = background, y = fung, color = sp_age)) +
  geom_hline(yintercept = 0, color = "grey") +
  stat_pointinterval(fatten_point = 5,
                     shape = 95,
                     point_interval = mean_hdci,
                     .width = c(0.95, 1),
                     position = position_dodge(0.5)) +
  scale_color_manual(values = col_pal, name = "Focal group") +
  scale_x_discrete(labels = c(~ atop(NA, atop(NA, textstyle(italic("M. vimineum")))),
                              ~ atop(NA, atop(textstyle("first-year"), 
                                              textstyle(italic("E. virginicus")))),
                              ~ atop(NA, atop(textstyle("adult"), 
                                              textstyle(italic("E. virginicus")))))) +
  labs(x = "Background group",
       y = "Fungicide effect on disease\nseverity response (log-odds)") +
  fig_theme +
  theme(legend.text = element_markdown(),
        axis.text.x = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 5))


#### combined figure ####

sev_com_fig <- sev_time_fig / (sev_fung_fig + sev_comp_fig + theme(legend.position = "none")) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        plot.margin = margin(5.5, 5.5, -1, 5.5))

ggsave("output/combined_severity_figure_2018_2019_density_exp.png",
       sev_com_fig,
       width = 6, height = 7)


#### values for text ####

# fungicide effects without neighbors
SevDraws %>%
  group_by(sp_age) %>%
  mean_hdci(resp_diff)

compSevDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(ctrl)

compSevDraws %>%
  group_by(sp_age, background) %>%
  mean_hdci(fung)

evSSevDraws %>%
  mutate(prob_none = 100 * plogis(b_Intercept + b_fungicide),
         prob_EvS = 100 * plogis(b_Intercept + b_fungicide + 
                                  `b_background_density:backgroundEvseedling` * 10 + 
                                  `b_fungicide:background_density:backgroundEvseedling` * 10),
         prob_EvA = 100 * plogis(b_Intercept + b_fungicide + 
                                   `b_background_density:backgroundEvadult` * 10 + 
                                   `b_fungicide:background_density:backgroundEvadult` * 10)) %>%
  transmute(EvSeffect = prob_EvS - prob_none,
            EvAeffect = prob_EvA - prob_none) %>%
  mean_hdci()
