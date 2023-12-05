##### outputs ####

# intermediate-data/mv_germination_disease_2018_density_exp.csv

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
mvGermD1Dat1 <- read_csv("data/mv_germination_disease_set_1_2018_density_exp.csv")
mvGermD1Dat2 <- read_csv("data/mv_germination_disease_set_2_2018_density_exp.csv")


#### edit data ####

# notes
unique(mvGermD1Dat1$notes_check_1) # seedlings with lesions in notes
unique(mvGermD1Dat1$notes_check_2) # may be a contaminate on plate
unique(mvGermD1Dat1$notes_check_3)
unique(mvGermD1Dat1$notes_germination_check_1) # some plates were put into fridge during one day
unique(mvGermD1Dat1$notes_germination_final) 
unique(mvGermD1Dat2$notes) # may be a contaminate on plate

# Mv data
# average across trials
mvGermD1Dat <- mvGermD1Dat1 %>%
  mutate(germination_final = ifelse(is.na(germination_final), germination_check_1, germination_final),
         seeds_dark_check_1 = rowSums(cbind(seeds_dark_check_1, seeds_pink_check_1, seeds_red_check_1, seeds_green_check_1), na.rm = T),
         seeds_dark_check_2 = rowSums(cbind(seeds_dark_check_2, seeds_pink_check_2, seeds_red_check_2, seeds_green_check_2), na.rm = T),
         seeds_seeds_dark_check_3 = rowSums(cbind(seeds_dark_check_3, seeds_red_check_3, seeds_green_check_3), na.rm = T),
         seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, seeds_dark_check_3, na.rm = T),
         seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, seeds_light_check_3, na.rm = T)) %>%
  select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light) %>%
  full_join(mvGermD1Dat2 %>%
              mutate(seeds_dark = pmax(seeds_dark_check_1, seeds_dark_check_2, na.rm = T),
                     seeds_light = pmax(seeds_light_check_1, seeds_light_check_2, na.rm = T)) %>%
              select(site_plot, trial, seeds, germination_final, seeds_dark, seeds_light)) %>%
  mutate(site = gsub(" .*$", "", site_plot),
         plot = gsub(".* ","", site_plot) %>% 
           gsub("[^[:digit:]]", "", .) %>% 
           as.numeric(),
         treatment = gsub(".* ","", site_plot) %>% 
           gsub("[^[:alpha:]]", "", .) %>% 
           as.factor() %>%
           dplyr::recode("F" = "fungicide", "W" = "water"),
         site = ifelse(site == "P1", "D1", site)) %>%
  select(site, plot, treatment, trial, seeds, germination_final, seeds_dark, seeds_light)

# save
write_csv(mvGermD1Dat, "intermediate-data/mv_germination_disease_2018_density_exp.csv")
