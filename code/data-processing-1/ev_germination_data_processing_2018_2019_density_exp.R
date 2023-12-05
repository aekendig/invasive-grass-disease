##### outputs ####

# intermediate-data/ev_germination_disease_2018_2019_density_exp.csv

#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)

# import data
evGermDat <- read_csv("data/ev_germination_2018_2019_density_exp.csv")


#### edit data ####

# look at notes
unique(evGermDat$notes) # may be a contaminate on plate

# correct reduction in emergents between week 3 and 4
# correct the increase in cut-tops
# correct repair of cut tops between weeks 3 and 4
evGermDat2 <- evGermDat %>%
  filter(seeds_planted > 0) %>%
  mutate(week_4_emerg = case_when(week_4_emerg < week_3_emerg ~ week_3_emerg,
                                  TRUE ~ week_4_emerg),
         week_3_cut_tops = case_when(week_4_cut_tops > week_3_cut_tops & week_4_cut_tops <= week_2_emerg ~ week_4_cut_tops,
                                     TRUE ~ week_3_cut_tops),
         week_4_cut_tops = case_when(week_4_cut_tops > week_2_emerg ~ week_3_cut_tops,
                                     week_4_cut_tops < week_3_cut_tops ~ week_3_cut_tops,
                                     TRUE ~ week_4_cut_tops),
         week_3_new_emerg = week_3_emerg - week_3_cut_tops,
         week_4_new_emerg = week_4_emerg - week_4_cut_tops,
         germinants = week_2_emerg + week_4_new_emerg + week_4_soil_germ) %>%
  select(year, site, plot, treatment, age, seeds_planted, germinants)

# save
write_csv(evGermDat2, "intermediate-data/ev_germination_2018_2019_density_exp.csv")
