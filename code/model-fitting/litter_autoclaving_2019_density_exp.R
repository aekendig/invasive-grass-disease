#### set-up ####

# clear environment
rm(list = ls())

# load packages
library(tidyverse)

# import data
litter_weights <- read_csv("data/litter_autoclave_effects_apr_2018_litter_exp.csv")


#### summarize ####

# calcualte difference
litter_weights2 <- litter_weights %>% 
  mutate(difference_g = post_autoclave_weight_g - pre_autoclave_weight_g,
         perc_change = 100 * difference_g / pre_autoclave_weight_g)

# mean change
mean(litter_weights2$difference_g)
mean(litter_weights2$perc_change)


#### t-test ####

t.test(litter_weights$pre_autoclave_weight_g,
       litter_weights$post_autoclave_weight_g,
       paired = T)
