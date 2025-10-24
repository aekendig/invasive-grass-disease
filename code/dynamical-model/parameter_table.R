#### set-up ####

# clear environment
rm(list = ls())

# import parameters (loads tidyverse)
source("code/dynamical-model/parameters.R")

# load other packages
library(janitor)


#### parameters from simulations ####

# parameter combinations
params_iters <- 15000 # all draws

# parameters with fungicide effect on germination
params_gfung <- params_fun(iters = params_iters, gA_type = "fungicide")

# parameters with infection effect on germination
params_ginf <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws)

# parameters with infection effect and priors
params_infP <- params_fun(iters = params_iters, gA_type = "infection",
                          draws = params_gfung[["healthy"]]$draws,
                          priors = "Stricker")


#### summary function ####

params_sum_fun <- function(input_list){
  
  map(names(input_list), 
      ~ {treatment <- .x
      map(input_list[[treatment]], 
          ~ summarise_all(.x, mean) %>%  
            pivot_longer(everything(), 
                         names_to = "name", values_to = "mean")) %>%
        list_rbind() %>%
        mutate(treatment = treatment)
      }
  ) %>% 
    list_rbind() %>% 
    pivot_wider(names_from = "treatment", values_from = "mean")
}


#### summarize and combine ####

# apply function
param_tab <- params_sum_fun(params_gfung) %>% 
  mutate(set1 = "x") %>% 
  full_join(params_sum_fun(params_ginf) %>% 
              mutate(set2 = "x")) %>% 
  full_join(params_sum_fun(params_infP) %>% 
              mutate(set3 = "x")) %>% 
  mutate(across(.cols = c("healthy", "disease"),
                ~ round_half_up(.x, 2)),
         across(starts_with("set"),
                ~ replace_na(.x, ""))) %>% 
  relocate(disease, .after = "name") %>% 
  arrange(name, desc(set1), desc(set2), desc(set3)) %>% 
  filter(name != ".draw")

# save
write_csv(param_tab, "output/parameter_value_mean_table.csv")

