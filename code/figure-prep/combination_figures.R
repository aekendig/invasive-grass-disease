#### outputs #####

# 


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(tidybayes)
library(ggtext)
library(patchwork)

# load figures
load("output/germination_fungicide_figure_2018_2019_density_exp.rda")
load("output/establishment_fungicide_figure_2018_2019_density_exp.rda")
load("output/survival_fungicide_figure_2018_2019_density_exp.rda")
load("output/seed_fungicide_figure_2018_2019_density_exp.rda")

load("output/establishment_competition_figure_2018_2019_density_exp.rda")
load("output/survival_competition_figure_2018_2019_density_exp.rda")
load("output/seed_competition_figure_2018_2019_density_exp.rda")


#### combine ####

fung_fig <- germ_fung_fig + est_fung_fig + surv_fung_fig + seed_fung_fig +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        plot.margin = margin(5.5, 5.5, -4, 5.5))

ggsave("output/combined_fungicide_figure_2018_2019_density_exp.png",
       fung_fig,
       width = 6, height = 6)

comp_fig <- est_comp_fig + theme(legend.position = "none", 
                                 axis.title.x = element_blank()) +
  surv_comp_fig + theme(legend.position = "none", 
                        axis.title.x = element_blank()) + 
  seed_comp_fig +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave("output/combined_competition_figure_2018_2019_density_exp.png",
       comp_fig,
       width = 6, height = 8)
