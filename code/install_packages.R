# packages used in this repository
list.of.packages <- c("tidyverse", "brms", "broom.mixed", "GGally", "betareg", "MASS", "lubridate", "caTools", "weathermetrics", "cowplot", "tidybayes", "car")

# packages in that list that you do not have
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

# install missing packages
if(length(new.packages)) install.packages(new.packages)