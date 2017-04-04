library(knockout)
library(data.table)
library(tidyverse)

fl   <- "~/PROJECTS/sysvirdrug_project/data/svd/integrated_data_files/rnai_screen_normalized.rds"
rnai <- readRDS(fl)
h    <- hyper.statistic(rnai)
