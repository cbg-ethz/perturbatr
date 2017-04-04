library(knockout)
library(data.table)
library(tidyverse)

fl   <- "~/PROJECTS/sysvirdrug_project/data/svd/integrated_data_files/rnai_screen_normalized.rds"
rnai <- readRDS(fl) %>% filter(Virus != "Virus")
h    <- hyper.statistic(rnai)
lmm    <- lmm(rnai)

padjust <- "BH"
summ.method <- "mean"
level <- "sirna"
do.summarization <- F
