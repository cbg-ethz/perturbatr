library(dtplyr)
library(dplyr)
library(tidyr)
library(knockout)
library(ggplot2)
library(uuid)


# output directories
dirs <- c("/Users/simondi/PROJECTS/sysvirdrug_project/results/plots",
          "/Users/simondi/PROJECTS/sysvirdrug_project/src/package/analysis/plots",
          "/Users/simondi/PROJECTS/sysvirdrug_project/docs/sysvirdrug_modelling_paper/plots"
)

path <- "/Users/simondi/PROJECTS/sysvirdrug_project/src/package/analysis/"
diff.path <- paste(path, "data/diffusion_primary_data_23_5.rds", sep="/")
diff.model <- readRDS(diff.path)

