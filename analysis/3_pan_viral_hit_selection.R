library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(knockout)
library(xtable)


wd <- "/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline/hit_selection/all_pathogen_hit_selection/"

.random.effects.model <- function(dat, bootstrap=10)
{
  model.form <- "Readout ~ Virus + (1 | GeneSymbol) + (1 | Virus:GeneSymbol) + (1 | ScreenType) + (1 | Virus:ScreenType)"
  lmm.fit     <- knockout::lmm(obj=dat,
                               bootstrap.cnt=100,
                               drop=T,
                               weights=list(pooled=1.5, single=1))
  lmm.fit
}

random.effects.model <- function(rnai.screen, file.name, boo=0)
{
  lmm.fit     <- .random.effects.model(rnai.screen, boo)
  re.file  <- paste0(wd,"/", file.name)
  message(paste("Writing random effect model hits to", re.file))


  gene.hits <- lmm.fit@.gene.hits %>%
    dplyr::select(GeneSymbol, Effect, Qval)
  gps <- lmm.fit@.gene.pathogen.effects %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)

  full.table <- dplyr::left_join(gene.hits, gps, by="GeneSymbol") %>%
    .[order(-abs(Effect))]

  saveRDS(lmm.fit, paste0(re.file, ".rds"))
  write.csv(full.table, paste0(re.file, ".csv"), quote=F, row.names=F)

  lmm.hits
}

do.pan.viral.analysis.on.different.data.sets <- function()
{

  setwd(wd)
  rnai.screen <- readRDS("integrated_data_files/rnai_screen_normalized.rds") %>% filter(Virus != "CVB")
  validation.screen <- readRDS("integrated_data_files/validation_screen_norm.rds")

  ## analysis on the primary data
  primary.lmm <- random.effects.model(rnai.screen, "random_effects_model_hits_primary_data")

  ## analysis on the combined data-sets
  dat.combined <- rbindlist(list(rnai.screen, validation.screen))
  lmm.dat.combined <- knockout::model.data.lmm(dat.combined, weights=list("pooled"=1.5, "single"=1))
  full.lmm <- random.effects.model(lmm.dat.combined, "random_effects_model_hits_full_data", 0)

}
