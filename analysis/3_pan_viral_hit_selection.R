#!/usr/bin/env Rscript

library(data.table)
library(dtplyr)
library(dplyr)
library(ggplot2)
library(knockdown)
library(xtable)

.random.effects.model <- function(dat, bootstrap=10)
{
  lmm.fit     <- knockdown::lmm(obj=dat,
                               bootstrap.cnt=bootstrap,
                               drop=T,
                               weights=list(pooled=1.5, single=1))
  lmm.fit
}

random.effects.model <- function(rnai.screen, boo=0)
{
  lmm.fit     <- .random.effects.model(rnai.screen, boo)

  gene.hits <- lmm.fit@.gene.hits %>%
    dplyr::select(GeneSymbol, Effect, Qval)
  gps <- lmm.fit@.gene.pathogen.effects %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)

  full.table <- dplyr::left_join(gene.hits, gps, by="GeneSymbol") %>%
    .[order(-abs(Effect))]

  list(fit=lmm.fit, full.table=full.table)
}


run <- function()
{
  path        <- "./"
  rnai.file   <- paste(path, "data/rnai_screen_normalized.rds", sep="/")
  lmm.outfile  <- paste0(path, "/", "lmm_fit", ".rds")

  rnai.screen <- readRDS(rnai.file)
  fit <- random.effects.model(rnai.screen)

  saveRDS(fit, lmm.outfile)
}
