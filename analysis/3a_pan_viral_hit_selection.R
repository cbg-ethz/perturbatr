#!/usr/bin/env Rscript

library(tibble)
library(dplyr)
library(ggplot2)
library(perturbatr)
library(xtable)


.random.effects.model <- function(dat, bootstrap=10)
{
  weights <- rep(1, nrow(dat))
  weights[dat$Design == "pooled"] <- 1.5

  lmm.fit <- perturbatr::hm(
      obj=as(dat, "PerturbationData"),
      "Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType)",
      bootstrap.cnt=0,
      weights=weights,
      drop=T)

  (lmm.fit)
}

random.effects.model <- function(rnai.screen, boo=0)
{
  lmm.fit     <- .random.effects.model(rnai.screen, boo)

  gene.hits <- geneHits(lmm.fit) %>%
    dplyr::select(GeneSymbol, Effect, Qval)
  gps <- nestedGeneHits(lmm.fit) %>%
    dplyr::select(Condition, GeneSymbol, Effect) %>%
    tidyr::spread(Condition, Effect)

  full.table <- dplyr::left_join(gene.hits, gps, by="GeneSymbol") %>%
    dplyr::arrange(desc(abs(Effect)))

  list(fit=lmm.fit, full.table=full.table)
}


run <- function()
{
  path        <- "./data"
  rnai.file   <- paste(path,  "/rnai_screen_normalized_2.rds", sep="/")
  lmm.outfile <- paste0(path, "/", "lmm_fit", ".rds")

  rnai.screen <- readRDS(rnai.file)
  fit <- random.effects.model(rnai.screen)

  saveRDS(fit, lmm.outfile)
}

run()
