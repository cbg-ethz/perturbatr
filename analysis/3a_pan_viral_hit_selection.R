#!/usr/bin/env Rscript

# perturbatr: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbatr
#
# perturbatr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbatr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbatr. If not, see <http://www.gnu.org/licenses/>.


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
  data.dir    <- "./data"
  rnai.file   <- paste(data.dir,  "/rnai_screen_normalized_2.rds", sep="/")
  lmm.outfile <- paste0(data.dir, "/", "lmm_fit_2", ".rds")

  rnai.screen <- readRDS(rnai.file)
  fit <- random.effects.model(rnai.screen)

  saveRDS(fit, lmm.outfile)
}

run()
