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


library(dplyr)
library(tibble)
library(perturbatr)
library(purrr)


network.propagation <- function(lmm.fit, graph.file)
{
  ret <- purrr:::map_dfr(seq(0.2, 0.5, by=0.1), function(.)
  {
    diffusion <- perturbatr::diffuse(
      lmm.fit, path=graph.file, delete.nodes.on.degree=1, r=.,
      take.largest.component=TRUE)

    ge <- geneEffects(diffusion)
    ge$restart <- as.character(.)

    ge
  })

  ret
}


run <- function()
{
  data.dir   <- out.dir <- "./data"
  graph.file <- paste(data.dir, "graph_file.tsv", sep="/")

  lmm.fit <- readRDS(paste(data.dir, "lmm_fit.rds", sep="/"))
  lmm.fit.full <- lmm.fit$fit
  lmm.fit.full@geneHits <- lmm.fit.full@geneEffects

  diff.full <- network.propagation(lmm.fit.full, graph.file)

  out.file <- paste0(out.dir, "/", "diffusion_2.rds")
  saveRDS(diff.full, out.file)
}

run()
