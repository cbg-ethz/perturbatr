#!/usr/bin/env Rscript

library(dplyr)
library(dtplyr)
library(data.table)
library(perturbatr)

network.propagation <- function(lmm, graph.file)
{

  diffusion <- perturbatr::diffuse(
    lmm, path=graph.file, delete.nodes.on.degree=1, r=0.20)

  diffusion
}

run <- function()
{
  path       <- "./data"
  out.dir    <- "./data"
  graph.file <- paste(path, "graph_file.tsv", sep="/")

  lmm.fit <- readRDS(paste(path, "lmm_fit.rds", sep="/"))
  diff    <- network.propagation(lmm.fit$fit, graph.file)

  out.file <- paste0(out.dir, "/", "diffusion.rds")
  saveRDS(diff, out.file)
}

run()
