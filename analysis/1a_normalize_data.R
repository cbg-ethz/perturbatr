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
library(readr)
library(dplyr)
library(perturbatr)

source("_preprocess.R")


normalize.chikv.kinome <- function(rnai.screen.raw)
{
  chikv.raw  <- dplyr::filter(rnai.screen.raw, Condition == "CHIKV")
  chikv.norm <- preprocess(
    chikv.raw,
    normalize=c("log",  "background", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    background.column=12,
    drop=T)

  chikv.norm
}


normalize.sars.kinome <- function(rnai.screen.raw)
{
  sars.raw  <- dplyr::filter(rnai.screen.raw, Condition == "SARS")
  sars.norm <- preprocess(
    sars.raw,
    normalize=c("log", "background", "robust-z.score"),
    summarization="mean",
    background.column=12,
    z.score.level="plate",
    drop=T)

  sars.norm
}


normalize.hcv.genome <- function(rnai.screen.raw)
{
  hcv.g.raw  <- dplyr::filter(rnai.screen.raw, Condition=="HCV",
                              Screen=="Genome")
  hcv.g.norm <- preprocess(
    hcv.g.raw,
    normalize=c( "log", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T
  )

  hcv.g.norm
}


normalize.denv.genome <- function(rnai.screen.raw)
{
  denv.g.raw  <- dplyr::filter(rnai.screen.raw, Condition=="DENV",
                               Screen=="Genome")
  denv.g.norm <- preprocess(
    denv.g.raw,
    normalize=c("log", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T)

  denv.g.norm
}


normalize.hcv.kinome <-  function(rnai.screen.raw)
{
  hcv.k.raw  <- dplyr::filter(rnai.screen.raw, Condition=="HCV",
                              Screen=="Kinome")
  hcv.k.norm <- preprocess(
    hcv.k.raw,
    normalize=c("log","loess", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells=c(.05, .95))

  hcv.k.norm
}


normalize.denv.kinome <- function(rnai.screen.raw)
{
  denv.k.raw  <- dplyr::filter(rnai.screen.raw,
                               Condition=="DENV", Screen=="Kinome")
  denv.k.norm <- preprocess(
    denv.k.raw,
    normalize=c("log", "loess", "b.score","robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells=c(.05, .95))

  denv.k.norm
}


normalize <- function(rnai.screen.raw)
{
  chikv.kinome.norm  <- normalize.chikv.kinome(rnai.screen.raw)
  sars.kinome.norm   <- normalize.sars.kinome(rnai.screen.raw)

  hcv.genome.norm    <- normalize.hcv.genome(rnai.screen.raw)
  denv.genome.norm   <- normalize.denv.genome(rnai.screen.raw)

  hcv.kinome.norm    <- normalize.hcv.kinome(rnai.screen.raw)
  denv.kinome.norm   <- normalize.denv.kinome(rnai.screen.raw)

  rnai.screen <- rbind(sars.kinome.norm,
                       chikv.kinome.norm,
                       hcv.genome.norm,
                       denv.genome.norm,
                       hcv.kinome.norm,
                       denv.kinome.norm)

  rnai.screen
}

run <- function()
{

  file.dir <- "./data"
  infile   <- paste0(file.dir, "/", "rnai_screen_raw", ".tsv")
  outfile.rds  <- paste0(file.dir, "/", "rnai_screen_normalized_2", ".rds")
  outfile.tsv  <- paste0(file.dir, "/", "rnai_screen_normalized_2", ".tsv")

  rnai.screen.raw <- readr::read_tsv(infile)
  rnai.screen     <- normalize(rnai.screen.raw)

  saveRDS(rnai.screen, outfile.rds)
  readr::write_tsv(x=rnai.screen, path=outfile.tsv)
}

run()
