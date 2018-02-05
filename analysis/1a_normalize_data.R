#!/usr/bin/env Rscript

library(dtplyr)
library(dplyr)
library(data.table)
library(knockdown)

normalize.chikv.kinome <- function(rnai.screen.raw)
{
  chikv.raw  <- knockdown::filter(rnai.screen.raw, Virus=="CHIKV")
  chikv.norm <- knockdown::preprocess(
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
  sars.raw  <- filter(rnai.screen.raw, Virus=="SARS")
  sars.norm <- knockdown::preprocess(
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
  hcv.g.raw        <- knockdown::filter(rnai.screen.raw, Virus=="HCV", Screen=="Genome")
  hcv.g.norm       <- knockdown::preprocess(
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
  denv.g.raw        <- knockdown::filter(rnai.screen.raw, Virus=="DENV", Screen=="Genome")
  denv.g.norm       <- knockdown::preprocess(
    denv.g.raw,
    normalize=c( "log",  "b.score",  "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T)

  denv.g.norm
}

normalize.hcv.kinome <-  function(rnai.screen.raw)
{
  hcv.k.raw        <- knockdown::filter(rnai.screen.raw, Virus=="HCV", Screen=="Kinome")
  hcv.k.norm       <- knockdown::preprocess(
    hcv.k.raw,
    normalize=c("log","loess", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells="q",
    outlier.well.range=c(.05, .95))

  hcv.k.norm
}

normalize.denv.kinome <- function(rnai.screen.raw)
{
  denv.k.raw        <- knockdown::filter(rnai.screen.raw, Virus=="DENV", Screen=="Kinome")
  denv.k.norm       <- knockdown::preprocess(
    denv.k.raw,
    normalize=c("log", "loess", "b.score","robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells="q",
    outlier.well.range=c(.05, .95))

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

  path    <- "./"
  outdir <- "./data"
  outfile <- paste0(outdir, "/", "rnai_screen_normalized", ".rds")

  rnai.screen.raw <- readRDS(
    paste(path, "data/rnai_screen_raw.rds", sep="/"))

  rnai.screen <- normalize(rnai.screen.raw)
  saveRDS(rnai.screen, outfile)
}

run()
