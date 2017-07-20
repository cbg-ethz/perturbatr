library(data.table)
library(dtplyr)
library(dplyr)
library(knockout)

normalize.chikv.kinome <- function(rnai.screen.raw)
{
  chikv.raw  <- filter(rnai.screen.raw, Virus=="CHIKV")
  chikv.norm <- knockout::preprocess(
    chikv.raw,
    normalize=c("log",  "background", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    background.column=12,
    drop=T
    #rm.cytotoxic="Scrambled", normalize.viability=T
  )
  chikv.norm
}

normalize.sars.kinome <- function(rnai.screen.raw)
{
  sars.raw  <- filter(rnai.screen.raw, Virus=="SARS")
  sars.norm <- knockout::preprocess(
    sars.raw,
    normalize=c("log", "background", "robust-z.score"),
    summarization="mean",
    background.column=12,
    z.score.level="plate",
    drop=T
    #rm.cytotoxic="Scrambled", normalize.viability=T
  )

  sars.norm
}

normalize.hcv.genome <- function(rnai.screen.raw)
{
  hcv.g.raw        <- filter(rnai.screen.raw, Virus=="HCV", Screen=="Genome")
  hcv.g.norm       <- knockout::preprocess(
    hcv.g.raw,
    normalize=c( "log", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T
    #rm.cytotoxic = "Scrambled", normalize.viability=T
  )
  hcv.g.norm
}

normalize.denv.genome <- function(rnai.screen.raw)
{
  denv.g.raw        <- filter(rnai.screen.raw, Virus=="DENV", Screen=="Genome")
  denv.g.norm       <- knockout::preprocess(
    denv.g.raw,
    normalize=c( "log",  "b.score",  "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T
    #rm.cytotoxic = "Scrambled", normalize.viability=T
  )
  denv.g.norm
}

normalize.hcv.kinome <-  function(rnai.screen.raw)
{
  hcv.k.raw        <- filter(rnai.screen.raw, Virus=="HCV", Screen=="Kinome")
  hcv.k.norm       <- knockout::preprocess(
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
  denv.k.raw        <- filter(rnai.screen.raw, Virus=="DENV", Screen=="Kinome")
  denv.k.norm       <- knockout::preprocess(
    denv.k.raw,
    normalize=c("log", "loess", "b.score","robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells="q",
    outlier.well.range=c(.05, .95))
  denv.k.norm
}

normalize.cvb.genome <- function(rnai.screen.raw)
{
  cvb.raw        <- filter(rnai.screen.raw, Virus=="CVB")
  cvb.norm       <- knockout::preprocess(
    cvb.raw, normalize=c("log", "loess", "b.score", "robust-z.score"),
    summarization="mean",
    z.score.level="plate",
    drop=T,
    rm.outlier.wells="q",
    outlier.well.range=c(.05, .95))
  cvb.norm
}

normalize <- function(rnai.screen.raw)
{


  chikv.kinome.norm  <- normalize.chikv.kinome(rnai.screen.raw)
  sars.kinome.norm   <- normalize.sars.kinome(rnai.screen.raw)
  hcv.genome.norm    <- normalize.hcv.genome(rnai.screen.raw)
  denv.genome.norm   <- normalize.denv.genome(rnai.screen.raw)
  hcv.kinome.norm    <- normalize.hcv.kinome(rnai.screen.raw)
  denv.kinome.norm   <- normalize.denv.kinome(rnai.screen.raw)
  cvb.genome.norm    <- normalize.cvb.genome(rnai.screen.raw)


  rnai.screen <- rbind(sars.kinome.norm,
                       chikv.kinome.norm,
                       hcv.genome.norm,
                       denv.genome.norm,
                       hcv.kinome.norm,
                       denv.kinome.norm,
                       cvb.genome.norm)

  out.file <- "/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline/integrated_data_files/rnai_screen_normalized.rds"
  message(paste0("Wrinting normalized data to ", out.file))
  saveRDS(rnai.screen, out.file)

  rnai.screen
}

.load.raw <- function()
{
  setwd("/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline")
  rnai.screen.raw <- readRDS("/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline/integrated_data_files/rnai_screen_raw.rds")
}

.load.normalized <- function()
{
  setwd("/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline")
  rnai.screen <- readRDS("~/PROJECTS/sysvirdrug_project/data/svd/integrated_data_files/rnai_screen_normalized.rds")
}

.load.validation <- function()
{
  setwd("/Users/simondi/PROJECTS/sysvirdrug_project/src/util/knockout_svd_pipeline")
  validation.screen <- readRDS("integrated_data_files/validation_screen_norm.rds")
}
