#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(lme4)
library(argparse)
library(perturbatr)
library(ggplot2)


.create.noiseless.data <- function(rep.cnt, virus.cnt, mean)
{
  genes.cnt     <- 100
  screens.cnt   <- 5

  viruses       <- paste0("V", 1:virus.cnt)
  virus.effects <- rnorm(virus.cnt, 0, 1)
  if (virus.cnt == 2) { virus.effects[1] <- virus.effects[2] * -1 }
  virus.table <- tibble(Virus=viruses, VirusEffect=virus.effects)

  genes <- c(paste0("G", 1:(genes.cnt-10)),
             paste0("G-", 1:5),
             paste0("G+", 1:5))
  gene.effects.non.hits <- rnorm(genes.cnt - 10, 0, 1)
  gene.effects.neg.hits <- rnorm(5, -mean, 1)
  gene.effects.pos.hits <- rnorm(5, mean, 1)

  gene.effects  <- c(gene.effects.non.hits, gene.effects.neg.hits, gene.effects.pos.hits)
  gene.table    <- tibble(GeneSymbol=genes, GeneEffect=gene.effects)

  screens         <- paste0("S", 1:screens.cnt)
  screen.effects  <- rnorm(screens.cnt, 0, 1)
  if (screens.cnt == 2) screen.effects[1] <- screen.effects[2] * -1
  screen.table    <- tibble(ScreenType=screens, ScreenEffect=screen.effects)

  VG  <- paste(sep=":", viruses, rep(genes, each=virus.cnt))
  vg.effects <- rnorm(genes.cnt*virus.cnt, 0, 1)
  vg.table   <-  tibble(VG=VG, VirusGeneEffect=vg.effects)

  VS  <- paste(sep=":", viruses, rep(screens, each=virus.cnt))
  vs.effects <- rnorm(virus.cnt*screens.cnt, 0, 1)
  vs.table   <- tibble(VS=VS, VirusScreenEffect=vs.effects)

  # combine effect names and add the respective sum of effects
  effect.data <- as.tibble(expand.grid(genes, screens, viruses)) %>%
    dplyr::rename(GeneSymbol=Var1, ScreenType=Var2, Virus=Var3) %>%
    dplyr::group_by(GeneSymbol, ScreenType, Virus) %>%
    dplyr::mutate(VG=paste(Virus, GeneSymbol, sep=":"),
                  VS=paste(Virus, ScreenType, sep=":")) %>%
    ungroup()

  effect.data <- dplyr::left_join(effect.data, virus.table, by="Virus")
  effect.data <- dplyr::left_join(effect.data, gene.table, by="GeneSymbol")
  effect.data <- dplyr::left_join(effect.data, screen.table, by="ScreenType")
  effect.data <- dplyr::left_join(effect.data, vg.table, by="VG")
  effect.data <- dplyr::left_join(effect.data, vs.table, by="VS")

  effect.data <- effect.data %>%
    dplyr::group_by(VS, VG, ScreenType, GeneSymbol, Virus) %>%
    dplyr::mutate(Effect=sum(VirusEffect, GeneEffect,
                             ScreenEffect, VirusGeneEffect,
                             VirusScreenEffect)) %>%
    ungroup()

  # repeat every line x times to have replicates
  effect.data <- effect.data[rep(1:nrow(effect.data), each=rep.cnt),]
  effect.data$Weight <- 1

  effect.data
}


.cv.sets <- function(model.data)
{
  dat <-
    model.data %>%
    dplyr::group_by(Condition, ScreenType, GeneSymbol) %>%
    { dplyr::mutate(ungroup(.), cnt=n(), grp=group_indices(.)) }

  max.cv <- min(dat$cnt)
  if (max.cv > 10) max.cv <- 10
  grps <- unique(dat$grp)

  res <- purrr::map_dfr(
    grps,
    function (g)
    {
        grp.dat           <- dplyr::filter(dat, grp==g)
        grp.dat$cvtestset <- rep(1:max.cv, length.out=nrow(grp.dat))

        grp.dat
      }
  )

  res            <- res %>% dplyr::select(-cnt, -grp)
  res$GeneSymbol <- as.character(res$GeneSymbol)

  res
}


.analyse <- function(md)
{
  pmm.fit     <- lme4::lmer(Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol),
                            data = md, weights = md$Weight, verbose = F)
  lmm.fit     <- lme4::lmer(Readout ~ Condition + (1 | GeneSymbol) + (1 | Condition:GeneSymbol) + (1 | ScreenType) + (1 | Condition:ScreenType),
                            data = md, weights = md$Weight, verbose = F)

  list(lmm.fit=lmm.fit, pmm.fit=pmm.fit)
}


mse <- function(cv.test, fits)
{
  .mse <- function(fit, data)
  {
    pred <- predict(fit, data) %>% as.tibble()
    mseq <- sum((pred - data$Readout)**2) / nrow(data)
    mseq
  }

  lmm.mse <- .mse(fits$lmm.fit, cv.test)
  pmm.mse <- .mse(fits$pmm.fit, cv.test)

  list(lmm=lmm.mse, pmm=pmm.mse)
}


.ranefs <- function(fits)
{
  list(lmm.fit.effects=lme4::ranef(fits$lmm.fit)[["GeneSymbol"]],
       pmm.fit.effects=lme4::ranef(fits$pmm.fit)[["GeneSymbol"]])
}


.cv.validate <- function(model.data)
{
  cv.sets <- .cv.sets(model.data)
  max.cv  <- max(cv.sets$cvtestset)

  li     <- list()
  models <- list()
  for (i in seq(max.cv))
  {
    cv.train    <- dplyr::filter(cv.sets, cvtestset != i)
    cv.test     <- dplyr::filter(cv.sets, cvtestset == i)
    fits        <- .analyse(cv.train)
    mses        <- mse(cv.test, fits)
    li[[i]]     <- tibble(Sampling="CV", LMM_MSE=mses$lmm, PMM_MSE=mses$pmm)
    models[[i]] <- .ranefs(fits)
  }

  df <- purrr::map_dfr(li, function(e) e) %>%
    dplyr::rename(LMM=LMM_MSE, PMM=PMM_MSE) %>%
    tidyr::gather(Model, MSE, 2:3)

  list(data=df, models=models)
}


benchmark.syn.predictability <- function(out.path, virs.cnt, var, rep.cnt)
{
  cat("Doing synthetic benchmark!\n")

  for (m in c(0.5, 1, 2))
  {
    noiseless.data <- .create.noiseless.data(
      rep.cnt=rep.cnt, virus.cnt=virs.cnt, mean=m)
    bench.list     <- list()
    noisy.data <- dplyr::mutate(
      noiseless.data, Readout = Effect+rnorm(nrow(noiseless.data), 0, var)) %>%
      dplyr::rename(Condition = Virus)

    for (vir.cnt in seq(2, virs.cnt))
    {
      viruses     <- paste0("V", 1:vir.cnt)
      model.data  <- dplyr::filter(noisy.data, Condition %in% viruses)

      el <- paste0("var:", var, ",r:",rep.cnt, ",vir:",vir.cnt)
      if (rep.cnt >= 7)
      {
        df.cv       <- .cv.validate(model.data)
        sse         <- rbind(df.cv$dat)
        bench.list[[ paste0(el,"bt") ]] <-
          list(Vir=vir.cnt, Rep=rep.cnt, Var=var, sse=sse,
               cv.model=df.cv$models, data=model.data)
        bench.list[[ paste0(el,"full") ]] <-
          list(Vir=vir.cnt, Rep=rep.cnt, Var=var,
               models=.ranefs(.analyse(model.data)), data=model.data)
      }
      else
      {
        bench.list[[el]] <- list(Vir=vir.cnt, Rep=rep.cnt, Var=var,
                                 models=.ranefs(analyse(model.data)),
                                 data=model.data)
      }
    }

    out.file <- paste0(out.path, "/", "lmm_predictability",
                        "_synthetic_data",
                        "_viruscnt_", virs.cnt,
                        "_repcnt_", rep.cnt,
                        "_var_", var,
                        "_genemean_", m,
                        ".rds")

    saveRDS(list(benchmark=bench.list, data=noisy.data), out.file)
  }
}


benchmark.bio.predictability <- function(model.data, out.path)
{
  cat("Doing full bio!\n")
  df.cv       <- .cv.validate(model.data)
  sse         <- rbind(df.cv$data)
  bench.list  <- list(full=list(sse=sse, cv.gene.effects=df.cv$models))

  vrs <- c("HCV", "DENV", "CHIKV", "SARS")
  for (idx in seq(2, length(vrs)))
  {
    cat(paste("Doing bio on", idx, "\n"))
    rnai.screen.sample <- dplyr::filter(model.data, Condition %in% vrs[1:idx])
    s <- paste0(vrs[1:idx], collapse="_")
    df.cv           <- .cv.validate(rnai.screen.sample)
    sse             <- rbind(df.cv$data)
    bench.list[[s]] <- list(Virus=s,
                            sse=sse,
                            cv.model=df.cv$models,
                            data=model.data)
  }

  out.file <- paste0(out.path, "/",
                     "lmm_predictability_bio_data",
                     ".rds")

  saveRDS(bench.list, out.file)
}


run <- function()
{
  data.dir <- "./data"
  out.dir  <-  "./data"

  option_list <- list(
    make_option(c("-v", "--virus"), action="store",
                help="virus count", type="integer"),
    make_option(c("-r", "--replicate"), action="store",
                help="replicate count", type="integer"),
    make_option(c("-s", "--sig"), action="store",
                help="standard deviation", type="integer"))
  opt_parser <- ArgumentParser(option_list=option_list)
  opt        <- opt_parser$parse_args(opt_parser)

  if (is.null(opt$virus) || is.null(opt$sig) || is.null(opt$replicate))
  {
    stop("Please provide correct arguments")
  }

  rna.file    <- paste(data.dir, "rnai_screen_normalized.rds", sep="/")
  rnai.screen <- readRDS(rna.file)

  model.data  <- dplyr::filter(rnai.screen, Condition != "CVB")
  weights     <- rep(1, nrow(model.data))
  weights[model.data$Design == "pooled"] <- 1.5
  model.data <- methods::as(rnai.screen, "PerturbationData")
  model.data <- perturbatr:::setModelData(model.data, drop=T, weights=weights)

  benchmark.bio.predictability(model.data, out.dir)
  benchmark.syn.predictability(out.dir, opt$virus, opt$sig, opt$replicate)
}

run()
