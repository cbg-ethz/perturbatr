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
library(argparse)
library(igraph)
library(mvtnorm)


.create.noiseless.data <- function(rep.cnt, virus.cnt, genes.cnt, gene.vcov)
{

  screens.cnt    <- 5

  viruses      <- paste0("V", 1:virus.cnt)
  virus.effects <- rnorm(virus.cnt, 0, 1)
  if (virus.cnt == 2) { virus.effects[1] <- virus.effects[2] * -1 }
  virus.table <- tibble(Virus=as.character(viruses), VirusEffect=virus.effects)

  genes         <- rownames(gene.vcov)
  gene.effects  <- as.vector(mvtnorm::rmvnorm(1, sigma=gene.vcov))
  gene.table    <- tibble(GeneSymbol=as.character(genes), GeneEffect=gene.effects)

  screens         <- paste0("S", 1:screens.cnt)
  screen.effects  <- rnorm(screens.cnt, 0, 1)
  if (screens.cnt == 2) screen.effects[1] <- screen.effects[2] * -1
  screen.table    <- tibble(ScreenType=as.character(screens), ScreenEffect=screen.effects)

  VG  <- paste(sep=":", viruses, rep(genes, each=virus.cnt))
  vg.effects <- rnorm(genes.cnt*virus.cnt, 0, 1)
  vg.table <-  tibble(VG=as.character(VG), VirusGeneEffect=vg.effects)

  VS  <- paste(sep=":", viruses, rep(screens, each=virus.cnt))
  vs.effects <- rnorm(virus.cnt*screens.cnt, 0, 1)
  vs.table <-  tibble(VS=as.character(VS), VirusScreenEffect=vs.effects)

  # combine effect names and add the respective sum of effects
  effect.data <- as.tibble(expand.grid(genes, screens, viruses)) %>%
    dplyr::rename(GeneSymbol=Var1, ScreenType=Var2, Virus=Var3) %>%
    dplyr::mutate(GeneSymbol=as.character(GeneSymbol),
                  ScreenType=as.character(ScreenType),
                  Virus=as.character(Virus)) %>%
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
    dplyr::mutate(Effect=sum(VirusEffect, GeneEffect, ScreenEffect,
                             VirusGeneEffect, VirusScreenEffect)) %>%
    ungroup()

  # repeat every line x times to have replicates
  effect.data <- effect.data[rep(1:nrow(effect.data), each=rep.cnt), ]
  effect.data$Weight <- 1
  effect.data$Control <- 0

  effect.data
}


.effects <- function(ran)
{
  tibble(GeneSymbol=rownames(ran),Effect=ran[,1])
}


.graph <- function(n)
{
  repeat {
    tryCatch({
      adj <- sapply(1:n, function(e) rnorm(n*10, 0, 1))
      g   <- barabasi.game(n, directed=F)
      edj <- igraph::as_edgelist(g)
      for (i in nrow(edj))
      {
        adj[,edj[i, ]] <- adj[,edj[i, ]] + sign(rnorm(1)) * rbeta(1, 2, 1)
      }
      vcov <- cov(adj)
      invisible(mvtnorm::rmvnorm(1, sigma=vcov))
      break
    }, error=function(r){print(r)}, warning=function(r){print(r)})
  }

  vcov
}


boot <- function(dat, rep.cnt, vir.cnt, v)
{
  bench.list <- list()
  i <- 1
  repeat
  {
    s <- paste0("var:", v, ",vir:", vir.cnt, ",rep:", rep.cnt, ",bootstrap:", i)
    tryCatch({
      subs <- perturbatr:::bootstrap(dat)
      bench.list[[i]] <-  list(Var=v,
                               Rep=rep.cnt,
                               Vir=vir.cnt,
                               Bootstrap=i,
                               model=analyse(subs))

      i <- i + 1
    }, error=function(e)   { print(paste0("Didnt fit ", i, ": ", e)); i <- 10000 })
    if (i >= 10) break
  }

  bench.list
}


analyse <- function(md)
{

  pmm.fit     <- lme4::lmer(Readout ~ Condition +
                              (1 | GeneSymbol) +
                              (1 | Condition:GeneSymbol),
                            data = md,
                            weights = md$Weight,
                            verbose = F)
  lmm.fit     <- lme4::lmer(Readout ~ Condition +
                              (1 | GeneSymbol) +
                              (1 | Condition:GeneSymbol) +
                              (1 | ScreenType) +
                              (1 | Condition:ScreenType),
                            data = md,
                            weights = md$Weight,
                            verbose = F)

  pmm.effects <- .effects(lme4::ranef(pmm.fit)[["GeneSymbol"]])
  lmm.effects <- .effects(lme4::ranef(lmm.fit)[["GeneSymbol"]])

  list(lmm.fit=lmm.effects, pmm.fit=pmm.effects)
}


ranking.stability.sythetic <- function(output.path, virs.cnt,
                                       rep.cnt, var, gene.cnt=100)
{
  cat("Synthetic stability ranking\n")
  graph           <- .graph(gene.cnt)
  rownames(graph) <- colnames(graph) <- paste0("G", 1:gene.cnt)
  noiseless.data  <- .create.noiseless.data(
    rep.cnt, virs.cnt, gene.cnt, graph)  %>%
    dplyr::rename(Condition = Virus)

  noisy.data <- dplyr::mutate(
    noiseless.data,
    Readout= Effect + rnorm(nrow(noiseless.data), 0, var))
  bench.list <- list(
    graph=graph,
    data=list(Var=0, Rep.cnt=rep.cnt, Vir=virs.cnt, Bootstrap=0, model=noiseless.data))

  s <- paste0("var:", var, ",vir:", virs.cnt, ",rep:", rep.cnt, ",bootstrap:0")
  m <- analyse(noisy.data)

  bench.list[[s]] <- list(Var=var, Rep=rep.cnt, Vir=virs.cnt, Bootstrap=0, model=m)
  for (vir.cnt in seq(2, virs.cnt))
  {
    viruses     <- paste0("V", 1:vir.cnt)
    sample.data <- dplyr::filter(noisy.data, Condition %in% viruses)
    bench.list  <- c(bench.list, boot(dat=sample.data,
                                      rep.cnt=rep.cnt,vir.cnt=vir.cnt, v=var))
  }

  data.path   <- paste0(output.path,
                        "/lmm_stability_",
                        "synthetic_data",
                        "_viruscnt_", virs.cnt,
                        "_repcnt_", rep.cnt,
                        "_var_", var,
                        ".rds")
  saveRDS(bench.list, data.path)
}


ranking.stability.bio <- function(model.data, output.path)
{
  cat("Biological stability ranking\n")
  rank.all.data <- analyse(model.data)
  bench.list    <- list(full=rank.all.data)

  vrs <- c("HCV", "DENV", "CHIKV", "SARS")
  for (idx in seq(2, length(vrs)))
  {
    dat <- dplyr::filter(model.data, Condition %in% vrs[1:idx])
    i   <- 1
    run <- 1
    repeat
    {
      cat(paste0("Bootstrap bio: ", i, ",v: ", paste0(vrs[1:idx], collapse="_")), "\n")
      tryCatch({
        rnai.screen.sample <- perturbatr::bootstrap(dat)
        s <- paste0("bootstrap:", i, "virs:", paste0(vrs[1:idx], collapse="_"))
        bench.list[[s]] <- list(Bootstrap=i,
                                Virus= paste0(vrs[1:idx], collapse="_"),
                                analyse(rnai.screen.sample))
        i <- i + 1
      }, error=function(e){ print(paste0("Didnt fit ", i, ": ", e)); i <<- 10000 },
         warning=function(e) { print(paste0("Bootstrap bio warning: ", e));  })
      if (i >= 10) break
      run <- run + 1
    }
  }

  dat.pth <- paste0(output.path, "/",  "lmm_stability_bio_data", ".rds")
  saveRDS(bench.list, dat.pth)
}


run <- function()
{
  data.dir <- "./data"
  out.dir  <- "./data"

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

  rna.file    <- paste(data.dir, "rnai_screen_normalized_2.rds", sep="/")
  rnai.screen <- readRDS(rna.file)

  model.data  <- dplyr::filter(rnai.screen, Condition != "CVB")
  weights     <- rep(1, nrow(model.data))
  weights[model.data$Design == "pooled"] <- 1.5
  model.data <- methods::as(rnai.screen, "PerturbationData")
  model.data <- perturbatr:::setModelData(model.data, drop=T, weights=weights)

  ranking.stability.bio(model.data, out.dir)
  ranking.stability.sythetic(out.dir, opt$virus, opt$replicate, opt$sig)

}

run()
