#!/usr/bin/env Rscript

library(dplyr)
library(dtplyr)
library(tidyr)
library(perturbatr)
library(ggplot2)
library(uuid)
library(xtable)

show.xtable <- function(tabl, caption, label, out.file)
{
  align <- rep("l", ncol(tabl) + 1)
  bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
  print.xtable(xtable(tabl, align=align, digits=5,
                      label=label, caption=caption),
              booktabs=T,
              sanitize.colnames.function=bold,
              sanitize.subheadings.function = bold,
              include.rownames = FALSE,  comment=FALSE,
              file=out.file)
}

# this is the compilation of the data-sets together
print.merged.tables <- function(lmm, diffu, path, out.dir)
{

  lmm.p <- lmm$fit@.gene.effects
  if (all(is.na(lmm.p$Qval))) lmm.p$Qval <- 0.2

  diff.hits <- dplyr::select(diffu@.data, GeneSymbol, DiffusionEffect) %>%
    .[order(-DiffusionEffect)]
  diff.hits <- dplyr::left_join(diff.hits,
                                dplyr::select(lmm.p, GeneSymbol, Qval),
                                by="GeneSymbol")

  gps <- lmm$fit@.gene.pathogen.effects %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)
  full.table <- diff.hits %>%
    .[order(-abs(DiffusionEffect))] %>%
    dplyr::mutate(RankNetworkAnalysis=seq(.N)) %>%
    dplyr::select(-DiffusionEffect)

  lmm.not.in <-
    lmm.p[!(lmm.p$GeneSymbol %in% full.table$GeneSymbol)] %>%
    dplyr::select(GeneSymbol, Qval) %>% dplyr::mutate(RankNetworkAnalysis = NA_integer_)

  merged.table <- rbindlist(list(lmm.not.in, full.table)) %>%
    dplyr::mutate(o=seq(.N)) %>%
    dplyr::left_join(dplyr::select(lmm$fit@.gene.effects, GeneSymbol, Effect),
                                   by="GeneSymbol") %>%
    dplyr::left_join(gps, by="GeneSymbol") %>%
    .[order(o)] %>% dplyr::select(-o)

  setDT(merged.table)[Effect == 0, Effect := NA_real_]

  merged.first.25 <-
    merged.table %>%
    dplyr::filter(!is.na(RankNetworkAnalysis)) %>%
    dplyr::rename(Rank=RankNetworkAnalysis) %>%
    dplyr::select(-Qval) %>%
    .[1:25]

  out.file <- paste(out.dir, "two_stage_final_hit_list",  sep="/")
  show.xtable(merged.first.25, "caption", "label", paste0(out.file, ".txt"))
  saveRDS(merged.first.25, paste0(out.file, ".rds"))
}

plot.subgraphs <- function(diffu, out.dir)
{
  genes <- c("elk1", "dyrk1b", "pkn3", "ssx2ip")
  for (gene in genes)
  {
    idx <- which(V(diffu@.graph)$name == gene)
    nei <- as.integer(neighbors(diffu@.graph, gene))
    subgraph <- igraph::induced_subgraph(diffu@.graph, c(idx, nei))
    V(subgraph)$color <- "gray40"
    V(subgraph)[gene]$color <- "darkblue"
    V(subgraph)$size <- 1.25
    V(subgraph)[gene]$size <- 3.25

    setEPS()
    postscript(paste0(out.dir, "/graph_plot-", gene, ".eps"))
    plot(subgraph, layout=layout.kamada.kawai, edge.curved = -.05,
              edge.width=2, edge.color="grey",
              vertex.label.family = "Helvetica", vertex.label.font=2,
              vertex.label.color=V(subgraph)$color,
              vertex.label.cex=V(subgraph)$size,
              asp=0, rescale = TRUE, vertex.shape="none")
    dev.off()
  }
}

run <- function()
{
  path     <- "./data"
  out.dir  <- "./plots"
  lmm      <- readRDS(paste(path, "lmm_fit.rds", sep="/"))
  diffu    <- readRDS(paste(path, "diffusion.rds", sep="/"))

  plot.subgraphs(diffu, out.dir)
  print.merged.tables(lmm, diffu, path, out.dir)
}

run()
