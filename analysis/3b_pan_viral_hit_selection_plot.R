#!/usr/bin/env Rscript

library(dtplyr)
library(dplyr)
library(tidyr)
library(lme4)
library(optparse)
library(perturbatr)
library(ggplot2)
library(hashmap)
library(hrbrthemes)


effect.matrices <- function(obj)
{
  g <- obj@.gene.hits %>%
    dplyr::select(GeneSymbol, Effect) %>%
    .[order(-abs(Effect))]
  pg <- obj@.gene.pathogen.effects %>%
    dplyr::select(Virus, GeneSymbol, Effect) %>%
    tidyr::spread(Virus, Effect)

  list(gene.effects=g, gene.pathogen.effects=pg)
}


plot.gene.effects  <- function(x)
{
  if ("Virus" %in% colnames(x))
  {
    x <- dplyr::filter(x, Control == 0) %>%
      .[order(abs(Effect), decreasing=TRUE), .SD[1:25], by=Virus] %>%
      dplyr::filter(!is.na(GeneSymbol))
  }
  else
  {
    x <- x[order(abs(Effect), decreasing=TRUE), .SD[1:25]] %>%
      dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  }
  x.pos.range <- max(abs(x$Effect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  LDcolors     <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(ggplot2::aes(x=GeneSymbol, y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_gradient2(expression(paste("Gene effect ", gamma)),
      low=LDcolors[1],high=LDcolors[11], na.value=LDcolors[6]) +
    ggplot2::scale_x_discrete(labels = rev(sort(x$GeneSymbol)),
                              limits=rev(sort(x$GeneSymbol))) +
    ggplot2::theme_bw() +
    ggplot2::coord_flip() +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    ggplot2::theme(text            = ggplot2::element_text(size = 20),
                   strip.text=ggplot2::element_text(face=x$font),
                   axis.text.x=element_text(size=14),
                   axis.text.y=element_text(size=14),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.text=element_text(size=14),
                   plot.title = element_text(hjust = 0.5),
                   strip.text.x = element_text(size = 14, hjust=.1),
                   strip.text.y = element_text(size = 14),
                   panel.spacing.y = ggplot2::unit(2, "lines"))

  pl
}

plot.gene.virus.effects <- function(x)
{

  effect.matrices <- effect.matrices(x)
  ge <- effect.matrices$gene.effects %>%
    .[order(-abs(Effect))]  %>%
    .[1:25]
  gpe <-  effect.matrices$gene.pathogen.effects %>%
    dplyr::filter(GeneSymbol %in% ge$GeneSymbol) %>%
    tidyr::gather(GeneSymbol)
  colnames(gpe) <- c("GeneSymbol", "Pathogen", "Effect")
  gpe$GeneSymbol <- factor(gpe$GeneSymbol, levels=rev(unique(gpe$GeneSymbol)))
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))

  pl <-
    ggplot2::ggplot(gpe, ggplot2::aes(GeneSymbol, Pathogen)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low      = LDcolors[1],
                                  high     = LDcolors[11],
                                  na.value = LDcolors[6],
                                  name     = expression(paste("Gene virus\neffect"))) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    hrbrthemes::theme_ipsum_rc(base_family="Helvetica") +
    ggplot2::theme(text            = ggplot2::element_text(size = 20),
                   aspect.ratio = 2,
                   axis.text.x=element_blank(),
                   axis.text.y=element_text(size=14),
                   axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   legend.text=element_text(size=14),
                   plot.title = element_text(hjust = 0.5),
                   strip.text.x = element_text(size = 14, hjust=.1),
                   strip.text.y = element_text(size = 14),
                   panel.spacing.y = ggplot2::unit(2, "lines"))


    pl
}

run <- function()
{
  path <- "./"
  out.dir <-  "./plots"

  data.file <- paste(path, "data/lmm_fit.rds", sep="/")
  fit       <- readRDS(data.file)

  pl <- plot.gene.effects(fit$fit@.gene.effects)
  ggsave(
    filename = paste0(out.dir, "/", "gene_effects", ".eps"),
    plot = pl,
    width = 6,
    height = 8)

  pl <- plot.gene.virus.effects(fit$fit)
  ggsave(
    filename = paste0(out.dir, "/", "effect_matrix", ".eps"),
    plot = pl,
    width = 8,
    height = 8)

}

run()
