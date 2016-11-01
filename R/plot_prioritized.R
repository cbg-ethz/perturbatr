# knockout: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2015 - 2016 Simon Dirmeier
#
# This file is part of knockout
#
# knockout is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# knockout is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with knockout. If not, see <http://www.gnu.org/licenses/>.

#' @export
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
#' @method plot svd.prioritized.pmm
plot.svd.prioritized.pmm <- function(x, y, ...)
{
  gen.pat <- x$gene.pathogen.effect.hits
  pl <- .plot.svd.prioritized.pmm(x$gene.effect.hits, main="Gene effects")
  pl2 <-
    .plot.svd.prioritized.pmm(gen.pat, main="Gene-virus effects") +
    ggplot2::facet_wrap(. ~ Virus, ncol=length(unique(gen.pat$Virus))/2)
  pl3 <- .multiplot(plotlist=list(pl, pl2), cols=2)
  pl3
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
.plot.svd.prioritized.pmm <- function(x, main, ...)
{
  x <- x[order(abs(Effect), decreasing=T), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  x.pos.range <- max(abs(x$Effect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank()) +
    ggplot2::coord_polar() +
    ggplot2::ggtitle(main)
  pl
}

#' Plot the hit counts of a PMM analysis for every virus as barplots
#'
#' @export
#' @import data.table
#'
#' @param x  an svd.prioritized.pmm object
#' @param ...  additional parameters
plot.pathogen.hit.counts <- function(x, ...)
{
  UseMethod("plot.pathogen.hit.counts")
}

#' @export
#' @importFrom dplyr group_by summarize mutate
plot.pathogen.hit.counts.svd.prioritized.pmm <- function(x, ...)
{
  obj <- x$gene.pathogen.effect.hits
  single.res <-
    dplyr::group_by(obj, Virus, Sign=sign(Effect)) %>%
    dplyr::summarize(cnt=n()) %>%
    ungroup %>%
    dplyr::mutate(Count=cnt*Sign)

  pl <-
    ggplot2::ggplot(single.res, aes(x=Virus, y=Count, fill=Sign)) +
    ggplot2::geom_bar(stat="identity", position="dodge") +
    ggplot2::scale_fill_distiller(palette="Spectral") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none",
                   text = element_text(size=22),
                   axis.text.y=element_blank()) +
    ggplot2::geom_hline(yintercept=0) +
    ggplot2::ylab("Count hits") +
    ggplot2::geom_text(aes(x=Virus, y=ifelse(Count>0, Count+1, Count-1),
                           label=abs(Count)), size=5, colour="black")
  pl
}

#' Plot the effects matrix of a PMM analysis
#'
#' @export
#' @import data.table
#'
#' @param x  an svd.prioritized.pmm object
#' @param ...  additional parameters
plot.effect.matrix <- function(x, ...)
{
  UseMethod("plot.effect.matrix")
}

#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
#' @method plot svd.prioritized.pmm.single.gene.matrices
plot.effect.matrix.svd.prioritized.pmm <- function(x, ...)
{
  # TODO:
  # plot the cool matrix here from the PMM paper
  # and similar matrices
  # use two columns: lefft column result of gene, right column matrix of gene hits in colors of gene-pathogen hits
  params <- list(...)
  size <- ifelse(hasArg(size), params$size, 14)
  dat <- tidyr::gather(x$cpg.mat, Virus, Effect, 2:5)
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  or <- rev(x$GeneSymbol)
  pl <-
    ggplot2::ggplot(dat, aes(Virus, GeneSymbol)) +
    ggplot2::geom_tile(aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0),
                              limits=or,
                              breaks=or,
                              labels=or) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1],
                                  high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene-pathogen\neffect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 16, family = "Helvetica"),
                   aspect.ratio=2,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=10),
                   axis.text.y=element_text(size=10),
                   axis.title=element_blank(),
                   axis.ticks=element_blank())
  pl
}

#' @export
#' @import data.table
#' @method plot svd.prioritized.hyper
plot.svd.prioritized.hyper <- function(x, y, ...)
{
  .plot.svd.prioritized(x, ...)
}

#' @export
#' @import data.table
#' @method plot svd.prioritized.tt
plot.svd.prioritized.tt <- function(x, y, ...)
{
  .plot.svd.prioritized(x, ...)
}

#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarize mutate filter
.plot.svd.prioritized <- function(x, ...)
{
  x <- x[order(abs(MeanEffect), decreasing=T), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(MeanEffect))
  x.pos.range <- max(abs(x$MeanEffect))
  x.lim  <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(aes(x=GeneSymbol,y=abs(MeanEffect), fill=MeanEffect),
                      stat="identity") +
    ggplot2::scale_fill_distiller(palette="Spectral", limits=x.lim) +
    ylab("Effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank()) +
    ggplot2::coord_polar()
  pl
}
