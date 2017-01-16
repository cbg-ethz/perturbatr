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
  gen.pat <- x$gene.pathogen.hits
  pl <- .plot.svd.prioritized.pmm(x$gene.hits, main="Gene effects", ...)
  pl2 <-
    .plot.svd.prioritized.pmm(gen.pat, main="Gene-virus effects", ...) +
    ggplot2::facet_wrap( ~ Virus, ncol=ceiling(length(unique(gen.pat$Virus))/2))
  return(list(pl, pl2))
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.plot.svd.prioritized.pmm <- function(x, main, ...)
{
  pars <- list(...)
  size <- ifelse(methods::hasArg(size), pars$size, 10)
  if ("Virus" %in% colnames(x))
  {
    x <- x %>%
      .[order(abs(Effect), decreasing=T), .SD[1:25], by=Virus] %>%
      dplyr::filter(!is.na(GeneSymbol))
  }
  else
  {
    x <- x[order(abs(Effect), decreasing=T), .SD[1:25]] %>%
      dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  }
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
                   axis.ticks=element_blank(),
                   text = element_text(size = size, family = "Helvetica")) +
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
show.pathogen.hit.counts <- function(x, ...)
{
  UseMethod("show.pathogen.hit.counts")
}

#' @export
#' @importFrom dplyr group_by summarize mutate
#' @method show.pathogen.hit.counts svd.prioritized.pmm
show.pathogen.hit.counts.svd.prioritized.pmm <- function(x, ...)
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
show.effect.matrix <- function(x, ...)
{
  UseMethod("show.effect.matrix")
}

#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
#' @method show.effect.matrix svd.prioritized.pmm
show.effect.matrix.svd.prioritized.pmm <- function(x, ...)
{
  # TODO: change for other prioritization than abs (see: svd.prioritize.pmm)
  effect.matrices <- effect.matrices(x)
  v <- effect.matrices$gene.effects %>% .[order(GeneSymbol, decreasing=T)]
  v$GeneSymbol <- factor(v$GeneSymbol, levels=v$GeneSymbol)
  v$Pathogen <- "All pathogens"

  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pl1 <-
    ggplot2::ggplot(v, aes(Pathogen, GeneSymbol)) +
    ggplot2::geom_tile(aes(fill=Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene effect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 14, family = "Helvetica"),
                   aspect.ratio=4,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=10),
                   axis.text.y=element_text(size=10),
                   axis.title=element_blank(), axis.ticks=element_blank())
  v <-  effect.matrices$gene.pathogen.hits %>% gather(GeneSymbol)
  colnames(v) <- c("GeneSymbol", "Pathogen", "Effect")
  v$GeneSymbol <- factor(v$GeneSymbol, levels=rev(unique(v$GeneSymbol)))
  pl2 <-
    ggplot2::ggplot(v, aes(GeneSymbol, Pathogen)) +
    ggplot2::geom_tile(aes(fill = Effect), colour=LDcolors[1]) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Gene-pathogen\neffect") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 14, family = "Helvetica"),
                   aspect.ratio=.5,
                   axis.text.x=element_text(angle=45,  hjust = 1, size=10),
                   axis.text.y=element_text(size=10),
                   axis.title=element_blank(),
                   axis.ticks=element_blank())
  return(list(gene.effects.plot=pl1, gene.pathogen.effects.plot=pl2))
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
