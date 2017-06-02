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


#' Plot an \code{knockout.lmm.analysed} object
#'
#' The method returns three individual plots as a list.
#'
#' @export
#' @method plot knockout.lmm.analysed
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
#' @param x  the object to be plotted
#' @param ...  additional parameters
plot.knockout.lmm.analysed <- function(x, ...)
{
  pl <- .plot.knockout.lmm.analysed (x@.gene.hits, main="Gene effects", ...)
  pl2 <-
    .plot.knockout.lmm.analysed(x@.gene.pathogen.hits, main="") +
    ggplot2::facet_wrap(
      ~ Virus,
      ncol=ceiling(length(unique(x@.gene.pathogen.hits$Virus))/2))
  pl3 <- .plot.effect.matrices.knockout.analysed.lmm(x)
  return(list(pl, pl2, pl3))
}

#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.plot.knockout.lmm.analysed  <- function(x, main, ...)
{
  pars <- list(...)
  size <- ifelse(methods::hasArg(size), pars$size, 10)
  if ("Virus" %in% colnames(x))
  {
    x <- dplyr::filter(x, Control == 0) %>%
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
  LDcolors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pl <-
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(ggplot2::aes(x=GeneSymbol,y=abs(Effect), fill=Effect),
                      stat="identity") +
    ggplot2::scale_fill_gradient2(low=LDcolors[1], high=LDcolors[11],
                                  na.value=LDcolors[6],
                                  name="Effect") +
    ggplot2::ylab("Gene effect strength") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=element_blank(),
                   axis.ticks=element_blank(),
                   text = element_text(size = 20, family = "Helvetica"),
                   axis.text.x = element_text(angle=45, size = 15, family = "Helvetica"),
                   strip.text=element_text(face=x$font))+
    ggplot2::coord_polar() +
    ggplot2::ggtitle(main)

  pl
}

#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
.plot.effect.matrices.knockout.analysed.lmm <- function(x, ...)
{

  effect.matrices <- .effect.matrices(x)
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
                                  name     = "Gene virus\neffect") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 8, family = "Helvetica"),
                   aspect.ratio = 2,
                   axis.text.x  = ggplot2::element_text(angle=45,  hjust = 1, size=9),
                   axis.text.y  = ggplot2::element_text(size=9),
                   axis.title   = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank())
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
