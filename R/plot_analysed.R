# perturbR: analysis of high-throughput gene perturbation screens
#
# Copyright (C) 2018 Simon Dirmeier
#
# This file is part of perturbR
#
# perturbR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# perturbR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with perturbR. If not, see <http://www.gnu.org/licenses/>.


#' @include class_data.R
#' @include class_analysed.R
#' @include util_enums.R

#' @noRd
#' @export
#' @method plot perturbation.hm.analysed
#' @import ggplot2
#' @import data.table
#' @importFrom dplyr filter
#'
#' @param x  the object to be plotted
#' @param size  size of letters
#' @param ...  additional parameters
#'
#' @return returns a list of plots
#'
plot.perturbation.hm.analysed <- function(x, size=10, main="", ...)
{
  pl <- try({
    .plot.perturbation.hm.analysed(
      x@.gene.hits, main=main, size=size)
  })
  pl2 <- try({
    .plot.perturbation.hm.analysed(x@.nested.gene.hits, main="", size=size) +
      ggplot2::facet_wrap( ~ Condition,
        ncol=ceiling(length(unique(x@.nested.gene.hits$Condition))/2))
  })
  pl3 <- try({
    .plot.effect.matrices.perturbation.analysed.hm(x, size)
  })

  return(list(gene.effect.barplot        = pl,
              nested.gene.effect.barplot = pl2,
              nested.gene.effect.matrix  = pl3))
}


#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom methods hasArg
.plot.perturbation.hm.analysed  <- function(x, main, size, ...)
{
  pars <- list(...)
  if ("Condition" %in% colnames(x))
  {
    x <- dplyr::filter(x, Control == 0) %>%
      .[order(abs(Effect), decreasing=TRUE), .SD[1:25], by=Condition] %>%
      dplyr::filter(!is.na(GeneSymbol))
  }
  else
  {
    x <- x[order(abs(Effect), decreasing=TRUE), .SD[1:25]] %>%
      dplyr::filter(!is.na(GeneSymbol), !is.na(Effect))
  }

  x.pos.range  <- max(abs(x$Effect))
  x.lim        <- c(-x.pos.range, x.pos.range) + c(-x.pos.range, x.pos.range)/5
  x$GeneSymbol <- factor(x$GeneSymbol, levels=rev(unique(x$GeneSymbol)))

  pl <- .plot.bars(x, size, main, ...)
  pl
}


#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr gather
.plot.effect.matrices.perturbation.analysed.hm <- function(x, size, ...)
{

  effect.matrices <- .effect.matrices(x)
  ge <- effect.matrices$gene.effects %>%
    .[order(-abs(Effect))]  %>%
    .[1:25]

  gpe <-  effect.matrices$nested.gene.effects %>%
    dplyr::filter(GeneSymbol %in% ge$GeneSymbol) %>%
    tidyr::gather(GeneSymbol)

  colnames(gpe) <- c("GeneSymbol", "Condition", "Effect")
  gpe$GeneSymbol <- factor(gpe$GeneSymbol, levels=rev(unique(gpe$GeneSymbol)))

  pl <-
    ggplot2::ggplot(gpe, ggplot2::aes(GeneSymbol, Condition)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Effect), colour="black") +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_gradient2(low      = .colors()$blue,
                                  high     = .colors()$red,
                                  na.value = "white",
                                  name     = "Nested gene effect") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(text=ggplot2::element_text(size = size, family = "Helvetica"),
                   aspect.ratio = 2,
                   axis.text.x  = ggplot2::element_text(angle=45,
                                                        hjust = 1,
                                                        size=size),
                   axis.text.y  = ggplot2::element_text(size=size),
                   axis.title   = ggplot2::element_blank(),
                   axis.ticks   = ggplot2::element_blank())
  pl
}


#' Plot a \code{perturbation.hyper.analysed} object
#'
#' @export
#' @import data.table
#' @method plot perturbation.hyper.analysed
#'
#' @param x  the object to plot
#' @param size  size of the text
#' @param ...  additional parameters
#'
#' @return returns a plot object
plot.perturbation.hyper.analysed <- function(x, size=10, ...)
{
  .plot.perturbation.analysed(x, ...)
}


#' Plot a \code{perturbation.tstatistic.analysed} object
#'
#' @export
#' @import data.table
#' @method plot perturbation.tstatistic.analysed
#'
#' @param x  the object to plot
#' @param size  size of the text
#' @param ...  additional parameters
#'
#' @return returns a plot object
plot.perturbation.tstatistic.analysed <- function(x, size=10, ...)
{
  .plot.perturbation.analysed(x, size, ...)
}


#' @noRd
#' @import data.table
#' @import ggplot2
#' @importFrom dplyr group_by summarize mutate filter
.plot.perturbation.analysed <- function(x, size, main="", ...)
{
  df <- x@.gene.hits[order(abs(MeanEffect), decreasing=TRUE), .SD[1:25]] %>%
    dplyr::filter(!is.na(GeneSymbol), !is.na(MeanEffect)) %>%
    dplyr::rename(Effect=MeanEffect)

  pl <- .plot.bars(df, size, main, ...)
  pl
}


.plot.bars <- function(x, size, main, ...)
{
    ggplot2::ggplot(x) +
    ggplot2::geom_bar(ggplot2::aes(x=GeneSymbol, y=abs(Effect),
                                  fill=factor(sign(Effect))),
    									stat="identity") +
    ggplot2::scale_fill_manual("Trend",
      values = c(.colors()$red, "grey", .colors()$blue),
      limits = c("1",  "0", "-1"),
      labels=c("Positive", "None", "Negative")) +
    ggplot2::scale_x_discrete(labels = rev(sort(x$GeneSymbol)),
                              limits = rev(sort(x$GeneSymbol))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = size - 2,
                                                       family = "Helvetica"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   text = ggplot2::element_text(size = size, family = "Helvetica"),
                   axis.text.x = ggplot2::element_text(
                                              size = size - 2,
                                              family = "Helvetica"),
                   strip.text=ggplot2::element_text(face=x$font)) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(main)

}
